#!/usr/bin/env python3
"""
For each extracted PICI GBFF region:
  1. Run PFAM hmmscan (optional)
  2. Classify every CDS into a PICI module category
  3. Score element completeness against a canonical PICI definition
  4. Assign architecture class (complete_canonical, packaging_positive_partial, etc.)
  5. Draw a gene schematic with module colour coding
  6. Write per-group annotation TSVs
  7. Write a master summary TSV (one row per element) for downstream analysis
"""

import argparse
import glob
import os
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Patch

import pandas as pd
from Bio import SeqIO

# ──────────────────────────────────────────────────────────────────────────────
# CONSTANTS
# ──────────────────────────────────────────────────────────────────────────────

MAX_PLOT_WIDTH = 60_000  # bp; elements wider than this are truncated in the plot

# Canonical PICI module definitions used for completeness scoring.
# Each key maps to a list of (product_keyword, pfam_keyword) pairs.
# A module is present if ANY pair matches.
PICI_MODULES = {
    "integrase":           ["integrase", "recombinase", "tyrosine recombinase",
                             "serine recombinase", "site-specific recombinase",
                             "site-specific integrase"],
    "regulator":           ["alpa", "merr", "stl", "ci repressor", "cro", "cox",
                             "regulator", "repressor", "hth", "helix-turn-helix",
                             "transcriptional regulator"],
    "replication":         ["replication protein", "replicase", "primase",
                             "replication initiator", "replication factor", "repa",
                             "dna replication", "rep protein"],
    "putative_replication":["aaa"],          # pfam-only, checked separately
    "terminase":           ["terminase", "ters", "large subunit terminase",
                             "small subunit terminase", "dna packaging",
                             "packaging atpase", "packaging protein",
                             "headful packaging"],
    "capsid":              ["capsid", "major capsid", "head protein",
                             "prohead", "portal protein", "portal"],
    "mobilisation":        ["mobilization", "mobilisation", "conjugation",
                             "conjugative transfer", "relaxase",
                             "type iv secretion", "trb", "mob protein"],
}

# Colour scheme — consistent with your existing heatmap palette
COLOR_MAP = {
    "integrase":            "#66c2a5",
    "regulator":            "#fc8d62",
    "replication":          "#8da0cb",
    "putative_replication": "#e41a1c",
    "terminase":            "#a6d854",
    "capsid":               "#e78ac3",
    "mobilisation":         "#ffd92f",
    "other":                "#d9d9d9",
}

# Architecture classification rules (evaluated in order, first match wins).
# Each rule is (label, required_modules, forbidden_modules)
ARCHITECTURE_RULES = [
    ("complete_canonical",
     {"integrase", "regulator", "replication", "terminase", "capsid"},
     set()),
    ("packaging_positive_partial",
     {"integrase", "regulator", "terminase"},
     set()),
    ("packaging_positive_partial",
     {"integrase", "regulator", "capsid"},
     set()),
    ("backbone_replication_no_terminase",
     {"integrase", "regulator", "replication"},
     {"terminase", "capsid"}),
    ("backbone_replication_no_terminase",
     {"integrase", "regulator", "putative_replication"},
     {"terminase", "capsid"}),
    ("mobilisation_associated",
     {"integrase", "mobilisation"},
     set()),
    ("backbone_only",
     {"integrase", "regulator"},
     {"replication", "putative_replication", "terminase", "capsid", "mobilisation"}),
    ("minimal_or_relic",
     {"integrase"},
     {"regulator", "replication", "putative_replication", "terminase", "capsid"}),
    ("partial_unclear",
     set(),
     set()),   # catch-all
]

# ──────────────────────────────────────────────────────────────────────────────
# ARGUMENT PARSING
# ──────────────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("gbff_dir",
                   help="Directory containing group subdirs with *.gbff files. "
                        "Or a flat dir of *.gbff files (treated as one group).")
    p.add_argument("--pfam", default="Pfam-A.hmm",
                   help="Path to Pfam-A.hmm (default: Pfam-A.hmm in CWD)")
    p.add_argument("--no-pfam", action="store_true",
                   help="Skip PFAM scanning entirely (use product names only)")
    p.add_argument("--rerun-pfam", action="store_true",
                   help="Re-run PFAM even if .pfam output already exists")
    p.add_argument("--outdir", default=None,
                   help="Output directory for master summary (default: gbff_dir)")
    p.add_argument("--regraph", action="store_true",
                   help="Skip PFAM and classification. Read existing *_annotations.tsv "
                        "and master_completeness_summary.tsv and regenerate all plots only.")
    return p.parse_args()

# ──────────────────────────────────────────────────────────────────────────────
# PFAM
# ──────────────────────────────────────────────────────────────────────────────

def run_pfamscan(faa: str, domtbl: str, pfam_db: str) -> bool:
    cmd = ["hmmscan", "--domtblout", domtbl, "--noali", "-E", "1e-3", pfam_db, faa]
    try:
        subprocess.run(cmd, check=True, capture_output=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"  [PFAM WARN] hmmscan failed on {faa}: {e.stderr.decode()[:200]}")
        return False


def parse_pfamscan(domtbl: str) -> dict:
    """
    Parse hmmscan --domtblout file.
    Returns dict: protein_id -> (domain_name, i_evalue, description, bitscore)

    domtblout columns (space-separated, 0-indexed after splitting):
      0   target name (domain)
      1   target accession
      2   tlen
      3   query name  (protein id)
      4   query accession
      5   qlen
      6   full-seq E-value
      7   full-seq score
      8   full-seq bias
      9   domain #
      10  of (total domains)
      11  c-Evalue (conditional)
      12  i-Evalue (independent)  <-- use this
      13  domain score
      14  domain bias
      ...
      22+ description
    """
    hits = {}
    if not os.path.exists(domtbl):
        return hits
    with open(domtbl) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 23:
                continue
            domain   = fields[0]
            prot_id  = fields[3].split("|")[-1]
            bitscore = float(fields[13])   # domain bitscore (col 13)
            i_evalue = float(fields[12])   # i-Evalue (col 12)
            desc     = " ".join(fields[22:])
            if bitscore >= 20:
                if prot_id not in hits or bitscore > hits[prot_id][3]:
                    hits[prot_id] = (domain, i_evalue, desc, bitscore)
    return hits

# ──────────────────────────────────────────────────────────────────────────────
# FAA EXTRACTION
# ──────────────────────────────────────────────────────────────────────────────

def extract_faa(gbff: str, out_faa: str):
    with open(out_faa, "w") as fh:
        for rec in SeqIO.parse(gbff, "genbank"):
            for feat in rec.features:
                if feat.type == "CDS" and "translation" in feat.qualifiers:
                    pid = feat.qualifiers.get("protein_id", ["unknown"])[0].split("|")[-1]
                    seq = feat.qualifiers["translation"][0]
                    fh.write(f">{pid}\n{seq}\n")

# ──────────────────────────────────────────────────────────────────────────────
# GENE CLASSIFICATION
# ──────────────────────────────────────────────────────────────────────────────

def classify_gene(product: str, pfam_desc: str = None) -> str:
    """
    Assign a PICI module category to a CDS.

    Priority order (important — don't reorder without understanding consequences):
      1. Mobilisation   — explicit conjugation/relaxase keywords
      2. DNA repair exclusion — glycosylase etc. should never be called replication
      3. Replication    — strong product-name keywords
      4. Terminase      — before capsid to avoid portal being called capsid
      5. Capsid
      6. Integrase
      7. Regulator
      8. Putative replication — AAA domain (pfam only, broad — placed last to avoid FPs)
      9. other
    """
    pl = (product  or "").lower()
    fl = (pfam_desc or "").lower()
    combined = pl + " " + fl

    # 1. Mobilisation
    mob_kws = ["mobiliz", "mobilis", "conjugat", "relaxase",
               "type iv secretion", "trb", "mob protein", "tra "]
    if any(k in combined for k in mob_kws):
        return "mobilisation"

    # 2. Exclude DNA repair (broad AAA / replication keywords can FP here)
    repair_kws = ["glycosylase", "dna repair", "excision repair",
                  "endonuclease", "resolvase", "recombination protein"]
    if any(k in fl for k in repair_kws):
        return "other"

    # 3. Replication (product name takes priority over pfam)
    rep_kws = ["replication protein", "replicase", "primase",
               "replication initiator", "replication factor",
               "rep protein", "repa", "dna replication"]
    if any(k in pl for k in rep_kws):
        return "replication"

    # 4. Terminase
    ter_kws = ["terminase", "ters", "large subunit terminase",
               "small subunit terminase", "dna packaging",
               "packaging atpase", "headful packaging"]
    if any(k in combined for k in ter_kws):
        return "terminase"

    # 5. Capsid
    cap_kws = ["major capsid", "capsid", "head protein",
               "prohead", "portal protein", "portal"]
    if any(k in combined for k in cap_kws):
        return "capsid"

    # 6. Integrase
    int_kws = ["integrase", "recombinase", "tyrosine recombinase",
               "serine recombinase", "site-specific recombinase",
               "site-specific integrase"]
    if any(k in combined for k in int_kws):
        return "integrase"

    # 7. Regulator
    reg_kws = ["alpa", "merr", "stl", "ci repressor", "cro", "cox",
               "regulator", "repressor", "helix-turn-helix", "hth",
               "transcriptional regulator"]
    if any(k in combined for k in reg_kws):
        return "regulator"

    # 8. Putative replication — AAA domain (pfam only, broad)
    if "aaa" in fl:
        return "putative_replication"

    return "other"

# ──────────────────────────────────────────────────────────────────────────────
# COMPLETENESS SCORING
# ──────────────────────────────────────────────────────────────────────────────

CANONICAL_MODULES = ["integrase", "regulator", "replication", "terminase", "capsid"]
ALL_MODULES       = list(PICI_MODULES.keys())


def score_completeness(genes: list) -> dict:
    """
    Given a list of classified gene dicts, return a completeness report.
    Returns dict with:
      - modules_present: set of module names found
      - canonical_score: fraction of canonical modules present (0.0–1.0)
      - canonical_complete: bool (all 5 canonical modules present)
      - architecture: string label
      - module_counts: dict module -> count of CDS assigned to it
    """
    module_counts = {m: 0 for m in ALL_MODULES}
    module_counts["other"] = 0

    for g in genes:
        cls = g["class"]
        if cls in module_counts:
            module_counts[cls] += 1

    modules_present = {m for m in ALL_MODULES if module_counts.get(m, 0) > 0}

    canonical_score = sum(
        1 for m in CANONICAL_MODULES if m in modules_present
    ) / len(CANONICAL_MODULES)

    canonical_complete = canonical_score == 1.0

    # Architecture classification
    architecture = "partial_unclear"
    for label, required, forbidden in ARCHITECTURE_RULES:
        if required.issubset(modules_present) and not forbidden.intersection(modules_present):
            architecture = label
            break

    return {
        "modules_present":   modules_present,
        "canonical_score":   round(canonical_score, 2),
        "canonical_complete": canonical_complete,
        "architecture":      architecture,
        "module_counts":     module_counts,
    }

# ──────────────────────────────────────────────────────────────────────────────
# GBFF PARSING
# ──────────────────────────────────────────────────────────────────────────────

def parse_gbff(filepath: str, pfam_hits: dict) -> list:
    genes = []
    for rec in SeqIO.parse(filepath, "genbank"):
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            start   = int(feat.location.start)
            end     = int(feat.location.end)
            strand  = feat.location.strand or 1
            product = feat.qualifiers.get("product",    ["unknown"])[0]
            pid     = feat.qualifiers.get("protein_id", ["unknown"])[0].split("|")[-1]

            pfam_domain = pfam_evalue = pfam_desc = pfam_bitscore = None
            if pid in pfam_hits:
                pfam_domain, pfam_evalue, pfam_desc, pfam_bitscore = pfam_hits[pid]

            genes.append({
                "Protein_ID":    pid,
                "start":         start,
                "end":           end,
                "strand":        strand,
                "product":       product,
                "pfam_domain":   pfam_domain,
                "pfam_desc":     pfam_desc,
                "pfam_evalue":   pfam_evalue,
                "pfam_bitscore": pfam_bitscore,
                "class":         classify_gene(product, pfam_desc),
            })
    return sorted(genes, key=lambda x: x["start"])

# ──────────────────────────────────────────────────────────────────────────────
# PLOTTING
# ──────────────────────────────────────────────────────────────────────────────

def draw_gene(ax, start, end, y, strand, color):
    width = end - start
    if width <= 0:
        return
    head = min(width * 0.4, max(80, width * 0.15))
    y0, y1, ym = y - 0.35, y + 0.35, y
    if strand >= 0:
        pts = [(start, y0), (end - head, y0), (end, ym),
               (end - head, y1), (start, y1)]
    else:
        pts = [(end, y0), (start + head, y0), (start, ym),
               (start + head, y1), (end, y1)]
    ax.add_patch(Polygon(pts, closed=True,
                         facecolor=color, edgecolor="black", linewidth=0.5))


def _render_schematic(ax, all_records, names, completeness_list):
    """
    Core rendering: draw all elements onto an existing Axes.
    Returns max_len (bp) so the caller can set xlim.
    """
    n = len(all_records)
    max_len = 0
    max_label_chars = max((len(n) for n in names), default=20)

    for i, (genes, name, comp) in enumerate(
            zip(all_records, names, completeness_list)):
        y = n - i
        local_min = min(g["start"] for g in genes)

        for g in genes:
            s = g["start"] - local_min
            e = g["end"]   - local_min
            draw_gene(ax, s, e, y, g["strand"], COLOR_MAP[g["class"]])
            max_len = max(max_len, e)

        arch  = comp["architecture"]
        score = comp["canonical_score"]
        short_name = name if len(name) <= 45 else name[:42] + "…"
        label = f"{short_name}  [{arch} | {score:.0%}]"
        ax.text(-400, y, label, ha="right", va="center", fontsize=7)

    return max_len


def plot_single(genes: list, name: str, comp: dict, output: str):
    """One PNG per element — always readable regardless of dataset size."""
    fig, ax = plt.subplots(figsize=(16, 1.8))

    local_min = min(g["start"] for g in genes)
    max_len = 0
    for g in genes:
        s = g["start"] - local_min
        e = g["end"]   - local_min
        draw_gene(ax, s, e, 1, g["strand"], COLOR_MAP[g["class"]])
        max_len = max(max_len, e)

    arch  = comp["architecture"]
    score = comp["canonical_score"]
    ax.set_title(f"{name}\n{arch}  |  completeness {score:.0%}",
                 fontsize=8, loc="left", pad=4)

    legend_patches = [
        Patch(facecolor=COLOR_MAP[k], edgecolor="black", label=k)
        for k in COLOR_MAP if k != "other"
    ] + [Patch(facecolor=COLOR_MAP["other"], edgecolor="black", label="other")]

    ax.legend(handles=legend_patches, loc="upper right",
              ncol=4, fontsize=7, frameon=False)

    xlim = min(max_len, MAX_PLOT_WIDTH) + 300
    ax.set_xlim(-max(xlim * 0.05, 500), xlim)
    ax.set_ylim(0.4, 1.6)
    ax.axis("off")

    plt.tight_layout()
    plt.savefig(output, dpi=150, bbox_inches="tight")
    plt.close()


def plot_group(all_records: list, names: list, completeness_list: list,
               output: str, title: str = "", batch_size: int = 30):
    """
    Multi-element schematic. If there are more than batch_size elements,
    the plot is split into pages/batches saved as output_batch1.png, etc.
    For <=batch_size elements, saves to output directly.
    """
    n = len(all_records)
    batches = [
        (all_records[i:i+batch_size],
         names[i:i+batch_size],
         completeness_list[i:i+batch_size])
        for i in range(0, n, batch_size)
    ]

    base, ext = os.path.splitext(output)

    for b_idx, (recs, nms, comps) in enumerate(batches):
        bn = len(recs)
        fig_h = max(3.0, 0.5 * bn + 1.8)
        fig, ax = plt.subplots(figsize=(20, fig_h))

        max_len = _render_schematic(ax, recs, nms, comps)

        legend_patches = [
            Patch(facecolor=COLOR_MAP[k], edgecolor="black", label=k)
            for k in COLOR_MAP
        ]
        ax.legend(handles=legend_patches,
                  loc="upper center", bbox_to_anchor=(0.5, 1.0),
                  ncol=min(8, len(COLOR_MAP)), fontsize=8, frameon=False)

        xlim = min(max_len, MAX_PLOT_WIDTH) + 500
        ax.set_xlim(-max(xlim * 0.30, 4000), xlim)
        ax.set_ylim(0.2, bn + 1.5)
        ax.axis("off")

        batch_label = (f"{title} ({b_idx+1}/{len(batches)})"
                       if len(batches) > 1 else title)
        if batch_label:
            ax.set_title(batch_label, fontsize=10, pad=32)

        out = f"{base}_batch{b_idx+1}{ext}" if len(batches) > 1 else output
        plt.tight_layout(rect=[0, 0, 1, 0.92])
        plt.savefig(out, dpi=150, bbox_inches="tight")
        plt.close()
        yield out   # yield each saved path so the caller can report it

# ──────────────────────────────────────────────────────────────────────────────
# OUTPUT WRITING
# ──────────────────────────────────────────────────────────────────────────────

def save_annotation_table(out_path: str, all_records: list, names: list):
    rows = []
    for genes, genome in zip(all_records, names):
        for g in genes:
            rows.append({
                "Genome":        genome,
                "Protein_ID":    g["Protein_ID"],
                "start":         g["start"],
                "end":           g["end"],
                "strand":        g["strand"],
                "product":       g["product"],
                "pfam_domain":   g["pfam_domain"],
                "pfam_desc":     g["pfam_desc"],
                "pfam_evalue":   g["pfam_evalue"],
                "pfam_bitscore": g["pfam_bitscore"],
                "class":         g["class"],
            })
    pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)


def save_completeness_table(out_path: str, names: list, completeness_list: list):
    rows = []
    for name, comp in zip(names, completeness_list):
        row = {
            "Genome":             name,
            "architecture":       comp["architecture"],
            "canonical_score":    comp["canonical_score"],
            "canonical_complete": comp["canonical_complete"],
            "modules_present":    ";".join(sorted(comp["modules_present"])),
        }
        for m in ALL_MODULES:
            row[f"n_{m}"] = comp["module_counts"].get(m, 0)
        row["n_other"] = comp["module_counts"].get("other", 0)
        rows.append(row)
    pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)

# ──────────────────────────────────────────────────────────────────────────────
# GROUP DISCOVERY
# ──────────────────────────────────────────────────────────────────────────────

def discover_groups(base_dir: str) -> list[tuple[str, str, list]]:
    """
    Return list of (group_name, group_path, [gbff_paths]).
    Handles both:
      - base_dir/group_*/  *.gbff      (grouped layout)
      - base_dir/          *.gbff      (flat layout → one synthetic group "all")
    """
    subdirs = [
        d for d in os.listdir(base_dir)
        if os.path.isdir(os.path.join(base_dir, d))
    ]
    groups = []
    for sub in sorted(subdirs):
        path  = os.path.join(base_dir, sub)
        gbffs = sorted(glob.glob(os.path.join(path, "*.gbff")))
        if gbffs:
            groups.append((sub, path, gbffs))

    # flat layout fallback
    if not groups:
        gbffs = sorted(glob.glob(os.path.join(base_dir, "*.gbff")))
        if gbffs:
            groups.append(("all", base_dir, gbffs))

    return groups

# ──────────────────────────────────────────────────────────────────────────────
# PIPELINE
# ──────────────────────────────────────────────────────────────────────────────

def run_pipeline(args):
    base_dir = args.gbff_dir
    out_root = args.outdir or base_dir
    os.makedirs(out_root, exist_ok=True)

    groups = discover_groups(base_dir)
    if not groups:
        print(f"No GBFF files found under {base_dir}")
        sys.exit(1)

    master_completeness = []
    master_names        = []

    for group_name, group_path, gbffs in groups:
        print(f"\n── {group_name} ({len(gbffs)} elements) ──")

        all_records       = []
        all_names         = []
        all_completeness  = []

        for gbff in gbffs:
            base = os.path.splitext(gbff)[0]
            name = os.path.basename(base)

            # ── PFAM ──────────────────────────────────────────────────────────
            pfam_hits = {}
            if not args.no_pfam and os.path.exists(args.pfam):
                faa    = base + ".faa"
                domtbl = base + ".pfam"

                if not os.path.exists(faa):
                    extract_faa(gbff, faa)

                if not os.path.exists(domtbl) or args.rerun_pfam:
                    print(f"  hmmscan {name} …", end=" ", flush=True)
                    ok = run_pfamscan(faa, domtbl, args.pfam)
                    print("OK" if ok else "FAILED")

                pfam_hits = parse_pfamscan(domtbl)
            elif not args.no_pfam and not os.path.exists(args.pfam):
                if group_name == groups[0][0] and gbff == gbffs[0]:
                    print(f"  [WARN] Pfam-A.hmm not found at {args.pfam} — "
                          f"running product-name classification only. "
                          f"Use --no-pfam to suppress this warning.")

            # ── Parse & classify ──────────────────────────────────────────────
            genes = parse_gbff(gbff, pfam_hits)
            if not genes:
                print(f"  [SKIP] {name} — no CDS found")
                continue

            comp = score_completeness(genes)

            all_records.append(genes)
            all_names.append(name)
            all_completeness.append(comp)

            master_names.append(f"{group_name}/{name}")
            master_completeness.append(comp)

            # Quick per-element summary to stdout
            present = ", ".join(sorted(comp["modules_present"])) or "none"
            print(f"  {name:50s}  {comp['architecture']:35s}  "
                  f"score={comp['canonical_score']:.0%}  modules=[{present}]")

        if not all_records:
            continue

        # ── Per-group TSV outputs ─────────────────────────────────────────────
        ann_tsv  = os.path.join(group_path, f"{group_name}_annotations.tsv")
        comp_tsv = os.path.join(group_path, f"{group_name}_completeness.tsv")
        save_annotation_table(ann_tsv, all_records, all_names)
        save_completeness_table(comp_tsv, all_names, all_completeness)
        print(f"  → {ann_tsv}")
        print(f"  → {comp_tsv}")

        # ── Per-element single PNG schematics ─────────────────────────────────
        singles_dir = os.path.join(group_path, "schematics")
        os.makedirs(singles_dir, exist_ok=True)
        for genes, name, comp in zip(all_records, all_names, all_completeness):
            out_png = os.path.join(singles_dir, f"{name}_schematic.png")
            plot_single(genes, name, comp, out_png)
        print(f"  → {singles_dir}/ ({len(all_records)} individual PNGs)")

        # ── Batched overview plots (30 elements per page) ─────────────────────
        overview_png = os.path.join(group_path, f"{group_name}_overview.png")
        overview_pdf = os.path.join(group_path, f"{group_name}_overview.pdf")
        saved_pngs = list(plot_group(all_records, all_names, all_completeness,
                                     overview_png, title=group_name, batch_size=30))
        list(plot_group(all_records, all_names, all_completeness,
                        overview_pdf, title=group_name, batch_size=30))
        for p in saved_pngs:
            print(f"  → {p}")

    # ── Master summary across all groups ─────────────────────────────────────
    master_tsv = os.path.join(out_root, "master_completeness_summary.tsv")
    save_completeness_table(master_tsv, master_names, master_completeness)
    print(f"\n── Master summary → {master_tsv}")

    # ── Architecture distribution ─────────────────────────────────────────────
    arch_counts = {}
    for c in master_completeness:
        arch_counts[c["architecture"]] = arch_counts.get(c["architecture"], 0) + 1

    print("\nArchitecture distribution:")
    for arch, n in sorted(arch_counts.items(), key=lambda x: -x[1]):
        print(f"  {arch:40s} {n:4d}")

    # ── Plots grouped by architecture ─────────────────────────────────────────
    # Build a lookup: genome name -> (genes, comp)
    # master_names are "group/element_name" — strip the group prefix for file lookup
    print("\n── Generating architecture-grouped schematics ──")

    arch_dir = os.path.join(out_root, "by_architecture")
    os.makedirs(arch_dir, exist_ok=True)

    # Index all parsed data by the bare element name (after group/ prefix)
    name_to_genes = {}
    name_to_comp  = {}
    for group_name, group_path, gbffs in groups:
        # Re-read the annotation TSV for this group to avoid re-running PFAM
        ann_tsv = os.path.join(group_path, f"{group_name}_annotations.tsv")
        if not os.path.exists(ann_tsv):
            continue
        ann_df = pd.read_csv(ann_tsv, sep="\t")
        for element_name, edf in ann_df.groupby("Genome"):
            genes = edf.to_dict("records")
            # Restore field names expected by plot functions
            for g in genes:
                g["start"]  = int(g["start"])
                g["end"]    = int(g["end"])
                g["strand"] = int(g["strand"])
            name_to_genes[element_name] = genes

    for full_name, comp in zip(master_names, master_completeness):
        bare = full_name.split("/", 1)[-1]
        name_to_comp[bare] = comp

    # Group and plot by architecture
    from collections import defaultdict
    arch_groups = defaultdict(list)
    for full_name, comp in zip(master_names, master_completeness):
        bare = full_name.split("/", 1)[-1]
        arch_groups[comp["architecture"]].append(bare)

    for arch, element_names in sorted(arch_groups.items()):
        # Filter to elements we have gene data for
        valid = [n for n in element_names if n in name_to_genes]
        if not valid:
            print(f"  [SKIP] {arch} — no gene data found")
            continue

        recs  = [name_to_genes[n] for n in valid]
        comps = [name_to_comp[n]  for n in valid]

        out_png = os.path.join(arch_dir, f"{arch}_schematic.png")
        out_pdf = os.path.join(arch_dir, f"{arch}_schematic.pdf")

        saved = list(plot_group(recs, valid, comps, out_png,
                                title=f"{arch}  (n={len(valid)})",
                                batch_size=30))
        list(plot_group(recs, valid, comps, out_pdf,
                        title=f"{arch}  (n={len(valid)})",
                        batch_size=30))

        for p in saved:
            print(f"  → {p}  [{len(valid)} elements]")


# ──────────────────────────────────────────────────────────────────────────────
# REGRAPH PIPELINE  (--regraph flag)
# ──────────────────────────────────────────────────────────────────────────────

def regraph_pipeline(args):
    """
    Read existing *_annotations.tsv and master_completeness_summary.tsv,
    regenerate all plots without touching PFAM or re-classifying anything.

    Expects the same directory layout produced by a previous full run:
      <gbff_dir>/
        [group_name]/
          <group_name>_annotations.tsv
        master_completeness_summary.tsv   (in outdir or gbff_dir)
    """
    base_dir = args.gbff_dir
    out_root = args.outdir or base_dir

    # ── Load master completeness ──────────────────────────────────────────────
    master_tsv = os.path.join(out_root, "master_completeness_summary.tsv")
    if not os.path.exists(master_tsv):
        print(f"[ERROR] master_completeness_summary.tsv not found at {master_tsv}")
        print("        Run without --regraph first to generate it.")
        sys.exit(1)

    master_df = pd.read_csv(master_tsv, sep="\t")
    print(f"Loaded {len(master_df)} elements from {master_tsv}")

    # Rebuild comp dicts keyed by bare element name
    name_to_comp  = {}
    for _, row in master_df.iterrows():
        bare = row["Genome"].split("/", 1)[-1]
        modules_present = set(
            m for m in str(row.get("modules_present") or "").split(";") if m
        )
        # Rebuild module_counts from n_* columns
        module_counts = {}
        for col in master_df.columns:
            if col.startswith("n_"):
                module_counts[col[2:]] = int(row[col])
        name_to_comp[bare] = {
            "architecture":       row["architecture"],
            "canonical_score":    float(row["canonical_score"]),
            "canonical_complete": bool(row["canonical_complete"]),
            "modules_present":    modules_present,
            "module_counts":      module_counts,
        }

    # ── Load gene data from annotation TSVs ──────────────────────────────────
    name_to_genes = {}
    groups = discover_groups(base_dir)

    for group_name, group_path, _ in groups:
        ann_tsv = os.path.join(group_path, f"{group_name}_annotations.tsv")
        if not os.path.exists(ann_tsv):
            print(f"  [WARN] Missing {ann_tsv} — skipping group {group_name}")
            continue
        ann_df = pd.read_csv(ann_tsv, sep="\t")
        for element_name, edf in ann_df.groupby("Genome"):
            genes = edf.to_dict("records")
            for g in genes:
                g["start"]  = int(g["start"])
                g["end"]    = int(g["end"])
                g["strand"] = int(g["strand"])
                # ensure 'class' key exists (column may be named differently)
                if "class" not in g:
                    g["class"] = g.get("gene_class", "other")
            name_to_genes[element_name] = genes

    # Flat layout fallback — annotation TSV at root
    if not name_to_genes:
        ann_tsv = os.path.join(base_dir, "all_annotations.tsv")
        if os.path.exists(ann_tsv):
            ann_df = pd.read_csv(ann_tsv, sep="\t")
            for element_name, edf in ann_df.groupby("Genome"):
                genes = edf.to_dict("records")
                for g in genes:
                    g["start"]  = int(g["start"])
                    g["end"]    = int(g["end"])
                    g["strand"] = int(g["strand"])
                    if "class" not in g:
                        g["class"] = g.get("gene_class", "other")
                name_to_genes[element_name] = genes

    if not name_to_genes:
        print("[ERROR] No *_annotations.tsv files found. Cannot regraph.")
        sys.exit(1)

    print(f"Loaded gene data for {len(name_to_genes)} elements")

    # ── Per-element single PNGs ───────────────────────────────────────────────
    singles_dir = os.path.join(out_root, "schematics")
    os.makedirs(singles_dir, exist_ok=True)
    n_drawn = 0
    for name, genes in name_to_genes.items():
        comp = name_to_comp.get(name)
        if comp is None:
            print(f"  [WARN] No completeness data for {name} — skipping")
            continue
        out_png = os.path.join(singles_dir, f"{name}_schematic.png")
        plot_single(genes, name, comp, out_png)
        n_drawn += 1
    print(f"\nWrote {n_drawn} individual schematics → {singles_dir}/")

    # ── Architecture-grouped overview plots ───────────────────────────────────
    from collections import defaultdict
    arch_groups = defaultdict(list)
    for name in name_to_genes:
        comp = name_to_comp.get(name)
        if comp:
            arch_groups[comp["architecture"]].append(name)

    arch_dir = os.path.join(out_root, "by_architecture")
    os.makedirs(arch_dir, exist_ok=True)
    print("\nGenerating architecture-grouped plots:")

    for arch, element_names in sorted(arch_groups.items()):
        valid = [n for n in element_names if n in name_to_genes]
        if not valid:
            continue
        recs  = [name_to_genes[n] for n in valid]
        comps = [name_to_comp[n]  for n in valid]

        out_png = os.path.join(arch_dir, f"{arch}_schematic.png")
        out_pdf = os.path.join(arch_dir, f"{arch}_schematic.pdf")

        saved = list(plot_group(recs, valid, comps, out_png,
                                title=f"{arch}  (n={len(valid)})",
                                batch_size=30))
        list(plot_group(recs, valid, comps, out_pdf,
                        title=f"{arch}  (n={len(valid)})",
                        batch_size=30))
        for p in saved:
            print(f"  → {p}  [{len(valid)} elements]")

    print("\nDone (regraph only — no PFAM or classification was re-run).")


# ──────────────────────────────────────────────────────────────────────────────
# ENTRY POINT
# ──────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    args = parse_args()
    if args.regraph:
        regraph_pipeline(args)
    else:
        run_pipeline(args)
