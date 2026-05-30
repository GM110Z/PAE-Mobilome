#!/usr/bin/env python3
"""
att_find_attB_first.py

Unbiased att site detection using an attB-first approach.

Unlike direct repeat searching (which assumes you know where the att site is),
this script:

  1. Extracts flanking protein sequences from each PICI element
  2. Uses DIAMOND to find the same chromosomal locus in PICI-negative genomes
  3. Extracts the intergenic sequence there — this IS the attB (empty site)
  4. Aligns attB against the full element sequence
  5. attB appears twice in the element — once as attL (left boundary)
     and once as attR (right boundary)
  6. The overlap between attB/attL/attR is the att site

This works for ALL elements regardless of:
  - Architecture class
  - Integrase type
  - Whether direct repeats were detectable before
  - Whether the element is complete or partial

Usage:
    # Step 1: Build negative genome DIAMOND database (once)
    python att_find_attB_first.py build \
        --pa_matrix pici_pa_matrix.tsv \
        --neg_gbff_dir /path/to/full/genomes \
        --outdir att_attB \
        --threads 8

    # Step 2: Find att sites for all elements
    python att_find_attB_first.py find \
        --summary extraction_summary.tsv \
        --gbff_dir pici_gbff_NEWregions \
        --db_dir att_attB \
        --outdir att_attB \
        --threads 8

Dependencies:
    diamond, biopython, pandas, numpy
"""

import argparse
import os
import sys
import subprocess
import tempfile
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# ──────────────────────────────────────────────────────────────────────────────
# CONSTANTS
# ──────────────────────────────────────────────────────────────────────────────

N_FLANK        = 5      # flanking genes each side
MAX_LOCUS_GAP  = 100000 # max bp between left/right anchor hits
MIN_ATT_LEN    = 8      # minimum att site length to report
ATT_MATCH_ID   = 0.85   # minimum identity for attB vs attL/attR alignment
DIAMOND_ID     = 40     # minimum % identity for flanking gene DIAMOND hits
DIAMOND_COV    = 70     # minimum query coverage for DIAMOND hits

# ──────────────────────────────────────────────────────────────────────────────
# ARGUMENT PARSING
# ──────────────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command")

    # Build subcommand
    build = sub.add_parser("build",
        help="Build DIAMOND database from PICI-negative genomes")
    build.add_argument("--pa_matrix",    required=True)
    build.add_argument("--neg_gbff_dir", required=True)
    build.add_argument("--outdir",       default="att_attB")
    build.add_argument("--threads",      type=int, default=4)

    # Find subcommand
    find = sub.add_parser("find",
        help="Find att sites using attB-first approach")
    find.add_argument("--summary",  required=True,
                      help="extraction_summary.tsv")
    find.add_argument("--gbff_dir", required=True,
                      help="Directory with extracted PICI GBFFs")
    find.add_argument("--db_dir",   required=True,
                      help="Directory with DIAMOND database (from build step)")
    find.add_argument("--outdir",   default="att_attB")
    find.add_argument("--threads",  type=int, default=4)
    find.add_argument("--neg_gbff_dir", required=True,
                      help="Directory with full genome GBFFs for negative genomes")
    find.add_argument("--completeness", default=None,
                      help="master_completeness_summary.tsv for architecture labels")

    return p.parse_args()

# ──────────────────────────────────────────────────────────────────────────────
# FILE UTILITIES
# ──────────────────────────────────────────────────────────────────────────────

def find_gbff(directory: str, name: str) -> Path | None:
    root = Path(directory)
    exact = root / f"{name}.gbff"
    if exact.exists():
        return exact
    for candidate in root.rglob(f"{name}.gbff"):
        return candidate
    for candidate in root.rglob(f"*{name}*.gbff"):
        return candidate
    return None


def find_neg_gbff(directory: str, genome_id: str) -> Path | None:
    root  = Path(directory)
    clean = genome_id
    for strip in [".g.fasta", ".fasta", ".fa", ".fna"]:
        if clean.endswith(strip):
            clean = clean[:-len(strip)]
            break
    for ext in [".gbff", ".gbk", ".gb", "_genomic.gbff"]:
        for candidate in root.rglob(f"*{clean}*{ext}"):
            return candidate
    return None


def get_cds_list(record: SeqRecord) -> list:
    return sorted(
        [f for f in record.features if f.type == "CDS"],
        key=lambda f: int(f.location.start)
    )


def get_trna_info(record: SeqRecord,
                   pici_start: int,
                   pici_end: int,
                   window: int = 500) -> tuple[bool, str | None]:
    """
    Check if any tRNA feature falls within window bp of the PICI boundaries.
    Returns (trna_at_boundary, trna_product_string).
    """
    trna_products = []
    for feat in record.features:
        if feat.type != "tRNA":
            continue
        fstart = int(feat.location.start)
        fend   = int(feat.location.end)
        # Check proximity to either boundary
        near_left  = abs(fstart - pici_start) < window or abs(fend - pici_start) < window
        near_right = abs(fstart - pici_end)   < window or abs(fend - pici_end)   < window
        if near_left or near_right:
            prod = feat.qualifiers.get("product", ["unknown"])[0]
            trna_products.append(prod)

    if trna_products:
        return True, "; ".join(sorted(set(trna_products)))
    return False, None

# ──────────────────────────────────────────────────────────────────────────────
# BUILD STEP
# ──────────────────────────────────────────────────────────────────────────────

def cmd_build(args):
    """
    Extract proteins from all PICI-negative genomes, build DIAMOND db,
    and save a protein location index.
    """
    os.makedirs(args.outdir, exist_ok=True)

    pa_df = pd.read_csv(args.pa_matrix, sep=r'\s+', header=None,
                         names=["genome", "pici_present", "n_systems"],
                         engine="python")
    pa_df = pa_df[~pa_df["pici_present"].str.contains("[A-Za-z]", na=False)]
    pa_df["pici_present"] = pa_df["pici_present"].astype(int)
    neg_genomes = pa_df[pa_df["pici_present"] == 0]["genome"].tolist()
    print(f"Building database from {len(neg_genomes)} PICI-negative genomes")

    db_faa  = os.path.join(args.outdir, "neg_genomes.faa")
    db_dmnd = os.path.join(args.outdir, "neg_genomes.dmnd")
    idx_tsv = os.path.join(args.outdir, "neg_protein_locations.tsv")

    loc_rows  = []
    n_written = 0

    with open(db_faa, "w") as fh:
        for i, genome in enumerate(neg_genomes):
            if (i + 1) % 500 == 0:
                print(f"  [{i+1}/{len(neg_genomes)}] proteins written: {n_written}")
            gbff_path = find_neg_gbff(args.neg_gbff_dir, genome)
            if not gbff_path:
                continue
            try:
                for rec in SeqIO.parse(str(gbff_path), "genbank"):
                    for feat in rec.features:
                        if feat.type != "CDS":
                            continue
                        pid = feat.qualifiers.get(
                            "protein_id", [""])[0].split("|")[-1]
                        seq = feat.qualifiers.get("translation", [""])[0]
                        if not pid or not seq:
                            continue
                        uid = f"{genome}||{pid}"
                        fh.write(f">{uid}\n{seq}\n")
                        loc_rows.append({
                            "uid":     uid,
                            "genome":  genome,
                            "contig":  rec.id,
                            "start":   int(feat.location.start),
                            "end":     int(feat.location.end),
                            "strand":  feat.location.strand,
                        })
                        n_written += 1
            except Exception:
                continue

    print(f"Total proteins: {n_written}")

    # Save location index
    loc_df = pd.DataFrame(loc_rows)
    loc_df.to_csv(idx_tsv, sep="\t", index=False)
    print(f"Location index saved: {idx_tsv}")

    # Build DIAMOND database
    print("Running diamond makedb...")
    subprocess.run([
        "diamond", "makedb",
        "--in",      db_faa,
        "--db",      db_dmnd,
        "--threads", str(args.threads),
        "--quiet"
    ], check=True)
    print(f"Database built: {db_dmnd}")

# ──────────────────────────────────────────────────────────────────────────────
# EXTRACT FLANKING PROTEINS
# ──────────────────────────────────────────────────────────────────────────────

def extract_flanking_proteins(gbff_path: Path,
                               first_hit: str,
                               last_hit: str) -> tuple[list, list, str, str]:
    """
    Extract flanking protein sequences and full element sequence.
    Returns (left_proteins, right_proteins, left_boundary_seq, right_boundary_seq)
    where each protein list is [(uid, sequence), ...] and boundary seqs are
    the 500bp immediately flanking the element.
    """
    try:
        records = list(SeqIO.parse(str(gbff_path), "genbank"))
    except Exception:
        return [], [], "", ""

    if not records:
        return [], [], "", ""

    record  = records[0]
    cds     = get_cds_list(record)
    seq_str = str(record.seq)

    # Find PICI boundary
    first_idx = last_idx = -1
    for i, feat in enumerate(cds):
        pid = feat.qualifiers.get("protein_id", [""])[0].split("|")[-1]
        if pid == first_hit and first_idx == -1:
            first_idx = i
        if pid == last_hit:
            last_idx = i

    if first_idx < 0 or last_idx < 0:
        # Positional fallback
        if len(cds) > N_FLANK * 2:
            first_idx = N_FLANK
            last_idx  = len(cds) - N_FLANK - 1
        else:
            return [], [], "", ""

    left_cds  = cds[max(0, first_idx - N_FLANK):first_idx]
    right_cds = cds[last_idx + 1:last_idx + 1 + N_FLANK]

    pici_start = int(cds[first_idx].location.start)
    pici_end   = int(cds[last_idx].location.end)

    # Boundary sequences (500bp flanking the element)
    left_boundary  = seq_str[max(0, pici_start - 500):pici_start]
    right_boundary = seq_str[pici_end:min(len(seq_str), pici_end + 500)]
    element_seq    = seq_str[pici_start:pici_end]

    # tRNA detection
    trna_at_boundary, trna_product = get_trna_info(record, pici_start, pici_end)

    def get_proteins(cds_list, side):
        proteins = []
        for j, feat in enumerate(cds_list):
            pid = feat.qualifiers.get("protein_id", [""])[0].split("|")[-1]
            seq = feat.qualifiers.get("translation", [""])[0]
            if pid and seq and pid != "unknown":
                proteins.append((f"{side}_{j}|{pid}", seq))
        return proteins

    return (get_proteins(left_cds, "LEFT"),
            get_proteins(right_cds, "RIGHT"),
            left_boundary,
            right_boundary,
            element_seq,
            trna_at_boundary,
            trna_product)

# ──────────────────────────────────────────────────────────────────────────────
# DIAMOND SEARCH
# ──────────────────────────────────────────────────────────────────────────────

def diamond_search(query_faa: str, db_dmnd: str, threads: int) -> pd.DataFrame:
    with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as tmp:
        out_tsv = tmp.name

    cmd = [
        "diamond", "blastp",
        "--query",        query_faa,
        "--db",           db_dmnd,
        "--out",          out_tsv,
        "--outfmt", "6",
        "qseqid", "sseqid", "pident", "length",
        "qstart",  "qend",  "sstart", "send",
        "evalue",  "bitscore",
        "--evalue",          "1e-5",
        "--id",              str(DIAMOND_ID),
        "--query-cover",     str(DIAMOND_COV),
        "--threads",         str(threads),
        "--quiet",
        "--max-target-seqs", "50",
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True)
    except subprocess.CalledProcessError:
        return pd.DataFrame()

    if not os.path.exists(out_tsv) or os.path.getsize(out_tsv) == 0:
        os.unlink(out_tsv)
        return pd.DataFrame()

    df = pd.read_csv(out_tsv, sep="\t", header=None,
                     names=["qseqid","sseqid","pident","length",
                             "qstart","qend","sstart","send",
                             "evalue","bitscore"])
    os.unlink(out_tsv)
    return df

# ──────────────────────────────────────────────────────────────────────────────
# FIND attB LOCUS
# ──────────────────────────────────────────────────────────────────────────────

def find_attB_locus(left_hits: pd.DataFrame,
                     right_hits: pd.DataFrame,
                     loc_index: dict) -> list[dict]:
    """
    Given DIAMOND hits for left and right flanking proteins,
    find negative genome loci where both sides hit on the same contig
    within MAX_LOCUS_GAP bp.

    Returns list of candidate loci:
        [{"genome", "contig", "attB_start", "attB_end"}, ...]
    """
    if left_hits.empty or right_hits.empty:
        return []

    # Group hits by genome+contig
    def hits_by_contig(hits, side):
        result = defaultdict(list)
        for _, row in hits.iterrows():
            uid = row["sseqid"]
            if uid not in loc_index:
                continue
            genome, contig, start, end, strand = loc_index[uid]
            result[(genome, contig)].append({
                "side":   side,
                "start":  int(start),
                "end":    int(end),
                "pident": row["pident"],
            })
        return result

    left_by_contig  = hits_by_contig(left_hits,  "LEFT")
    right_by_contig = hits_by_contig(right_hits, "RIGHT")

    shared = set(left_by_contig.keys()) & set(right_by_contig.keys())
    loci   = []

    for genome, contig in shared:
        left_positions  = left_by_contig[(genome, contig)]
        right_positions = right_by_contig[(genome, contig)]

        # Find best left anchor (closest to right hits)
        for lp in left_positions:
            for rp in right_positions:
                gap = abs(rp["start"] - lp["end"])
                if gap > MAX_LOCUS_GAP:
                    continue

                # attB is the sequence between these two anchors
                attB_start = min(lp["end"],   rp["start"])
                attB_end   = max(lp["start"], rp["end"])

                loci.append({
                    "genome":     genome,
                    "contig":     contig,
                    "attB_start": attB_start,
                    "attB_end":   attB_end,
                    "gap_bp":     gap,
                    "left_pident":  lp["pident"],
                    "right_pident": rp["pident"],
                })

    # Sort by gap (smallest = most likely same locus)
    loci.sort(key=lambda x: x["gap_bp"])
    return loci[:5]  # return top 5 candidate loci

# ──────────────────────────────────────────────────────────────────────────────
# EXTRACT attB AND FIND att SITE
# ──────────────────────────────────────────────────────────────────────────────

def extract_attB(locus: dict, neg_gbff_dir: str) -> str | None:
    """Extract the attB sequence from a PICI-negative genome."""
    gbff_path = find_neg_gbff(neg_gbff_dir, locus["genome"])
    if not gbff_path:
        return None
    try:
        for rec in SeqIO.parse(str(gbff_path), "genbank"):
            if rec.id != locus["contig"]:
                continue
            try:
                seq = str(rec.seq).upper()
            except Exception:
                return None
            # Extract with a window around the locus
            start = max(0, locus["attB_start"] - 200)
            end   = min(len(seq), locus["attB_end"] + 200)
            return seq[start:end]
    except Exception:
        return None
    return None


def find_att_site(attB: str,
                   left_boundary: str,
                   right_boundary: str,
                   element_seq: str,
                   min_len: int = MIN_ATT_LEN,
                   min_id: float = ATT_MATCH_ID) -> dict | None:
    """
    Find the att site by aligning attB against the element boundaries.

    attB appears as:
      - attL: near the left boundary of the element
      - attR: near the right boundary of the element

    The overlap between attB/attL/attR is the att site.

    Strategy:
      For each k-mer in attB (length min_len to len(attB)):
        Search for it in left_boundary AND right_boundary
        If found in both → candidate att site
        Score by length and identity
    """
    attB_upper  = attB.upper()
    left_upper  = left_boundary.upper()
    right_upper = right_boundary.upper()

    best_hit = None

    for length in range(min(len(attB_upper), 50), min_len - 1, -1):
        for i in range(len(attB_upper) - length + 1):
            kmer = attB_upper[i:i + length]

            # Skip low complexity
            if len(set(kmer)) < 3:
                continue

            # Search in left boundary
            left_pos = left_upper.find(kmer)
            if left_pos == -1:
                # Try with allowed mismatches
                left_pos = fuzzy_find(kmer, left_upper, min_id)

            if left_pos == -1:
                continue

            # Search in right boundary
            right_pos = right_upper.find(kmer)
            if right_pos == -1:
                right_pos = fuzzy_find(kmer, right_upper, min_id)

            if right_pos == -1:
                continue

            # Found in both boundaries — this is the att site
            gc = (kmer.count('G') + kmer.count('C')) / length * 100
            hit = {
                "att_sequence":       kmer,
                "att_length":         length,
                "attB_position":      i,
                "attL_position":      left_pos,
                "attR_position":      right_pos,
                "att_gc_pct":         round(gc, 1),
                "attB_source_genome": None,  # filled in later
            }

            if best_hit is None or length > best_hit["att_length"]:
                best_hit = hit

        if best_hit is not None:
            break  # Found best length — stop

    return best_hit


def fuzzy_find(kmer: str, seq: str, min_id: float) -> int:
    """Find kmer in seq with allowed mismatches. Returns position or -1."""
    n = len(kmer)
    for i in range(len(seq) - n + 1):
        candidate = seq[i:i + n]
        matches   = sum(a == b for a, b in zip(kmer, candidate))
        if matches / n >= min_id:
            return i
    return -1

# ──────────────────────────────────────────────────────────────────────────────
# FIND STEP
# ──────────────────────────────────────────────────────────────────────────────

def cmd_find(args):
    os.makedirs(args.outdir, exist_ok=True)
    align_dir = os.path.join(args.outdir, "att_alignments")
    os.makedirs(align_dir, exist_ok=True)

    # ── Load inputs ───────────────────────────────────────────────────────────
    db_dmnd = os.path.join(args.db_dir, "neg_genomes.dmnd")
    idx_tsv = os.path.join(args.db_dir, "neg_protein_locations.tsv")

    if not os.path.exists(db_dmnd):
        print(f"[ERROR] DIAMOND database not found: {db_dmnd}")
        print("Run the 'build' step first.")
        sys.exit(1)

    print(f"Loading protein location index...")
    loc_df    = pd.read_csv(idx_tsv, sep="\t")
    loc_index = dict(zip(
        loc_df["uid"],
        zip(loc_df["genome"], loc_df["contig"],
            loc_df["start"],  loc_df["end"], loc_df["strand"])
    ))
    print(f"  {len(loc_index)} proteins indexed")

    summary = pd.read_csv(args.summary, sep="\t")
    summary = summary[summary["status"] == "OK"].copy()
    summary["element"] = summary["Genome"] + "_" + summary["sys_id"]
    print(f"Processing {len(summary)} elements")

    # Architecture map
    arch_map = {}
    if args.completeness and os.path.exists(args.completeness):
        comp_df = pd.read_csv(args.completeness, sep="\t")
        comp_df["bare"] = comp_df["Genome"].apply(
            lambda x: x.split("/", 1)[-1])
        arch_map = dict(zip(comp_df["bare"], comp_df["architecture"]))
        print(f"Architecture data for {len(arch_map)} elements")

    # ── Process each element ──────────────────────────────────────────────────
    results = []

    for idx, row in summary.iterrows():
        element   = row["element"]
        first_hit = row["first_hit"]
        last_hit  = row["last_hit"]

        if (idx + 1) % 100 == 0 or (idx + 1) == len(summary):
            print(f"  [{idx+1}/{len(summary)}] {element}")

        result = {
            "element":           element,
            "genome":            row["Genome"],
            "sys_id":            row["sys_id"],
            "architecture":      arch_map.get(element, "unknown"),
            "element_len_bp":    row.get("bp_length",
                                  row.get("bp_end", 0) - row.get("bp_start", 0)),
            "att_found":         False,
            "att_sequence":      None,
            "att_length":        None,
            "att_gc_pct":        None,
            "attB_genome":       None,
            "attL_position":     None,
            "attR_position":     None,
            "trna_at_boundary":  False,
            "trna_product":      None,
            "n_loci_found":      0,
            "notes":             "",
        }

        # ── Find extracted GBFF ───────────────────────────────────────────────
        gbff_path = find_gbff(args.gbff_dir, element)
        if not gbff_path:
            result["notes"] = "GBFF not found"
            results.append(result)
            continue

        # ── Extract flanking proteins and boundary sequences ──────────────────
        out = extract_flanking_proteins(gbff_path, first_hit, last_hit)
        if len(out) == 3:
            left_proteins, right_proteins, left_boundary = out
            right_boundary   = ""
            element_seq      = ""
            trna_at_boundary = False
            trna_product     = None
        elif len(out) == 5:
            left_proteins, right_proteins, left_boundary, right_boundary, element_seq = out
            trna_at_boundary = False
            trna_product     = None
        else:
            left_proteins, right_proteins, left_boundary, right_boundary, element_seq, trna_at_boundary, trna_product = out

        result["trna_at_boundary"] = trna_at_boundary
        result["trna_product"]     = trna_product

        if not left_proteins or not right_proteins:
            result["notes"] = "no flanking proteins extracted"
            results.append(result)
            continue

        # ── Write query FAA ───────────────────────────────────────────────────
        query_faa = os.path.join(args.outdir, f"_q_{element}.faa")
        with open(query_faa, "w") as fh:
            for uid, seq in left_proteins:
                fh.write(f">LEFT|{uid}\n{seq}\n")
            for uid, seq in right_proteins:
                fh.write(f">RIGHT|{uid}\n{seq}\n")

        # ── DIAMOND search ────────────────────────────────────────────────────
        hits = diamond_search(query_faa, db_dmnd, args.threads)
        os.unlink(query_faa)

        if hits.empty:
            result["notes"] = "no DIAMOND hits in negative genomes"
            results.append(result)
            continue

        left_hits  = hits[hits["qseqid"].str.startswith("LEFT|")]
        right_hits = hits[hits["qseqid"].str.startswith("RIGHT|")]

        # ── Find attB loci ────────────────────────────────────────────────────
        loci = find_attB_locus(left_hits, right_hits, loc_index)
        result["n_loci_found"] = len(loci)

        if not loci:
            result["notes"] = "DIAMOND hits found but no co-located loci"
            results.append(result)
            continue

        # ── Try each candidate locus ──────────────────────────────────────────
        att_hit = None
        for locus in loci:
            attB = extract_attB(locus, args.db_dir.replace(
                "att_attB", "") or ".")

            # Try neg_gbff_dir — need to pass it through
            # Use a hacky approach: check if db_dir/../ has the GBFFs
            # Better: check neg_gbff_dir argument
            # For now search common locations
            attB = None
            gbff = find_neg_gbff(args.neg_gbff_dir, locus["genome"])
            if gbff:
                try:
                    for rec in SeqIO.parse(str(gbff), "genbank"):
                        if rec.id != locus["contig"]:
                            continue
                        try:
                            seq = str(rec.seq).upper()
                        except Exception:
                            break
                        start = max(0, locus["attB_start"] - 200)
                        end   = min(len(seq), locus["attB_end"] + 200)
                        attB  = seq[start:end]
                        break
                except Exception:
                    pass

            if not attB:
                continue

            # ── Find att site by comparing attB to element boundaries ─────────
            att_hit = find_att_site(
                attB, left_boundary, right_boundary, element_seq)

            if att_hit:
                att_hit["attB_source_genome"] = locus["genome"]
                att_hit["locus_gap_bp"]        = locus["gap_bp"]
                break

        if att_hit:
            result.update({
                "att_found":    True,
                "att_sequence": att_hit["att_sequence"],
                "att_length":   att_hit["att_length"],
                "att_gc_pct":   att_hit["att_gc_pct"],
                "attB_genome":  att_hit["attB_source_genome"],
                "attL_position": att_hit["attL_position"],
                "attR_position": att_hit["attR_position"],
            })

            # Save alignment
            align_path = os.path.join(align_dir, f"{element}.txt")
            with open(align_path, "w") as f:
                f.write(f"Element:      {element}\n")
                f.write(f"Architecture: {result['architecture']}\n")
                f.write(f"att sequence: {att_hit['att_sequence']}\n")
                f.write(f"att length:   {att_hit['att_length']} bp\n")
                f.write(f"att GC:       {att_hit['att_gc_pct']:.1f}%\n")
                f.write(f"attB genome:  {att_hit['attB_source_genome']}\n\n")
                att = att_hit["att_sequence"]
                f.write(f"attB (empty site, att marked):\n")
                f.write(attB.replace(att, f"[{att}]") + "\n\n")
                f.write(f"Left boundary (attL, att marked):\n")
                f.write(left_boundary.upper().replace(att, f"[{att}]") + "\n\n")
                f.write(f"Right boundary (attR, att marked):\n")
                f.write(right_boundary.upper().replace(att, f"[{att}]") + "\n")

        results.append(result)

    # ── Save outputs ──────────────────────────────────────────────────────────
    df = pd.DataFrame(results).sort_values(
        ["att_found", "att_length"], ascending=[False, False])

    out_tsv = os.path.join(args.outdir, "att_attB_results.tsv")
    df.to_csv(out_tsv, sep="\t", index=False)

    found_tsv = os.path.join(args.outdir, "att_attB_found.tsv")
    df[df["att_found"]].to_csv(found_tsv, sep="\t", index=False)

    # Summary stats
    n_total    = len(df)
    n_hits     = int(df["att_found"].sum())
    n_loci     = int((df["n_loci_found"] > 0).sum())

    print(f"\n── Summary ──────────────────────────────────────────")
    print(f"  Elements processed:          {n_total}")
    print(f"  Elements with attB locus:    {n_loci} ({n_loci/n_total*100:.1f}%)")
    print(f"  att sites found:             {n_hits} ({n_hits/n_total*100:.1f}%)")

    # Per-architecture
    if arch_map:
        arch_summary = df.groupby("architecture").agg(
            n_elements   = ("att_found", "count"),
            n_loci_found = ("n_loci_found", lambda x: (x > 0).sum()),
            n_att_found  = ("att_found",    "sum"),
            mean_att_len = ("att_length",   "mean"),
        ).round(1)
        arch_summary["pct_loci"] = (
            arch_summary["n_loci_found"] /
            arch_summary["n_elements"] * 100).round(1)
        arch_summary["pct_att"] = (
            arch_summary["n_att_found"] /
            arch_summary["n_elements"] * 100).round(1)
        arch_summary = arch_summary.sort_values("pct_att", ascending=False)
        print(f"\n  Per-architecture:")
        print(arch_summary.to_string())
        arch_summary.to_csv(
            os.path.join(args.outdir, "att_by_architecture.tsv"), sep="\t")

    # Plot
    if n_hits > 0:
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        found = df[df["att_found"]]

        found["att_length"].astype(int).hist(
            bins=range(MIN_ATT_LEN, 55), ax=axes[0],
            color="#185FA5", edgecolor="black")
        axes[0].set_xlabel("att site length (bp)")
        axes[0].set_ylabel("Count")
        axes[0].set_title(f"att site length (n={n_hits})")

        if "architecture" in df.columns:
            arch_pct = arch_summary["pct_att"].sort_values(ascending=False)
            arch_pct.plot(kind="barh", ax=axes[1],
                          color="#fc8d62", edgecolor="black")
            axes[1].set_xlabel("% elements with att site")
            axes[1].set_title("att site rate by architecture")
            axes[1].invert_yaxis()

        plt.tight_layout()
        plt.savefig(os.path.join(args.outdir, "att_summary.png"),
                    dpi=150, bbox_inches="tight")
        plt.close()

    print(f"\n→ {out_tsv}")
    print(f"→ {found_tsv}")
    print(f"→ {align_dir}/")
    print("\nDone.")

# ──────────────────────────────────────────────────────────────────────────────
# ENTRY POINT
# ──────────────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    if args.command == "build":
        cmd_build(args)
    elif args.command == "find":
        cmd_find(args)
    else:
        print("Specify a command: build or find")
        sys.exit(1)


if __name__ == "__main__":
    main()
