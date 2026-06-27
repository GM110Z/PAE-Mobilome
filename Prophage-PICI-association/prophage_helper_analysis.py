
#!/usr/bin/env python3
"""
Inputs:
  --prophages       geNomad all_prophages.tsv
  --taxmyphage      taxMyPhage Summary_taxonomy.tsv
  --completeness    all_completeness.tsv
  --overlap_summary optional element_inovirus_overlap_summary.tsv
"""

import argparse
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact


UNKNOWN = {
    "", "unknown", "new_genus", "new genus", "not defined yet",
    "none", "nan", "unclassified", "unclassified viruses",
    "unclassified virus"
}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--prophages", required=True)
    p.add_argument("--taxmyphage", required=True)
    p.add_argument("--completeness", required=True)
    p.add_argument("--overlap_summary", default=None)
    p.add_argument("--outdir", default="combined_prophage_helper_results")
    p.add_argument("--architecture_wide_threshold", type=float, default=50.0)
    p.add_argument("--subgroup_threshold", type=float, default=10.0)
    return p.parse_args()


def clean_genome(x):
    x = str(x).strip().split("/")[-1]
    for suffix in [".g.fasta", ".fasta", ".fna", ".fa", ".gbff", ".faa"]:
        if x.endswith(suffix):
            x = x[:-len(suffix)]
            break
    return x


def genome_from_element(x):
    x = clean_genome(x)
    if "_UserReplicon" in x:
        x = x.rsplit("_UserReplicon", 1)[0]
    return x


def normalise_prophage_id(x):
    x = str(x).strip().replace("|", "_")
    x = re.sub(r"\s+", "", x)
    return x


def is_classified(x):
    return str(x).strip().lower() not in UNKNOWN


def bh_fdr(pvalues):
    p = np.asarray(pvalues, dtype=float)
    if len(p) == 0:
        return np.array([])
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    out = np.empty(len(p))
    out[order] = adjusted
    return out


def odds_ratio_ci(a, b, c, d):
    cells = np.array([a, b, c, d], dtype=float)
    if np.any(cells == 0):
        cells += 0.5
    a2, b2, c2, d2 = cells
    OR = (a2 * d2) / (b2 * c2)
    se = math.sqrt(1/a2 + 1/b2 + 1/c2 + 1/d2)
    z = 1.959963984540054
    low = math.exp(math.log(OR) - z * se)
    high = math.exp(math.log(OR) + z * se)
    return OR, low, high


def load_architectures(path):
    df = pd.read_csv(path, sep="\t")
    element_col = "Genome" if "Genome" in df.columns else df.columns[0]
    arch_col = "architecture" if "architecture" in df.columns else df.columns[1]

    out = df[[element_col, arch_col]].copy()
    out.columns = ["element", "architecture"]
    out["genome"] = out["element"].map(genome_from_element)

    genome_arch = (
        out.groupby("genome")["architecture"]
        .apply(lambda x: set(str(v) for v in x if pd.notna(v)))
        .to_dict()
    )

    sizes = (
        out.groupby("architecture")
        .agg(
            n_elements=("element", "size"),
            n_genomes=("genome", "nunique")
        )
        .reset_index()
        .sort_values("n_elements", ascending=False)
    )

    return out, genome_arch, sizes


def build_presence_matrix(all_genomes, genome_to_taxa):
    taxa = sorted({
        taxon for values in genome_to_taxa.values() for taxon in values
    })
    matrix = pd.DataFrame(0, index=sorted(all_genomes), columns=taxa, dtype=int)
    for genome, taxa_here in genome_to_taxa.items():
        if genome in matrix.index:
            matrix.loc[genome, list(taxa_here)] = 1
    matrix.index.name = "genome"
    return matrix


def association_table(matrix, genome_arch, source, rank):
    architectures = sorted({
        arch for values in genome_arch.values() for arch in values
    })
    all_genomes = set(matrix.index)
    rows = []

    for arch in architectures:
        focus = {g for g in all_genomes if arch in genome_arch.get(g, set())}
        control = all_genomes - focus

        for taxon in matrix.columns:
            a = int(matrix.loc[list(focus), taxon].sum()) if focus else 0
            c = int(matrix.loc[list(control), taxon].sum()) if control else 0
            b = len(focus) - a
            d = len(control) - c

            OR_fisher, p = fisher_exact([[a, b], [c, d]])
            OR_ci, low, high = odds_ratio_ci(a, b, c, d)
            pct_focus = 100 * a / len(focus) if focus else np.nan
            pct_control = 100 * c / len(control) if control else np.nan

            rows.append({
                "source": source,
                "rank": rank,
                "taxon": taxon,
                "architecture": arch,
                "n_focus_present": a,
                "n_focus_total": len(focus),
                "pct_focus": pct_focus,
                "n_control_present": c,
                "n_control_total": len(control),
                "pct_control": pct_control,
                "prevalence_difference_pp": pct_focus - pct_control,
                "odds_ratio_fisher": OR_fisher,
                "odds_ratio_ci_estimate": OR_ci,
                "ci95_low": low,
                "ci95_high": high,
                "p_value": p,
            })

    out = pd.DataFrame(rows)

    # FDR is applied within source/rank/architecture, because each architecture
    # is a separate biological question and taxon panel.
    out["p_fdr"] = np.nan
    for keys, idx in out.groupby(["source", "rank", "architecture"]).groups.items():
        out.loc[idx, "p_fdr"] = bh_fdr(out.loc[idx, "p_value"].values)

    return out


def classify_signal(df, wide_thr, subgroup_thr):
    def classify(row):
        if row["p_fdr"] >= 0.05:
            return "no_significant_association"
        if row["odds_ratio_ci_estimate"] < 1:
            return "depleted"
        if row["pct_focus"] >= wide_thr:
            return "architecture-wide_candidate"
        if row["pct_focus"] >= subgroup_thr:
            return "subgroup_association"
        return "rare_association"

    df = df.copy()
    df["interpretation"] = df.apply(classify, axis=1)
    return df


def genomad_inoviridae_matrix(prophages_path, all_genomes):
    df = pd.read_csv(prophages_path, sep="\t")
    if "genome" not in df.columns or "taxonomy" not in df.columns:
        raise ValueError("geNomad table must contain genome and taxonomy columns")

    df["host_genome"] = df["genome"].map(clean_genome)
    tax = df["taxonomy"].fillna("").astype(str)

    df["is_inoviridae"] = (
        tax.str.contains("Hofneiviricota", case=False, na=False) |
        tax.str.contains("Faserviricetes", case=False, na=False) |
        tax.str.contains("Tubulavirales", case=False, na=False) |
        tax.str.contains("Inoviridae", case=False, na=False)
    )

    positive = set(df.loc[df["is_inoviridae"], "host_genome"])
    genome_to_taxa = {
        g: {"Inoviridae"} if g in positive else set()
        for g in all_genomes
    }
    return build_presence_matrix(all_genomes, genome_to_taxa)


def taxmyphage_matrices(tax_path, prophages_path, all_genomes):
    tax = pd.read_csv(tax_path, sep="\t")
    proph = pd.read_csv(prophages_path, sep="\t")

    if "Genome" not in tax.columns:
        raise ValueError("taxMyPhage table must contain Genome column")
    if "seq_name" not in proph.columns or "genome" not in proph.columns:
        raise ValueError("geNomad table must contain seq_name and genome columns")

    tax = tax.copy()
    proph = proph.copy()

    tax["prophage_key"] = tax["Genome"].map(normalise_prophage_id)
    proph["prophage_key"] = proph["seq_name"].map(normalise_prophage_id)
    proph["host_genome"] = proph["genome"].map(clean_genome)

    joined = tax.merge(
        proph[["prophage_key", "host_genome"]].drop_duplicates(),
        on="prophage_key",
        how="left"
    )

    genus_map = {}
    family_map = {}

    for genome, grp in joined[joined["host_genome"].notna()].groupby("host_genome"):
        genus_map[genome] = set(
            str(x) for x in grp["Genus"] if is_classified(x)
        )
        family_map[genome] = set(
            str(x) for x in grp["Family"] if is_classified(x)
        )

    genus_matrix = build_presence_matrix(all_genomes, genus_map)
    family_matrix = build_presence_matrix(all_genomes, family_map)

    match_stats = {
        "taxmyphage_rows": len(tax),
        "matched_unique_ids": joined.loc[
            joined["host_genome"].notna(), "prophage_key"
        ].nunique(),
        "unmatched_unique_ids": joined.loc[
            joined["host_genome"].isna(), "prophage_key"
        ].nunique(),
        "classified_genus_rows": int(
            joined["Genus"].map(is_classified).fillna(False).sum()
        ),
        "classified_family_rows": int(
            joined["Family"].map(is_classified).fillna(False).sum()
        ),
    }

    return joined, genus_matrix, family_matrix, match_stats


def add_overlap_context(df, overlap_path):
    df = df.copy()
    df["coordinate_overlap_context"] = "not_assessed"

    if not overlap_path:
        return df

    path = Path(overlap_path)
    if not path.exists():
        return df

    ov = pd.read_csv(path, sep="\t")
    if "architecture" not in ov.columns:
        return df

    arch_overlap = (
        ov.groupby("architecture")
        .agg(
            n_elements_overlap_checked=("element_id", "size"),
            n_any_overlap=("has_any_overlap", "sum"),
            n_suspected_misclassification=("suspected_misclassification", "sum")
        )
        .reset_index()
    )

    df = df.merge(arch_overlap, on="architecture", how="left")
    is_inovirus = (
        (df["source"] == "geNomad") &
        (df["taxon"] == "Inoviridae")
    )

    df.loc[is_inovirus, "coordinate_overlap_context"] = (
        df.loc[is_inovirus, "n_suspected_misclassification"]
        .fillna(0)
        .astype(int)
        .astype(str)
        + " suspected overlapping elements"
    )
    return df


def plot_prevalence_comparison(sig, outpath):
    if sig.empty:
        return

    plot_df = sig.copy()
    plot_df["label"] = (
        plot_df["architecture"] + " — " +
        plot_df["taxon"] + " (" + plot_df["source"] + ")"
    )
    plot_df = plot_df.sort_values(
        ["interpretation", "pct_focus"],
        ascending=[True, False]
    )

    y = np.arange(len(plot_df))
    h = 0.38

    fig, ax = plt.subplots(figsize=(13, max(6, len(plot_df) * 0.48)))
    ax.barh(y - h/2, plot_df["pct_focus"], height=h, label="Architecture genomes")
    ax.barh(y + h/2, plot_df["pct_control"], height=h, label="Control genomes")
    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["label"], fontsize=8)
    ax.set_xlabel("Genomes carrying taxon (%)")
    ax.set_title("Absolute prevalence of significant architecture–phage associations")
    ax.legend()
    ax.invert_yaxis()

    for i, row in plot_df.reset_index(drop=True).iterrows():
        ax.text(
            max(row["pct_focus"], row["pct_control"]) + 0.7,
            i,
            row["interpretation"].replace("_", " "),
            va="center",
            fontsize=7,
        )

    plt.tight_layout()
    plt.savefig(outpath, dpi=220, bbox_inches="tight")
    plt.close()


def plot_forest(sig, outpath):
    if sig.empty:
        return

    plot_df = sig[
        np.isfinite(sig["odds_ratio_ci_estimate"]) &
        np.isfinite(sig["ci95_low"]) &
        np.isfinite(sig["ci95_high"])
    ].copy()

    plot_df["label"] = (
        plot_df["architecture"] + " — " +
        plot_df["taxon"] + " (" + plot_df["source"] + ")"
    )
    plot_df = plot_df.sort_values("odds_ratio_ci_estimate")
    y = np.arange(len(plot_df))

    OR = plot_df["odds_ratio_ci_estimate"].to_numpy()
    low = plot_df["ci95_low"].to_numpy()
    high = plot_df["ci95_high"].to_numpy()
    xerr = np.vstack([
        np.maximum(0, OR - low),
        np.maximum(0, high - OR)
    ])

    fig, ax = plt.subplots(figsize=(12, max(6, len(plot_df) * 0.48)))
    ax.errorbar(OR, y, xerr=xerr, fmt="o", capsize=3)
    ax.axvline(1, linestyle="--", linewidth=1)
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["label"], fontsize=8)
    ax.set_xlabel("Odds ratio with 95% confidence interval")
    ax.set_title("Statistical enrichment/depletion with absolute prevalence shown in labels")

    for i, row in plot_df.reset_index(drop=True).iterrows():
        ax.text(
            row["ci95_high"] * 1.08,
            i,
            f"{int(row['n_focus_present'])}/{int(row['n_focus_total'])} "
            f"({row['pct_focus']:.1f}%)",
            va="center",
            fontsize=7,
        )

    plt.tight_layout()
    plt.savefig(outpath, dpi=220, bbox_inches="tight")
    plt.close()


def plot_heatmap(all_assoc, outpath):
    sig_taxa = all_assoc.loc[
        all_assoc["p_fdr"] < 0.05,
        ["source", "rank", "taxon"]
    ].drop_duplicates()

    if sig_taxa.empty:
        return

    keys = {
        (r.source, r.rank, r.taxon)
        for r in sig_taxa.itertuples()
    }

    subset = all_assoc[
        all_assoc.apply(
            lambda r: (r["source"], r["rank"], r["taxon"]) in keys,
            axis=1
        )
    ].copy()

    subset["taxon_label"] = (
        subset["taxon"] + " [" + subset["source"] + "]"
    )

    matrix = subset.pivot_table(
        index="architecture",
        columns="taxon_label",
        values="prevalence_difference_pp",
        aggfunc="first"
    ).fillna(0)

    fig, ax = plt.subplots(
        figsize=(max(9, matrix.shape[1] * 1.2),
                 max(6, matrix.shape[0] * 0.5))
    )
    im = ax.imshow(matrix.values, aspect="auto", cmap="coolwarm")
    ax.set_xticks(np.arange(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(np.arange(matrix.shape[0]))
    ax.set_yticklabels(matrix.index, fontsize=8)
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Prevalence difference (percentage points)")
    ax.set_title("Architecture–phage prevalence differences")
    plt.tight_layout()
    plt.savefig(outpath, dpi=220, bbox_inches="tight")
    plt.close()


def architecture_interpretation(sig, sizes):
    rows = []
    for arch in sizes["architecture"]:
        grp = sig[sig["architecture"] == arch].copy()
        positive = grp[
            grp["interpretation"].isin([
                "architecture-wide_candidate",
                "subgroup_association",
                "rare_association"
            ])
        ]
        depleted = grp[grp["interpretation"] == "depleted"]

        if positive.empty:
            overall = "no_positive_candidate_detected"
        elif (positive["interpretation"] == "architecture-wide_candidate").any():
            overall = "architecture-wide_candidate_present"
        elif (positive["interpretation"] == "subgroup_association").any():
            overall = "subgroup_candidate_only"
        else:
            overall = "rare_candidate_only"

        best = ""
        if not positive.empty:
            best_row = positive.sort_values(
                ["pct_focus", "p_fdr"],
                ascending=[False, True]
            ).iloc[0]
            best = (
                f"{best_row['taxon']} ({best_row['source']}): "
                f"{int(best_row['n_focus_present'])}/"
                f"{int(best_row['n_focus_total'])} "
                f"({best_row['pct_focus']:.1f}%), "
                f"OR={best_row['odds_ratio_ci_estimate']:.2f}, "
                f"FDR={best_row['p_fdr']:.3g}"
            )

        rows.append({
            "architecture": arch,
            "overall_interpretation": overall,
            "best_positive_association": best,
            "n_positive_significant_taxa": len(positive),
            "n_depleted_significant_taxa": len(depleted),
        })

    return sizes.merge(pd.DataFrame(rows), on="architecture", how="left")


def write_summary(outdir, sizes, match_stats, sig, arch_interp):
    lines = [
        "COMBINED geNomad + taxMyPhage HELPER-PHAGE ANALYSIS",
        "===================================================",
        f"MGE-positive genomes tested: {sizes['n_genomes'].sum()} "
        "(sum across overlapping architecture groups)",
        f"taxMyPhage rows: {match_stats['taxmyphage_rows']}",
        f"Matched taxMyPhage IDs: {match_stats['matched_unique_ids']}",
        f"Unmatched taxMyPhage IDs: {match_stats['unmatched_unique_ids']}",
        f"Classified genus rows: {match_stats['classified_genus_rows']}",
        f"Classified family rows: {match_stats['classified_family_rows']}",
        "",
        "Interpretation rule:",
        "  Significant enrichment is not automatically architecture-wide.",
        "  >=50% prevalence = architecture-wide candidate",
        "  10–<50% prevalence = subgroup association",
        "  <10% prevalence = rare association",
        "",
        "Architecture-level interpretation:",
    ]

    for _, row in arch_interp.iterrows():
        lines.append(
            f"  {row['architecture']} "
            f"(elements={int(row['n_elements'])}, genomes={int(row['n_genomes'])}): "
            f"{row['overall_interpretation']}"
        )
        if row["best_positive_association"]:
            lines.append(f"    best: {row['best_positive_association']}")

    lines.extend(["", "All significant associations:"])

    if sig.empty:
        lines.append("  None")
    else:
        for _, r in sig.sort_values(
            ["architecture", "interpretation", "pct_focus"],
            ascending=[True, True, False]
        ).iterrows():
            lines.append(
                f"  {r['architecture']} | {r['source']} | "
                f"{r['rank']} | {r['taxon']} | "
                f"{int(r['n_focus_present'])}/{int(r['n_focus_total'])} "
                f"({r['pct_focus']:.1f}%) vs "
                f"{int(r['n_control_present'])}/{int(r['n_control_total'])} "
                f"({r['pct_control']:.1f}%) | "
                f"Δ={r['prevalence_difference_pp']:.1f} pp | "
                f"OR={r['odds_ratio_ci_estimate']:.2f} "
                f"(95% CI {r['ci95_low']:.2f}–{r['ci95_high']:.2f}) | "
                f"FDR={r['p_fdr']:.3g} | "
                f"{r['interpretation']}"
            )

    with open(Path(outdir) / "analysis_summary.txt", "w") as fh:
        fh.write("\n".join(lines) + "\n")


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    arch_df, genome_arch, sizes = load_architectures(args.completeness)
    sizes.to_csv(outdir / "architecture_group_sizes.tsv", sep="\t", index=False)
    all_genomes = set(genome_arch)

    # geNomad Inoviridae
    inv_matrix = genomad_inoviridae_matrix(args.prophages, all_genomes)
    inv_assoc = association_table(
        inv_matrix, genome_arch, "geNomad", "family"
    )
    inv_assoc.to_csv(
        outdir / "genomad_inoviridae_associations.tsv",
        sep="\t", index=False
    )

    # taxMyPhage family/genus
    joined, genus_matrix, family_matrix, match_stats = taxmyphage_matrices(
        args.taxmyphage, args.prophages, all_genomes
    )
    joined.to_csv(
        outdir / "taxmyphage_genomad_joined.tsv",
        sep="\t", index=False
    )

    fam_assoc = association_table(
        family_matrix, genome_arch, "taxMyPhage", "family"
    )
    gen_assoc = association_table(
        genus_matrix, genome_arch, "taxMyPhage", "genus"
    )

    fam_assoc.to_csv(
        outdir / "taxmyphage_family_associations.tsv",
        sep="\t", index=False
    )
    gen_assoc.to_csv(
        outdir / "taxmyphage_genus_associations.tsv",
        sep="\t", index=False
    )

    combined = pd.concat(
        [inv_assoc, fam_assoc, gen_assoc],
        ignore_index=True
    )
    combined = classify_signal(
        combined,
        args.architecture_wide_threshold,
        args.subgroup_threshold
    )
    combined = add_overlap_context(combined, args.overlap_summary)

    combined.to_csv(
        outdir / "combined_all_associations.tsv",
        sep="\t", index=False
    )

    sig = combined[combined["p_fdr"] < 0.05].copy()
    sig.to_csv(
        outdir / "combined_significant_associations.tsv",
        sep="\t", index=False
    )

    arch_interp = architecture_interpretation(sig, sizes)
    arch_interp.to_csv(
        outdir / "architecture_interpretation.tsv",
        sep="\t", index=False
    )

    plot_prevalence_comparison(
        sig,
        outdir / "prevalence_comparison.png"
    )
    plot_forest(
        sig,
        outdir / "odds_ratio_forest.png"
    )
    plot_heatmap(
        combined,
        outdir / "prevalence_difference_heatmap.png"
    )

    write_summary(
        outdir,
        sizes,
        match_stats,
        sig,
        arch_interp
    )

    print(f"[DONE] Results written to {outdir}")
    print(f"[DONE] Summary: {outdir / 'analysis_summary.txt'}")
    print(f"[DONE] Significant associations: "
          f"{outdir / 'combined_significant_associations.tsv'}")


if __name__ == "__main__":
    main()
