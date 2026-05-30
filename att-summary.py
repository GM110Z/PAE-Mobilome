#!/usr/bin/env python3
"""
att_site_summary.py

Summarise and cluster att site sequences from att site detection output.
Compatible with both att_site_finder.py and att_find_attB_first.py outputs.

Groups att sites by:
  1. Sequence identity (clusters identical/near-identical att sites)
  2. Architecture class
  3. Integration locus (tRNA type)

Shows how many genomes share each att site sequence and which
architecture classes use which integration loci.

Usage:
    # From att_find_attB_first.py (recommended)
    python att_site_summary.py \
        --input att_attB/att_attB_found.tsv \
        --outdir att_site_summary

    # From att_site_finder.py
    python att_site_summary.py \
        --input att_sites_confirmed/att_sites_found.tsv \
        --outdir att_site_summary
"""

import argparse
import os
from collections import defaultdict

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input",  required=True,
                   help="att_sites_found.tsv from att_site_finder.py")
    p.add_argument("--outdir", default="att_site_summary")
    p.add_argument("--min_att_len", type=int, default=15,
                   help="Minimum att site length to consider high-confidence "
                        "(default: 15)")
    return p.parse_args()


def cluster_att_sequences(sequences: list[str],
                           min_id: float = 0.85) -> dict[str, str]:
    """
    Cluster att sequences by identity.
    Returns dict: sequence -> cluster_representative
    """
    clusters = {}
    reps     = []

    for seq in sequences:
        seq = seq.upper()
        assigned = False
        for rep in reps:
            if len(seq) != len(rep):
                continue
            matches  = sum(a == b for a, b in zip(seq, rep))
            identity = matches / len(seq)
            if identity >= min_id:
                clusters[seq] = rep
                assigned = True
                break
        if not assigned:
            reps.append(seq)
            clusters[seq] = seq

    return clusters


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.input, sep="\t")
    print(f"Loaded {len(df)} elements with att sites")
    print(f"Architecture classes: {df['architecture'].nunique()}")

    # Normalise column names — handle both att_site_finder.py and
    # att_find_attB_first.py output formats
    if "trna_at_boundary" not in df.columns:
        df["trna_at_boundary"] = False
    else:
        # Ensure boolean — pandas may read True/False as strings from TSV
        df["trna_at_boundary"] = df["trna_at_boundary"].map(
            lambda x: str(x).strip().lower() in ("true", "1", "yes")
        )
    if "trna_product" not in df.columns:
        df["trna_product"] = None
    if "empty_site_confirmed" not in df.columns:
        if "attB_genome" in df.columns:
            df["empty_site_confirmed"] = df["attB_genome"].notna()
        else:
            df["empty_site_confirmed"] = False
    else:
        df["empty_site_confirmed"] = df["empty_site_confirmed"].map(
            lambda x: str(x).strip().lower() in ("true", "1", "yes")
        )
    if "att_identity" not in df.columns:
        df["att_identity"] = None
    if "genome" not in df.columns:
        df["genome"] = df["element"].str.rsplit("_UserReplicon", n=1).str[0]

    # ── Overall stats ─────────────────────────────────────────────────────────
    print(f"\natt site length distribution:")
    print(df["att_length"].value_counts().sort_index().to_string())

    print(f"\natt sites ≥{args.min_att_len} bp: "
          f"{(df['att_length'] >= args.min_att_len).sum()}")
    print(f"att sites with tRNA at boundary: {df['trna_at_boundary'].sum()}")
    print(f"att sites confirmed (empty site): {df['empty_site_confirmed'].sum()}")

    # ── Cluster att sequences ─────────────────────────────────────────────────
    seqs    = df["att_sequence"].dropna().str.upper().tolist()
    cluster_map = cluster_att_sequences(seqs)
    df["att_cluster"] = df["att_sequence"].str.upper().map(cluster_map)

    # Count genomes per cluster
    cluster_counts = df.groupby("att_cluster").agg(
        n_genomes       = ("genome",          "nunique"),
        n_elements      = ("element",         "count"),
        architectures   = ("architecture",    lambda x: "; ".join(sorted(set(x)))),
        trna_types      = ("trna_product",    lambda x: "; ".join(
                               sorted(set(str(v) for v in x if pd.notna(v))))),
        mean_length     = ("att_length",      "mean"),
        mean_identity   = ("att_identity",    "mean"),
        n_confirmed     = ("empty_site_confirmed", "sum"),
        n_tRNA_boundary = ("trna_at_boundary", "sum"),
    ).round(2).sort_values("n_genomes", ascending=False)

    cluster_tsv = os.path.join(args.outdir, "att_clusters.tsv")
    cluster_counts.to_csv(cluster_tsv, sep="\t")
    print(f"\n── att site clusters ────────────────────────────────")
    print(f"Total unique clusters: {len(cluster_counts)}")
    print(f"\nTop 20 clusters by genome count:")
    print(cluster_counts.head(20).to_string())
    print(f"\n→ {cluster_tsv}")

    # ── Per-architecture breakdown ────────────────────────────────────────────
    arch_breakdown = df.groupby(["architecture", "att_cluster"]).agg(
        n_genomes = ("genome", "nunique"),
        trna      = ("trna_product", lambda x: "; ".join(
                        sorted(set(str(v) for v in x if pd.notna(v))))),
        confirmed = ("empty_site_confirmed", "sum"),
    ).reset_index().sort_values(["architecture", "n_genomes"],
                                 ascending=[True, False])

    arch_tsv = os.path.join(args.outdir, "att_by_architecture.tsv")
    arch_breakdown.to_csv(arch_tsv, sep="\t", index=False)
    print(f"\n── Per-architecture att site breakdown ──────────────")
    for arch, grp in arch_breakdown.groupby("architecture"):
        print(f"\n  {arch} ({grp['n_genomes'].sum()} elements with att sites):")
        for _, row in grp.iterrows():
            print(f"    [{row['n_genomes']:3d} genomes] "
                  f"{row['att_cluster'][:30]:30s} "
                  f"tRNA={row['trna'] or 'none':15s} "
                  f"confirmed={int(row['confirmed'])}")
    print(f"\n→ {arch_tsv}")

    # ── Integration locus summary ─────────────────────────────────────────────
    trna_summary = df[df["trna_at_boundary"]].groupby(
        ["trna_product", "architecture"]
    ).agg(
        n_elements  = ("element",  "count"),
        n_att_seqs  = ("att_cluster", "nunique"),
    ).reset_index().sort_values("n_elements", ascending=False)

    trna_tsv = os.path.join(args.outdir, "integration_loci.tsv")
    trna_summary.to_csv(trna_tsv, sep="\t", index=False)
    print(f"\n── Integration loci (tRNA) ──────────────────────────")
    print(trna_summary.to_string(index=False))
    print(f"\n→ {trna_tsv}")

    # ── Plots ─────────────────────────────────────────────────────────────────

    # 1. Heatmap: architecture vs att cluster (top 20 clusters)
    top_clusters = cluster_counts.head(20).index.tolist()
    archs        = sorted(df["architecture"].unique())
    matrix       = pd.DataFrame(0, index=archs, columns=top_clusters)

    for _, row in df.iterrows():
        if row["att_cluster"] in top_clusters:
            matrix.loc[row["architecture"], row["att_cluster"]] += 1

    fig, ax = plt.subplots(figsize=(max(14, len(top_clusters) * 0.7),
                                    max(6,  len(archs) * 0.5)))
    im = ax.imshow(matrix.values, aspect="auto", cmap="Blues")
    ax.set_xticks(range(len(top_clusters)))
    ax.set_xticklabels([s[:20] for s in top_clusters],
                        rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(archs)))
    ax.set_yticklabels(archs, fontsize=8)
    plt.colorbar(im, ax=ax, label="Number of elements")
    ax.set_title("att site clusters by architecture class\n(top 20 clusters)")
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "att_cluster_heatmap.png"),
                dpi=150, bbox_inches="tight")
    plt.close()

    # 2. Bar chart: top att clusters coloured by architecture
    fig, ax = plt.subplots(figsize=(14, 5))
    top20 = cluster_counts.head(20)
    colors = plt.cm.tab20(np.linspace(0, 1, len(archs)))
    arch_color = dict(zip(archs, colors))

    bottoms = np.zeros(len(top20))
    for arch in archs:
        vals = []
        for cluster in top20.index:
            n = len(df[(df["att_cluster"] == cluster) &
                       (df["architecture"] == arch)])
            vals.append(n)
        vals = np.array(vals)
        if vals.sum() > 0:
            ax.bar(range(len(top20)), vals, bottom=bottoms,
                   label=arch, color=arch_color[arch])
            bottoms += vals

    ax.set_xticks(range(len(top20)))
    ax.set_xticklabels([s[:20] for s in top20.index],
                        rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("Number of elements")
    ax.set_title("Top 20 att site clusters — composition by architecture")
    ax.legend(fontsize=7, bbox_to_anchor=(1.01, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "att_cluster_by_arch.png"),
                dpi=150, bbox_inches="tight")
    plt.close()

    # 3. Integration locus pie chart
    if df["trna_at_boundary"].sum() > 0:
        trna_counts = df[df["trna_at_boundary"]]["trna_product"].value_counts()
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.pie(trna_counts.values,
               labels=trna_counts.index,
               autopct="%1.1f%%",
               startangle=90)
        ax.set_title("Integration loci (tRNA type)")
        plt.tight_layout()
        plt.savefig(os.path.join(args.outdir, "integration_loci_pie.png"),
                    dpi=150, bbox_inches="tight")
        plt.close()

    print(f"\nPlots saved to {args.outdir}/")
    print("Done.")


if __name__ == "__main__":
    main()
