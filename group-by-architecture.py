#!/usr/bin/env python3
"""
group_by_architecture.py

Organise extracted GBFF files into subdirectories by architecture class
using the master_completeness_summary.tsv produced by PICIasso.

Creates:
    <outdir>/
        complete_canonical/
            GCF_XXXXX_PICI_1.gbff
            GCF_XXXXX_PICI_1.faa      (if present)
            GCF_XXXXX_PICI_1.pfam     (if present)
        replication_mobilisation/
            ...
        backbone_replication_no_terminase/
            ...
        ...

Usage:
    python group_by_architecture.py \
        --summary CLaudEanalysis/master_completeness_summary.tsv \
        --gbff_dir CLaudEanalysis \
        --outdir CLaudEanalysis_grouped
"""

import argparse
import os
import shutil
import pandas as pd
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--summary",
                   default="CLaudEanalysis/master_completeness_summary.tsv",
                   help="master_completeness_summary.tsv from PICIasso")
    p.add_argument("--gbff_dir",
                   default="CLaudEanalysis",
                   help="Directory containing the GBFF files")
    p.add_argument("--outdir",
                   default="CLaudEanalysis_grouped",
                   help="Output directory for grouped files")
    p.add_argument("--copy", action="store_true",
                   help="Copy files instead of creating symlinks (default: symlinks)")
    p.add_argument("--extensions", default="gbff,faa,pfam",
                   help="Comma-separated file extensions to group (default: gbff,faa,pfam)")
    return p.parse_args()


def main():
    args = parse_args()
    extensions = [f".{e.lstrip('.')}" for e in args.extensions.split(",")]

    # Load summary
    df = pd.read_csv(args.summary, sep="\t")
    print(f"Loaded {len(df)} elements from {args.summary}")

    # Strip group prefix from Genome column if present (e.g. "all/GCF_XXX" -> "GCF_XXX")
    df["bare_name"] = df["Genome"].apply(lambda x: x.split("/", 1)[-1])

    # Count per architecture
    arch_counts = df["architecture"].value_counts()
    print("\nArchitecture distribution:")
    for arch, n in arch_counts.items():
        print(f"  {arch:45s} {n:4d}")

    # Create output directories
    os.makedirs(args.outdir, exist_ok=True)
    for arch in df["architecture"].unique():
        os.makedirs(os.path.join(args.outdir, arch), exist_ok=True)

    # Group files
    moved   = 0
    missing = 0
    counts  = defaultdict(int)

    for _, row in df.iterrows():
        arch = row["architecture"]
        name = row["bare_name"]

        for ext in extensions:
            src = os.path.join(args.gbff_dir, f"{name}{ext}")
            if not os.path.exists(src):
                # Try without extension variants
                continue

            dst = os.path.join(args.outdir, arch, f"{name}{ext}")

            if args.copy:
                shutil.copy2(src, dst)
            else:
                # Symlink with absolute path so links work from any directory
                src_abs = os.path.abspath(src)
                if os.path.exists(dst) or os.path.islink(dst):
                    os.remove(dst)
                os.symlink(src_abs, dst)

            if ext == ".gbff":
                moved += 1
                counts[arch] += 1

    print(f"\n{'Copied' if args.copy else 'Linked'} files into {args.outdir}/")
    print(f"  Total GBFFs placed: {moved}")
    print(f"  Missing GBFFs:      {missing}")
    print("\nPer-architecture GBFF counts:")
    for arch, n in sorted(counts.items(), key=lambda x: -x[1]):
        print(f"  {arch:45s} {n:4d}")

    # Write a summary per architecture folder
    for arch in df["architecture"].unique():
        arch_df = df[df["architecture"] == arch]
        out_tsv = os.path.join(args.outdir, arch, "elements.tsv")
        arch_df.to_csv(out_tsv, sep="\t", index=False)

    print(f"\nDone. Open any folder in Artemis or use with PICIasso --regraph.")
    print(f"To load a group schematic in Artemis:")
    print(f"  File -> Open -> select any .gbff in {args.outdir}/<architecture>/")


if __name__ == "__main__":
    main()
