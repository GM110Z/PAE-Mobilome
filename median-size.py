#!/usr/bin/env python3

import os
import sys
import glob
import pandas as pd
from Bio import SeqIO

if len(sys.argv) != 2:
    sys.exit("Usage: python pici_median_size_by_group.py Grouped")

base = sys.argv[1]

rows = []

# loop through each group folder
for group_dir in sorted(glob.glob(os.path.join(base, "*"))):

    if not os.path.isdir(group_dir):
        continue

    group = os.path.basename(group_dir)

    gbff_files = sorted(glob.glob(os.path.join(group_dir, "*.gbff")))

    for gbff in gbff_files:

        total_length = 0
        n_records = 0

        try:
            for record in SeqIO.parse(gbff, "genbank"):
                total_length += len(record.seq)
                n_records += 1

        except Exception as e:
            print(f"WARNING: could not parse {gbff}: {e}")
            continue

        if n_records == 0:
            print(f"WARNING: no GenBank records found in {gbff}")
            continue

        rows.append({
            "Group": group,
            "File": os.path.basename(gbff),
            "PICI_length_bp": total_length,
            "n_records_in_gbff": n_records
        })

df = pd.DataFrame(rows)

if df.empty:
    sys.exit("No GBFF files found or parsed.")

# save per-PICI sizes
per_locus_out = os.path.join(base, "pici_sizes_per_locus.tsv")
df.to_csv(per_locus_out, sep="\t", index=False)

# summary per group
summary = (
    df.groupby("Group")["PICI_length_bp"]
    .agg(["count", "median", "mean", "min", "max"])
    .reset_index()
    .rename(columns={
        "count": "n_PICIs",
        "median": "median_size_bp",
        "mean": "mean_size_bp",
        "min": "min_size_bp",
        "max": "max_size_bp"
    })
)

summary_out = os.path.join(base, "pici_size_summary_by_group.tsv")
summary.to_csv(summary_out, sep="\t", index=False)

print(f"Saved per-locus sizes: {per_locus_out}")
print(f"Saved group summary:   {summary_out}")
print()
print(summary.to_string(index=False))
