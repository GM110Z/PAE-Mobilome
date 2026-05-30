#!/usr/bin/env python3

import os
import sys
import glob
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO

if len(sys.argv) != 2:
    sys.exit("Usage: python pici_size_box_violin_from_grouped_gbff.py /path/to/Grouped")

base = sys.argv[1]

rows = []

for group_dir in sorted(glob.glob(os.path.join(base, "*"))):
    if not os.path.isdir(group_dir):
        continue

    group = os.path.basename(group_dir)

    for gbff in sorted(glob.glob(os.path.join(group_dir, "*.gbff"))):
        records = list(SeqIO.parse(gbff, "genbank"))

        if len(records) == 0:
            print(f"WARNING: no records found in {gbff}")
            continue

        length = sum(len(record.seq) for record in records)

        rows.append({
            "Group": group,
            "File": os.path.basename(gbff),
            "PICI_size_bp": length,
            "n_records": len(records)
        })

df = pd.DataFrame(rows)

if df.empty:
    sys.exit("No .gbff files found. Check the path to your Grouped folder.")

# save per-locus table
per_locus_out = os.path.join(base, "pici_sizes_per_locus.tsv")
df.to_csv(per_locus_out, sep="\t", index=False)

# summary
summary = (
    df.groupby("Group")["PICI_size_bp"]
    .agg(["count", "median", "mean", "min", "max"])
    .reset_index()
    .rename(columns={
        "count": "n_PICIs",
        "median": "median_size_bp",
        "mean": "mean_size_bp",
        "min": "min_size_bp",
        "max": "max_size_bp"
    })
    .sort_values("median_size_bp")
)

summary_out = os.path.join(base, "pici_size_summary_by_group.tsv")
summary.to_csv(summary_out, sep="\t", index=False)

group_order = summary["Group"].tolist()

plot_data = [
    df.loc[df["Group"] == group, "PICI_size_bp"].values
    for group in group_order
]

# -------------------------
# BOXPLOT
# -------------------------
fig, ax = plt.subplots(figsize=(12, 6))

ax.boxplot(plot_data, labels=group_order, showfliers=True)

ax.set_xlabel("PICI architecture group")
ax.set_ylabel("PICI size (bp)")
ax.set_title("PICI size distribution by architecture group")

plt.xticks(rotation=45, ha="right")
plt.tight_layout()

box_png = os.path.join(base, "pici_size_boxplot_by_group.png")
box_pdf = os.path.join(base, "pici_size_boxplot_by_group.pdf")

plt.savefig(box_png, dpi=300)
plt.savefig(box_pdf)
plt.close()

# -------------------------
# VIOLIN PLOT
# -------------------------
fig, ax = plt.subplots(figsize=(12, 6))

ax.violinplot(
    plot_data,
    showmeans=False,
    showmedians=True,
    showextrema=True
)

ax.set_xticks(range(1, len(group_order) + 1))
ax.set_xticklabels(group_order, rotation=45, ha="right")

ax.set_xlabel("PICI architecture group")
ax.set_ylabel("PICI size (bp)")
ax.set_title("PICI size distribution by architecture group")

plt.tight_layout()

violin_png = os.path.join(base, "pici_size_violin_by_group.png")
violin_pdf = os.path.join(base, "pici_size_violin_by_group.pdf")

plt.savefig(violin_png, dpi=300)
plt.savefig(violin_pdf)
plt.close()

print(f"Saved per-locus table: {per_locus_out}")
print(f"Saved summary table:   {summary_out}")
print(f"Saved boxplot PNG:     {box_png}")
print(f"Saved boxplot PDF:     {box_pdf}")
print(f"Saved violin PNG:      {violin_png}")
print(f"Saved violin PDF:      {violin_pdf}")
print()
print(summary.to_string(index=False))
