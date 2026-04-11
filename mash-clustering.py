#!/usr/bin/env python3
import sys
import os

all_genomes_file = sys.argv[1]
mash_edges_file = sys.argv[2]
mash_threshold = float(sys.argv[3])
output_file = sys.argv[4]

parent = {}
rank = {}
genomes = []

def norm(x):
    x = x.strip()
    x = os.path.basename(x)
    x = x.replace(".fasta", "").replace(".fna", "")
    return x

def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x

def union(a, b):
    ra = find(a)
    rb = find(b)
    if ra == rb:
        return
    if rank[ra] < rank[rb]:
        parent[ra] = rb
    elif rank[ra] > rank[rb]:
        parent[rb] = ra
    else:
        parent[rb] = ra
        rank[ra] += 1

# Initialize every genome so none are lost
with open(all_genomes_file) as f:
    for line in f:
        g = norm(line)
        if not g:
            continue
        genomes.append(g)
        parent[g] = g
        rank[g] = 0

# Merge only close-enough pairs
with open(mash_edges_file) as f:
    for i, line in enumerate(f, 1):
        g1, g2, dist, pval, shared = line.rstrip("\n").split("\t")
        if float(dist) <= mash_threshold:
            union(norm(g1), norm(g2))
        if i % 10000000 == 0:
            print(f"Processed {i:,} lines", flush=True)

cluster_roots = {}
cluster_id = 0

with open(output_file, "w") as out:
    out.write("genome\tcluster_id\n")
    for g in sorted(genomes):
        root = find(g)
        if root not in cluster_roots:
            cluster_id += 1
            cluster_roots[root] = cluster_id
        out.write(f"{g}\t{cluster_roots[root]}\n")

print(f"Done. {len(genomes):,} genomes assigned to {cluster_id:,} clusters.")
