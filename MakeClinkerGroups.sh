#!/usr/bin/env bash

set -euo pipefail

if [ "$#" -lt 4 ]; then
    echo "Usage: $0 REPS MASH OUTDIR MAX_DIST"
    exit 1
fi

REPS="$1"
MASH="$2"
OUTDIR="$3"
MAX_DIST="$4"

# secondary filter to break weak chaining
MIN_SHARED=3

echo "Using:"
echo "  REPS       = $REPS"
echo "  MASH       = $MASH"
echo "  OUTDIR     = $OUTDIR"
echo "  MAX_DIST   = $MAX_DIST"
echo "  MIN_SHARED = $MIN_SHARED"

mkdir -p "$OUTDIR"

python3 <<PY
import os, shutil
from collections import defaultdict

reps_file = "$REPS"
mash_file = "$MASH"
outdir = "$OUTDIR"
max_dist = float("$MAX_DIST")
min_shared = int("$MIN_SHARED")

# ---------- load representatives ----------
reps=set()

with open(reps_file) as f:
    next(f)
    for line in f:
        line=line.strip()
        if not line: continue
        name=line.split("\\t")[0]
        base=os.path.basename(name).replace(".fna","")
        reps.add(base)

print(f"Loaded {len(reps)} reps")

# ---------- build adjacency ----------
adj=defaultdict(set)

def clean(x):
    return os.path.basename(x).replace(".fna","").replace(".msh","").replace(".gbff","")

def shared(x):
    try:
        return int(x.split("/")[0])
    except:
        return 0

with open(mash_file) as f:
    for line in f:
        if not line.strip(): continue
        parts=line.strip().split("\\t")
        if len(parts)<5: continue

        a,b,dist,pval,sh=parts[:5]

        a=clean(a)
        b=clean(b)

        if a not in reps or b not in reps: continue
        if a==b: continue

        dist=float(dist)

        # IMPORTANT: dual filter
        if dist <= max_dist and shared(sh) >= min_shared:
            adj[a].add(b)
            adj[b].add(a)

# ---------- greedy clustering (NO chaining explosion) ----------
unassigned=set(reps)
groups=[]

while unassigned:
    seed=unassigned.pop()
    group={seed}

    # only keep nodes directly connected to seed OR strongly connected
    neighbors=adj[seed].copy()

    for n in list(neighbors):
        # require mutual connectivity (stronger than simple chaining)
        if seed in adj[n]:
            group.add(n)

    # remove assigned
    for x in group:
        unassigned.discard(x)

    groups.append(sorted(group))

print(f"Built {len(groups)} groups")

# ---------- write output ----------
for i,comp in enumerate(groups,1):
    gdir=os.path.join(outdir,f"group_{i}")
    os.makedirs(gdir,exist_ok=True)

    for locus in comp:
        gbff=f"{locus}.gbff"
        if os.path.exists(gbff):
            shutil.copy(gbff,gdir)

print("Done.")
PY
