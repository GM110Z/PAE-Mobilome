#!/usr/bin/env bash
set -euo pipefail

################ CONFIG ################
THRESH=0.01
THREADS=1
PARALLEL_JOBS=24
WORKDIR=${1:-.}
BATCH_SIZE=1000
SORT_THREADS=8
SORT_MEM=6G
TMPDIR=${TMPDIR:-/tmp}
#######################################

cd "$WORKDIR"

mkdir -p batch_lists batch_sketches tmp_pairs tmp_pairs/jobs NR_genomes

find . -maxdepth 1 -type f -name "*.fna" | sort > fasta_list.txt
TOTAL=$(wc -l < fasta_list.txt)
echo "[*] Found $TOTAL genomes"

if [[ "$TOTAL" -eq 0 ]]; then
    echo "[!] No .fna files found"
    exit 1
fi

n_sketch=$(find sketches -maxdepth 1 -name "*.msh" | wc -l || true)
if [[ "$n_sketch" -ne "$TOTAL" ]]; then
    echo "[!] sketches/*.msh count ($n_sketch) does not match genomes ($TOTAL)"
    echo "    Recreate missing sketches first."
    exit 1
fi

echo "[*] Writing genome order map..."
awk '{print NR "\t" $0}' fasta_list.txt > genome_order.tsv

echo "[*] Splitting genomes into batches of $BATCH_SIZE ..."
rm -f batch_lists/list_*
split -d -l "$BATCH_SIZE" fasta_list.txt batch_lists/list_
mapfile -t LISTS < <(find batch_lists -maxdepth 1 -type f -name 'list_*' | sort)
NBATCH=${#LISTS[@]}
echo "[*] $NBATCH batches"

echo "[*] Building batch sketch files..."
for lst in "${LISTS[@]}"; do
    out="batch_sketches/$(basename "$lst")"
    if [[ ! -f "${out}.msh" ]]; then
        mapfile -t files < "$lst"
        msh_files=()
        for g in "${files[@]}"; do
            base=$(basename "$g")
            root=${base%.fna}
            msh_files+=("sketches/${root}.msh")
        done
        mash paste "$out" "${msh_files[@]}" >/dev/null 2>&1
    fi
done
echo "[*] Batch sketches ready"

echo "[*] Preparing batch-pair job list..."
rm -f tmp_pairs/pair_jobs.tsv
for ((i=0; i<NBATCH; i++)); do
    for ((j=i; j<NBATCH; j++)); do
        printf "%s\t%s\n" \
            "batch_sketches/$(basename "${LISTS[$i]}").msh" \
            "batch_sketches/$(basename "${LISTS[$j]}").msh" \
            >> tmp_pairs/pair_jobs.tsv
    done
done

pair_total=$(wc -l < tmp_pairs/pair_jobs.tsv)
echo "[*] $pair_total batch-pair jobs"

rm -f tmp_pairs/jobs/*.tsv

echo "[*] Running mash dist in parallel..."
export THRESH THREADS
awk -F '\t' '{print $1; print $2}' tmp_pairs/pair_jobs.tsv | \
xargs -P "$PARALLEL_JOBS" -n 2 bash -c '
    a="$1"
    b="$2"
    out="tmp_pairs/jobs/$(basename "$a")__$(basename "$b").tsv"
    mash dist -p "$THREADS" -d "$THRESH" "$a" "$b" > "$out"
' _

echo "[*] Concatenating pair outputs..."
cat tmp_pairs/jobs/*.tsv > tmp_pairs/pairs_raw.tsv

echo "[*] Converting close pairs to ordered indices..."
python3 << 'PY1'
path_to_idx = {}
with open("genome_order.tsv") as f:
    for line in f:
        idx, path = line.rstrip("\n").split("\t", 1)
        path_to_idx[path] = int(idx)

kept = 0
with open("tmp_pairs/pairs_raw.tsv") as fin, open("tmp_pairs/pairs_indexed.tsv", "w") as out:
    for line in fin:
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        a, b, d = parts[0], parts[1], float(parts[2])
        if a == b:
            continue
        ia = path_to_idx.get(a)
        ib = path_to_idx.get(b)
        if ia is None or ib is None:
            continue
        later, earlier = (ia, ib) if ia > ib else (ib, ia)
        out.write(f"{later}\t{earlier}\t{d:.6f}\n")
        kept += 1

print(f"[*] kept {kept} close pairs")
PY1

echo "[*] Sorting close pairs on disk..."
mkdir -p "$TMPDIR"
LC_ALL=C sort -T "$TMPDIR" --parallel="$SORT_THREADS" -S "$SORT_MEM" \
    -k1,1n -k2,2n -k3,3g tmp_pairs/pairs_indexed.tsv > close_pairs.sorted.tsv

echo "[*] Greedy clustering..."
python3 << 'PY2'
with open("fasta_list.txt") as f:
    genomes = [line.strip() for line in f if line.strip()]

n = len(genomes)
pair_fh = open("close_pairs.sorted.tsv")

def next_pair():
    line = pair_fh.readline()
    if not line:
        return None
    later, earlier, dist = line.rstrip("\n").split("\t")
    return int(later), int(earlier), float(dist)

current = next_pair()
rep_indices = set()
rep_to_cluster = {}
reps = []

with open("clusters.tsv", "w") as out:
    out.write("cluster_id\trepresentative\tmember\tdistance_to_rep\n")

    for i in range(1, n + 1):
        if i % 1000 == 0:
            print(f"[*] processed {i}/{n}")
        g = genomes[i - 1]
        neighbors = []

        while current is not None and current[0] == i:
            _, earlier, dist = current
            neighbors.append((earlier, dist))
            current = next_pair()

        if i == 1:
            rep_indices.add(i)
            rep_to_cluster[i] = 1
            reps.append(i)
            out.write(f"1\t{g}\t{g}\t0.00000\n")
            continue

        best_rep = None
        best_dist = None
        for earlier, dist in neighbors:
            if earlier in rep_indices:
                if best_dist is None or dist < best_dist:
                    best_dist = dist
                    best_rep = earlier

        if best_rep is None:
            cid = len(reps) + 1
            rep_indices.add(i)
            rep_to_cluster[i] = cid
            reps.append(i)
            out.write(f"{cid}\t{g}\t{g}\t0.00000\n")
        else:
            cid = rep_to_cluster[best_rep]
            rep_g = genomes[best_rep - 1]
            out.write(f"{cid}\t{rep_g}\t{g}\t{best_dist:.5f}\n")

with open("representatives.txt", "w") as out:
    for idx in reps:
        out.write(genomes[idx - 1] + "\n")

pair_fh.close()
print(f"[*] {n} genomes -> {len(reps)} clusters")
PY2

echo "[*] Copying representatives..."
rep_total=$(wc -l < representatives.txt)
i=0
while read -r g; do
    cp "$g" "NR_genomes/$(basename "$g")"
    i=$((i+1))
    printf "\rCopying reps: %d/%d" "$i" "$rep_total"
done < representatives.txt
echo
echo "[✓] Done"
