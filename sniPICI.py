#!/usr/bin/env python3
"""

Usage:
    python3 02_extract_pici_regions.py \
        --summary pici_systems_summary.tsv \
        --gbff_dir /path/to/gbff_files \
        --outdir pici_gbff_regions \
        --flank 5
"""

import argparse
import csv
import os
import sys
from pathlib import Path
from Bio import SeqIO

# ── argument parsing ──────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--summary",  default="pici_systems_summary.tsv",
                   help="Output of 01_build_pici_summary.sh")
    p.add_argument("--gbff_dir", default=".",
                   help="Directory containing .gbff files (searched recursively)")
    p.add_argument("--outdir",   default="pici_gbff_regions")
    p.add_argument("--flank",    type=int, default=5,
                   help="CDS units to add upstream/downstream (default: 5)")
    return p.parse_args()

# ── helpers ───────────────────────────────────────────────────────────────────

def find_gbff(gbff_dir: str, genome: str) -> Path | None:
    """Search for a GBFF file matching the genome accession."""
    root = Path(gbff_dir)
    # Try common naming patterns
    for suffix in [".gbff", ".gbk", ".gb", "_genomic.gbff"]:
        for candidate in root.rglob(f"*{genome}*{suffix}"):
            return candidate
    return None


def build_protein_index(records):
    """
    Return dict: protein_id -> (record, cds_list, cds_index)
    Indexes ALL CDS across ALL contigs in one pass.
    """
    index = {}
    for rec in records:
        cds = [f for f in rec.features if f.type == "CDS"]
        for i, f in enumerate(cds):
            pid = f.qualifiers.get("protein_id", [""])[0]
            if pid:
                index[pid] = (rec, cds, i)
    return index


def extract_region(rec, cds_list, idx_start, idx_end, flank):
    """
    Slice a SeqRecord around CDS[idx_start:idx_end], adding flank CDS
    on each side. Uses SeqRecord slicing to avoid private _shift calls.
    Returns (sub_record, bp_start, bp_end, n_cds).
    """
    s = max(0, idx_start - flank)
    e = min(len(cds_list) - 1, idx_end + flank)
    chosen = cds_list[s : e + 1]

    bp_start = min(int(f.location.start) for f in chosen)
    bp_end   = max(int(f.location.end)   for f in chosen)

    # SeqRecord slicing correctly shifts all feature coordinates
    sub = rec[bp_start:bp_end]
    sub.annotations["molecule_type"] = "DNA"
    sub.id          = rec.id
    sub.name        = rec.name
    sub.description = rec.description

    return sub, bp_start, bp_end, len(chosen)


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    log_path = os.path.join(args.outdir, "extraction_summary.tsv")
    log = open(log_path, "w")
    log.write(
        "Genome\tsys_id\tstatus\tcontig\t"
        "first_hit\tlast_hit\t"
        "cds_start_idx\tcds_end_idx\t"
        "bp_start\tbp_end\tbp_length\tn_cds\toutfile\n"
    )

    ok_count = fail_count = split_count = 0

    with open(args.summary) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            genome    = row["Genome"]
            sys_id    = row["sys_id"]
            first_hit = row["first_hit"]
            last_hit  = row["last_hit"]

            # ── locate GBFF ──────────────────────────────────────────────────
            gbff_path = find_gbff(args.gbff_dir, genome)
            if gbff_path is None:
                print(f"[MISSING GBFF] {genome}", file=sys.stderr)
                log.write(f"{genome}\t{sys_id}\tMISSING_GBFF\tNA\t"
                          f"{first_hit}\t{last_hit}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n")
                fail_count += 1
                continue

            try:
                records = list(SeqIO.parse(str(gbff_path), "genbank"))
            except Exception as exc:
                print(f"[PARSE ERROR] {gbff_path}: {exc}", file=sys.stderr)
                fail_count += 1
                continue

            # ── build protein index ──────────────────────────────────────────
            prot_index = build_protein_index(records)

            hit1 = prot_index.get(first_hit)
            hit2 = prot_index.get(last_hit)

            if hit1 is None or hit2 is None:
                missing = [h for h, v in [(first_hit, hit1), (last_hit, hit2)] if v is None]
                print(f"[PROTEIN NOT FOUND] {genome} {sys_id}: {missing}", file=sys.stderr)
                log.write(f"{genome}\t{sys_id}\tPROTEIN_NOT_FOUND\tNA\t"
                          f"{first_hit}\t{last_hit}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n")
                fail_count += 1
                continue

            rec1, cds1, idx1 = hit1
            rec2, cds2, idx2 = hit2

            # ── same contig: normal extraction ───────────────────────────────
            if rec1.id == rec2.id:
                if idx1 > idx2:
                    idx1, idx2 = idx2, idx1
                    first_hit, last_hit = last_hit, first_hit

                sub, bp_start, bp_end, n_cds = extract_region(
                    rec1, cds1, idx1, idx2, args.flank
                )
                sub.description += f" | {sys_id} predicted PICI region"

                outfile = os.path.join(args.outdir, f"{genome}_{sys_id}.gbff")
                SeqIO.write(sub, outfile, "genbank")

                log.write(
                    f"{genome}\t{sys_id}\tOK\t{rec1.id}\t"
                    f"{first_hit}\t{last_hit}\t"
                    f"{idx1+1}\t{idx2+1}\t"
                    f"{bp_start+1}\t{bp_end}\t{bp_end - bp_start}\t{n_cds}\t{outfile}\n"
                )
                ok_count += 1

            # ── different contigs: extract each separately ───────────────────
            else:
                print(f"[SPLIT CONTIGS] {genome} {sys_id}: "
                      f"{rec1.id} / {rec2.id}", file=sys.stderr)
                split_count += 1

                for label, (rec, cds, idx) in [
                    ("partA", (rec1, cds1, idx1)),
                    ("partB", (rec2, cds2, idx2)),
                ]:
                    sub, bp_start, bp_end, n_cds = extract_region(
                        rec, cds, idx, idx, args.flank
                    )
                    sub.description += (
                        f" | {sys_id} predicted PICI region "
                        f"[SPLIT {label}]"
                    )
                    outfile = os.path.join(
                        args.outdir, f"{genome}_{sys_id}_{label}.gbff"
                    )
                    SeqIO.write(sub, outfile, "genbank")

                    log.write(
                        f"{genome}\t{sys_id}\tSPLIT_{label}\t{rec.id}\t"
                        f"{first_hit}\t{last_hit}\t"
                        f"{idx+1}\t{idx+1}\t"
                        f"{bp_start+1}\t{bp_end}\t{bp_end - bp_start}\t{n_cds}\t{outfile}\n"
                    )

    log.close()

    print(f"\nDone.")
    print(f"  OK (same contig):   {ok_count}")
    print(f"  Split (diff contig):{split_count}")
    print(f"  Failed:             {fail_count}")
    print(f"  Log: {log_path}")


if __name__ == "__main__":
    main()
