#!/usr/bin/env bash

# ==========================================================
# sniPICI.sh
#
# Extract predicted PICI loci as mini-GBFF files
# Supports:
#   - single-record complete genomes
#   - multi-contig draft genomes
#   - warns if first_hit / last_hit are on different contigs
#
# Input:
#   pici_systems_summary.tsv
#   matching *.gbff
#
# Output:
#   pici_gbff_regions/
#      *.gbff
#      extraction_summary.tsv
# ==========================================================

SUMMARY="pici_systems_summary.tsv"
OUTDIR="pici_gbff_regions"
FLANK=3

mkdir -p "$OUTDIR"

python3 <<'PY'
import csv, os
from Bio import SeqIO

summary = "pici_systems_summary.tsv"
outdir  = "pici_gbff_regions"
flank   = 3

os.makedirs(outdir, exist_ok=True)

sumfile = open(f"{outdir}/extraction_summary.tsv", "w")
sumfile.write(
    "Genome\tsys_id\tstatus\tcontig\tfirst_hit\tlast_hit\t"
    "cds_start\tcds_end\tbp_start\tbp_end\tn_cds\toutfile\n"
)

with open(summary) as fh:
    reader = csv.DictReader(fh, delimiter="\t")

    for row in reader:

        genome = row["Genome"].replace(".g.fasta","")
        sysid  = row["sys_id"]

        gbff = genome + ".gbff"

        if not os.path.exists(gbff):
            print("Missing:", gbff)
            continue

        first_hit = row["first_hit"]
        last_hit  = row["last_hit"]

        try:
            records = list(SeqIO.parse(gbff, "genbank"))
        except Exception as e:
            print("Cannot read", gbff, e)
            continue

        found_first = None
        found_last  = None

        # search every contig
        for rec in records:

            cds = [f for f in rec.features if f.type == "CDS"]

            for i, f in enumerate(cds):
                pid = f.qualifiers.get("protein_id", [""])[0]

                if pid == first_hit and found_first is None:
                    found_first = (rec, cds, i)

                if pid == last_hit and found_last is None:
                    found_last = (rec, cds, i)

        if found_first is None or found_last is None:
            print("Protein match failed:", genome, first_hit, last_hit)
            continue

        rec1, cds1, idx1 = found_first
        rec2, cds2, idx2 = found_last

        # hits on different contigs
        if rec1.id != rec2.id:
            print("Different contigs:", genome, first_hit, last_hit)
            sumfile.write(
                f"{genome}\t{sysid}\tDIFFERENT_CONTIGS\tNA\t"
                f"{first_hit}\t{last_hit}\tNA\tNA\tNA\tNA\tNA\tNA\n"
            )
            continue

        rec = rec1
        cds = cds1

        if idx1 > idx2:
            idx1, idx2 = idx2, idx1
            first_hit, last_hit = last_hit, first_hit

        s = max(0, idx1 - flank)
        e = min(len(cds)-1, idx2 + flank)

        chosen = cds[s:e+1]

        bp_start = min(int(f.location.start) for f in chosen)
        bp_end   = max(int(f.location.end) for f in chosen)

        sub = rec[bp_start:bp_end]

        new_feats = []
        for f in chosen:
            nf = f._shift(-bp_start)
            new_feats.append(nf)

        sub.features = new_feats
        sub.annotations["molecule_type"] = "DNA"
        sub.id = rec.id
        sub.name = rec.name
        sub.description = f"{rec.description} | {sysid} predicted PICI region"

        outfile = f"{outdir}/{genome}_{sysid}.gbff"

        SeqIO.write(sub, outfile, "genbank")

        sumfile.write(
            f"{genome}\t{sysid}\tOK\t{rec.id}\t{first_hit}\t{last_hit}\t"
            f"{s+1}\t{e+1}\t{bp_start+1}\t{bp_end}\t{len(chosen)}\t{outfile}\n"
        )

sumfile.close()
print("Done.")
PY
