#!/usr/bin/env bash


SYSTEMS="pici_systems_summary.tsv"
PA="pici_pa_matrix.tsv"

# MacSyFinder v2 best_solution.tsv columns (1-indexed):
#  1  replicon
#  2  hit_id         (protein_id — use this for GBFF lookup)
#  3  gene_name      (e.g. integrase, terminase, ...)
#  4  hit_pos        (ordinal CDS position on replicon)
#  5  hit_status     (mandatory / accessory / neutral)
#  6  hit_seq_len
#  7  hit_i_eval
#  8  hit_score
#  9  hit_profile_cov
# 10  hit_seq_cov
# 11  hit_begin_match
# 12  hit_end_match
# 13  counterpart     (may be empty — shifts remaining cols by 1 in older versions)
# 14  sys_id
# 15  sys_locus_id
# 16  sys_wholeness
# 17  sys_score
# 18  sys_occ_nb
# 19  model_fqn
# 20  hit_gene_ref
# 21  hit_system_ref
# 22  used_in


echo -e "Genome\tsys_id\tfirst_hit\tlast_hit\tfirst_pos\tlast_pos\tn_markers\tsys_score\tsys_wholeness" \
    > "$SYSTEMS"

for d in sat_*/; do
    [ -d "$d" ] || continue
    file="${d}best_solution.tsv"
    [ -s "$file" ] || continue

    # Strip known prefixes to get clean genome accession
    genome=$(basename "$d")
    genome="${genome#sat_PICI_}"
    genome="${genome#sat_cfPICI_}"
    genome="${genome#sat_P4_}"
    genome="${genome#sat_PLE_}"

    awk -F '\t' -v g="$genome" '
    BEGIN { OFS="\t"; sys_col=-1; score_col=-1; whole_col=-1 }

    # Skip comment lines
    /^#/ { next }

    # Parse header to find column positions dynamically
    # (guards against minor version differences)
    NR==1 || $1=="replicon" {
        for (i=1; i<=NF; i++) {
            if ($i == "sys_id")        sys_col   = i
            if ($i == "sys_score")     score_col = i
            if ($i == "sys_wholeness") whole_col = i
            if ($i == "hit_id")        hit_col   = i
            if ($i == "hit_pos")       pos_col   = i
        }
        next
    }

    # Skip rows without enough fields
    NF < 10 { next }

    {
        if (sys_col < 0) {
            # Fallback to v2 defaults if header was absent
            sys_col=14; score_col=17; whole_col=16; hit_col=2; pos_col=4
        }
        sys   = $sys_col
        hit   = $hit_col
        pos   = $pos_col + 0
        score = $score_col
        whole = $whole_col

        n[sys]++
        if (!(sys in minpos) || pos < minpos[sys]) {
            minpos[sys]   = pos
            first_hit[sys] = hit
        }
        if (!(sys in maxpos) || pos > maxpos[sys]) {
            maxpos[sys]   = pos
            last_hit[sys]  = hit
        }
        sys_score[sys] = score
        sys_whole[sys] = whole
    }

    END {
        for (sys in n) {
            print g, sys, first_hit[sys], last_hit[sys], \
                  minpos[sys], maxpos[sys], n[sys], \
                  sys_score[sys], sys_whole[sys]
        }
    }' "$file"

done | sort >> "$SYSTEMS"

echo "Written: $SYSTEMS"

# ---------- presence/absence matrix ----------
echo -e "Genome\tPICI_present\tn_systems" > "$PA"

awk -F '\t' '
NR==1 { next }
{
    genome = $1
    seen[genome][$2] = 1
}
END {
    for (g in seen) {
        c = 0
        for (s in seen[g]) c++
        print g "\t1\t" c
    }
}' "$SYSTEMS" | sort >> "$PA"

# Add genomes with zero systems
for d in sat_*/; do
    [ -d "$d" ] || continue
    genome=$(basename "$d")
    genome="${genome#sat_PICI_}"
    genome="${genome#sat_cfPICI_}"
    genome="${genome#sat_P4_}"
    genome="${genome#sat_PLE_}"
    grep -q "^${genome}[[:space:]]" "$PA" || \
        echo -e "${genome}\t0\t0" >> "$PA"
done

sort -o "$PA" "$PA"

echo "Written: $PA"
echo "Done."
