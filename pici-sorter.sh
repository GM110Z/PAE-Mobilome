#!/usr/bin/env bash


SYSTEMS="pici_systems_summary.tsv"
PA="pici_pa_matrix.tsv"

# ---------- detailed systems table ----------

echo -e "Genome\tsys_id\tfirst_hit\tlast_hit\tfirst_pos\tlast_pos\tn_markers\tsys_score\tsys_wholeness" > "$SYSTEMS"

for d in sat_*; do
    [ -d "$d" ] || continue

    file="$d/best_solution.tsv"
    [ -s "$file" ] || continue

    genome=$(basename "$d")
    genome=${genome#sat_PICI_}
    genome=${genome#sat_cfPICI_}
    genome=${genome#sat_P4_}
    genome=${genome#sat_PLE_}

    awk -F '\t' -v g="$genome" '
    BEGIN{OFS="\t"}

    /^#/ {next}
    $1=="replicon" {next}
    NF<10 {next}

    {
        sys=$6
        hit=$2
        pos=$4+0
        whole=$9
        score=$10

        n[sys]++

        if(!(sys in min) || pos < min[sys]){
            min[sys]=pos
            first_hit[sys]=hit
        }

        if(!(sys in max) || pos > max[sys]){
            max[sys]=pos
            last_hit[sys]=hit
        }

        sys_score[sys]=score
        sys_whole[sys]=whole
    }

    END{
        for(sys in n){
            print g, sys, first_hit[sys], last_hit[sys], \
                  min[sys], max[sys], n[sys], \
                  sys_score[sys], sys_whole[sys]
        }
    }' "$file"

done | sort >> "$SYSTEMS"

# ---------- presence absence matrix ----------

echo -e "Genome\tPICI_present\tn_systems" > "$PA"

awk -F '\t' '
NR==1 {next}

{
    genome=$1
    seen[genome][$2]=1
}

END{
    for(g in seen){
        c=0
        for(s in seen[g]) c++
        print g "\t1\t" c
    }
}' "$SYSTEMS" | sort >> "$PA"

# ---------- add genomes with zero systems ----------

for d in sat_*; do
    [ -d "$d" ] || continue

    genome=$(basename "$d")
    genome=${genome#sat_PICI_}
    genome=${genome#sat_cfPICI_}
    genome=${genome#sat_P4_}
    genome=${genome#sat_PLE_}

    grep -q "^${genome}[[:space:]]" "$PA" || \
    echo -e "${genome}\t0\t0" >> "$PA"

done

sort -o "$PA" "$PA"

echo "Done."
echo "Created:"
echo "  $SYSTEMS"
echo "  $PA"
