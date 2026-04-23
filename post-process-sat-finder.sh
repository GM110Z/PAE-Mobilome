#!/usr/bin/env bash

echo -e "Genome\tPICI_present\tfirst_hit\tlast_hit\tfirst_pos\tlast_pos\tn_markers"

for d in sat_*; do
    [ -d "$d" ] || continue

    genome="${d#sat_}"
    genome="${genome#*_}"

    file="$d/best_solution.tsv"

    if [ ! -s "$file" ]; then
        echo -e "${genome}\t0\tNA\tNA\tNA\tNA\t0"
        continue
    fi

    awk -F '\t' -v g="$genome" '
    BEGIN{
        min=""
        max=""
        n=0
    }

    /^#/ {next}
    $1=="replicon" {next}
    NF<5 {next}

    {
        pos=$4+0
        id=$2

        n++

        if(min=="" || pos<min){
            min=pos
            min_id=id
        }

        if(max=="" || pos>max){
            max=pos
            max_id=id
        }
    }

    END{
        if(n>0)
            print g "\t1\t" min_id "\t" max_id "\t" min "\t" max "\t" n
        else
            print g "\t0\tNA\tNA\tNA\tNA\t0"
    }' "$file"

done | sort
