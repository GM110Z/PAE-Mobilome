#!/bin/bash

out="results.tsv"

# header
echo -e "accession\torganism" > "$out"

while read id; do
    # skip empty lines
    [ -z "$id" ] && continue

    printf "%s\t" "$id" >> "$out"

    res=$(esummary -db nuccore -id "$id" 2>/dev/null | \
          xtract -pattern DocumentSummary -element Organism Title)

    if [ -z "$res" ]; then
        echo "NA" >> "$out"
    else
        echo "$res" >> "$out"
    fi

    sleep 0.3
done < $1
