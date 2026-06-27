#  Prophage Helper Analysis

`prophage_helper_analysis.py` combines **geNomad**, **taxMyPhage**, and SatelliteFinder/PICIasso architecture results to test whether particular phage taxa are enriched or depleted among genomes carrying each predicted mobile-element architecture.

The script reports both statistical enrichment and the absolute prevalence of each taxon.

## Analysis overview

For every architecture–taxon combination, the script compares:

- **Focal genomes:** genomes carrying at least one element assigned to the architecture being tested.
- **Control genomes:** all other genomes represented in the architecture input table.

It calculates:

- number and percentage of focal genomes carrying the taxon;
- number and percentage of control genomes carrying the taxon;
- absolute prevalence difference in percentage points;
- Fisher's exact test odds ratio and P value;
- an odds ratio with a 95% confidence interval;
- Benjamini–Hochberg false-discovery rate;
- a prevalence-based interpretation category.

The script performs three sets of association tests:

1. **geNomad Inoviridae-like prophages**
2. **taxMyPhage family assignments**
3. **taxMyPhage genus assignments**

## Interpretation categories

By default, statistically significant positive associations are divided into:

| Category | Definition |
|---|---|
| `architecture-wide_candidate` | FDR < 0.05, odds ratio > 1, and taxon prevalence in focal genomes ≥50% |
| `subgroup_association` | FDR < 0.05, odds ratio > 1, and prevalence from 10% to <50% |
| `rare_association` | FDR < 0.05, odds ratio > 1, and prevalence <10% |
| `depleted` | FDR < 0.05 and odds ratio <1 |
| `no_significant_association` | FDR ≥0.05 |

The 50% and 10% thresholds are descriptive labels rather than additional statistical filters. They can be changed from the command line.

## Requirements
- `numpy`
- `pandas`
- `matplotlib`
- `scipy`


## Input files

### 1. geNomad prophage table

Provided with:

```text
--prophages all_prophages.tsv
```

Required columns:

| Column | Description |
|---|---|
| `genome` | Host genome name or path |
| `seq_name` | Prophage sequence identifier used for matching to taxMyPhage |
| `taxonomy` | geNomad taxonomic annotation |


### 2. taxMyPhage taxonomy summary

Provided with:

```text
--taxmyphage Summary_taxonomy.tsv
```

Required columns:

| Column | Description |
|---|---|
| `Genome` | Prophage identifier to match against geNomad `seq_name` |
| `Family` | Assigned phage family |
| `Genus` | Assigned phage genus |

Unclassified or placeholder assignments are excluded. These include values such as `unknown`, `new_genus`, `not defined yet`, `unclassified`, and empty values.

Before joining, prophage identifiers are normalised by:

- removing surrounding whitespace;
- replacing `|` with `_`;
- removing internal whitespace.

The script writes the complete join table so unmatched identifiers can be inspected.

### 3. Architecture/completeness table

Provided with:

```text
--completeness all_completeness.tsv
```

Expected columns:

| Column | Description |
|---|---|
| `Genome` | SatelliteFinder element identifier |
| `architecture` | Architecture assigned to that element |


Association tests are performed at the **genome level**, not the element level. A genome is positive for an architecture when it contains at least one element assigned to that architecture.

A genome carrying multiple architecture classes is included in each relevant focal group. Consequently, architecture groups are not necessarily mutually exclusive.

### 4. Optional coordinate-overlap summary

Provided with:

```text
--overlap_summary element_inovirus_overlap_summary.tsv
```

This file is optional. When present, it adds architecture-level context about possible overlap between predicted elements and Inoviridae-like prophages.

Expected columns include:

| Column | Description |
|---|---|
| `architecture` | Element architecture |
| `element_id` | Element identifier |
| `has_any_overlap` | Whether an overlap was detected |
| `suspected_misclassification` | Whether the overlap suggests possible misclassification |

The overlap data do not alter the association tests. They are added as interpretive context to the output table.

## Basic usage

```bash
python3 combined_prophage_helper_analysis.py \
  --prophages all_prophages.tsv \
  --taxmyphage Summary_taxonomy.tsv \
  --completeness all_completeness.tsv
```



