**mash-dereplication.sh** : A script that uses mash distances to dereplicate a dataset of nucleotide fasta files. It needs pre-computed sketches in a folder (sketches/*.sh).In addition it needs the original *.fna files with the same basenames

**pici-sorter.sh**: Combines the ouput of satellite finder (https://research.pasteur.fr/en/software/satellitefinder/) into a summary file, providing info on first and last gene identified, as well as a presence/absence matrix of PICIs/genome. 

**sniPICI.py** :Uses summary from pici-sorter.sh to extract PICI regions, with 3 genes UP/DW. Excludes hits on different contigs

Runs as : python sniPICI.py      --summary pici_systems_summary.tsv  --gbff_dir /path/to/gbff_files  --outdir pici_gbff_regions --flank 5

**MakeClinkerGroups.sh**. Uses dereplicated picis output from *mash-derep.sh script* and the mash matrix to group gbff files in groups that have a set similarity so they look decent on Clinker
./makeClinkGroups.sh summaryclust.tsv pici_mash_all_vs_all.tsv outfolder mash_distance_of_choice


**fetch-strain-name.sh** Uses a list of nuccore ID to fetch strain name.Run as *fetch-strain-name.sh id_file.txt*

----------------------------------------------------------------------------------------------------------------------------------------------------------
**PICIasso**

A  pipeline for the  extraction, classification, and
visualisation of Pathogenicity Island-like Chromosomal Islands (PICIs) predicted
by SatelliteFinder in bacterial genomes.
---

## Overview

PICIasso takes the  output of sniPICI.sh and produces:

- A clean per-system summary table with genomic coordinates and gene product
- PFAM domain annotation of every CDS within each element
- Module completeness scoring against a canonical PICI definition
- Architecture classification (e.g. `complete_canonical`, `backbone_only`)
- Gene schematic plots coloured by module identity
- Per-architecture overview plots for cross-element comparison


## Usage

#### Full run (PFAM + classification + plots)
```bash
python PICIasso.py pici_gbff_regions \
    --pfam /path/to/Pfam-A.hmm
```

#### Without PFAM (product names only, faster, less accurate)
```bash
python pici_annotate_and_classify.py pici_gbff_regions --no-pfam
```

#### Regraph only (plots only, reuses existing TSV outputs)
Use this after a full run if you want to regenerate plots without re-running
PFAM or reclassifying. Requires `master_completeness_summary.tsv` and
`*_annotations.tsv` files from a previous run.

```bash
python pici_annotate_and_classify.py pici_gbff_regions --regraph
```

If your output files are in a different directory than your GBFFs:
```bash
python pici_annotate_and_classify.py pici_gbff_regions \
    --outdir /path/to/previous/output \
    --regraph
```

| Argument | Default | Description |
|----------|---------|-------------|
| `--pfam` | `Pfam-A.hmm` | Path to pressed Pfam-A.hmm database |
| `--no-pfam` | — | Skip PFAM, use product name keywords only |
| `--rerun-pfam` | — | Force re-run PFAM even if `.pfam` files exist |
| `--outdir` | same as `gbff_dir` | Output directory for master summary and plots |
| `--regraph` | — | Skip all computation, regenerate plots from existing TSVs |

---

## Architecture Classification

Elements are classified in the following priority order:

| Architecture | Required modules | Notes |
|---|---|---|
| `complete_canonical` | integrase + regulator + replication + terminase + capsid | Full PICI |
| `packaging_positive_partial` | integrase + regulator + terminase or capsid | Incomplete packaging |
| `backbone_replication_no_terminase` | integrase + regulator + replication (or putative) | No packaging genes |
| `mobilisation_associated` | integrase + mobilisation | IME-like |
| `backbone_only` | integrase + regulator only | Minimal backbone |
| `minimal_or_relic` | integrase only | Likely degraded |
| `partial_unclear` | anything else | Ambiguous |

The **canonical completeness score** (0–100%) reflects how many of the five
canonical modules (integrase, regulator, replication, terminase, capsid) are
present, regardless of architecture class.
