**mash-dereplication.sh** : A script that uses mash distances to dereplicate a dataset of nucleotide fasta files. It needs pre-computed sketches in a folder (sketches/*.sh).In addition it needs the original *.fna files with the same basenames

**post-process-sat-finder.sh**: Combines the ouput of satellite finder (https://research.pasteur.fr/en/software/satellitefinder/) into a summary file, providing info on first and last gene identified, as well as a presence/absence matrix of PICIs/genome. Runs as:post-process-sat-finder.sh>>summary.tsv
