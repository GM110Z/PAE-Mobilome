Run mash sketch as normal
Run mash dist as normal
Then run LC_ALL=C awk -F'\t' '$1!=$2 && $3<=0.05' mash_all_out/all_vs_all.tsv > mash_all_out/all_vs_all_le_005.tsv
*This will remove all the self-comparisons and it will further pick out unique genomes that are at a mash dist >0.1 from antything else*
Run : find . -maxdepth 1 -name "*.fna" | sort > all_genomes.txt 

Then run the 
python mash-clustering.py <list-of-all-genomes> <mash-matrix-cleaned-with-awk> <your-threshold-for-clustering> <outfile-name>
