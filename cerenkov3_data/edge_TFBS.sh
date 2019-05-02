#!/bin/bash

version="3.5.171"
biogrid_dir="./edge/snp-gene/TFBS/BioGRID/"

wget "https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-${version}/BIOGRID-ALL-${version}.tab2.zip"
unzip -q BIOGRID-ALL-${version}.tab2.zip -d .  # => BIOGRID-ALL-${version}.tab2.txt

# save the header
head -n1 BIOGRID-ALL-${version}.tab2.txt | cut -f2,3,8,9,16,17 > _temp.tsv
# append human records (Taxonomy ID: 9606) 
cut -f2,3,8,9,16,17 BIOGRID-ALL-${version}.tab2.txt | grep -P "9606\t9606" >> _temp.tsv
# cut off the Taxonmy ID columns
cut -f1,2,3,4 _temp.tsv > _temp2.tsv
# remove duplicate lines; see https://unix.stackexchange.com/a/30178
awk '!seen[$0]++' _temp2.tsv > BIOGRID-Human-${version}.tab2.tsv

# delete the temporary file, the original txt file and the original zip file
rm _temp.tsv
rm _temp2.tsv
rm BIOGRID-ALL-${version}.tab2.txt
rm BIOGRID-ALL-${version}.tab2.zip

mv BIOGRID-Human-${version}.tab2.tsv ${biogrid_dir}