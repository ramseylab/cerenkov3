#!/bin/bash

VERSION="3.5.171"

DIRECTORY="./resource/BioGRID"
FILENAME="BIOGRID-Human.tab2.tsv"

if [ -f "${DIRECTORY}/${FILENAME}" ]; then
    # Control will enter here if that file exists.
    echo "[BioGRID] file '${DIRECTORY}/${FILENAME}' exists; skip downloading"
else
    wget "https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-${VERSION}/BIOGRID-ALL-${VERSION}.tab2.zip"
    unzip -q BIOGRID-ALL-${VERSION}.tab2.zip -d .  # => BIOGRID-ALL-${VERSION}.tab2.txt

    # save the header
    head -n1 BIOGRID-ALL-${VERSION}.tab2.txt | cut -f2,3,8,9,16,17 > _temp.tsv
    # save human genomic records (Taxonomy ID: 9606) 
    cut -f2,3,8,9,16,17 BIOGRID-ALL-${VERSION}.tab2.txt | grep -P "9606\t9606" >> _temp.tsv
    # cut off the Taxonmy ID columns
    cut -f1,2,3,4 _temp.tsv > _temp2.tsv
    # remove duplicate lines; see https://unix.stackexchange.com/a/30178
    awk '!seen[$0]++' _temp2.tsv > ${FILENAME}

    # delete the temporary file, the original txt file and the original zip file
    rm _temp.tsv
    rm _temp2.tsv
    rm BIOGRID-ALL-${VERSION}.tab2.txt
    rm BIOGRID-ALL-${VERSION}.tab2.zip

    # save the edgelist to its destination
    mv ${FILENAME} ${DIRECTORY}
fi

python3 edge_biogrid.py

echo "[BioGRID] done!"