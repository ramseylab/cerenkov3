#!/bin/bash

DIRECTORY="./resource/Util_DEI"
FILENAME="discontinued_entrez_id.tsv"

if [ -f "${DIRECTORY}/${FILENAME}" ]; then
    # Control will enter here if that file exists.
    echo "[Util_DEI] file '${DIRECTORY}/${FILENAME}' exists; skip downloading"
else
    curl "ftp://ftp.ncbi.nih.gov/gene/DATA/gene_history.gz" | gunzip -c | cut -f 3,4 > ${FILENAME}

    mv ${FILENAME} ${DIRECTORY}
fi

echo "[Util_DEI] done!"