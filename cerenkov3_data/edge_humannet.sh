#!/bin/bash

DIRECTORY="./resource/HumanNet"
FILENAME="HumanNet-XN.tsv"

if [ -f "${DIRECTORY}/${FILENAME}" ]; then
    # Control will enter here if that file exists.
    echo "[HumanNet] file '${DIRECTORY}/${FILENAME}' exists; skip downloading"
else
    wget -O ${FILENAME} "https://www.inetbio.org/humannet/networks/HumanNet-XN.tsv"
    mv ${FILENAME} ${DIRECTORY}
fi

python3 edge_humannet.py

echo "[HumanNet] done!"