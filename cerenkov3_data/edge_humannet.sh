#!/bin/bash

DIRECTORY="./resource/HumanNet"

if [ -f "${DIRECTORY}/HumanNet-XN.tsv" ]; then
    # Control will enter here if that file exists.
    echo "[HumanNet] file '${DIRECTORY}/HumanNet-XN.tsv' exists; no downloading"
else
    wget "https://www.inetbio.org/humannet/networks/HumanNet-XN.tsv"
    mv HumanNet-XN.tsv ${DIRECTORY}
fi

python3 edge_humannet.py

echo "[HumanNet] done!"