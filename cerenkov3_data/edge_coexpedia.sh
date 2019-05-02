#!/bin/bash

DIRECTORY="./resource/Coexpedia"

if [ -d "${DIRECTORY}/Hsa" ]; then
    # Control will enter here if that folder exists.
    echo "[Coexpedia] folder '${DIRECTORY}/Hsa' exists; no downloading"
else
    wget "http://www.coexpedia.org/dump/human_cx_net.zip"
    unzip -q human_cx_net.zip -d ${DIRECTORY}
    rm human_cx_net.zip
fi

python3 edge_coexpedia.py

echo "[Coexpedia] done!"