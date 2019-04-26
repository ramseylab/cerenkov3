#!/bin/bash

wget "https://www.inetbio.org/humannet/networks/HumanNet-XN.tsv"

mv HumanNet-XN.tsv ./edge/gene-gene/co-expression/HumanNet

python3 edge_humannet.py