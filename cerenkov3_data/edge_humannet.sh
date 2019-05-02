#!/bin/bash

wget "https://www.inetbio.org/humannet/networks/HumanNet-XN.tsv"

mv HumanNet-XN.tsv ./resource/HumanNet

python3 edge_humannet.py