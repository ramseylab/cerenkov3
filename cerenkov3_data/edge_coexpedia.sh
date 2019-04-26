#!/bin/bash

wget "http://www.coexpedia.org/dump/human_cx_net.zip"
unzip -q human_cx_net.zip -d ./edge/gene-gene/co-expression/Coexpedia/
mv human_cx_net.zip ./edge/gene-gene/co-expression/Coexpedia/

python3 edge_coexpedia.py