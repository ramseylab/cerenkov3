#!/bin/bash

wget "http://www.coexpedia.org/dump/human_cx_net.zip"
unzip -q human_cx_net.zip -d ./resource/Coexpedia/
rm human_cx_net.zip

python3 edge_coexpedia.py