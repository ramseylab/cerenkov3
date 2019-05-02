#!/bin/bash

DIRECTORY="./resource/Ensembl"
FILENAME="mart_export.txt"

if [ -f "${DIRECTORY}/${FILENAME}" ]; then
    # Control will enter here if that file exists.
    echo "[Ensembl] file '${DIRECTORY}/${FILENAME}' exists; no downloading"
else
    wget -O ${FILENAME} 'http://grch37.ensembl.org/biomart/martservice/results?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="0" count="" datasetConfigVersion="0.6"><Dataset name="hsapiens_gene_ensembl" interface="default"><Filter name="chromosome_name" value="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y"/><Attribute name="chromosome_name"/><Attribute name="transcript_start"/><Attribute name="transcript_end"/><Attribute name="transcription_start_site"/><Attribute name="strand"/><Attribute name="ensembl_gene_id"/><Attribute name="external_gene_name"/><Attribute name="start_position"/><Attribute name="end_position"/><Attribute name="ensembl_transcript_id"/><Attribute name="external_transcript_name"/></Dataset></Query>'

    mv ${FILENAME} ${DIRECTORY}
fi

Rscript vertex_ensembl.R

echo "[Ensembl] done!"