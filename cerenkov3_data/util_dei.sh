curl "ftp://ftp.ncbi.nih.gov/gene/DATA/gene_history.gz" | gunzip -c | cut -f 3,4 > discontinued_entrez_id.tsv

mv discontinued_entrez_id.tsv ./util/discontinued_entrez_id