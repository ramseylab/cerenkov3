# gene2ensembl

wget 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz'
zcat gene2ensembl.gz | grep ENSG0 | cut -f 1,2,3 > gene2ensembl_9606.tsv

rm gene2ensembl.gz
mv gene2ensembl_9606.tsv ./vertex/gene/Entrez/

# BioMart

wget -O GRCh37_p13_mart_export.txt 'http://grch37.ensembl.org/biomart/martservice/results?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query> <Query  virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="0" count="" datasetConfigVersion="0.6"><Dataset name="hsapiens_gene_ensembl" interface="default"><Attribute name="ensembl_gene_id"/><Attribute name="ensembl_transcript_id"/><Attribute name="hgnc_id"/></Dataset></Query>'

wget -O GRCh38_p12_mart_export.txt 'http://apr2019.archive.ensembl.org/biomart/martservice/results?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="0" count="" datasetConfigVersion="0.6"><Dataset name="hsapiens_gene_ensembl" interface="default"><Attribute name="ensembl_gene_id" /><Attribute name="ensembl_transcript_id" /><Attribute name="hgnc_id" /></Dataset></Query>'

mv GRCh37_p13_mart_export.txt ./vertex/gene/Entrez/
mv GRCh38_p12_mart_export.txt ./vertex/gene/Entrez/

# Gencode

wget 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh37_mapping/gencode.v29lift37.metadata.EntrezGene.gz'
gunzip gencode.v29lift37.metadata.EntrezGene.gz

mv gencode.v29lift37.metadata.EntrezGene ./vertex/gene/Entrez/

wget 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.metadata.EntrezGene.gz'
gunzip gencode.v29.metadata.EntrezGene.gz

mv gencode.v29.metadata.EntrezGene ./vertex/gene/Entrez/

# HGNC

wget -O HGNC_custom.txt 'https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_pub_chrom_map&col=gd_pub_refseq_ids&col=gd_pub_eg_id&col=gd_pub_ensembl_id&col=family.id&col=family.name&col=md_eg_id&col=md_ensembl_id&col=md_refseq_id&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit'

mv HGNC_custom.txt ./vertex/gene/Entrez/

# Mapping

python3 vertex_entrez.py