## 1. Mapping Ensembl Gene IDs to Entrez

We are going to build a SNP-gene network and currently the gene-gene edges from HumanNet and Coexpedia are Entrez-only; therefore if we want to connect SNP-gene edges whose gene node is a Ensembl ID to the known gene-gene networks, we need to transform Ensembl Gene IDs to Entrez.

The mapping is performed on 3 datasets:

- `gene2ensembl`
- `GRCh37-Biomart`
- `GRCh38-Biomart`

### 1.1 `gene2ensembl`

NCBI reports matches between NCBI and Ensembl annotation based on comparison of rna and protein features. This dataset can be download from [ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz](ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz), and relevant information is contained in [its README file](ftp://ftp.ncbi.nih.gov/gene/DATA/README).

Here I quote the section on `gene2ensembl` from that README file:

> This file reports matches between NCBI and Ensembl annotation based on comparison of rna and protein features.
> 
> Matches are collected as follows. For a protein to be identified as a match between RefSeq and Ensembl, there must be at least 80% overlap between the two. Furthermore, splice site matches must meet certain conditions: either 60% or more of the splice sites must match, or there may be at most one splice site mismatch.
> 
> For rna features, the best match between RefSeq and Ensembl is selected based on splice site and overlap comparisons. For coding transcripts, there is no minimum threshold for reporting other than the protein comparison criteria above. For non-coding transcripts, the splice site criteria are the same as for protein matching, but the overlap threshold is reduced to 50%.
> 
> Furthermore, both the rna and the protein features must meet these minimum matching criteria to be considered a good match. In addition, only the best matches will be reported in this file. Other matches that satisified the matching criteria but were not the best matches will not be reported in this file.
> 
> A summary report of species that have been compared is contained in another FTP file, README_ensembl (see next item).
> 
> More notes about this file:
> 
> tab-delimited： one line per match between RefSeq and Ensembl rna/protein Column header line is the first line in the file.

The columns of `gene2ensembl` are explained below:

| colname                          | remark                                                                            | 
|----------------------------------|-----------------------------------------------------------------------------------| 
| tax_id                           | the unique identifier provided by NCBI Taxonomy for the species or strain/isolate | 
| GeneID                           | the unique identifier for a gene                                                  | 
| Ensembl_gene_identifier          | the matching Ensembl identifier for the gene                                      | 
| RNA nucleotide accession.version | the identifier for the matching RefSeq rna                                        | 
|                                  | will be null (-) if only the protein matched                                      | 
| Ensembl_rna_identifier           | the identifier for the matching Ensembl rna                                       | 
|                                  | may include a version number                                                      | 
|                                  | will be null (-) if only the protein matched                                      | 
| protein accession.version        | the identifier for the matching RefSeq protein                                    | 
|                                  | will be null (-) if only the mRNA matched                                         | 
| Ensembl_protein_identifier       | the identifier for the matching Ensembl protein                                   | 
|                                  | may include a version number                                                      | 
|                                  | will be null (-) if only the mRNA matched                                         | 

We further extract the mapping of human genes only:

```bash
zcat gene2ensembl.gz | grep ENSG0 > gene2ensembl_9606.tsv
```

and then the columns `GeneID` (Entrez ID) and `Ensembl_gene_identifier` are what we need

### 1.2 `GRCh37-Biomart`

We fetch the full table of genes with 3 columns, 

- `Gene stable ID` (Ensembl ID), 
- `Transcript stable ID`, and 
- `HGNC ID`, 

from [GRCH37 Ensembl](https://grch37.ensembl.org/index.html) using Biomart. Then we adopt the following ways for mapping:

| method                         | remark                                                   | 
|--------------------------------|----------------------------------------------------------| 
| joining `Gencode` table        | mapping each `Transcript stable ID` to an Entrez Gene ID | 
| joining `HGNC` table           | mapping each `HGNC ID` to an Entrez Gene ID              | 
| querying `MyGene.info` API     | mapping each `Gene stable ID` to an Entrez Gene ID       | 

Therefore for any Ensembl Gene ID $\epsilon$, we can find 3 possible Entrez Gene IDs:

- $\xi_{gc}$ from Gencode
- $\xi_{hgnc}$ from HGNC (some of the HGNC entries have no Ensembl IDs)
- $\xi_{mg}$ from MyGene.info

$\xi_{gc}$, $\xi_{mg}$ and $\xi_{hgnc}$ can be null or inconsistent with each other.

To determine a final Entrez Gene ID $\xi$, we check $\xi_{hgnc}$, ~~$\xi_{mg}$,~~ $\xi_{gc}$ and $\xi_{bm}$ in order and pick the first non-null value. 

```python
def first_non_nan(row):
    ordered_keys = ["Entrez_Gene_ID_hgnc_2", "Entrez_Gene_ID_hgnc_1", "Entrez_Gene_ID_mg", "Entrez_Gene_ID_gc"]

    # next(iterator, default): Retrieve the next item from the iterator by calling its __next__() method. If default is given, it is returned if the iterator is exhausted, otherwise StopIteration is raised.
    return next((row[key] for key in ordered_keys if not pd.isnull(row[key])), None)
```

Discontinued Entrez Gene IDs will be eliminated.

### 1.3 `GRCh38-Biomart`

Similar to `GRCh37-Biomart` but here we fetch the initial gene table from GRCh38 Ensembl (at the time of writing, we used [Ensembl 96: Apr 2019](http://apr2019.archive.ensembl.org/)).

Since the mapping does not involve gene locations but only gene IDs, we think it's acceptable to use GRCh38 data here although our SNP data are from GRCh37.

## 2. `MyGene.info` API

The version of `MyGene.info` API we used is 3.0.0:

```bash
ramseylab:~$ pip3 show mygene
Name: mygene
Version: 3.0.0
Summary: Python Client for MyGene.Info services.
Home-page: https://github.com/suLab/mygene.py
Author: Chunlei Wu, Cyrus Afrasiabi, Sebastien Lelong
Author-email: cwu@scripps.edu
License: BSD
Location: /usr/local/lib/python3.5/dist-packages
Requires: requests
Required-by:
```

## 3. Data

- `gene2ensembl`
    - [gene2ensembl.gz](ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz)
- `GRCh37-Biomart`
    - [Initial table from Biomart](http://grch37.ensembl.org/biomart/martservice/results?query=%3C?xml%20version=%221.0%22%20encoding=%22UTF-8%22?%3E%3C!DOCTYPE%20Query%3E%20%3CQuery%20virtualSchemaName=%22default%22%20formatter=%22TSV%22%20header=%221%22%20uniqueRows=%220%22%20count=%22%22%20datasetConfigVersion=%220.6%22%3E%3CDataset%20name=%22hsapiens_gene_ensembl%22%20interface=%22default%22%3E%3CAttribute%20name=%22ensembl_gene_id%22/%3E%3CAttribute%20name=%22ensembl_transcript_id%22/%3E%3CAttribute%20name=%22hgnc_id%22/%3E%3C/Dataset%3E%3C/Query%3E)
        - Filters: None
        - Attributes: `Gene stable ID`, `Transcript stable ID`, `HGNC ID`
    - `Gencode`
        - Website: [**Gencode** Release 29 (GRCh37)](https://www.gencodegenes.org/human/release_29lift37.html)
            - => Metadata files 
                - => Entrez gene ids
        - [Direct Download Link](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh37_mapping/gencode.v29lift37.metadata.EntrezGene.gz)
    - `HGNC`
        - Downloaded from [**HGNC** Custom Downloads](https://www.genenames.org/download/custom/)
            - [Direct URL](https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_pub_chrom_map&col=gd_pub_refseq_ids&col=gd_pub_eg_id&col=gd_pub_ensembl_id&col=family.id&col=family.name&col=md_eg_id&col=md_ensembl_id&col=md_refseq_id&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit)
        - Query criteria:
            - Columns:
                * HGNC ID	
                * Approved symbol	
                * Approved name	
                * Status	
                * Chromosome	
                * RefSeq IDs	
                * NCBI Gene ID	
                * Ensembl gene ID	
                * Gene group ID	
                * Gene group name	
                * NCBI Gene ID(supplied by NCBI)
                    - *NB*. this column could be supplementary to "NCBI Gene ID" when some entries in that column are NULL
                * Ensembl ID(supplied by Ensembl)	
                * RefSeq(supplied by NCBI)
            - Status: Approved only
            - Chromosomes:
                - All, including
                    - 1-22, X, Y
                    - reserved loci
                    - mitochondrial
                    - pseudoautosomal
            - Order by HGNC ID
- `GRCh38-Biomart`
    - [Initial table from Biomart](http://apr2019.archive.ensembl.org/biomart/martservice/results?query=%3C?xml%20version=%221.0%22%20encoding=%22UTF-8%22?%3E%3C!DOCTYPE%20Query%3E%3CQuery%20virtualSchemaName=%22default%22%20formatter=%22TSV%22%20header=%221%22%20uniqueRows=%220%22%20count=%22%22%20datasetConfigVersion=%220.6%22%3E%3CDataset%20name=%22hsapiens_gene_ensembl%22%20interface=%22default%22%3E%3CAttribute%20name=%22ensembl_gene_id%22%20/%3E%3CAttribute%20name=%22ensembl_transcript_id%22%20/%3E%3CAttribute%20name=%22hgnc_id%22%20/%3E%3C/Dataset%3E%3C/Query%3E)
        - Filters: None
        - Attributes: `Gene stable ID`, `Transcript stable ID`, `HGNC ID`
    - `Gencode`
        - Website: [**Gencode** Release 29 (GRCh38.p12)](https://www.gencodegenes.org/human/release_29.html)
            - => Metadata files 
                - => Entrez gene ids
        - [Direct Download Link](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.metadata.EntrezGene.gz)
    - `HGNC`
        - Same as `GRCh37-Biomart`
        - [Direct URL](https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_pub_chrom_map&col=gd_pub_refseq_ids&col=gd_pub_eg_id&col=gd_pub_ensembl_id&col=family.id&col=family.name&col=md_eg_id&col=md_ensembl_id&col=md_refseq_id&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit)
