## `gene2ensembl.gz`

- [Download URL](ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz)
- Date Modified: 14/03/2019, 18:35:00
- Date Downloaded: 2019-03-15 14:34:32

===========================================================================  
gene2ensembl                                    recalculated daily
---------------------------------------------------------------------------  
           This file reports matches between NCBI and Ensembl annotation
           based on comparison of rna and protein features.

           Matches are collected as follows.
           For a protein to be identified as a match between RefSeq and
           Ensembl, there must be at least 80% overlap between the two.
           Furthermore, splice site matches must meet certain conditions:
           either 60% or more of the splice sites must match, or there may 
           be at most one splice site mismatch.

           For rna features, the best match between RefSeq and Ensembl is 
           selected based on splice site and overlap comparisons. For coding 
           transcripts, there is no minimum threshold for reporting other than 
           the protein comparison criteria above. For non-coding transcripts, 
           the splice site criteria are the same as for protein matching, but 
           the overlap threshold is reduced to 50%.

           Furthermore, both the rna and the protein features must meet these 
           minimum matching criteria to be considered a good match.  In 
           addition, only the best matches will be reported in this file.  
           Other matches that satisified the matching criteria but were
           not the best matches will not be reported in this file.

           A summary report of species that have been compared is contained
           in another FTP file, README_ensembl (see next item).

           More notes about this file:

           tab-delimited
           one line per match between RefSeq and Ensembl rna/protein
           Column header line is the first line in the file.


---------------------------------------------------------------------------

tax_id:
           the unique identifier provided by NCBI Taxonomy
           for the species or strain/isolate

GeneID:
           the unique identifier for a gene

Ensembl_gene_identifier:
           the matching Ensembl identifier for the gene

RNA nucleotide accession.version:
           the identifier for the matching RefSeq rna
           will be null (-) if only the protein matched

Ensembl_rna_identifier:
           the identifier for the matching Ensembl rna
           may include a version number
           will be null (-) if only the protein matched

protein accession.version:
           the identifier for the matching RefSeq protein
           will be null (-) if only the mRNA matched

Ensembl_protein_identifier:
           the identifier for the matching Ensembl protein
           may include a version number
           will be null (-) if only the mRNA matched

---------------------------------------------------------------------------

Transformation:

```bash
zcat gene2ensembl.gz | grep ENSG0 > gene2ensembl_9606.tsv
```