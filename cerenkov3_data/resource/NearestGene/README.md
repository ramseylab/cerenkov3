## Data Processing

For each SNP we have, we will draw an edge to its nearest gene. The distance will be measured between the SNP's start position and the candidate gene's TSS start, same strandedness not required. 

Note that in our BED files, for any record, `chromStart` (2nd column) is always lower than `chromEnd` (3rd column). Therefore,

- if a gene is on the forward strand, its TSS start will be its `chromStart`;
- otherwise `chromEnd`.

Instead of making new BED files, we will reuse 

- `cerenkov3/cerenkov3_data/vertex/SNP/OSU18/osu18_SNP.bed`, 
- `cerenkov3/cerenkov3_data/vertex/gene/Ensembl/p2_Ensembl_transcript.bed` and
- `cerenkov3/cerenkov3_data/vertex/gene/Entrez/p3_Entrez_transcript.bed`

calculate their start positions in Python and call the Python Bedtools API to run `bedtools closest`.

For accurate distance calculation, we will on purpose

- set `chromStart` == `chromEnd` == SNP's start position, and
- set `chromStart` == `chromEnd` == gens's TSS start.
