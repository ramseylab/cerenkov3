## Connect SNPs to Genes via TFBSs

First we will query the `wgEncodeRegTfbsClusteredV3` table of UCSC Genome Browser via MySQl to get all associated (SNP, TFBS) pairs. Then for each TFBS, we will query `MyGene.info` API to get associated (TFBS, gene) pairs. Finally we will merge these 2 set of pairs to output (SNP, gene) edges.

## Issues with Unofficial Symbols

The TFBS symbols in tabel `wgEncodeRegTfbsClusteredV3` of UCSC Genome Browser are mostly official, with a few exceptions: `FAM48A`, `GRp20`, `KAP1`, `RDBP`, `RPC155`, `SIN3AK20` and `SREBP1`. `MyGene.info` API cannot find the associated gene IDs for these 7 unofficial symbols, so we will revise them as indicated below:

- Rename to new symbols that do not exist in the `wgEncodeRegTfbsClusteredV3`
    - `FAM48A` => `SUPT20H`; see https://www.ncbi.nlm.nih.gov/gene/55578
    - `RDBP` => `NELFE`; see https://www.ncbi.nlm.nih.gov/gene/7936
    - `RPC155` => `POLR3A`; see https://www.ncbi.nlm.nih.gov/gene/11128
    - `SREBP1` => `SREBF1`; see https://www.ncbi.nlm.nih.gov/gene/6720
- Merge with the existing symbols in `wgEncodeRegTfbsClusteredV3`
    - `KAP1` => `TRIM28`; see https://www.ncbi.nlm.nih.gov/gene/10155
    - `GRp20` => `NR3C1`; see https://www.ncbi.nlm.nih.gov/gene/2908
        - Because `GRp20` == `GR (P-20)`, we see `GRp20` as a special `GR` and `GR` is synonymous to `NR3C1`
    - `SIN3AK20` => `SIN3A`; see https://www.ncbi.nlm.nih.gov/gene/25942
        - Because `SIN3AK20` == `SIN3A (K-20)`

