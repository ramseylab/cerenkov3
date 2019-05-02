## HumanNet

HumanNet (v2) is a human gene networks for disease research.

- [Project site](https://www.inetbio.org/humannet/)
    - [Download page](https://www.inetbio.org/humannet/download.php)
        - [Fully extended functional gene network - HumanNet-XN.tsv](https://www.inetbio.org/humannet/networks/HumanNet-XN.tsv)
            - Column _LLS_ means _log-likelihood score_

Previously we used its [v1 version](http://www.functionalnet.org/humannet/about.html):

- [HumanNet.v1.join.txt](http://www.functionalnet.org/humannet/HumanNet.v1.join.txt)
    - File format explained: [HumanNet.v1.evidence_code.txt](http://www.functionalnet.org/humannet/HumanNet.v1.evidence_code.txt)

## Data Processing

Our script aims to:

- Remove the LLS column from HumanNet tsv; rename columns to _Gene_A_ and _Gene_B_
    - The resulted dataset contains only two columns, with each row representing 2 co-expressed genes
