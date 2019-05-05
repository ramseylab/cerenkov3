## BIOGRID Data

We will use BIOGRID version `3.5.171` here. The file to download is [BIOGRID-ALL-3.5.171.tab2.zip](https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-3.5.171/BIOGRID-ALL-3.5.171.tab2.zip). Note that it's recommended to use the `tab2` file instead of the `tab` format. The header of this `tab2` file can be found at [BioGRID TAB 2.0 format](https://wiki.thebiogrid.org/doku.php/biogrid_tab_version_2.0)

BIOGRID reports interactions between TFs and essentially between genes. We will simply extract gene-gene edges from this table.

