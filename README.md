Repository for data and code for Hamm 2017 paper on the haplotype of *Danaus plexippis* and the evolution of chromosomes in the Danainae.

The original raw images used in this research have been deposited on FigShare:
Hamm, Christopher (2017): Monarch images. figshare.
https://doi.org/10.6084/m9.figshare.4633390.v1

The repo state at the time of the pre-print has been archived with Zenodo and can be cited:
[![DOI](https://zenodo.org/badge/81117728.svg)](https://zenodo.org/badge/latestdoi/81117728)

2017-02-06 - Initial commit

2017-02-08 - Added Markdown for the R code

Directory organization:

- `/Dan_chrom/` - input and output files for the `chromEvol` analyses. `CR` = Constant Rate, `LR` =  Linear Rate, `LRCR` = full model with all four parameters.

- `/Danainae_files/` -  figures made from `RMarkdown` needed for html.

- `/Data/`
  - `/acc/` - accession numbers for species by locus:  elongation factor 1 alpha (EF1a), mitochondrial DNA (mtDNA), and wingless (wg).
  - `Bayes` - all files output from `MrBayes`
  - `/fasta/` - Unaligned fasta files for each gene:
  - `/nexus/` - Aligned Nexus files for individual loci output from `Geneious v6.1.8` and concatenated file used for input to `MrBayes`.
  - `Danainae_karyotypes.csv` - file contain references to literature and counts.
  - `Danainae_karyotypes2.csv` - file contains only species names and haplotype counts.
  - `Danainae_phy.notes.txt` - Notes taken while working on the data.
  

- `/ms_images/` - Images generated for use in the manuscript.
