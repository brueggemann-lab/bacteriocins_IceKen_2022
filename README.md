# bacteriocins_IceKen_2022

Investigating the distribution of bacteriocin biosynthetic gene clusters in *Streptococcus pneumoniae* genomes recovered from Iceland and Kenya. 

~Citation here~ 

Correspondance: Prof Angela Brueggemann (angela.brueggemann@ndph.ox.ac.uk) or Dr Madeleine Butler (m.butler17@imperial.ac.uk). 

## Bacteriocin processing and analysis 

### Data files
* data/bigsdb_annotated_export.csv - a full export of genomes with bacteriocin gene annotations and associated metadata used by the processing and analysis code will be added at the time of publication when genomic data are made available. 
* data/continuity_cat_outputs/cluster_cont_Kenya_2500_2500.csv and cluster_cont_VICE_2500_2500.csv - outputs from contiguity_cat.py (see below) describing which bacteriocin genes are found as contiguous clusters within each genome from each dataset. 

### Code 
Code is provided as a Jupyter notebook, where all processing and analysis functions are called and outputs are saved, and also as two text files where the functions are defined. 

* code/bacteriocins/processing_and_visualisation.ipynb - Jupyter notebook calling all processing and analysis functions. 
* code/bacteriocins/processing.py and analysis.py - text files defining functions for processing and analysing bacteriocin data, and for generating visualisations and summaries. 

## Bacteriocin cluster contiguity assessment 

### Generating annotated sequence files 
* code/bigs_genbankerator/bigs_genbankerator.py - a command line tool for fetching annotated sequence data from the private BIGSdb database in which whole genome sequence data, annotations and metadata are stored. 

### Checking bacteriocin cluster contiguity 
* code/contiguity_cat/contiguity_cat.py - a command line tool that takes BIGS_genbankerator.py output files and assesses which bacteriocin gene clusters are contiguous according to user-defined thresholds. 
* code/contiguity_cat/dummy_input.gbk - example output ganbank file from bigs_genbankerator.py for input in contiguity_cat.py. 
