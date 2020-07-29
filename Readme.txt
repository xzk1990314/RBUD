1. MDGM database introduction

Mdgm database includes species data set and function data set. The functional data set contains gene annotation information from different databases. The folder is called annotation database, including information from ARDB, CARD, CAZY, COG, eggNOG and KEGG. In addition, we provide a file called Linkdb_gene, the users can link it with gene ID, GI ID and protein to find the comments. the method of establishing the database and the running code was gave as follows.


2. The code of MDGM database construction

#Preparation

You should download all microorganism data from NCBI, the data should include the suffixes of fna,gbk,ffn and gff. In this process,all microorganism contain archaea,bacteria,fungi,virus.

# 2bwt-builder-species.sh

it is for creating binary of reference genome for species dataset.

# species_annotation.pl 

it is for extracting the species annotation.

# 2bwt-builder-CDS.sh

it is for creating binary of reference genome for functional dataset.

# gene annotation.pl
 
it is for extracting the funcitonal annotation of microorganism.

# According to gene ID, GI ID and protein ID to find the annotation in databases of ARDB,CARD,CAZy, COG,KEGG. If you can't find the annotation,it is the best way to align DNA sequence to these databases for functional annotation using Blast software.

#Finally, based on the above biological information to construct MDGM database.

3. The software of RBUD method

# Preparation

Before you use RBUD method,you should build your development environment in Liux. In addition, SOAPaligner/soap2 should be installed.

## soap.sh

Align query to genome and store output. 

## length.py

Based on fasta format file, this script could caculate the length of reference sequence that aligned by sequecing data.


# Runing RBUD method

## RBUD.py

A script for RBUD method. length_file is the length of every items of reference sequence that aligned by sequencing data.input_file is the comparison results by SOAPaligner/soap2. output_file is the abundance of species or genes.In addition,you could run "python RBUD.py -h" for searching help.  

Comments: python RBUD.py -l length_file -i input_file -o output_file

## According the results of RBUD method, you could use other statistical tools/softwares to compare the difference in your metagenomic research.


