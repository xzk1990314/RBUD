1. MDGM database introduction

MDGM database has two datasets, including Species dataset and Functional dataset. Each dataset consists of two floders of Sequence and Annotation. In functional dataset, we added gene annotation information of different database in the folder named Annotation Database, containg ARDB, CARD,CAZy,COG,eggNOG,KEGG. Moreover, we provided a folder named linkdb_genes to find the annotation of the above database using gene id ,gi id and protein id.However,due to the size of this database is too large,we provided the method of building it.


2. The code of MDGM database construction

#Preparation

it should download the data of all microorganism with the suffixes of fna,gbk,ffn and gff.In this process,all microorganism contain archaea,bacteria,fungi,virus and so on.

# 2bwt-builder-species.sh

it was used to create binary of reference genome for species dataset.

# species_annotation.pl 

it was used to extract the species annotation.

# 2bwt-builder-CDS.sh

it was used to create binary of reference genome for functional dataset.

# gene annotation.pl
 
it was used to extract the funcitonal annotation of microorganism.

# According to gene id, gi id and protein id to find the annotation in the databases of ARDB,CARD,CAZy, COG,KEGG.If you can't find the annotation,it is the best way to align DNA sequence to these databases for functional annotation using Blast software.

#Finally, based on the above biological information to construct MDGM database.

3. The software of RBUD method

# Preparation

Before you use RBUD method,you should build your development environment in Liux. In addition,SOAPaligner/soap2 should be installed.

## soap.sh

Align query to genome and store output. 

## length.py

Based on fasta format file, this script could caculate the length of reference sequence that aligned by sequecing data.


# Runing the software of RBUD method

## RBUD.py

A script for RBUD method. length_file is the length of every items of reference sequence that aligned by sequencing data.input_file is the comparison results by SOAPaligner/soap2. output_file is the abundance of species or genes.In addition,you could run "python RBUD.py -h" for searching help.  

Comments: python RBUD.py -l length_file -i input_file -o output_file

## According the results of the software of RBUD method, you could use other statistical tools to compare the difference of different groups in your metagenomic research.


