**Welcome to the RBUD 1.0 tutorial**
===========================================
RBUD - a new functional potential analysis approach for whole microbial genome shotgun sequencing

Zhikai Xing edited this page on 8 August, 2020. first version

**Website: https://github.com/DMsiast/RBUD.git**

Welcome to the RBUD 1.0 tutorial, which provides
software, documentation, and tutorial for method for microbial
community profiling developed by Zhikai Xing. Most tools are supported
as individual software packages (typically `Python`,`Perl` or `R`). 

------------------------------------------------------------------------
## **1. Requirements**

- Linux
- SOAPaligner/soap2
- BLAST/Diamond
- MDGM database

Note: MDGM database is essential for RBUD method. For the purpose to study bacteria, viruses and fungus separately, we provided codes to establish a separated database that can save time and improve degree of accuracy. Befoe you use RBUD method, you should build your development envrionment in Linux. In addition, SOAPaligner/soap2 should be installed. 

## **2. MDGM database** 

### **2.1 MDGM database introduction**

MDGM database includes Species dataset and Functional dataset. The functional dataset contains gene annotation information from different databases. The folder is called Annotation Database, including information from Antibiotic Resistance Gene Database (ARDB), The Comprehensive Antibiotic Resistance Database(CARD), Carbohydrate-Active enZymes Database (CAZy)，eggNOG, Kyoto Encyclopedia of Genes and Genomes (KEGG), Universal Protein (UniProt), Metabolic Pathways From all Domains of Life (MetaCyc). In additon, we provide a file called linkdb_genes, the users can link it with gene ID, GI ID and protein ID to find the comments.The method of establishing the database and the running code gave as follows.

### **2.2 The code of MDGM database construction**

### **2.2.1 download microbial data**

you should download all microorganism data from NCBI,the data should includ the suffixes of .fna, .gbk, .ffn, .gff. In this process, all microorganism contain archaea, bacteria, fungi, virus.

### **2.2.2 create binary of reference genome for species dataset**

it is for creating binary of reference genome for species dataset.(2bwt-builder-species.sh)

	> cd /MDGM database/Species dataset/Sequence/NCBI_Species
	> tar -xvf all.fna.tar
	> cat * >all.fna
	> 2bwt-builder all.fna

### **2.2.3 extract the species taxonomy annotation**

it is for extracting the species annotation.(species_annotation.pl)

	> perl species_annotation.pl
Note: In this process, you could obtain the species taxnomy annotation information of all microorganisms, including Phylum, Class, Order, Family, Genus, Species.(species annotation) 

### **2.2.4 create binary of reference genome for functional dataset**

it is for creating binary of reference genome for functional dataset.(2bwt-builder-CDS.sh)

	> cd /MDGM database/Functional dataset/Sequence/NCBI_CDS
	> tar -xvf all.ffn.tar
	> cat * >all.ffn
	> 2bwt-builder all.ffn
	
### **2.2.5 extract the functional annotation**

it is for extracting the functional annotation of microorganisms.(gene_annotation.pl)

	> perl gene_annotation.pl
Note: In this process, you could obtain the functional annotation information of microbial genes, including gene ID, GI ID, protein ID, COG functions, genetically coded function, location of species origin.(CDS.gff/gene.gff)

### **2.2.6 link different functional database with gene ID, GI ID and protein ID to find the comments**

#### 1. transform different name in gene2accession

1> Download https://ftp.ncbi.nih.gov/gene/DATA/gene2accession file

2> According to gene ID and protein ID in gene functional annotation file (CDS.gff or gene.gff), you can extract protein_accession (the six column), protein_gi (the seven column), genomic_nucleotide_accession (the eight column), genomic_nucleotide_gi (the nine column) in the gene2accession file. And then you can transform different name in gene2accession.

#### 2. functional annotation of ARDB database

1> genomeblast.tab, ar_genes.tab, class2info.tab, resistance_profile.tab are opened in ARDB folder of linkdb_genes folder.

2> For genomeblast.tab file, according to protein ID (the second column), you can obtained special ID of ARDB database (the third column).

3> For ar_genes.tab file, you can transfer special ID (the first column) to antibiotic resistance gene name (the second column).

4> Based on antibiotic resistance gene name,  types of antibotic and functions of these genes are performed in resistance_profile.tab and classinfo.tab files separately.

another way1:

	> perl ardbAnno.pl
	
another way2:

	> cd /MDGM database/Functional dataset/Annotation/Annotation Database/ARDB/ardbAnno1.0/blastdb
	> formatdb -i resisGenes.pfasta -p T
	> blastall -i /MDGM database/Functional dataset/Sequence/all.ffn -d resisGenes.pfasta -o /MDGM database/Functional dataset/Annotation/Annotation Database/ARDB/resisGenes -p blastx -m 8 -e 1e-10

5> Based on the results above, you could establish ARDB annotation file.

#### 2. functional annotation of CARD database

.

2> 

3>Based on protein accession, you can obtain ARO name (the third column) and ARO accession (the fourth column) from aro_index file. And then, aro and CARD-aro_categories_index files provide antibiotic name (the second column)/description (the third column) and ARO Category Name (the third column) separately. 

### 3. functional annotation of CAZy database

1> 




Note: According to gene ID, GI ID and protein ID to find the annotation in databases of ARDB, CARD, CAZy, COG, KEGG, UniProt, MetaCyc. If you can't find the annotation, it is the best way to align DNA sequence to these databases for functional annotation using blast software. Finally, based on the above biological information to construct MDGM database.

## **Downstream analysis and statistics**

The methods in this section generally provide quantitative models for
interpreting microbial community profiles as generated by the methods
above. For example, this includes identifying significant associations
of sample metadata (phenotype, environment, health status, etc.) with
microbial taxonomic or functional composition. Please click on the links
below for detailed tutorials:

[![HAllA](https://github.com/biobakery/biobakery/blob/master/images/2517623131-HAllA.png)](https://github.com/biobakery/biobakery/wiki/halla) [![ARepA](https://github.com/biobakery/biobakery/blob/master/images/4142907121-ARepA.png)](http://huttenhower.sph.harvard.edu/arepa/tutorial) [![CCREPE](https://github.com/biobakery/biobakery/blob/master/images/3539496555-CCREPE.png)](https://github.com/biobakery/biobakery/wiki/ccrepe) [![LEfSe](https://github.com/biobakery/biobakery/blob/master/images/2196154061-LEfSe.png)](https://github.com/biobakery/biobakery/wiki/lefse) [![MaAsLin](https://github.com/biobakery/biobakery/blob/master/images/2350879162-MaAsLin.png)](https://github.com/biobakery/biobakery/wiki/maaslin2) [![MMUPHin](https://github.com/biobakery/biobakery/blob/master/images/2357044001-MMUPHin_alt.png)](https://bioconductor.org/packages/release/bioc/vignettes/MMUPHin/inst/doc/MMUPHin.html) [![microPITA](https://github.com/biobakery/biobakery/blob/master/images/255233476-MicroPITA.png)](https://github.com/biobakery/biobakery/wiki/micropita) [![SparseDOSSA](https://github.com/biobakery/biobakery/blob/master/images/3488299857-SparseDOSSA.png)](https://github.com/biobakery/biobakery/wiki/SparseDOSSA) [![BAN](https://github.com/biobakery/biobakery/blob/master/images/1871766089-BAnOCC.png)](https://github.com/biobakery/biobakery/wiki/banocc)

## **Infrastructure and utilities**

[![GraPhlAn](https://github.com/biobakery/biobakery/blob/master/images/3212034723-GraPhlAn.png)](https://github.com/biobakery/biobakery/wiki/graphlan) [![KneadData](https://github.com/biobakery/biobakery/blob/master/images/3968267398-KneadData.png)](https://github.com/biobakery/biobakery/wiki/kneaddata) [![AnADAMA](https://github.com/biobakery/biobakery/blob/master/images/1668386270-AnADAMA.png)](https://github.com/biobakery/biobakery/wiki/anadama2) [![workflows](https://github.com/biobakery/biobakery/blob/master/images/3531676205-workflows.png)](https://github.com/biobakery/biobakery/wiki/biobakery_workflows)

------------------------------------------------------------------------

Notes
-----

-   In case of any issue with any tool, please feel free post it
    in [bioBakery Support Forum](http://forum.biobakery.org/).
-   [Source code repository](https://github.com/biobakery)
-   [Galaxy server](http://huttenhower.sph.harvard.edu/galaxy)
-   [Lab website](https://huttenhower.sph.harvard.edu)

