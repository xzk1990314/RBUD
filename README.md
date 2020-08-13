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

**Note:** MDGM database is essential for RBUD method. For the purpose to study bacteria, viruses and fungus separately, we provided codes to establish a separated database that can save time and improve degree of accuracy. Befoe you use RBUD method, you should build your development envrionment in Linux. In addition, SOAPaligner/soap2 should be installed. 

## **2. MDGM database** 

### **2.1 MDGM database introduction**

MDGM database includes Species dataset and Functional dataset. The functional dataset contains gene annotation information from different databases. The folder is called Annotation Database, including information from Antibiotic Resistance Gene Database (ARDB), The Comprehensive Antibiotic Resistance Database(CARD), Carbohydrate-Active enZymes Database (CAZy)，eggNOG, Kyoto Encyclopedia of Genes and Genomes (KEGG), Universal Protein (UniProt), Metabolic Pathways From all Domains of Life (MetaCyc). In additon, we provide a file called linkdb_genes, the users can link it with gene ID, GI ID and protein ID to find the comments.The method of establishing the database and the running code gave as follows.

### **2.2 The code of MDGM database construction**

### **2.2.1 Download microbial data**

you should download all microorganism data from NCBI,the data should includ the suffixes of .fna, .gbk, .ffn, .gff. In this process, all microorganism contain archaea, bacteria, fungi, virus.

### **2.2.2 Create binary of reference genome for species dataset**

it is for creating binary of reference genome for species dataset.(2bwt-builder-species.sh)

	> cd /MDGM database/Species dataset/Sequence/NCBI_Species
	> tar -xvf all.fna.tar
	> cat * >all.fna
	> 2bwt-builder all.fna
	
**Note:** In this step, you can divided all.fna into four groups to improve the speed of alignment.

### **2.2.3 Extract the species taxonomy annotation**

it is for extracting the species annotation.(species_annotation.pl)

	> cd /MDGM database/Species dataset/Annotation/NCBI species annotation
	> tar -xvf all.gbk.tar.gz
	> perl /RBUD code/MDGM database code/species_annotation.pl
	
**Note:** In this process, you could obtain the species taxnomy annotation information of all microorganisms, including Phylum, Class, Order, Family, Genus, Species.(species annotation) 

### **2.2.4 Create binary of reference genome for functional dataset**

it is for creating binary of reference genome for functional dataset.(2bwt-builder-CDS.sh)

	> cd /MDGM database/Functional dataset/Sequence/NCBI_CDS
	> tar -xvf all.ffn.tar
	> cat * >all.ffn
	> 2bwt-builder all.ffn

**Note:** In this step, you can divided all.ffn into four groups to improve the speed of alignment.

### **2.2.5 Extract the functional annotation**

it is for extracting the functional annotation of microorganisms.(gene_annotation.pl)

	> cd /MDGM database/Functional dataset/Annotation/NCBI functional annotation
	> tar -xvf all.gff.tar.gz
	> cat * > all.gff
	> perl /RBUD code/MDGM database code/gene_annotation.pl
	
**Note:** In this process, you could obtain the functional annotation information of microbial genes, including gene ID, GI ID, protein ID, COG functions, genetically coded function, location of species origin.(CDS.gff/gene.gff)

### **2.2.6 Link different functional database with gene ID, GI ID and protein ID to find the comments**

#### 2.2.6.1 Transform different name in gene2accession

1> Download https://ftp.ncbi.nih.gov/gene/DATA/gene2accession file

2> According to gene ID and protein ID in gene functional annotation file (CDS.gff or gene.gff), you can extract protein_accession (the six column),     protein_gi (the seven column), genomic_nucleotide_accession (the eight column), genomic_nucleotide_gi (the nine column) in the gene2accession file. And then you can transform different name in gene2accession.

#### 2.2.6.2 Functional annotation of ARDB database

1> genomeblast.tab, ar_genes.tab, class2info.tab, resistance_profile.tab are opened in ARDB folder of linkdb_genes folder. Protein ID and Gene ID are obtained from CDS.gff file.

2> For genomeblast.tab file, according to protein ID (the second column), you can obtained special ID of ARDB database (the third column).

3> For ar_genes.tab file, you can transfer special ID (the first column) to antibiotic resistance gene name (the second column).

4> Based on antibiotic resistance gene name,  types of antibotic and functions of these genes are performed in resistance_profile.tab and classinfo.tab files separately.

**another way1:**

	> perl /RBUD code/MDGM database code/ardbAnno.pl
	
**another way2:**

	> cd /MDGM database/Functional dataset/Annotation/Annotation Database/ARDB/ardbAnno1.0/blastdb
	> formatdb -i resisGenes.pfasta -p T
	> blastall -i /MDGM database/Functional dataset/Sequence/all.ffn -d resisGenes.pfasta -o /MDGM database/Functional dataset/Annotation/Annotation Database/ARDB/resisGenes -p blastx -m 8 -e 1e-10

5> Based on the above results, you could establish ARDB annotation file.

#### 2.2.6.3 Functional annotation of CARD database

1> According to gene ID and protein ID (CDS.gff/gene.gff), you can obtain protein accession in gene2accession file.

2> Based on protein accession, you can obtain ARO name (the third column) and ARO accession (the fourth column) from aro_index file.

3> aro and CARD-aro_categories_index files provide antibiotic name (the second column)/description (the third column) and ARO Category Name (the third column) separately.

**another way:**

	> cd /MDGM database/Functional dataset/Annotation/Annotation Database/CARD
	> formatdb -i ARmeta-genes.fa -p T
	> blastall -i /MDGM database/Functional dataset/Sequence/all.ffn -d ARmeta-genes.fa -o /MDGM database/Functional dataset/Annotation/Annotation Database/CARD/CARD-ARmeta-genes -p blastx -m 8 -e 1e-10
	
4> Based on the above results, you could establish CARD annotation file.

#### 2.2.6.4 Functional annotation of CAZy database

1> According to gene ID and protein ID (CDS.gff/gene.gff), you can obtain genebank name in gene2accession file.

2>Based on genebank name, you can obatin family of CAzy database （the second column）, EC number (the third column) and protein name (the fourth column) in CAZyDB-ec-info.txt.07-20-2017 file.

3> Category of protein (the third column) in CAZy database is obtained from CAZyDB-ec-info.txt file. 

4> FamInfo.txt file details all information of each family of CAZy database, including ncbi-cdd, cazy-class, cazy-note and cazy-activities.

**another way:**

1> download all microorganism data with the suffixes of .faa file from NCBI.

2> download all.hmm.ps.len, dbCAN-fam-HMMs.txt, hsmmscan-parser.sh form dbCAN website.

	wget http://csbl.bmb.uga.edu/dbCAN/download/all,hmm.ps.len
	wget http://csbl.bmb.uga.edu/dbCAN/download/dbCAN-fam-HMMs.txt
	wget http://csbl.bmb.uga.edu/dbCAN/download/hmmscan-parser.sh

3> download hmmer software from ftp://selab.janelia.org/pub/software/hmmer3/3.0/hmmer-3.0.tar.gz

4> annotation for target protein (running dbCAN)
	> cd /MDGM database/Functional dataset/Annotation/Annotation Database/CAZy
	> hmmpress dbCAN-fam-HMMs.txt
	> hmmscan dbCAN-fam-HMMs.txt species_protein.faa >CAZyme.dbCAN
	> hmmscan -parser.sh CAZyme.dbCAN >CAZyme.annot

5> Based on the above results, you could establish CAZy annotation file.

**Note:** you can use dimond or blast software to align species_protein.faa with CAZyDB.fa. Moreover, you also run dbcan software to find the annotation of target protein.

#### 2.2.6.5 Functional annotation of COG database (eggNOG database)

**The Best Way**

you can obtain COG number directly in CDS.gff file. In COG database, cognames2003-2014.tab, fun2003-2014.tab, prot2003-2014.tab are used to find COG functions and protein name.
	
**another way1**

1> According to gene ID and protein ID, you can extract the GI ID from gene2accession file.

2> Based on GI ID (the first column), you can extract COG number from COG-gi-number.csv file.

3> cognames2003-2014.tab, fun2003-2014.tab, prot2003-2014.tab are used to find COG functions and protein name in COG database.

**another way2**

	> cd /MDGM database/Functional dataset/Annotation/Annotation Database/COG
	> formatdb -i prot2003-2014.fa -p T
	> blastall -i /MDGM database/Functional dataset/Sequence/all.ffn -d pro2003-2014.fa -o /MDGM database/Functional dataset/Annotation/Annotation Database/COG/COG_annotation -p blastx -m 8 -e 1e-10
	
4> Based on the above results, you could establish COG annotation file.

#### 2.2.6.6 Functional annotation of KEGG database

1> According to gene ID and GI ID in CDS.gff file, you can extract KEGG name from KEGG-geneid.list or KEGG-gi.list files.

2> Based on KEGG name, you can obtain KO number, uniprot number, pathway from genes_ko.list, genes_uniprot.lsit and genes_pathway,list files in KEGG database.

**another way**

	> cd /MDGM database/Functional dataset/Annotation/Annotation Database/KEGG
	> formatdb -i genes.nuc -p T
	> blastall -i /MDGM database/Functional dataset/Sequence/all.ffn -d genes.nuc -o /MDGM database/Functional dataset/Annotation/Annotation Database/KEGG/KEGG_annotation -p blastx -m 8 -e 1e-10
	
4> Based on the above results, you could establish COG annotation file.

#### 2.2.6.7 Functional annotation of UniProt database and MetaCyc database

1> download ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/idmapping/idmapping_selected.tab file. Here, the tab-delimited table which includes the following mappings delimited by tab: UniProtKB-AC, UniProtKB-ID, GeneID (EntrezGene), RefSeq, GI, PDB, GO, UniRef100, UniRef90, UniRef50, UniParc, PIR, NCBI-taxon, MIM, UniGene, PubMed, EMBL, EMBL-CDS, Ensembl, Ensembl_TRS, Ensembl_PRO, Additional PubMed.

2> According to Gene ID and  GI ID (CDS.gff), you can extract the name of UniRef100, UniRef90, UniRef50.

3> Based on the name of UniRef100, UniRef90, UniRef50, you can obtain the annotation of protein in UniProt database.

**another way**

1> download ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50 file.

2> alignment

	> cd /MDGM database/Functional dataset/Annotation/Annotation Database/UniProt/Uniref
	> diamond makedb --in uniref50.fasta -d uniref50.fasta -p 32
	> diamond blastx -c 1 --db uniref50.fasta.dmnd -t /tmp -p 34 -q /MDGM database/Functional dataset/Sequence/all.ffn -o uniref_annotation --outfmt 6 

3> Based on the above results, you could establish UniProt annotation file.

**Note:** UniRef90 and UniRef100 are the same operation with Uniref50.

#### MetaCyc

Based on the above Uniref name, you can obtain metacyc name by metacyc_reactions_level4ec_only.uniref, metacyc_pathways_structured, metacyc_pathway. And then the annotation file of MetaCyc is built. 

### **2.2.7 The complete MDGM database**

According to gene ID, GI ID and protein ID to find the annotation in databases of ARDB, CARD, CAZy, COG, KEGG, UniProt, MetaCyc. If you can't find the annotation, it is the best way to align DNA sequence to these databases for functional annotation using blast software. Finally, based on the above biological information to construct MDGM database.

------------------------------------------------------------------------

## **3. Runing RBUD Method** 

### **3.1 Calculate the length of microbial genome and CDS sequence**

Based on above sequence of all.fna and all.ffn, a script of length.py is used to extract the length of every microbial genome and CDS sequence separately.

	> python /RBUD code/RBUD running code/length.py
	
**Note:** In this step, species_length and function_length are obtained for the length of every microbial gneome and CDS sequence separately using this script.

### **3.2 Alignment next generation sequencing data with MDGM database**

Align query to genome












