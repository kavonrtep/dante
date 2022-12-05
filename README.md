# Domain based annotation of transposable elements  - DANTE #

### Authors 
 Nina Hostakova, Petr Novak, Pavel Neumann, Jiri Macas
 Biology Centre CAS, Czech Republic
 
 
### Introduction

* Protein Domains Finder [dante.py]
	* Script performs scanning of given DNA sequence(s) in (multi)fasta format in order to discover protein domains using our protein domains database.
	* Domains searching is accomplished engaging LASTAL alignment tool.
	* Domains are subsequently annotated and classified - in case certain domain has multiple annotations assigned, classifation is derived from the common classification level of all of them. 	
			
* Proteins Domains Filter [dante_gff_output_filtering.py]
	* filters GFF3 output from previous step to obtain certain kind of domain and/or allows to adjust quality filtering  
        
### DEPENDENCIES ###

* python3.4 or higher with packages:	
	* numpy
	* biopython
* [lastal](http://last.cbrc.jp/doc/last.html) 744 or higher
* ProfRep/DANTE modules:
	* configuration.py 

### Installation using conda ####
[![Anaconda-Server Badge](https://anaconda.org/petrnovak/dante/badges/version.svg)](https://anaconda.org/petrnovak/dante)

     conda install -c conda-forge -c bioconda -c petrnovak dante


### Protein Domains Finder ###

This tool provides **preliminary** output of all domains types which are not filtered for quality.

#### INPUTS ####

* DNA sequence [multiFasta]
		
#### OUTPUTS ####
		
* **All protein domains GFF3** - individual domains are reported per line as regions (start-end) on the original DNA sequence including the seq ID and strand orientation. The last "Attributes" column contains several comma-separated information related to the domain annotation, alignment and its quality. This file can undergo further filtering using Protein Domain Filter tool.		

#### USAGE ####

		usage: dante.py [-h] -q QUERY -pdb PROTEIN_DATABASE -cs
								  CLASSIFICATION [-oug DOMAIN_GFF] [-nld NEW_LDB]
								  [-dir OUTPUT_DIR] [-thsc THRESHOLD_SCORE]
								  [-wd WIN_DOM] [-od OVERLAP_DOM]
								  
		optional arguments:
		  -h, --help            show this help message and exit
		  -oug DOMAIN_GFF, --domain_gff DOMAIN_GFF
								output domains gff format (default: None)
		  -nld NEW_LDB, --new_ldb NEW_LDB
								create indexed database files for lastal in case of
								working with new protein db (default: False)
		  -dir OUTPUT_DIR, --output_dir OUTPUT_DIR
								specify if you want to change the output directory
								(default: None)
		  -thsc THRESHOLD_SCORE, --threshold_score THRESHOLD_SCORE
								percentage of the best score in the cluster to be
								tolerated when assigning annotations per base
								(default: 80)
		  -wd WIN_DOM, --win_dom WIN_DOM
								window to process large input sequences sequentially
								(default: 10000000)
		  -od OVERLAP_DOM, --overlap_dom OVERLAP_DOM
								overlap of sequences in two consecutive windows
								(default: 10000)

		required named arguments:
		  -q QUERY, --query QUERY
								input DNA sequence to search for protein domains in a
								fasta format. Multifasta format allowed. (default:
								None)
		  -pdb PROTEIN_DATABASE, --protein_database PROTEIN_DATABASE
								protein domains database file (default: None)
		  -cs CLASSIFICATION, --classification CLASSIFICATION
								protein domains classification file (default: None)


		
#### HOW TO RUN EXAMPLE ####
		./protein_domains.py -q PATH_TO_INPUT_SEQ -pdb PATH_TO_PROTEIN_DB -cs PATH_TO_CLASSIFICATION_FILE
		
	 When running for the first time with a new database use -nld option allowing lastal to create indexed database files:

         -nld True

	use other arguments if you wish to rename your outputs or they will be created automatically with standard names 
	
### Protein Domains Filter ###
		
The script performs Protein Domains Finder output filtering for quality and/or extracting specific type of protein domain or mobile elements of origin. For the filtered domains it reports their translated protein sequence of original DNA.

WHEN NO PARAMETERS GIVEN, IT PERFORMS QUALITY FILTERING USING THE DEFAULT PARAMETRES (optimized for Viridiplantae species)

#### INPUTS ####
* GFF3 file produced by protein_domains.py OR already filtered GFF3
	
#### Filtering options ####
* QUALITY: 
	- Min relative length of alignemnt to the protein domain from DB (without gaps)
	- Identity 
	- Similarity (scoring matrix: BLOSUM80)
	- Interruption in the reading frame (frameshifts + stop codons) per every starting 100 AA
	- Max alignment proportion to the original length of database domain sequence
* DOMAIN TYPE: 'Name' attribute in GFF - see choices bellow
Records for ambiguous domain type (e.g. INT/RH) are filtered out automatically

* MOBILE ELEMENT TYPE:
arbitrary substring of the element classification ('Final_Classification' attribute in GFF)
		
#### OUTPUTS ####
* filtered GFF3 file
* fasta file of translated protein sequences for the aligned domains that match the filtering criteria 
	! as it is taken from the best hit alignment reported by LAST, it does not neccessary cover the whole region reported as domain in GFF
	
#### USAGE ####		

		usage: dante_gff_output_filtering.py [-h] -dg DOM_GFF [-ouf DOMAINS_FILTERED]
                            [-dps DOMAINS_PROT_SEQ]
                            [-thl {float range 0.0..1.0}]
                            [-thi {float range 0.0..1.0}]
                            [-ths {float range 0.0..1.0}] [-ir INTERRUPTIONS]
                            [-mlen MAX_LEN_PROPORTION]
                            [-sd {All,GAG,INT,PROT,RH,RT,aRH,CHDCR,CHDII,TPase,YR,HEL1,HEL2,ENDO}]
                            [-el ELEMENT_TYPE] [-dir OUTPUT_DIR]



		optional arguments:
		  -h, --help            show this help message and exit
		  -ouf DOMAINS_FILTERED, --domains_filtered DOMAINS_FILTERED
								output filtered domains gff file (default: None)
		  -dps DOMAINS_PROT_SEQ, --domains_prot_seq DOMAINS_PROT_SEQ
								output file containg domains protein sequences
								(default: None)
		  -thl {float range 0.0..1.0}, --th_length {float range 0.0..1.0}
								proportion of alignment length threshold (default:
								0.8)
		  -thi {float range 0.0..1.0}, --th_identity {float range 0.0..1.0}
								proportion of alignment identity threshold (default:
								0.35)
		  -ths {float range 0.0..1.0}, --th_similarity {float range 0.0..1.0}
								threshold for alignment proportional similarity
								(default: 0.45)
		  -ir INTERRUPTIONS, --interruptions INTERRUPTIONS
								interruptions (frameshifts + stop codons) tolerance
								threshold per 100 AA (default: 3)
		  -mlen MAX_LEN_PROPORTION, --max_len_proportion MAX_LEN_PROPORTION
								maximal proportion of alignment length to the original
								length of protein domain from database (default: 1.2)
		  -sd {All,GAG,INT,PROT,RH,RT,aRH,CHDCR,CHDII,TPase,YR,HEL1,HEL2,ENDO}, --selected_dom {All,GAG,INT,PROT,RH,RT,aRH,CHDCR,CHDII,TPase,YR,HEL1,HEL2,ENDO}
								filter output domains based on the domain type
								(default: All)
		  -el ELEMENT_TYPE, --element_type ELEMENT_TYPE
								filter output domains by typing substring from
								classification (default: )
		  -dir OUTPUT_DIR, --output_dir OUTPUT_DIR
								specify if you want to change the output directory
								(default: None)

		required named arguments:
		  -dg DOM_GFF, --dom_gff DOM_GFF
								basic unfiltered gff file of all domains (default:
								None)



#### HOW TO RUN EXAMPLE ####
e.g. getting quality filtered integrase(INT) domains of all gypsy transposable elements:
	
	./domains_filtering.py -dom_gff PATH_TO_INPUT_GFF -pdb PATH_TO_PROTEIN_DB -cs PATH_TO_CLASSIFICATION_FILE --selected_dom INT --element_type Ty3/gypsy 


### Extract Domains Nucleotide Sequences ###

This tool extracts nucleotide sequences of protein domains from reference DNA based on DANTE's output. It can be used e.g. for deriving phylogenetic relations of individual mobile elements classes within a species. 

#### INPUTS ####

* original DNA sequence in multifasta format to extract the domains from 
* GFF3 file of protein domains (**DANTE's output** - preferably filtered for quality and specific domain type)
* Domains database classification table (to check the classification level)

#### OUTPUTS ####

* fasta files of domains nucleotide sequences for individual transposons lineages
* txt file of domains counts extracted for individual lineages

**- For GALAXY usage all concatenated in a single fasta file**

#### USAGE ####	
		usage: dante_gff_to_dna.py [-h] -i INPUT_DNA -d DOMAINS_GFF -cs
			CLASSIFICATION [-out OUT_DIR] [-ex EXTENDED]

		optional arguments:
		  -h, --help            show this help message and exit
		  -i INPUT_DNA, --input_dna INPUT_DNA
								path to input DNA sequence
		  -d DOMAINS_GFF, --domains_gff DOMAINS_GFF
								GFF file of protein domains
		  -cs CLASSIFICATION, --classification CLASSIFICATION
								protein domains classification file
		  -out OUT_DIR, --out_dir OUT_DIR
								output directory
		  -ex EXTENDED, --extended EXTENDED
								extend the domains edges if not the whole datatabase
								sequence was aligned

#### HOW TO RUN EXAMPLE ####
	./extract_domains_seqs.py --domains_gff PATH_PROTEIN_DOMAINS_GFF --input_dna PATH_TO_INPUT_DNA  --classification PROTEIN_DOMAINS_DB_CLASS_TBL --extended True

	







