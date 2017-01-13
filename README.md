# TOOLS TO FIND REPETITIONS AND PROTEIN DOMAINS IN A SEQUENCE #


## PROFREP ##

This tool uses former clustering output (.cls files) from Repeat Explorer (list of all clusters and belonging reads), annotation table that assigns repetitive class to every cluster and list of all sequencing reads to create repetitive profile of a sequence as well as to report protein domains. You can either provide your own files or you can exploit our prepared datasets for certain species listed above. At first Blast+ run similarity search to find hits in a given sequence from the reads database which are subsequently filtered to gain only that ones which possess high similarity and appropriate length. Other parameters of blast search are also adjustable. Using window parameter can speed up the analysis especially in case of having large input data as it allows paralell processing of the sequence with a certain overlap. Be careful to set the overlap between windows at least of the length of reads so that hits on the border of windows can be found and of course the window size must be greater than overlap itself. 
When using our annotation datasets copy numbers are reported by default, if you wish to convert hits to copy numbers in your custom annotation data you will be asked to provide genome size of the species so that the sequencing coverage can be calculated. 

Protein domains search utilizes LAST tool to find hits between input sequence and our database of protein domains - this is accomplished using protein_domains_pd.py [SEE BELLOW] which also serves as a module for PROFREP.

OUTPUT:		
	
* table reporting copy numbers/hits at every position for all repeats occuring in the sequence
* GFF file reporting N regions in the sequence
* GFF file reporting repetitive regions of a certain length and above copy numbers/hits threshold
* GFF file reporting protein domains, classification of domain, chain orientation and alignment sequence
* HTML reporting domains and profiles summary visualization that can be further explored using a link to JBrowse

### Dependencies ###

* python 3.4 or higher 
* python packages:
	* numpy
	* matplotlib
* [BLAST 2.2.28+] (https://www.ncbi.nlm.nih.gov/books/NBK279690/) or higher
* [JBrowse 1.12.1] (http://jbrowse.org/) or higher
* protein_domains_pd.py
* gff.py
* visualization.py
* configuration.py 

## PROTEIN DOMAINS TOOLS ##
Two tools are available to explore protein domains in your sequence:
* Protein Domains Finder [protein_domains_pd.py]
	* makes overall scan of a sequence and finds all domains using our protein domains database 
	* classifies each domain with a certain repetitive annotation
		
		OUTPUT:		
		
		* GFF3 format reporting positions of domains, type, repeat class and lineage, chain orientation, protein sequence and alignment with quality of alignment parameters
		* summary txt of all domains amounts	
* Proteins Domains Filter [domains_filtering.py]
	* filters GFF3 output from previous step to obtain certain kind of domain and allows to adjust quality filtering  
	
		OUTPUT:
	
		* filtered GFF3 format 
		* list of protein sequences of the filtered domains in fasta format

Both are implemented on galaxy web enviroment or can be used as standalone python
scripts from the command line   
        
### Dependencies ###

* python 3.4 or higher 
* python packages:
	* numpy
* [lastal] (http://last.cbrc.jp/doc/last.html) 744
* configuration.py 


### Running the scripts ###


* Protein Domains Finder
	
		usage: protein_domains_pd.py [-h] -q QUERY [-pdb PROTEIN_DATABASE]
                             [-cs CLASSIFICATION] [-oug DOMAIN_GFF]
                             [-npd NEW_PDB] [-sum SUMMARY_FILE]
                             [-dir OUTPUT_DIR]

		optional arguments:
		  -h, --help            show this help message and exit
		  -q QUERY, --query QUERY
								reference sequence to find protein domains
		  -pdb PROTEIN_DATABASE, --protein_database PROTEIN_DATABASE
								protein domains database file
		  -cs CLASSIFICATION, --classification CLASSIFICATION
								protein domains classification file
		  -oug DOMAIN_GFF, --domain_gff DOMAIN_GFF
								output domains gff format
		  -npd NEW_PDB, --new_pdb NEW_PDB
								create new protein database for lastal
		  -sum SUMMARY_FILE, --summary_file SUMMARY_FILE
								output summary file containing overview of amount of
								domains in individual seqs
		  -dir OUTPUT_DIR, --output_dir OUTPUT_DIR
								specify if you want to change the output directory

	HOW TO RUN EXAMPLE:

		./protein_domains_pd.py -q PATH_TO_INPUT_SEQ -pdb PATH_TO_PROTEIN_DB -cs PATH_TO_CLASSIFICATION_FILE
	when running for the first time use -npd option to create a new lastal database

		-npd True
	use other arguments if you wish to rename your outputs or they will be created automatically with standard names 
	
* Protein Domains Filter

		usage: domains_filtering.py [-h] [-dom_gff DOMAIN_GFF] [-ouf DOMAINS_FILTERED]
                            [-dps DOMAINS_PROT_SEQ]
                            [-thl {float range 0.0..1.0}]
                            [-thi {float range 0.0..1.0}] [-fr FRAMESHIFTS]
                            [-sd {All,Ty1-GAG,Ty1-INT,Ty1-PROT,Ty1-RH,Ty1-RT,Ty3-GAG,Ty3-INT,Ty3-PROT,Ty3-RT,Ty3-gRH,Ty3-aRH,Ty3-CHDII,Ty3-CHDCR}]
                            [-dir OUTPUT_DIR]

		optional arguments:
		  -h, --help            show this help message and exit
		  -dom_gff DOMAIN_GFF, --domain_gff DOMAIN_GFF
								basic unfiltered gff file of all domains
		  -ouf DOMAINS_FILTERED, --domains_filtered DOMAINS_FILTERED
								output filtered domains gff file
		  -dps DOMAINS_PROT_SEQ, --domains_prot_seq DOMAINS_PROT_SEQ
								output file containg domains protein sequences
		  -thl {float range 0.0..1.0}, --th_length {float range 0.0..1.0}
								proportion of alignment length threshold
		  -thi {float range 0.0..1.0}, --th_identity {float range 0.0..1.0}
								proportion of alignment identity threshold
		  -fr FRAMESHIFTS, --frameshifts FRAMESHIFTS
								frameshifts tolerance threshold per 100 bp
		  -sd {All,Ty1-GAG,Ty1-INT,Ty1-PROT,Ty1-RH,Ty1-RT,Ty3-GAG,Ty3-INT,Ty3-PROT,Ty3-RT,Ty3-gRH,Ty3-aRH,Ty3-CHDII,Ty3-CHDCR}, --selected_dom {All,Ty1-GAG,Ty1-INT,Ty1-PROT,Ty1-RH,Ty1-RT,Ty3-GAG,Ty3-INT,Ty3-PROT,Ty3-RT,Ty3-gRH,Ty3-aRH,Ty3-CHDII,Ty3-CHDCR}
								filter output domains based on the domain type
		  -dir OUTPUT_DIR, --output_dir OUTPUT_DIR
								specify if you want to change the output directory


	HOW TO RUN EXAMPLE
					
		./domains_filtering.py [-dom_gff PATH_TO_DOM_GFF]
	
	NOTE! - When you run protein_domain_pd.py script with default parameters as written in EXAMPLE (meaning you do not specify your own outputs) then you do not need to enter path to GFF file (-dom_gff), it will load it automatically from previous domain search. Otherwise you have to specify the whole path to the GFF file that you wish to filter.

