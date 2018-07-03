# REPEATS ANNOTATION TOOLS FOR ASSEMBLIES  #


## 1. PROFREP ##
*- **PROF**iles of **REP**eats -* 


The ProfRep main tool engages outputs of RepeatExplorer for repeats annotation in DNA sequences (typically assemblies but not necessarily). Moreover, it provides repetitive profiles of the sequence, pointing out quantitative representation of individual repeats along the sequence as well as the overall repetitiveness.

### Dependencies ###


* python 3.4 or higher with packages:
	* numpy
	* matplotlib
	* biopython
* [BLAST 2.2.28+](https://www.ncbi.nlm.nih.gov/books/NBK279690/) or higher
* [wigToBigWig](http://hgdownload.cse.ucsc.edu/admin/exe/)
* [cd-hit](http://weizhongli-lab.org/cd-hit/)
* [JBrowse](http://jbrowse.org/install/) - **Only bin needed, does not have to be installed under a web server**
  
* ProfRep Modules:
	* gff.py
	* visualization.py
	* configuration.py 
	* protein_domains.py
	* domains_filtering.py
	

### Running ProfRep ###

	usage: profrep.py [-h] -q QUERY [-d DATABASE] [-a ANNOTATION_TBL] [-c CLS]
					  [-tbl DATASETS_TBL] [-id DB_ID] [-i IDENTICAL]
					  [-l ALIGN_LENGTH] [-m MAX_ALIGNMENTS] [-e E_VALUE]
					  [-df DUST_FILTER] [-ws WORD_SIZE] [-t TASK] [-n NEW_DB]
					  [-w WINDOW] [-o OVERLAP] [-pd PROTEIN_DOMAINS]
					  [-pdb PROTEIN_DATABASE] [-cs CLASSIFICATION] [-wd WIN_DOM]
					  [-od OVERLAP_DOM] [-thsc THRESHOLD_SCORE] [-lg LOG_FILE]
					  [-ouf OUTPUT_GFF] [-oug DOMAIN_GFF] [-oun N_GFF]
					  [-hf HTML_FILE] [-hp HTML_PATH] [-cn COPY_NUMBERS]
					  [-gs GENOME_SIZE] [-thr THRESHOLD_REPEAT]
					  [-ths THRESHOLD_SEGMENT] [-td TOOL_DIR] [-jb JBROWSE_BIN]


	optional arguments:
	  -h, --help            show this help message and exit

	required arguments:
	  -q QUERY, --query QUERY
							input DNA sequence in (multi)fasta format (default:
							None)
	  -d DATABASE, --database DATABASE
							blast database of all sequencing reads (default: None)
	  -a ANNOTATION_TBL, --annotation_tbl ANNOTATION_TBL
							clusters annotation table, tab-separated number of
							cluster and its classification (default: None)
	  -c CLS, --cls CLS     cls file containing reads assigned to clusters
							(hitsort.cls) (default: None)

	alternative required arguments - prepared datasets:
	  -tbl DATASETS_TBL, --datasets_tbl DATASETS_TBL
							table with prepared annotation datasets (default:
							None)
	  -id DB_ID, --db_id DB_ID
							annotation dataset ID (first column of datasets table)
							(default: None)

	optional arguments - BLAST Search:
	  -i IDENTICAL, --identical IDENTICAL
							blast filtering option: percentage indentity threshold
							between query and mapped read from db (default: 95)
	  -l ALIGN_LENGTH, --align_length ALIGN_LENGTH
							blast filtering option: minimal alignment length
							threshold in bp (default: 40)
	  -m MAX_ALIGNMENTS, --max_alignments MAX_ALIGNMENTS
							blast filtering option: maximal number of alignments
							in the output (default: 10000000)
	  -e E_VALUE, --e_value E_VALUE
							blast setting option: e-value (default: 1e-15)
	  -df DUST_FILTER, --dust_filter DUST_FILTER
							dust filters low-complexity regions during BLAST
							search (default: '20 64 1')
	  -ws WORD_SIZE, --word_size WORD_SIZE
							blast search option: initial word size for alignment
							(default: 11)
	  -t TASK, --task TASK  type of blast to be triggered (default: blastn)
	  -n NEW_DB, --new_db NEW_DB
							create a new blast database (default: False)

	optional arguments - Parallel Processing:
	  -w WINDOW, --window WINDOW
							sliding window size for parallel processing (default:
							5000)
	  -o OVERLAP, --overlap OVERLAP
							overlap for parallely processed regions, set greater
							than a read size (default: 150)

	optional arguments - Protein Domains:
	  -pd PROTEIN_DOMAINS, --protein_domains PROTEIN_DOMAINS
							use module for protein domains (default: True)
	  -pdb PROTEIN_DATABASE, --protein_database PROTEIN_DATABASE
							protein domains database (default: None)
	  -cs CLASSIFICATION, --classification CLASSIFICATION
							protein domains classification file (default: None)
	  -wd WIN_DOM, --win_dom WIN_DOM
							protein domains module: sliding window to process
							large input sequences sequentially (default: 10000000)
	  -od OVERLAP_DOM, --overlap_dom OVERLAP_DOM
							protein domains module: overlap of sequences in two
							consecutive windows (default: 10000)
	  -thsc THRESHOLD_SCORE, --threshold_score THRESHOLD_SCORE
							protein domains module: percentage of the best score
							within the cluster to significant domains (default:
							80)

	optional arguments - Output Paths:
	  -lg LOG_FILE, --log_file LOG_FILE
							path to log file (default: log.txt)
	  -ouf OUTPUT_GFF, --output_gff OUTPUT_GFF
							path to output gff of repetitive regions (default:
							output_repeats.gff)
	  -oug DOMAIN_GFF, --domain_gff DOMAIN_GFF
							path to output gff of protein domains (default:
							output_domains.gff)
	  -oun N_GFF, --n_gff N_GFF
							path to output gff of N regions (default:
							N_regions.gff)
	  -hf HTML_FILE, --html_file HTML_FILE
							path to output html file (default: output.html)
	  -hp HTML_PATH, --html_path HTML_PATH
							path to html extra files (default: output_dir)

	optional arguments - Copy Numbers/Hits :
	  -cn COPY_NUMBERS, --copy_numbers COPY_NUMBERS
							convert hits to copy numbers (default: False)
	  -gs GENOME_SIZE, --genome_size GENOME_SIZE
							genome size is required when converting hits to copy
							numbers and you use custom data (default: None)
	  -thr THRESHOLD_REPEAT, --threshold_repeat THRESHOLD_REPEAT
							threshold for hits/copy numbers per position to be
							considered repetitive (default: 5)
	  -ths THRESHOLD_SEGMENT, --threshold_segment THRESHOLD_SEGMENT
							threshold for the length of repetitive segment to be
							reported (default: 80)

	optional arguments - Enviroment Variables:
	  -td TOOL_DIR, --tool_dir TOOL_DIR
							Galaxy tool data directory in galaxy (default: None)
	  -jb JBROWSE_BIN, --jbrowse_bin JBROWSE_BIN
							path to JBrowse bin directory (default: None)
							
### ProfRep Data Preparation ###

In case of using custom input datasets these tools can be used for easy obtaining the correct files and to prepare the reduced datasets to speed up the main ProfRep analysis:

* Extract Data For ProfRep (extract_data_for_profrep.py)
* ProfRep DB Reducing (profrep_db_reducing.py)

### ProfRep Supplementary Tools ###

These additional tools can be used for further work with the ProfRep outputs: 

* ProfRep Refiner (profrep_refining.py) 
* ProfRep Masker (profrep_masking.py)
* GFF Region Selector (gff_selection.py)


## 2. DANTE ##
*- **D**omain based **AN**notation of **T**ransposable **E**lements -* 


* Protein Domains Finder [protein_domains.py]
	* Script performs scanning of given DNA sequence(s) in (multi)fasta format in order to discover protein domains using our protein domains database.
	* Domains searching is accomplished engaging LASTAL alignment tool.
	* Domains are subsequently annotated and classified - in case certain domain has multiple annotations assigned, classifation is derived from the common classification level of all of them. 	
			
* Proteins Domains Filter [domains_filtering.py]
	* filters GFF3 output from previous step to obtain certain kind of domain and/or allows to adjust quality filtering  

        
### Dependencies ###

* python3.4 or higher with packages:	
	* numpy
* [lastal](http://last.cbrc.jp/doc/last.html) 744 or higher
* ProfRep/DANTE modules:
	* configuration.py 


### Protein Domains Finder ###

		usage: protein_domains.py [-h] -q QUERY -pdb PROTEIN_DATABASE -cs
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

#### DEPENDENCIES ####
* python 3.4 or higher
* ProfRep modules:
	- configuration.py 

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
	
			
		usage: domains_filtering.py [-h] -dg DOM_GFF [-ouf DOMAINS_FILTERED]
                            [-dps DOMAINS_PROT_SEQ]
                            [-thl {float range 0.0..1.0}]
                            [-thi {float range 0.0..1.0}]
                            [-ths {float range 0.0..1.0}] [-ir INTERRUPTIONS]
                            [-mlen MAX_LEN_PERCENT]
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
			  -mlen MAX_LEN_PERCENT, --max_len_percent MAX_LEN_PERCENT
									maximal percentage of alignment length to the original
									length of protein domain from database (default: 120)
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


### GALAXY implementation ###

#### Dependencies ####

* python3.4 or higher with packages:	
	* numpy
	* matplotlib
	* biopython
* [BLAST 2.2.28+](https://www.ncbi.nlm.nih.gov/books/NBK279671/) or higher
* [LAST](http://last.cbrc.jp/doc/last.html) 744 or higher:
	* [download](http://last.cbrc.jp/)
	* [install](http://last.cbrc.jp/doc/last.html)
* [wigToBigWig](http://hgdownload.cse.ucsc.edu/admin/exe/)
* [cd-hit](http://weizhongli-lab.org/cd-hit/)
* [JBrowse](http://jbrowse.org/install/) - **Only bin needed, does not have to be installed under a web server**

#### Source ######

https://nina_h@bitbucket.org/nina_h/profrep.git

branch "cerit"

#### Configuration #####

Add tools

	<section id="experimental" name="Experimental Tools" >
		<label id="profrep_prepare" text="ProfRep Data Preparation" />
		  <tool file="profrep/extract_data_for_profrep.xml" />
		  <tool file="profrep/db_reducing.xml" />
		<label id="profrep_main" text="Profrep" />
		  <tool file="profrep/profrep.xml" />
		<label id="profrep_supplementary" text="Profrep Supplementary" />
		  <tool file="profrep/profrep_refine.xml" />
		  <tool file="profrep/profrep_masking.xml" />
		  <tool file="profrep/gff_select_region.xml" />
	</section>
	
to 

	$__root_dir__/config/tool_conf.xml
	
------------------------------------------------------------------------

Place PROFREP_DB files to

	$__tool_data_path__/profrep


Place DANTE_DB files to

	$__tool_data_path__/protein_domains

------------------------------------------------------------------------	
	
Create
	
	$__root_dir__/database/dependencies/profrep/1.0.0/env.sh
	
containing:
	
	export JBROWSE_BIN=PATH_TO_JBROWSE_DIR/bin

------------------------------------------------------------------------	
	
Link the following files into galaxy tool-data dir

	ln -s $__tool_directory__/profrep/domains_data/select_domain.txt $__tool_data_path__
	ln -s $__tool_directory__/profrep/profrep_data/prepared_datasets.txt $__tool_data_path__
	
	







