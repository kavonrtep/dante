# Domain based annotation of transposable elements - DANTE #
Package for domain based annotation of transposable elements in DNA sequences.
## Table of contents
* [Introduction](#introduction)
* [Citation](#citation)
* [DEPENDENCIES](#dependencies)
* [Installation using conda](#installation-using-conda)
* [Usage](#usage)
  * [Annotation of protein domains using `dante`](#annotation-of-protein-domains-using-dante)
  * [Protein Domains Filter](#protein-domains-filter)
  * [Extract Domains Nucleotide Sequences](#extract-domains-nucleotide-sequences)
  * [Identification of protein domains with `dante.py`](#identification-of-protein-domains-with-dantepy)


### Introduction

Main programs in DANTE package are:
* `dante` is a tool for domain based annotation of transposable elements in DNA sequences.
	* Script performs scanning of given DNA sequence(s) in (multi)fasta format in order to discover protein domains using our protein domains database.
	* Domains searching is accomplished engaging LASTAL alignment tool.
	* Domains are subsequently annotated and classified - in case certain domain has multiple annotations assigned, classification is derived from the common classification level of all of them. 	
			
* `dante_gff_output_filtering.py` is a tool for filtering of the domains output from previous step.
	* filters GFF3 output from previous step to obtain certain kind of domain and/or allows to adjust quality filtering  

### Citation
Novak, P., Hostakova, N., Neumann, P., Macas, J. (2024) – DANTE and DANTE_LTR: lineage-centric annotation pipelines for long terminal repeat retrotransposons in plant genomes. NAR Genomics and Bioinformatics 6:113. [https://doi.org/10.1093/nargab/lqae113]        

### DEPENDENCIES
See `requirements.txt`

### Installation using conda ####
[![Anaconda-Server Badge](https://anaconda.org/petrnovak/dante/badges/version.svg)](https://anaconda.org/petrnovak/dante)
```bash
conda install -c conda-forge -c bioconda -c petrnovak dante
```

### Usage ###

#### Annotation of protein domains using `dante` 

Script `dante` is used to annotate protein domains in DNA sequences. It requires a DNA sequence in (multi)fasta format. The script uses LASTAL to align the DNA sequence to the protein domains database. The script outputs a GFF3 file with the annotated protein domains. Note that there is also script `dante.py` which is kept for compatibility reasons - see details below.

To annotate protein domains in a DNA sequence, run the following command:

```bash
# download example data:
wget https://raw.githubusercontent.com/kavonrtep/dante_ltr/main/test_data/sample_genome.fasta
dante -q sample_genome.fasta -o DANTE_output.gff3 -c 10 run dante
```
Output from `dante` is GFF3 file with annotated protein domains. Individual domains are reported per line as regions (start-end) on the original DNA sequence including the seq ID and strand orientation. The last "Attributes" column contains several comma-separated information related to the domain annotation, alignment and its quality. This ouput can be used as input for =dante_ltr= (https://github.com/kavonrtep/dante_ltr) to annotate LTR retrotransposons in the genome. If you want to use protein domains for subsequent analysis, it is recommended that you perform filtering of the domains using `dante_gff_output_filtering.py` script. This script allows you to filter the domains based on their quality and/or extract specific types of protein domains or mobile elements of origin. See **Protein Domains Filter** section for more details.

##### Atributes in GFF3 file #####
- **Name**: The attibute correspond to name of protein 
  domain (RT. RH, PROT, ...)
- **Best_Hit**: Information about the best match of the protein domain to a known database entry.
- **Best_Hit_DB_Pos**: Position of the best hit within the database sequence.
- **DB_Seq**: The database sequence that corresponds to the best hit.
- **Region_Seq**: The sequence of the region where the query sequence was aligned to the database sequence.
- **Query_Seq**: The sequence of the query used to find the best hit.
- **Identity**: The percentage identity of the best hit match.
- **Similarity**: The similarity score of the best hit match.
- **Relat_Length**: The relative length of the match compared to the database sequence.
- **Relat_Interruptions**: Indicates the relative number of interruptions in the 
  domain sequence.Interuption could be either stop codon or frameshift.
- **Hit_to_DB_Length**: The length of the hit compared to the database sequence length.

##### Detailed description of `dante` options

```text

usage: dante [-h] -q QUERY [-D {Viridiplantae_v4.0,Viridiplantae_v3.0,Metazoa_v3.1,Viridiplantae_v2.2,Metazoa_v3.0}] -o DOMAIN_GFF [-dir OUTPUT_DIR] [-M {BL80,BL62,MIQS}] [-thsc THRESHOLD_SCORE] [-wd WIN_DOM] [-od OVERLAP_DOM]
             [-c CPU] [-e EXTRA_DATABASE]

Script performs similarity search on given DNA sequence(s) in (
multi)fasta against our protein domains database of all Transposable element for 
certain group of organisms (Viridiplantae or Metazoans). Domains are subsequently 
annotated and classified - in case certain domain has multiple annotations 
assigned, classification is derived from the common classification level of all of 
them. Domains search is accomplished engaging LASTAL alignment tool.
    

options:
  -h, --help            show this help message and exit
  -D {Viridiplantae_v4.0,Viridiplantae_v3.0,Metazoa_v3.1,Viridiplantae_v2.2,Metazoa_v3.0}, --database {Viridiplantae_v4.0,Viridiplantae_v3.0,Metazoa_v3.1,Viridiplantae_v2.2,Metazoa_v3.0}
                        Select version of RExDB database to use for similarity search (default: Viridiplantae_v4.0)
  -o DOMAIN_GFF, --domain_gff DOMAIN_GFF
                        output domains gff format (default: None)
  -dir OUTPUT_DIR, --output_dir OUTPUT_DIR
                        specify if you want to change the output directory (default: .)
  -M {BL80,BL62,MIQS}, --scoring_matrix {BL80,BL62,MIQS}
                        specify scoring matrix to use for similarity search (BL80, BL62, MIQS) (default: BL80)
  -thsc THRESHOLD_SCORE, --threshold_score THRESHOLD_SCORE
                        percentage of the best score in the cluster to be tolerated when assigning annotations per base (default: 80)
  -wd WIN_DOM, --win_dom WIN_DOM
                        window to process large input sequences sequentially (default: 10000000)
  -od OVERLAP_DOM, --overlap_dom OVERLAP_DOM
                        overlap of sequences in two consecutive windows (default: 10000)
  -c CPU, --cpu CPU     number of threads to use (default: 1)
  -e EXTRA_DATABASE, --extra_database EXTRA_DATABASE
                        extra database to use for similarity search (default: None)

required named arguments:
  -q QUERY, --query QUERY
                        input DNA sequence to search for protein domains in a fasta format. Multifasta format allowed. (default: None)

    Extra database format:
    Extra database is FASTA file with protein domains sequences. This file is appended 
    to selected REXdb database. Header of sequences must contain information about 
    classification compatible with REXdb classification system and also protein domain 
    type. Example of FASTA header:
    
       >MNCI01000001.1:152848-153282 RH Class_I|LTR|Ty3/gypsy|non-chromovirus|OTA|Tat
       CQEALDNIMRELAQVSTVYSSQNDKSFYIYLTISDisissllcQKLDDGVELsvyylsha
       litYET*YIEVEKFFLALVVSFKK*rnylfrshINVICKDKVLRDITTNIYKNSRIA**K
       DILDEFGfhyisqa*TKGQVIATQLT
    
    where:
    
    MNCI01000001.1:152848-153282 is unique identifier of sequence in database.
    
    RH is type of protein domain.
    
    Class_I|LTR|Ty3/gypsy|non-chromovirus|OTA|Tat is classification of protein domain.

```


	
#### Protein Domains Filter
		
The script performs Protein Domains Finder output filtering for quality and/or extracting specific type of protein domain or mobile elements of origin. For the filtered domains it reports their translated protein sequence of original DNA.

When no parameters are given, the script performs quality filtering using the default parameters optimized for Viridiplantae species.

##### INPUTS 
* GFF3 file produced by `dante`.
	
##### Filtering options 
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
		
##### OUTPUTS 
* filtered GFF3 file
* fasta file of translated protein sequences for the aligned domains that match the filtering criteria 
	! as it is taken from the best hit alignment reported by LAST, it does not neccessary cover the whole region reported as domain in GFF
	
##### USAGE 		

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



##### HOW TO RUN EXAMPLE ####
e.g. getting quality filtered integrase(INT) domains of all gypsy transposable elements:
	
	./domains_filtering.py -dom_gff PATH_TO_INPUT_GFF -pdb PATH_TO_PROTEIN_DB -cs PATH_TO_CLASSIFICATION_FILE --selected_dom INT --element_type Ty3/gypsy 


#### Extract Domains Nucleotide Sequences 

This tool extracts nucleotide sequences of protein domains from reference DNA based on DANTE's output. It can be used e.g. for deriving phylogenetic relations of individual mobile elements classes within a species. 

##### INPUTS ####

* original DNA sequence in multifasta format to extract the domains from 
* GFF3 file of protein domains (**DANTE's output** - preferably filtered for quality and specific domain type)
* Domains database classification table (to check the classification level)

##### OUTPUTS ####

* fasta files of domains nucleotide sequences for individual transposons lineages
* txt file of domains counts extracted for individual lineages
* GFF3 file of domains nucleotide sequences for individual transposons lineages


##### USAGE ####	
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

##### HOW TO RUN EXAMPLE ####
	./extract_domains_seqs.py --domains_gff PATH_PROTEIN_DOMAINS_GFF --input_dna PATH_TO_INPUT_DNA  --classification PROTEIN_DOMAINS_DB_CLASS_TBL --extended True

	

#### Identification of protein domains with `dante.py` ###

This tool is kept for compatibility reasons. It performs same function as `dante` but with different interface.

##### INPUTS ####

* DNA sequence [multiFasta]
		
##### OUTPUTS ####
		
* **All protein domains GFF3** - individual domains are reported per line as regions (start-end) on the original DNA sequence including the seq ID and strand orientation. The last "Attributes" column contains several comma-separated information related to the domain annotation, alignment and its quality. This file can undergo further filtering using Protein Domain Filter tool.		

##### USAGE ####

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


		
##### HOW TO RUN EXAMPLE ####
		./protein_domains.py -q PATH_TO_INPUT_SEQ -pdb PATH_TO_PROTEIN_DB -cs PATH_TO_CLASSIFICATION_FILE
		
	 When running for the first time with a new database use -nld option allowing lastal to create indexed database files:

         -nld True

	use other arguments if you wish to rename your outputs or they will be created automatically with standard names 






