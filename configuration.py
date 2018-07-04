#!/usr/bin/env python3
import os
''' configuration file to set up the paths and constants '''

######## PROFREP #######################################################
## Constansts 
N_segment = 50
MAX_FILES_SUBPROFILES = 1000
MAX_PIC_NUM = 50
IMAGE_RES = 300
FASTA_LINE = 60
SEQ_LEN_VIZ = 200000
FORBIDDEN_CHARS = "\\/"
HTML_STR = '''
	<!DOCTYPE html>
	<html>
	<body>
		<h2>PROFREP OUTPUT</h2>
		<h4> Sequences processed: </h4>
		{}
		<h4> Total length: </h4>
		<pre> {} bp </pre>
		<h4> Database: </h4>
		<pre> {} </pre>
		<hr>
		<h3> Repetitive profile(s)</h3> </br>
		{} <br/>
		<h4>References: </h4>
		{}
		</h6>
	</body>
	</html>
	'''
## IO
DOMAINS_GFF = "output_domains.gff"
N_GFF = "N_regions.gff"
REPEATS_GFF = "output_repeats.gff"
HTML = "output.html"
LOG_FILE = "log.txt"
PROFREP_DATA = "profrep_data"
PROFREP_TBL = "prepared_datasets.txt"
PROFREP_OUTPUT_DIR = "profrep_output_dir"
## JBrowse and Tracks Conf
jbrowse_data_dir = "data"
JSON_CONF_R = """{"hooks" : {"modify": "function( track, f, fdiv ) {fdiv.style.backgroundColor = '#278ECF'}"}}"""
JSON_CONF_N = """{"hooks" : {"modify": "function( track, f, fdiv ) {fdiv.style.background = '#474747'}"}}"""
COLORS_HEX = ["#7F7F7F", "#00FF00", "#0000FF", "#FF0000", "#01FFFE", "#FFA6FE", "#FFDB66", "#006401", "#010067", "#95003A", "#007DB5", "#FF00F6", "#774D00", "#90FB92", "#0076FF", "#D5FF00", "#FF937E", "#6A826C", "#FF029D", "#FE8900", "#7A4782", "#7E2DD2", "#85A900", "#FF0056", "#A42400", "#00AE7E", "#683D3B", "#BDC6FF", "#263400", "#BDD393", "#00B917", "#9E008E", "#001544", "#C28C9F", "#FF74A3", "#01D0FF", "#004754", "#E56FFE", "#788231", "#0E4CA1", "#91D0CB", "#BE9970", "#968AE8", "#BB8800", "#43002C", "#DEFF74", "#00FFC6", "#FFE502", "#620E00", "#008F9C", "#98FF52", "#7544B1", "#B500FF", "#00FF78", "#FF6E41", "#005F39", "#6B6882", "#5FAD4E", "#A75740", "#A5FFD2", "#FFB167", "#009BFF", "#E85EBE" ]
COLORS_RGB = ["127,127,127", "0,255,0","0,0,255","255,0,0","1,255,254","255,166,254","255,219,102","0,100,1","1,0,103","149,0,58","0,125,181","255,0,246","119,77,0","144,251,146","0,118,255","213,255,0","255,147,126","106,130,108","255,2,157","254,137,0","122,71,130","126,45,210","133,169,0","255,0,86","164,36,0","0,174,126","104,61,59","189,198,255","38,52,0","189,211,147","0,185,23","158,0,142","0,21,68","194,140,159","255,116,163","1,208,255","0,71,84","229,111,254","120,130,49","14,76,161","145,208,203","190,153,112","150,138,232","187,136,0","67,0,44","222,255,116","0,255,198","255,229,2","98,14,0","0,143,156","152,255,82","117,68,177","181,0,255","0,255,120","255,110,65","0,95,57","107,104,130","95,173,78","167,87,64","165,255,210","255,177,103","0,155,255","232,94,190"]
TRACK_LIST = '''
	\t,{}\n
	\t"storeClass" : "JBrowse/Store/SeqFeature/BigWig",
	\t"urlTemplate" : "{}",
	\t"type" : "JBrowse/View/Track/Wiggle/XYPlot",
	\t"label" : "{}",
	\t"key" : "{}",
	\t"style": {}
	\t\t"pos_color": "{}"
	\t {},
	\t"scale" : "log"
	\t{}\n
	'''
## GFF tracks
HEADER_GFF = "##gff-version 3"
SOURCE_PROFREP = "profrep"
SOURCE_DANTE = "dante"
PHASE = "."
DOMAINS_FEATURE = "protein_domain"
REPEATS_FEATURE = "repeat"
N_NAME = "N"
N_FEATURE = "N_region"
HEADER_WIG = "variableStep\tchrom=" 
GFF_EMPTY = "."


######### BIG WIG ######################################################
CHROM_SIZES_FILE = "chrom_sizes.txt"

######### EXTRACT_DATA_DOR_PROFREP #####################################
HITSORT_CLS = "seqclust/clustering/hitsort.cls"
READS_ALL = "seqclust/sequences/sequences.fasta"
ANNOTATION = "PROFREP_CLASSIFICATION_TEMPLATE.csv"

######### PROFREP_DB_REDUCING ##########################################
MEM_LIM = 1500 # MB
CLS_REDUCED = "hitsort_reduced.cls"
READS_ALL_REDUCED = "reads_all_reduced"

######### PROFREP_REFINING #############################################
WITH_DOMAINS = "mobile_element"
QUALITY_DIFF_TO_REMOVE = 0.05 # 5% tolerance of PID

######### DANTE ##############################################
MAIN_GIT_DIR = os.path.dirname(os.path.realpath(__file__)) 
DOMAINS_DATA = os.path.join(MAIN_GIT_DIR, "domains_data")
TMP = "tmp"
SC_MATRIX = os.path.join(DOMAINS_DATA, "blosum80.txt")
AMBIGUOUS_TAG = "Ambiguous_domain"
## IO
CLASS_FILE = "ALL.classification-new"
LAST_DB_FILE = "ALL_protein-domains_05.fasta"
DOM_PROT_SEQ = "dom_prot_seq.fa" 
FILT_DOM_GFF = "domains_filtered.gff"
EXTRACT_DOM_STAT = "domains_counts.txt"
EXTRACT_OUT_DIR = "extracted_domains"











