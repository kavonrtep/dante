#!/usr/bin/env python3
import os
""" configuration file to set up the paths and constants """

# PATHS
MAIN_GIT_DIR = os.path.dirname(os.path.realpath(__file__)) 
######## PROFREP
JBROWSE_BIN = "/home/ninah/programs/JBrowse-1.12.3/bin"
jbrowse_data_dir = "data"
jbrowse_link_to_galaxy = "data/galaxy_profrep_files"
jbrowse_link_to_profrep = "data/profrep_data"
LINK_PART1 = "http://nina/JBrowse-1.12.1/index.html?data="
LINK_PART1_PC = "http://nina/JBrowse-1.12.1/index.html?data="
PROFREP_DATA = os.path.join(MAIN_GIT_DIR, "profrep_data")
#TBL = os.path.join(PROFREP_DATA, "prepared_datasets.txt")
######### DOMAINS
MAIN_RESULTS_DIR = os.getcwd()
DOMAINS_DATA = os.path.join(MAIN_GIT_DIR, "domains_data")
TMP = os.path.join(MAIN_RESULTS_DIR, "tmp")
SC_MATRIX = os.path.join(DOMAINS_DATA, "blosum80.txt")
# FILE_NAMES
CLASS_FILE = "ALL.classification-new"
LAST_DB_FILE = "ALL_protein-domains_05.fasta"

# CONSTANTS 
N_segment = 50
MAX_PIC_NUM = 50
IMAGE_RES = 300
FASTA_LINE = 60 
COLORS_HEX = ["#7F7F7F", "#00FF00", "#0000FF", "#FF0000", "#01FFFE", "#FFA6FE", "#FFDB66", "#006401", "#010067", "#95003A", "#007DB5", "#FF00F6", "#FFEEE8", "#774D00", "#90FB92", "#0076FF", "#D5FF00", "#FF937E", "#6A826C", "#FF029D", "#FE8900", "#7A4782", "#7E2DD2", "#85A900", "#FF0056", "#A42400", "#00AE7E", "#683D3B", "#BDC6FF", "#263400", "#BDD393", "#00B917", "#9E008E", "#001544", "#C28C9F", "#FF74A3", "#01D0FF", "#004754", "#E56FFE", "#788231", "#0E4CA1", "#91D0CB", "#BE9970", "#968AE8", "#BB8800", "#43002C", "#DEFF74", "#00FFC6", "#FFE502", "#620E00", "#008F9C", "#98FF52", "#7544B1", "#B500FF", "#00FF78", "#FF6E41", "#005F39", "#6B6882", "#5FAD4E", "#A75740", "#A5FFD2", "#FFB167", "#009BFF", "#E85EBE" ]
COLORS_RGB = ["127,127,127","0,255,0","0,0,255","255,0,0","1,255,254","255,166,254","255,219,102","0,100,1","1,0,103","149,0,58","0,125,181","255,0,246","255,238,232","119,77,0","144,251,146","0,118,255","213,255,0","255,147,126","106,130,108","255,2,157","254,137,0","122,71,130","126,45,210","133,169,0","255,0,86","164,36,0","0,174,126","104,61,59","189,198,255","38,52,0","189,211,147","0,185,23","158,0,142","0,21,68","194,140,159","255,116,163","1,208,255","0,71,84","229,111,254","120,130,49","14,76,161","145,208,203","190,153,112","150,138,232","187,136,0","67,0,44","222,255,116","0,255,198","255,229,2","98,14,0","0,143,156","152,255,82","117,68,177","181,0,255","0,255,120","255,110,65","0,95,57","107,104,130","95,173,78","167,87,64","165,255,210","255,177,103","0,155,255","232,94,190"]
s_info_header = "seq_id\tseq_legth\tseq_count\tfile_start_pos\tfile_end_pos\n"
html_str_domains = '''
	<!DOCTYPE html>
	<html>
	<body>
		<h2>PROTEIN DOMAINS VISUALIZATION</h2>
		<a href="{}" target="_blank" >Link to JBrowse</a> </br>
		<hr>
	</body>
	</html>
	'''

#GFF tracks constants
HEADER_GFF = "##gff-version 3"
SOURCE = "profrep"
PHASE = "."
DOMAINS_FEATURE = "protein_domain"
REPEATS_FEATURE = "repeat"
R_SCORE = "."
R_STRAND = "." 
R_PHASE = "."
N_NAME = "N"
N_FEATURE = "N_region"
N_SCORE = "."
N_STRAND = "."
N_PHASE = "."

# JBROWSE TRACKS CONF
JSON_CONF_R = """{"hooks" : {"modify": "function( track, f, fdiv ) {fdiv.style.backgroundColor = '#278ECF'}"}}"""
JSON_CONF_N = """{"hooks" : {"modify": "function( track, f, fdiv ) {fdiv.style.background = '#474747'}"}}"""

# OUPUTS
ADJUSTED_GFF = "adjusted_domains.gff"
#DOMAINS_GFF = "{}/{}".format(TMP,"output_domains.gff")
DOMAINS_GFF = "output_domains.gff"
INPUT_DOMAINS_GFF = os.path.join(TMP, "output_domains.gff")
HTML = os.path.join(TMP,"output.html")
#HTML = "output.html"
HTML_DOM = "dom_vis.html"
REPEATS_GFF = os.path.join(TMP,"output_repeats.gff")
#REPEATS_GFF = "output_repeats.gff"
#N_REG = "N_regions.gff"
N_REG = os.path.join(TMP,"N_regions.gff")
REPEATS_TABLE = os.path.join(TMP,"output_table.csv")
#REPEATS_TABLE = "output_table.csv"
SEQ_INFO = "seq_info.csv"
#DOM_PROT_SEQ = os.path.join(TMP,"dom_prot_seq.txt") 
DOM_PROT_SEQ = "dom_prot_seq.txt" 
#FILT_DOM_GFF = os.path.join(TMP,"domains_filtered.gff")
FILT_DOM_GFF = "domains_filtered.gff"
#DOM_SUMMARY = os.path.join(TMP,"dom_summary.txt")
DOM_SUMMARY = "dom_summary.txt"
#ELEM_TABLE = os.path.join(TMP,"elements_table.txt")
ELEM_TABLE = "elements_table.txt"



