#!/usr/bin/env python3
''' configuration file to set up the paths and constants '''
import os

MAIN_GIT_DIR = os.path.dirname(os.path.realpath(__file__))
TOOL_DATA = os.path.join(MAIN_GIT_DIR, "tool-data")
TMP = "tmp"
SC_MATRIX = os.path.join(TOOL_DATA, "blosum80.txt.sample")
AMBIGUOUS_TAG = "Ambiguous_domain"
## IO
CLASS_FILE = "ALL.classification-new"
LAST_DB_FILE = "ALL_protein-domains_05.fasta"
DOM_PROT_SEQ = "dom_prot_seq.fa"
FILT_DOM_GFF = "domains_filtered.gff"
EXTRACT_DOM_STAT = "domains_counts.txt"
EXTRACT_OUT_DIR = "extracted_domains"
FASTA_LINE = 60
SOURCE_PROFREP = "profrep"
SOURCE_DANTE = "dante"
DOMAINS_FEATURE = "protein_domain"
PHASE = "."
HEADER_GFF = "##gff-version 3"
DOMAINS_GFF = "output_domains.gff"
