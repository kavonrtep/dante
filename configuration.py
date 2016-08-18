#!/usr/bin/env python3
import os
""" configuration file to set up the paths """

# PATHS
JBROWSE_BIN="/home/galaxy/bin/JBrowse-1.12.1/bin"
JBROWSE_BIN_PC = "/var/www/JBrowse-1.12.1/bin"
jbrowse_data_dir = "jbrowse"
jbrowse_link_to_galaxy = "data/galaxy_files"
LINK_PART1 = "http://nod6/JBrowse-1.12.1/index.html?data="
LAST_DB = "proteins_all"
MAIN_DIR = os.path.dirname(os.path.realpath(__file__))
DATA = "{}/{}".format(MAIN_DIR, "data")
TMP = "{}/{}".format(MAIN_DIR, "tmp")


# FILES
CLASSIFICATION = "{}/{}".format(DATA, "classification.csv")
LAST_DB = "{}/{}".format(DATA, "proteins_all")


# OUPUTS

DOMAINS_GFF = "{}/{}".format(TMP,"output_domains.gff")
#DOMAINS_PIC = "{}/{}".format(TMP,"output_domains.png")
HTML = "{}/{}".format(TMP,"output.html")
REPEATS_GFF = "{}/{}".format(TMP,"output_repeats.gff")
REPEATS_TABLE = "{}/{}".format(TMP,"output_table.csv")
