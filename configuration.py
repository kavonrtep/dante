#!/usr/bin/env python3
import os
""" configuration file to set up the paths """

# PATHS
JBROWSE_BIN = "/home/galaxy/bin/JBrowse-1.12.1/bin"
JBROWSE_BIN_PC = "/var/www/html/JBrowse-1.12.1/bin"
jbrowse_data_dir = "jbrowse"
jbrowse_link_to_galaxy = "data/galaxy_files"
LINK_PART1 = "http://nod6/JBrowse-1.12.1/index.html?data="
LINK_PART1_PC = "http://nina/JBrowse-1.12.1/index.html?data="
LAST_DB = "proteins_all"
MAIN_DIR = os.path.dirname(os.path.realpath(__file__))
DATA = "{}/{}".format(MAIN_DIR, "data")
TMP = "{}/{}".format(MAIN_DIR, "tmp")

# CONSTANTS
MAX_PIC_NUM = 5
COLORS = ["#00FF00", "#0000FF", "#FF0000", "#01FFFE", "#FFA6FE", "#FFDB66", "#006401", "#010067", "#95003A", "#007DB5", "#FF00F6", "#FFEEE8", "#774D00", "#90FB92", "#0076FF", "#D5FF00", "#FF937E", "#6A826C", "#FF029D", "#FE8900", "#7A4782", "#7E2DD2", "#85A900", "#FF0056", "#A42400", "#00AE7E", "#683D3B", "#BDC6FF", "#263400", "#BDD393", "#00B917", "#9E008E", "#001544", "#C28C9F", "#FF74A3", "#01D0FF", "#004754", "#E56FFE", "#788231", "#0E4CA1", "#91D0CB", "#BE9970", "#968AE8", "#BB8800", "#43002C", "#DEFF74", "#00FFC6", "#FFE502", "#620E00", "#008F9C", "#98FF52", "#7544B1", "#B500FF", "#00FF78", "#FF6E41", "#005F39", "#6B6882", "#5FAD4E", "#A75740", "#A5FFD2", "#FFB167", "#009BFF", "#E85EBE" ]


# FILES
CLASSIFICATION = "{}/{}".format(DATA, "classification.csv")
LAST_DB = "{}/{}".format(DATA, "proteins_all")


# OUPUTS

DOMAINS_GFF = "{}/{}".format(TMP,"output_domains.gff")
#DOMAINS_PIC = "{}/{}".format(TMP,"output_domains.png")
HTML = "{}/{}".format(TMP,"output.html")
REPEATS_GFF = "{}/{}".format(TMP,"output_repeats.gff")
REPEATS_TABLE = "{}/{}".format(TMP,"output_table.csv")
SEQ_INFO = "{}/{}".format(TMP,"seq_info.csv")
