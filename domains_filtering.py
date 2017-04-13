#!/usr/bin/env python3

import time
import configuration
import os
import textwrap


class Range():
    '''
    This class is used to check float range in argparse
    '''

    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __eq__(self, other):
        return self.start <= other <= self.end

    def __str__(self):
        return "float range {}..{}".format(self.start, self.end)

    def __repr__(self):
        return "float range {}..{}".format(self.start, self.end)

			
def filter_qual_dom(OUTPUT_DOMAIN, FILT_DOM_GFF, TH_IDENTITY, TH_SIMILARITY, TH_LENGTH, TH_FRAMESHIFTS, SELECTED_DOM):
	''' Filter gff output based on domain and quality of alignment '''
	if SELECTED_DOM != "All":
		selected_dom_type = SELECTED_DOM.split("-")[1]
		element_type = SELECTED_DOM.split("-")[0]
	else:
		selected_dom_type = None
		element_type = None
	with open(OUTPUT_DOMAIN, "r") as gff_all:
		next(gff_all)
		with open (FILT_DOM_GFF, "w") as gff_filtered:
			gff_filtered.write("##gff-version 3\n")
			for line in gff_all:
				attributes = line.rstrip().split("\t")[-1]
				classification = attributes.split(",")[1]
				if classification != "Classification=Ambiguous_domain":
					al_identity = float(attributes.split(",")[-4].split("=")[1])
					al_similarity = float(attributes.split(",")[-3].split("=")[1])
					al_length = float(attributes.split(",")[-2].split("=")[1])
					relat_frameshifts = float(attributes.split("\t")[-1].split(",")[-1].split("=")[1])
					dom_type = attributes.split(",")[0].split("=")[1]
					if al_identity >= TH_IDENTITY and al_similarity >= TH_SIMILARITY and al_length >= TH_LENGTH and relat_frameshifts <= TH_FRAMESHIFTS and ((dom_type == selected_dom_type and element_type in classification) or SELECTED_DOM == "All"):
						gff_filtered.writelines(line)
					
	
def get_domains_protseq(FILT_DOM_GFF, DOMAIN_PROT_SEQ):
	''' Get the translated protein sequence of original DNA seq for all the filtered domains regions '''
	with open(FILT_DOM_GFF, "r") as filt_gff:
		next(filt_gff)
		with open(DOMAIN_PROT_SEQ, "w") as dom_prot_file:
			for line in filt_gff: 
				attributes = line.rstrip().split("\t")[8]
				positions = attributes.split(",")[3].split("=")[1].split(":")[1].split("[")[0]
				dom = attributes.split(",")[0].split("=")[1]
				dom_class = "/".join(attributes.split(",")[3].split("=")[1].split(":")[0].split("/")[1:])
				seq_id = line.rstrip().split("\t")[0]
				prot_seq = line.rstrip().split("\t")[8].split(",")[5].split("=")[1]
				header_prot_seq = ">{}:{} {} {}".format(seq_id, positions, dom, dom_class)
				dom_prot_file.write("{}\n{}\n".format(header_prot_seq, textwrap.fill(prot_seq, configuration.FASTA_LINE)))


def main(args):
	
	t = time.time()
	
	OUTPUT_DOMAIN = args.domain_gff
	DOMAIN_PROT_SEQ = args.domains_prot_seq
	TH_IDENTITY = args.th_identity
	TH_LENGTH = args.th_length 
	TH_FRAMESHIFTS = args.frameshifts
	TH_SIMILARITY = args.th_similarity
	FILT_DOM_GFF = args.domains_filtered
	SELECTED_DOM = args.selected_dom
	OUTPUT_DIR = args.output_dir
	
	if DOMAIN_PROT_SEQ is None:
		DOMAIN_PROT_SEQ = configuration.DOM_PROT_SEQ
	if FILT_DOM_GFF is None:
		FILT_DOM_GFF = configuration.FILT_DOM_GFF
	

	if OUTPUT_DIR is None:
		OUTPUT_DIR = os.path.dirname(os.path.abspath(OUTPUT_DOMAIN))
	if not os.path.exists(OUTPUT_DIR):
		os.makedirs(OUTPUT_DIR)		
	FILT_DOM_GFF = os.path.join(OUTPUT_DIR, os.path.basename(FILT_DOM_GFF))
	DOMAIN_PROT_SEQ = os.path.join(OUTPUT_DIR, os.path.basename(DOMAIN_PROT_SEQ))
	
	#if SELECTED_DOM is not "All":
		#selected_dom_type = SELECTED_DOM.split("-")[1]
		#element_type = SELECTED_DOM.split("-")[0]
	#else:
		#selected_dom_type = None
		#element_type = None
	filter_qual_dom(OUTPUT_DOMAIN, FILT_DOM_GFF, TH_IDENTITY, TH_SIMILARITY, TH_LENGTH, TH_FRAMESHIFTS, SELECTED_DOM)
	get_domains_protseq(FILT_DOM_GFF, DOMAIN_PROT_SEQ)
	

	print("ELAPSED_TIME_DOMAINS = {} s".format(time.time() - t))
	
if __name__ == "__main__":
	import argparse
	from argparse import RawDescriptionHelpFormatter
	
	class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
		pass
    
	
	
	parser = argparse.ArgumentParser(
		description='''Script performs filtering of gff output which is result of protein_domains_pd.py and contains all types of domains 
		without quality filtering. The script enables to obtain results for specific kind of domain separately for individual types of repetitive elements and/or filter out domains that do not reach appropriate length, similarity or have more frameshifts per 100 bp than set by threshold. Records for ambiguous domain type (e.g. INT/RH) are filtered out automatically. Based on filtered gff file protein sequences are reported in separate file - these translations of original DNA sequence are taken from the LASTAL output . NOTE -  foundscanning of given DNA sequence(s) in (multi)fasta format in order to discover protein domains using our protein domains database. Domains are subsequently annotated and classified - in case certain domain has multiple annotations assigned, classifation is derived from the common classification level of all of them. Domains searching is accomplished engaging LASTAL search tool.
		
	DEPENDANCIES:
		- python 3.4 or higher
		- configuration.py module

	EXAMPLE OF USAGE:
		
		./protein_domains_pd.py -q PATH_TO_INPUT_SEQ -pdb PATH_TO_PROTEIN_DB -cs PATH_TO_CLASSIFICATION_FILE

	
		''',
		epilog="""""",
		formatter_class=CustomFormatter)
	requiredNamed = parser.add_argument_group('required named arguments')
	requiredNamed.add_argument("-dom_gff", "--domain_gff",type=str, required=True,
						help="basic unfiltered gff file of all domains")
	parser.add_argument("-ouf","--domains_filtered",type=str, 
						help="output filtered domains gff file") 
	parser.add_argument("-dps","--domains_prot_seq",type=str, 
						help="output file containg domains protein sequences")
	parser.add_argument("-thl","--th_length",type=float, choices=[Range(0.0, 1.0)],
						default= 0.8, help="proportion of alignment length threshold")
	parser.add_argument("-thi","--th_identity",type=float, choices=[Range(0.0, 1.0)],
						default= 0.35, help="proportion of alignment identity threshold")
	parser.add_argument("-ths","--th_similarity",type=float, choices=[Range(0.0, 1.0)],
						default= 0, help="threshold for alignment proportional similarity")
	parser.add_argument("-fr","--frameshifts",type=int, default=1,
						help="frameshifts tolerance threshold per 100 bp")
	parser.add_argument("-sd","--selected_dom",type=str, default="All", choices=[
						"All",
						"Ty1-GAG",
						"Ty1-INT",
						"Ty1-PROT",
						"Ty1-RH",
						"Ty1-RT",
						"Ty3-GAG",
						"Ty3-INT",
						"Ty3-PROT",
						"Ty3-RT",
						"Ty3-gRH",
						"Ty3-aRH",
						"Ty3-CHDII",
						"Ty3-CHDCR",
						"CACTA-TPase",
						"DIRS-RH",
						"DIRS-RT",
						"DIRS-YR",
						"Harbinger-TPase",
						"hAT-TPase",
						"Helitron-HEL1",
						"Helitron-HEL2",
						"Kolobok-TPase",
						"LINE-ENDO",
						"LINE-RH",
						"LINE-RT",
						"Mariner-TPase",
						"Merlin-TPase",
						"MuDR-TPase",
						"Novosib-TPase",
						"PARA-PROT",
						"PARA-RH",
						"PARA-RT",
						"Penelope-RT",
						"PiggyBac-TPase",
						"P-TPase",
						"Sola1-TPase",
						"Sola2-TPase"
						],
						help="filter output domains based on the domain type")
	parser.add_argument("-dir","--output_dir",type=str, default=None,
						help="specify if you want to change the output directory")
	args = parser.parse_args()
	main(args)
