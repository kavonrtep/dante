#!/usr/bin/env python3

import time
import configuration
import os
from tempfile import NamedTemporaryFile
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

def filter_qual(OUTPUT_DOMAIN, FILT_DOM_GFF, TH_IDENTITY, TH_LENGTH, TH_FRAMESHIFTS):
	''' Filter gff only based on quality of alignment without domain type considering '''
	with open(OUTPUT_DOMAIN, "r") as gff_all:
		next(gff_all)
		for line in gff_all:
			attributes = line.rstrip().split("\t")[-1]
			al_identity = float(attributes.split(",")[-3].split("=")[1])
			al_length = float(attributes.split(",")[-2].split("=")[1])
			relat_frameshifts = float(attributes.split("\t")[-1].split(",")[-1].split("=")[1])
			dom_type = "-".join([attributes.split(",")[1].split("=")[1].split("/")[0], attributes.split(",")[0].split("=")[1]])
			if al_identity >= TH_IDENTITY and al_length >= TH_LENGTH and relat_frameshifts <= TH_FRAMESHIFTS :
				with open(FILT_DOM_GFF, "a") as gff_filtered:
					gff_filtered.writelines(line)
					
			
def filter_qual_dom(OUTPUT_DOMAIN, FILT_DOM_GFF, TH_IDENTITY, TH_LENGTH, TH_FRAMESHIFTS, SELECTED_DOM):
	''' Filter gff output based on domain and quality of alignment '''
	with open (FILT_DOM_GFF, "a") as gff_filtered:
		with open(OUTPUT_DOMAIN, "r") as gff_all:
			next(gff_all)
			for line in gff_all:
				attributes = line.rstrip().split("\t")[-1]
				al_identity = float(attributes.split(",")[-3].split("=")[1])
				al_length = float(attributes.split(",")[-2].split("=")[1])
				relat_frameshifts = float(attributes.split("\t")[-1].split(",")[-1].split("=")[1])
				dom_type = "-".join([attributes.split(",")[1].split("=")[1].split("/")[0], attributes.split(",")[0].split("=")[1]])
				if al_identity >= TH_IDENTITY and al_length >= TH_LENGTH and relat_frameshifts <= TH_FRAMESHIFTS and dom_type == SELECTED_DOM:
					gff_filtered.writelines(line)
					
	
def get_domains_protseq(FILT_DOM_GFF, DOMAIN_PROT_SEQ):
	''' Get the original nucleic sequence and protein sequence of the reported domains regions '''
	with open(FILT_DOM_GFF, "r") as filt_gff:
		next(filt_gff)
		for line in filt_gff: 
			start = int(line.rstrip().split("\t")[3])
			end = int(line.rstrip().split("\t")[4])
			attributes = line.rstrip().split("\t")[8]
			dom = attributes.split(",")[0].split("=")[1]
			dom_class = "{}/{}".format(attributes.split(",")[1].split("=")[1], attributes.split(",")[2].split("=")[1])
			seq_id = line.rstrip().split("\t")[0]
			prot_seq = line.rstrip().split("\t")[8].split(",")[4].split("=")[1]
			header_prot_seq = ">{}:{}-{} {} {}".format(seq_id, start, end, dom, dom_class)
			with open(DOMAIN_PROT_SEQ, "a") as dom_prot_file:
				dom_prot_file.write("{}\n{}\n".format(header_prot_seq, textwrap.fill(prot_seq, configuration.FASTA_LINE)))

					
def filter_params(reference_seq, alignment_seq, protein_len):
	''' Calculate basic statistics of the quality of alignment '''
	num_ident = 0
	count_frm = 0
	alignment_len = 0
	for i,j in zip(reference_seq, alignment_seq):
		if i==j:
			num_ident += 1
		if j == "/" or j == "\\":
			count_frm += 1
		if i.isalpha():
			alignment_len += 1
	relat_align_len = round(alignment_len/protein_len, 3) 
	align_identity = round(num_ident/len(alignment_seq), 2)
	relat_frameshifts = round(count_frm/(len(alignment_seq)/100),2)
	return align_identity, relat_align_len, relat_frameshifts	


def main(args):
	
	t = time.time()
	
	OUTPUT_DOMAIN = args.domain_gff
	DOMAIN_PROT_SEQ = args.domains_prot_seq
	TH_IDENTITY = args.th_identity
	TH_LENGTH = args.th_length 
	TH_FRAMESHIFTS = args.frameshifts
	FILT_DOM_GFF = args.domains_filtered
	SELECTED_DOM = args.selected_dom
	OUTPUT_DIR = args.output_dir
	
	if not os.path.exists(OUTPUT_DIR) and not os.path.exists(FILT_DOM_GFF):
		os.makedirs(OUTPUT_DIR)
		FILT_DOM_GFF = os.path.join(OUTPUT_DIR, os.path.basename(FILT_DOM_GFF))
		DOMAIN_PROT_SEQ = os.path.join(OUTPUT_DIR, os.path.basename(DOMAIN_PROT_SEQ))
	elif os.path.exists(OUTPUT_DIR) and not os.path.exists(FILT_DOM_GFF):
		if not(os.path.dirname(FILT_DOM_GFF)):
			FILT_DOM_GFF = os.path.join(OUTPUT_DIR, os.path.basename(FILT_DOM_GFF))
		if not(os.path.dirname(DOMAIN_PROT_SEQ)):
			DOMAIN_PROT_SEQ = os.path.join(OUTPUT_DIR, os.path.basename(DOMAIN_PROT_SEQ))
			
	with open (FILT_DOM_GFF, "a") as gff_filtered:
		gff_filtered.write("##gff-version 3")
	
	if SELECTED_DOM != "All":
		filter_qual_dom(OUTPUT_DOMAIN, FILT_DOM_GFF, TH_IDENTITY, TH_LENGTH, TH_FRAMESHIFTS, SELECTED_DOM)
	else:
		filter_qual(OUTPUT_DOMAIN, FILT_DOM_GFF, TH_IDENTITY, TH_LENGTH, TH_FRAMESHIFTS)
	get_domains_protseq(FILT_DOM_GFF, DOMAIN_PROT_SEQ)
	

	print("ELAPSED_TIME_DOMAINS = {} s".format(time.time() - t))
	
if __name__ == "__main__":
	import argparse

	INPUT_DOMAINS_GFF = configuration.INPUT_DOMAINS_GFF
	DOM_PROT_SEQ = configuration.DOM_PROT_SEQ
	FILT_DOM_GFF = configuration.FILT_DOM_GFF

	parser = argparse.ArgumentParser()
	parser.add_argument("-dom_gff", "--domain_gff",type=str, default=INPUT_DOMAINS_GFF,
						help="basic unfiltered gff file of all domains")
	parser.add_argument("-ouf","--domains_filtered",type=str, default=FILT_DOM_GFF,
						help="output filtered domains gff file") 
	parser.add_argument("-dps","--domains_prot_seq",type=str, default=DOM_PROT_SEQ,
						help="output file containg domains protein sequences")
	parser.add_argument("-thl","--th_length",type=float, choices=[Range(0.0, 1.0)],
						default= 0.8, help="proportion of alignment length threshold")
	parser.add_argument("-thi","--th_identity",type=float, choices=[Range(0.0, 1.0)],
						default= 0.35, help="proportion of alignment identity threshold")
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
						"Ty3-CHDCR"
						],
						help="filter output domains based on the domain type")
	parser.add_argument("-dir","--output_dir",type=str, default=configuration.TMP,
						help="specify if you want to change the output directory")
	args = parser.parse_args()
	main(args)
