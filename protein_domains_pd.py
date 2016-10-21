#!/usr/bin/env python3

import numpy as np
import subprocess
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import sys
import csv
import time
from operator import itemgetter
import os
import configuration
from tempfile import NamedTemporaryFile


def domain_annotation(element, CLASSIFICATION):
	""" assign protein domain to each hit from protein database  """
	domain = []
	rep_type = []
	rep_lineage = []
	with open(CLASSIFICATION, "r") as cl_tbl:
		annotation = {}
		header_classification = cl_tbl.readline().strip().split("\t")
		for line in cl_tbl:
			record = line.rstrip().split("\t")
			annotation[record[0]] = [record[1],record[2]]		
	for i in range(len(element)):
		domain.append(element[i].split("__")[0].split("-")[1])
		element_name = element[i].split("__")[1]
		if element_name in annotation.keys():			
			rep_type.append(annotation[element_name][0])
			rep_lineage.append(annotation[element_name][1])
		else:
			rep_type.append(" ")
			rep_lineage.append(" ")
	return domain, rep_type, rep_lineage
	

def hits_processing(sequence_hits):
	""" gain hits intervals separately for forward and reverse strand """
	seq_length = sequence_hits[0,5]
	reverse_strand_idx = np.where(sequence_hits[:,4] == "-")[0]
	if not reverse_strand_idx.any():
		start_pos_plus = sequence_hits[:,2]
		end_pos_plus = sequence_hits[:,3]
		regions_plus = list(zip(start_pos_plus, end_pos_plus))
		regions_minus = []
	else:
		reverse_strand_idx = reverse_strand_idx[0]
		start_pos_plus = sequence_hits[0:reverse_strand_idx,2]
		end_pos_plus = sequence_hits[0:reverse_strand_idx,3]
		start_pos_minus = seq_length - sequence_hits[reverse_strand_idx:,3]
		end_pos_minus = seq_length - sequence_hits[reverse_strand_idx:,2]
		regions_plus= list(zip(start_pos_plus, end_pos_plus))
		regions_minus = list(zip(start_pos_minus, end_pos_minus))
	return reverse_strand_idx, regions_plus, regions_minus, seq_length


def overlapping_regions(input_data):
	""" join all overalapping intervals """
	if input_data: 
		sorted_idx, sorted_data = zip(*sorted([(index,data) for index,data in enumerate(input_data)], key=itemgetter(1)))
		merged_ends = input_data[sorted_idx[0]][1]
		intervals = []
		output_intervals = [] 
		for i in sorted_idx:
			if input_data[i][0] < merged_ends:
				merged_ends = max(input_data[i][1], merged_ends)
				intervals.append(i)
			else:
				output_intervals.append(intervals)
				intervals = []
				intervals.append(i)
				merged_ends = input_data[i][1]		
		output_intervals.append(intervals)
	else:
		output_intervals = []
	return output_intervals


def best_score(scores, indices):
	""" from overlapping intervals take the one with the highest score """
	best_scores = []
	best_idx = []
	for idx in indices:
		# if more hits have the same best score take only the first one
		best_idx.append(idx[np.where(scores[idx] == max(scores[idx]))[0][0]])
	return best_idx


def create_gff(sequence_hits, best_idx, seq_id, regions, OUTPUT_DOMAIN, DOMAINS_PROT_SEQ, CLASSIFICATION):
	""" track format of domains found in query sequence(s)
	- custom atributes reported:
	Domain Name, Repet.type, Repet.lineage, 
	Original protein sequence, Aligned sequence 
	"""
	t2 = time.time()
	
	SOURCE = "profrep"
	FEATURE = "protein_domain"
	PHASE = "."
	xminimal = []
	xmaximal = []
	scores = []
	strands = []
	domains = []
	count = 0
	[domain, rep_type, rep_lineage] = domain_annotation(sequence_hits[:,6][best_idx], CLASSIFICATION)
	for i in best_idx:
		alignment_start = regions[i][0]
		xminimal.append(alignment_start)
		alignment_end = regions[i][1]
		xmaximal.append(alignment_end)
		strand = sequence_hits[i,4]
		strands.append(strand)
		score = sequence_hits[i,0]
		scores.append(score)
		sequence = sequence_hits[i,7]
		alignment = sequence_hits[i,8]
		dom_annotation = list(map(lambda rtype, rlin: "{}/{}".format(rtype, rlin), rep_type, rep_lineage))
		with open(OUTPUT_DOMAIN, "a") as gff:	
			gff.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tName={},Rep_type={},Rep_lineage={},Sequence={},Alignment={}\n".format(seq_id, SOURCE, FEATURE, alignment_start, alignment_end, score, strand, PHASE, domain[count], rep_type[count], rep_lineage[count], sequence, alignment))	
		if DOMAINS_PROT_SEQ:
			with open(DOMAINS_PROT_SEQ, "a") as domains_prot_file:
				header_dom_seq = ">{}:{}-{} {} {}".format(seq_id, alignment_start, alignment_end, domain[count], dom_annotation[count])
				domains_prot_file.write(("{}\n{}\n".format(header_dom_seq, sequence)))
		count += 1	
	return xminimal, xmaximal, scores, strands, domain, dom_annotation

def multifasta(QUERY):
	''' Create single fasta temporary files to be processed sequentially '''
	PATTERN = ">"
	fasta_list = []
	with open(QUERY, "r") as fasta:
		reader = fasta.read()
		splitter = reader.split(PATTERN)[1:]
		if len(splitter) > 1:
			for fasta_num, part in enumerate(splitter):
				ntf = NamedTemporaryFile(delete=False)
				ntf.write("{}{}".format(PATTERN, part).encode("utf-8"))
				fasta_list.append(ntf.name)
				ntf.close()
			return fasta_list
		else:
			fasta_list.append(QUERY)
			return fasta_list


def fasta_read(subfasta):
	''' Read fasta, gain header and sequence without gaps '''
	sequence_lines = []
	with open(subfasta, "r") as fasta:
		header = fasta.readline().strip().split(" ")[0][1:]
		for line in fasta:
			clean_line = line.strip()			
			if clean_line:				
				sequence_lines.append(clean_line)
	sequence = "".join(sequence_lines)
	return header, sequence
	
def cut_domains_seq(QUERY, DOMAINS_SEQ, DOMAINS_PROT_SEQ, xminimal, xmaximal, domains_all, seq_ids, dom_annotation_all):
	fasta_list = multifasta(QUERY)
	for subfasta in fasta_list:
		[header, sequence] = fasta_read(subfasta)
		if header in seq_ids:
			seq_idx = seq_ids.index(header)
			for count in list(range(len(domains_all[seq_idx]))):
				seq_cut = sequence[xminimal[seq_idx][count] - 1 : xmaximal[seq_idx][count]]
				#seq_prot =alignment[seq_idx][count]
				dom = domains_all[seq_idx][count]
				dom_class = dom_annotation_all[seq_idx][count]
				header_dom_seq = ">{}:{}-{} {} {}".format(seq_ids[seq_idx], xminimal[seq_idx][count], xmaximal[seq_idx][count], dom, dom_class)
				with open(DOMAINS_SEQ, "a") as dom_seq_file:
					dom_seq_file.write("{}\n{}\n".format(header_dom_seq, seq_cut))		 
				#with open(DOMAINS_SEQ, "a") as dom_prot_file:
					#dom_prot_file.write("{}\n{}\n".format(header_dom_seq, alignment))		 



#def domains_protein_seq(DOMAINS_PROT_SEQ, ):
	#with open(DOMAINS_PROT_SEQ, "a") as dom_seq_file:
			#header_prot_seq = ">{}:{}-{} {} {}".format(seq_ids[seq_idx], xminimal[seq_idx][count], xmaximal[seq_idx][count], dom, dom_class)
			#dom_seq_file.write("{}\n{}\n".format(header_dom_seq, seq_cut))	

def get_identity(reference_seq, alignment_seq):
	num_ident = 0
	for i,j in zip(reference_seq, alignment_seq):
		if i==j:
			num_ident += 1
	align_len = len(reference_seq)
	return num_ident/align_len * 100, align_len	


def domain_search(QUERY, LAST_DB, CLASSIFICATION, OUTPUT_DOMAIN, DOMAIN_PROT_SEQ, TH_IDENTITY, TH_LENGTH):		
	""" search for protein domains using our protein database and external tool LAST,
	stdout is parsed in real time and hits for a single sequence undergo further processing
	- tabular format(TAB) to get info about position, score, orientation 
	- MAF format to gain alignment and original sequence
	"""	
	seq_ids = []
	xminimal_all = []
	xmaximal_all = []
	domains_all = []
	dom_annotation_all = []
	header_gff = "##gff-version 3"
	sequence_hits = np.empty((0,9))
	with open(QUERY, "r") as fasta:
		seq_id = fasta.readline().strip().split(" ")[0][1:]
	with open(OUTPUT_DOMAIN, "a") as gff:
		gff.write("{}\n".format(header_gff))
	tab = subprocess.Popen("lastal -F15 {} {} -f TAB ".format(LAST_DB, QUERY), stdout=subprocess.PIPE, shell=True)
	maf = subprocess.Popen("lastal -F15 {} {} -f MAF ".format(LAST_DB, QUERY), stdout=subprocess.PIPE, shell=True)
	maf.stdout.readline()
	for line_tab in tab.stdout:
		line_tab = line_tab.decode("utf-8")
		if not line_tab.startswith('#'):
			line_maf = [maf.stdout.readline() for line_count in range(4)]
			reference_seq = line_maf[1].decode("utf-8").rstrip().split(" ")[-1]
			alignment_seq = line_maf[2].decode("utf-8").rstrip().split(" ")[-1]
			#############################################################
			[percent_ident, align_len] = get_identity(reference_seq, alignment_seq)
			if percent_ident >= TH_IDENTITY and align_len >= TH_LENGTH:
			#############################################################
				line = line_tab.rstrip().split("\t")
				line_maf = []
				element_name = line[1]
				if np.all(sequence_hits==0):
					seq_id = line[6]
					seq_ids.append(seq_id)
				if line[6] != seq_id: 
					[reverse_strand_idx, regions_plus, regions_minus, seq_length] = hits_processing(sequence_hits)
					if reverse_strand_idx == []:
						positions = overlapping_regions(regions_plus)
						best_idx = best_score(sequence_hits[:,0], positions)
						[xminimal, xmaximal, scores, strands, domain, dom_annotation] = create_gff(sequence_hits, best_idx, seq_id, regions_plus, OUTPUT_DOMAIN, DOMAIN_PROT_SEQ, CLASSIFICATION)
					else:
						positions_plus = overlapping_regions(regions_plus)
						positions_minus = overlapping_regions(regions_minus)
						regions = regions_plus + regions_minus
						positions = positions_plus + [x + reverse_strand_idx for x in positions_minus]
						best_idx = best_score(sequence_hits[:,0], positions)
						[xminimal, xmaximal, scores, strands, domain, dom_annotation] = create_gff(sequence_hits, best_idx, seq_id, regions, OUTPUT_DOMAIN, DOMAIN_PROT_SEQ, CLASSIFICATION)
					sequence_hits = np.empty((0,9))
					seq_id = line[6]
					seq_ids.append(seq_id)
					xminimal_all.append(xminimal)
					xmaximal_all.append(xmaximal)
					domains_all.append(domain)
					dom_annotation_all.append(dom_annotation)
				line_parsed = np.array([int(line[0]), seq_id, int(line[7]), int(line[7]) + int(line[8]), line[9], int(line[10]), element_name, reference_seq, alignment_seq], dtype=object)
				sequence_hits = np.append(sequence_hits, [line_parsed], axis=0)
		else:
			maf.stdout.readline()
	if not np.all(sequence_hits==0):	
		[reverse_strand_idx, regions_plus, regions_minus, seq_length] = hits_processing(sequence_hits)
		if reverse_strand_idx == []:
			positions = overlapping_regions(regions_plus)
			best_idx = best_score(sequence_hits[:,0], positions)
			[xminimal, xmaximal, scores, strands, domain, dom_annotation] = create_gff(sequence_hits, best_idx, seq_id, regions_plus, OUTPUT_DOMAIN, DOMAIN_PROT_SEQ, CLASSIFICATION)
		else:
			positions_plus = overlapping_regions(regions_plus)
			positions_minus = overlapping_regions(regions_minus)
			regions = regions_plus + regions_minus
			positions = positions_plus + [x + reverse_strand_idx for x in positions_minus]
			best_idx = best_score(sequence_hits[:,0], positions)
			[xminimal, xmaximal, scores, strands, domain, dom_annotation] = create_gff(sequence_hits, best_idx, seq_id, regions, OUTPUT_DOMAIN, DOMAIN_PROT_SEQ, CLASSIFICATION)
		xminimal_all.append(xminimal)
		xmaximal_all.append(xmaximal)
		domains_all.append(domain)
		dom_annotation_all.append(dom_annotation)
	
	return xminimal_all, xmaximal_all, domains_all, seq_ids, dom_annotation_all
	
	
def main(args):
	QUERY = args.query
	LAST_DB = args.protein_database
	CLASSIFICATION = args.classification
	OUTPUT_DOMAIN = args.domain_gff
	NEW_PDB = args.new_pdb
	DOMAINS_SEQ = args.domains_seq
	DOMAIN_PROT_SEQ = args.domains_prot_seq
	TH_IDENTITY = args.th_identity
	TH_LENGTH = args.th_length 
	
	
	if NEW_PDB:
		subprocess.call("lastdb -p -cR01 {} {}".format(LAST_DB, LAST_DB), shell=True)
	
	if not os.path.exists(LAST_DB):
		CLASSIFICATION = os.path.join(configuration.PTOOL_DATA_DIR, CLASSIFICATION)
		LAST_DB = os.path.join(configuration.TOOL_DATA_DIR, LAST_DB) 
	
	[xminimal, xmaximal, domains_all, seq_ids, dom_annotation_all] = domain_search(QUERY, LAST_DB, CLASSIFICATION, OUTPUT_DOMAIN, DOMAIN_PROT_SEQ, TH_IDENTITY, TH_LENGTH)
	cut_domains_seq(QUERY, DOMAINS_SEQ, DOMAIN_PROT_SEQ, xminimal, xmaximal, domains_all, seq_ids, dom_annotation_all)

if __name__ == "__main__":
	import argparse
	
	LAST_DB = configuration.LAST_DB
	CLASSIFICATION = configuration.CLASSIFICATION
	DOMAINS_GFF = configuration.DOMAINS_GFF
	DOM_SEQ = "tmp/dom_seq.txt"
	DOM_PROT_SEQ = "tmp/dom_prot_seq.txt"

	parser = argparse.ArgumentParser()
	parser.add_argument("-q","--query",type=str, required=True,
						help="query sequence to find protein domains in")
	parser.add_argument('-pdb', '--protein_database', type=str, default=LAST_DB, 
                        help='protein domains database')
	parser.add_argument('-cs', '--classification', type=str, default=CLASSIFICATION, 
                        help='protein domains classification file')
	parser.add_argument("-oug", "--domain_gff",type=str, default=DOMAINS_GFF,
						help="output domains gff format")
	parser.add_argument("-npd","--new_pdb",type=bool, default=False,
						help="create new protein database for last")
	parser.add_argument("-ds","--domains_seq",type=str, default=DOM_SEQ,
						help="file containg domains nucleic acid sequences")
	parser.add_argument("-dps","--domains_prot_seq",type=str, default=DOM_PROT_SEQ,
						help="file containg domains protein sequences")
	parser.add_argument("-thl","--th_length",type=int, default=20,
						help="length threshold for alignment")
	parser.add_argument("-thi","--th_identity",type=int, default=20,
						help="identity threshold for alignment")
	
	args = parser.parse_args()
	main(args)
