#!/usr/bin/env python3

import numpy as np
import subprocess
import time
from operator import itemgetter
from collections import Counter
import os
import configuration
from tempfile import NamedTemporaryFile


def domain_annotation(element, CLASSIFICATION):
	''' Assign protein domain to each hit from protein database  '''
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
	''' Gain hits intervals separately for forward and reverse strand '''
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
	''' Join all overalapping intervals '''
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
	''' From overlapping intervals take the one with the highest score '''
	best_scores = []
	best_idx = []
	for idx in indices:
		# if more hits have the same best score take only the first one
		best_idx.append(idx[np.where(scores[idx] == max(scores[idx]))[0][0]])
	return best_idx


def create_gff(sequence_hits, best_idx, seq_id, regions, OUTPUT_DOMAIN, CLASSIFICATION):
	''' Create track format of domains found in the query sequence(s)
	- custom atributes reported:
	Domain Name, Repet.type, Repet.lineage, 
	Original protein sequence, Aligned sequence, %Identity, 
	Relat.alignment lenght, Frameshifts per 100bp 
	
	In case of running protein domains module itself, save also reported
	domains	sequences in fasta format 
	'''
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
		dom_len = sequence_hits[i,9]
		[percent_ident, relat_align_len, relat_frameshifts] = filter_params(sequence, alignment, dom_len)
		with open(OUTPUT_DOMAIN, "a") as gff:	
			gff.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tName={},Rep_type={},Rep_lineage={},Sequence={},Alignment={},Identity={},Relat_length={},Relat_frameshifts={}\n".format(seq_id, SOURCE, FEATURE, alignment_start, alignment_end, score, strand, PHASE, domain[count], rep_type[count], rep_lineage[count], sequence, alignment, percent_ident, relat_align_len, relat_frameshifts))
		count += 1	
	return xminimal, xmaximal, scores, strands, domain


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


def domains_stat(domains_all, seq_ids, SUMMARY):
	'''  Create a file containing amounts of individual domains types'''
	with open(SUMMARY, "w") as sumfile:
		count_seq = 0
		for seq_id in seq_ids:
			sumfile.write("{}\n".format(seq_id))
			dom_in_seq = Counter(domains_all[count_seq])
			[sumfile.write("\t{}:{}\n".format(k,v)) for k,v in dom_in_seq.items()]
			count_seq += 1


def domain_search(QUERY, LAST_DB, CLASSIFICATION, OUTPUT_DOMAIN):		
	''' Search for protein domains using our protein database and external tool LAST,
	stdout is parsed in real time and hits for a single sequence undergo further processing
	- tabular format(TAB) to get info about position, score, orientation 
	- MAF format to gain alignment and original sequence
	'''	
	seq_ids = []
	xminimal_all = []
	xmaximal_all = []
	domains_all = []
	header_gff = "##gff-version 3"
	sequence_hits = np.empty((0,10))
	with open(QUERY, "r") as fasta:
		seq_id = fasta.readline().strip().split(" ")[0][1:]
	with open(OUTPUT_DOMAIN, "a") as gff:
		gff.write("{}\n".format(header_gff))
	tab = subprocess.Popen("lastal -F15 {} {} -n 100 -f TAB".format(LAST_DB, QUERY), stdout=subprocess.PIPE, shell=True)
	maf = subprocess.Popen("lastal -F15 {} {} -n 100 -f MAF".format(LAST_DB, QUERY), stdout=subprocess.PIPE, shell=True)
	maf.stdout.readline()
	for line_tab in tab.stdout:
		line_tab = line_tab.decode("utf-8")
		if not line_tab.startswith('#'):
			line_maf = [maf.stdout.readline() for line_count in range(4)]
			reference_seq = line_maf[1].decode("utf-8").rstrip().split(" ")[-1]
			alignment_seq = line_maf[2].decode("utf-8").rstrip().split(" ")[-1]
			line = line_tab.rstrip().split("\t")
			line_maf = []
			element_name = line[1]
			dom_len = line[5]
			if np.all(sequence_hits==0):
				seq_id = line[6]
				seq_ids.append(seq_id)
			if line[6] != seq_id: 
				[reverse_strand_idx, regions_plus, regions_minus, seq_length] = hits_processing(sequence_hits)
				if reverse_strand_idx == []:
					positions = overlapping_regions(regions_plus)
					best_idx = best_score(sequence_hits[:,0], positions)
					[xminimal, xmaximal, scores, strands, domain] = create_gff(sequence_hits, best_idx, seq_id, regions_plus, OUTPUT_DOMAIN, CLASSIFICATION)
				else:
					positions_plus = overlapping_regions(regions_plus)
					positions_minus = overlapping_regions(regions_minus)
					regions = regions_plus + regions_minus
					positions = positions_plus + [x + reverse_strand_idx for x in positions_minus]
					best_idx = best_score(sequence_hits[:,0], positions)
					[xminimal, xmaximal, scores, strands, domain] = create_gff(sequence_hits, best_idx, seq_id, regions, OUTPUT_DOMAIN, CLASSIFICATION)
				sequence_hits = np.empty((0,10))
				seq_id = line[6]
				seq_ids.append(seq_id)
				xminimal_all.append(xminimal)
				xmaximal_all.append(xmaximal)
				domains_all.append(domain)
			line_parsed = np.array([int(line[0]), seq_id, int(line[7]), int(line[7]) + int(line[8]), line[9], int(line[10]), element_name, reference_seq, alignment_seq, int(dom_len)], dtype=object)
			sequence_hits = np.append(sequence_hits, [line_parsed], axis=0)
		else:
			maf.stdout.readline()
	# The last (or the only one) sequence domains search
	if not np.all(sequence_hits==0):
		[reverse_strand_idx, regions_plus, regions_minus, seq_length] = hits_processing(sequence_hits)
		if reverse_strand_idx == []:
			positions = overlapping_regions(regions_plus)
			best_idx = best_score(sequence_hits[:,0], positions)
			[xminimal, xmaximal, scores, strands, domain] = create_gff(sequence_hits, best_idx, seq_id, regions_plus, OUTPUT_DOMAIN, CLASSIFICATION)
		else:
			positions_plus = overlapping_regions(regions_plus)
			positions_minus = overlapping_regions(regions_minus)
			regions = regions_plus + regions_minus
			positions = positions_plus + [x + reverse_strand_idx for x in positions_minus]
			best_idx = best_score(sequence_hits[:,0], positions)
			[xminimal, xmaximal, scores, strands, domain] = create_gff(sequence_hits, best_idx, seq_id, regions, OUTPUT_DOMAIN, CLASSIFICATION)
		xminimal_all.append(xminimal)
		xmaximal_all.append(xmaximal)
		domains_all.append(domain)
	
	return xminimal_all, xmaximal_all, domains_all, seq_ids
	
	
def main(args):
	
	t = time.time()
	
	QUERY = args.query
	LAST_DB = args.protein_database
	CLASSIFICATION = args.classification
	OUTPUT_DOMAIN = args.domain_gff
	NEW_PDB = args.new_pdb
	SUMMARY = args.summary_file
	
	if NEW_PDB:
		subprocess.call("lastdb -p -cR01 {} {}".format(LAST_DB, LAST_DB), shell=True)
	
	if not os.path.exists(LAST_DB):
		CLASSIFICATION = os.path.join(configuration.TOOL_DATA_DIR, CLASSIFICATION)
		LAST_DB = os.path.join(configuration.TOOL_DATA_DIR, LAST_DB) 
	
	[xminimal, xmaximal, domains_all, seq_ids] = domain_search(QUERY, LAST_DB, CLASSIFICATION, OUTPUT_DOMAIN)
	domains_stat(domains_all, seq_ids, SUMMARY)
	
	
	print("ELAPSED_TIME_DOMAINS = {} s".format(time.time() - t))
	
if __name__ == "__main__":
	import argparse
	
	LAST_DB = configuration.LAST_DB
	CLASSIFICATION = configuration.CLASSIFICATION
	DOMAINS_GFF = configuration.DOMAINS_GFF
	DOM_SUMMARY = configuration.DOM_SUMMARY

	parser = argparse.ArgumentParser()
	parser.add_argument("-q","--query",type=str, required=True,
						help="query sequence to find protein domains in")
	parser.add_argument('-pdb', "--protein_database", type=str, default=LAST_DB, 
                        help='protein domains database')
	parser.add_argument('-cs', '--classification', type=str, default=CLASSIFICATION, 
                        help='protein domains classification file')
	parser.add_argument("-oug", "--domain_gff",type=str, default=DOMAINS_GFF,
						help="output domains gff format")
	parser.add_argument("-npd","--new_pdb",type=bool, default=False,
						help="create new protein database for last")
	parser.add_argument("-sum","--summary_file",type=str, default=DOM_SUMMARY,
						help="summary file containing overview of amount of domains in individual seqs")
	
	args = parser.parse_args()
	main(args)
