#!/usr/bin/env python3

import numpy as np
import subprocess
import matplotlib.pyplot as plt
import sys
import csv
import time
from operator import itemgetter
import os
import configuration


def domain_annotation(element, CLASSIFICATION):
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
		rep_type.append(annotation[element_name][0])
		rep_lineage.append(annotation[element_name][1])
	return domain, rep_type, rep_lineage
	

def hits_processing(sequence_hits):
	seq_length = sequence_hits[0,5]
	reverse_strand_idx = np.where(sequence_hits[:,4] == "-")[0]
	print(reverse_strand_idx)
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
	best_scores = []
	best_idx = []
	for idx in indices:
		print(np.where(scores[idx] == max(scores[idx])))
		best_idx.append(idx[np.where(scores[idx] == max(scores[idx]))[0][0]])
	return best_idx


def create_gff(sequence_hits, best_idx, seq_id, regions, OUTPUT, CLASSIFICATION):
	# Predefine the standard columns of GFF3 format
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
		with open(OUTPUT, "a") as gff:	
			gff.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tDomain={},Rep_type={},Rep_lineage={},Sequence={},Alignment={}\n".format(seq_id, SOURCE, FEATURE, alignment_start, alignment_end, score, strand, PHASE, domain[count], rep_type[count], rep_lineage[count], sequence, alignment))	
		count += 1	
	return xminimal, xmaximal, scores, strands, domain
	
	
def visualization(sequence_hits, seq_length, reverse_strand_idx, xminimal, xmaximal, scores, strands, domain, OUTPUT_PIC, fig, ax, profrep_module):
	if np.all(sequence_hits==0):
		seq_id = "sequence" 
	else:
		seq_id = sequence_hits[0,1]
	colors = [strand.replace("+", "red").replace("-", "blue") for strand in strands]
	if profrep_module == 0:
		plt.xlim([0, seq_length])
		ax.hlines(scores, xminimal, xmaximal, color=colors, lw=2)
		ax.set_title(seq_id)
		ax.set_xlabel('sequence bp')
		ax.set_ylabel('score')
		for count in range(len(domain)):
			plt.text(xminimal[count], scores[count], "{} {}".format(domain[count], strands[count]), size=9)
	else:
		y_upper_lim = ax.get_ylim()[1]
		ax.hlines([y_upper_lim + y_upper_lim/10]*len(xminimal), xminimal, xmaximal, color=colors, lw=2)
		print(y_upper_lim)
		for count in range(len(domain)):
			ax.text(xminimal[count], y_upper_lim + y_upper_lim/10, "{} {}".format(domain[count], strands[count]), size=9)
	if not os.path.exists(OUTPUT_PIC):
		os.makedirs(OUTPUT_PIC)
	output_pic_png = "{}/{}.png".format(OUTPUT_PIC, seq_id)
	fig.savefig(output_pic_png, bbox_inches="tight", format='png')
	return seq_id


def main(SEQUENCE, LAST_DB, CLASSIFICATION, OUTPUT, OUTPUT_PIC, fig_list, ax_list, profrep_module, NEW_PDB):		
	t2 = time.time()
	
	labels = []
	strand = []
	score = []
	rep_type = []
	rep_lineage = []
	start_pos = []
	end_pos = []
	line_counter = 1
	sequence_hits = np.empty((0,9))
	with open(SEQUENCE, "r") as fasta:
		seq_id = fasta.readline.strip()[1:]
	if NEW_PDB:
		subprocess.call("lastdb -p -cR01 {} {}".format(LAST_DB, LAST_DB), shell=True)
	tab = subprocess.Popen("lastal -F15 {} {} -f TAB ".format(LAST_DB, SEQUENCE), stdout=subprocess.PIPE, shell=True)
	maf = subprocess.Popen("lastal -F15 {} {} -f MAF ".format(LAST_DB, SEQUENCE), stdout=subprocess.PIPE, shell=True)
	maf.stdout.readline()
	seq_count = 0
	seq_names = []
	for line_tab in tab.stdout:
		line_tab = line_tab.decode("utf-8")
		if not line_tab.startswith('#'):
			line_maf = [maf.stdout.readline() for line_count in range(4)]
			reference_seq = line_maf[1].decode("utf-8").rstrip().split(" ")[-1]
			alignment_seq = line_maf[2].decode("utf-8").rstrip().split(" ")[-1]
			line = line_tab.rstrip().split("\t")
			line_maf = []
			element_name = line[1]
			#if line[6] != seq_id and seq_id != None:
			if line[6] != seq_id 
				[reverse_strand_idx, regions_plus, regions_minus, seq_length] = hits_processing(sequence_hits)
				print(seq_length)
				if reverse_strand_idx == []:
					positions = overlapping_regions(regions_plus)
					best_idx = best_score(sequence_hits[:,0], positions)
					[xminimal, xmaximal, scores, strands, domain] = create_gff(sequence_hits, best_idx, seq_id, regions_plus, OUTPUT, CLASSIFICATION)
				else:
					positions_plus = overlapping_regions(regions_plus)
					positions_minus = overlapping_regions(regions_minus)
					regions = regions_plus + regions_minus
					positions = positions_plus + [x + reverse_strand_idx for x in positions_minus]
					best_idx = best_score(sequence_hits[:,0], positions)
					[xminimal, xmaximal, scores, strands, domain] = create_gff(sequence_hits, best_idx, seq_id, regions, OUTPUT, CLASSIFICATION)
				seq_name = visualization(sequence_hits, seq_length, reverse_strand_idx, xminimal, xmaximal, scores, strands, domain, OUTPUT_PIC, fig_list[seq_count], ax_list[seq_count], profrep_module)
				seq_names.append(seq_name)
				seq_count += 1
				sequence_hits = np.empty((0,9))
			seq_id = line[6]
			line_parsed = np.array([int(line[0]), seq_id, int(line[7]), int(line[7]) + int(line[8]), line[9], int(line[10]), element_name, reference_seq, alignment_seq], dtype=object)
			sequence_hits = np.append(sequence_hits, [line_parsed], axis=0)
		else:
			maf.stdout.readline()
	
	if sequence_hits:		
		[reverse_strand_idx, regions_plus, regions_minus, seq_length] = hits_processing(sequence_hits)
		if reverse_strand_idx == []:
			positions = overlapping_regions(regions_plus)
			best_idx = best_score(sequence_hits[:,0], positions)
			[xminimal, xmaximal, scores, strands, domain] = create_gff(sequence_hits, best_idx, seq_id, regions_plus, OUTPUT, CLASSIFICATION)
		else:
			positions_plus = overlapping_regions(regions_plus)
			positions_minus = overlapping_regions(regions_minus)
			regions = regions_plus + regions_minus
			positions = positions_plus + [x + reverse_strand_idx for x in positions_minus]
			best_idx = best_score(sequence_hits[:,0], positions)
			[xminimal, xmaximal, scores, strands, domain] = create_gff(sequence_hits, best_idx, seq_id, regions, OUTPUT, CLASSIFICATION)
		seq_id = visualization(sequence_hits, seq_length, reverse_strand_idx, xminimal, xmaximal, scores, strands, domain, OUTPUT_PIC, fig_list[-1], ax_list[-1], profrep_module)
	seq_names.append(seq_id)
	
	print("ELAPSED_TIME_DOMAINS = {}".format(time.time() - t2))
	return seq_names

if __name__ == "__main__":
	import argparse
	
	LAST_DB = configuration.LAST_DB
	CLASSIFICATION = configuration.CLASSIFICATION
	DOMAINS_GFF = configuration.DOMAINS_GFF
	TMP = configuration.TMP
	print(LAST_DB)
	print(CLASSIFICATION)
	print(DOMAINS_GFF)
	
	parser = argparse.ArgumentParser()
	parser.add_argument("-q","--query",type=str, required=True,
						help="query sequence to find protein domains in")
	#parser.add_argument("-pdb","--protein_database",type=str, default="/mnt/raid/users/ninah/profrep/proteins_all",
						#help="protein domains database")					
	#parser.add_argument("-cs","--classification",type=str, default="/mnt/raid/users/ninah/profrep/classification.csv", 
						#help="protein domains classification file")	
	#parser.add_argument("-ouf","--output_gff", type=str,
						#help="output domains gff format")		
	parser.add_argument("-oup","--output_pic", type=str, default="/mnt/raid/users/ninah/profrep/",
						help="output domains picture")
	#parser.add_argument("-th","--threshold",type=int, default=100,
						#help="threshold for score")
	parser.add_argument("-npd","--new_pdb",type=str, default=False,
						help="create new protein database for last")
	
	args = parser.parse_args()
	SEQUENCE = args.query
	#LAST_DB = args.protein_database
	#CLASSIFICATION = args.classification
	#DOMAINS_GFF = args.output_gff
	#OUTPUT_PIC = args.output_pic
	#THRESHOLD = args.threshold
	NEW_PDB = args.new_pdb
	
	fig = plt.figure(figsize=(18, 8))
	ax = fig.add_subplot(111)		
	profrep_module = 0	
	print(DOMAINS_GFF)
	
	main(SEQUENCE, LAST_DB, CLASSIFICATION, DOMAINS_GFF, TMP, [fig], [ax], profrep_module, NEW_PDB)
