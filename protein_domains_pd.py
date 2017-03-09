#!/usr/bin/env python3

import numpy as np
import subprocess
import math
import time
from operator import itemgetter
from collections import Counter
from itertools import groupby
import os
import configuration
from tempfile import NamedTemporaryFile

np.set_printoptions(threshold=np.nan)


def domain_annotation(elements, CLASSIFICATION):
	''' Assign protein domain to each hit from protein database  '''
	domains = []
	annotations = []
	with open(CLASSIFICATION, "r") as cl_tbl:
		annotation = {}
		for line in cl_tbl:
			record = line.rstrip().split("\t")
			annotation[record[0]] = record[1:]	
	for i in range(len(elements)):
		domains.append(elements[i].split("__")[0].split("-")[1])
		element_name = "__".join(elements[i].split("__")[1:])
		if element_name in annotation.keys():			
			annotations.append("/".join([elements[i].split("__")[0].split("-")[1],("/".join(annotation[element_name]))]))
		else:
			rep_type.append(" ")
			rep_lineage.append(" ")
			annotations.append("unknown/unknown")
	return annotations
	

def hits_processing(seq_len, start, end, strand):
	''' Gain hits intervals separately for forward and reverse strand '''
	reverse_strand_idx = np.where(strand == "-")[0]
	if not reverse_strand_idx.any():
		start_pos_plus = start + 1
		end_pos_plus = end
		regions_plus = list(zip(start_pos_plus, end_pos_plus))
		regions_minus = []
	else:
		reverse_strand_idx = reverse_strand_idx[0]
		start_pos_plus = start[0:reverse_strand_idx] + 1
		end_pos_plus = end[0:reverse_strand_idx]
		start_pos_minus = seq_len[0] - end[reverse_strand_idx:] + 1
		end_pos_minus = seq_len[0] - start[reverse_strand_idx:]
		regions_plus= list(zip(start_pos_plus, end_pos_plus))
		regions_minus = list(zip(start_pos_minus, end_pos_minus))
	return reverse_strand_idx, regions_plus, regions_minus


def overlapping_regions(input_data):
	''' Join all overalapping intervals '''
	if input_data: 
		sorted_idx, sorted_data = zip(*sorted([(index,data) for index,data in enumerate(input_data)], key=itemgetter(1)))
		merged_ends = input_data[sorted_idx[0]][1]
		intervals = []
		data =[]
		output_intervals = [] 
		output_data = []
		for i,j in zip(sorted_idx, sorted_data):
			if input_data[i][0] < merged_ends:
				merged_ends = max(input_data[i][1], merged_ends)
				intervals.append(i)
				data.append(j)
			else:
				output_intervals.append(intervals)
				output_data.append(data)
				intervals = []
				data =[]
				intervals.append(i)
				data.append(j)
				merged_ends = input_data[i][1]		
		output_intervals.append(intervals)
		output_data.append(data)
		mins = [x[0][0] for x in output_data]
		maxs = [max(x, key=itemgetter(1))[1] for x in output_data]
	else:
		mins = []
		maxs = []
		output_intervals = []
		output_data = []
	return mins, maxs, output_data, output_intervals
	
def annotations_dict(annotations):
	# hash table based on unique annotations of the hits, each annotation has serial number assigned which indexes the row in the score_table
	classes_dict = {classes: idx for idx, classes in enumerate(set(annotations))}
	return classes_dict

def score_table(mins, maxs, data, annotations, scores, CLASSIFICATION):
	classes_dict = annotations_dict(annotations)
	score_matrix = np.zeros((len(classes_dict),maxs-mins + 1), dtype=int)
	count = 0
	for item in annotations:
		saved_scores = score_matrix[classes_dict[item], data[count][0]-mins : data[count][1]-mins + 1]
		new_scores = [scores[count]]*len(saved_scores)
		score_matrix[classes_dict[item], data[count][0]-mins : data[count][1]-mins + 1] = [max(*pos_score) for pos_score in zip(saved_scores, new_scores)]	
		count += 1
	return score_matrix, classes_dict

	
def score_matrix_evaluation(score_matrix, classes_dict):
	ann_per_reg = []
	max_scores_reg = []
	for position in score_matrix.T:
		max_scores_reg.append(max(position))
		# 80% of the best score
		threshold = max(position) * 0.8
		above_th = [idx for idx, score in enumerate(position) if position[idx] >= threshold]
		# select unique annotations in one position that are above threshold
		ann_per_pos = list(set([key for key, value in classes_dict.items() if value in above_th]))
		ann_per_reg.append(ann_per_pos)
	return ann_per_reg, max_scores_reg
	
def group_annot_regs(ann_per_reg, mins, maxs):
	# tranform list of lists (potential multiple annotations for every position ) to flat list of all annotations
	unique_annotations = list(set([item for sublist in ann_per_reg for item in sublist]))
	domain_type = list(set([annotation.split("/")[0] for annotation in unique_annotations]))
	#classification_list = [("/".join(annotation.split("/")[1:]) for annotation in unique_annotations]
	classification_list = [os.path.join(annotation,"") for annotation in unique_annotations]
	ann_substring = os.path.commonprefix(classification_list).rpartition("/")[0]
	#ann_substring = os.path.commonprefix(classification_list)
	domain_type = "/".join(domain_type) 
	return domain_type, ann_substring, unique_annotations
	
def best_score(scores, region):
	''' From overlapping intervals take the one with the highest score '''
	# if more hits have the same best score take only the first one
	best_idx = region[np.where(scores == max(scores))[0][0]]
	best_idx_reg = np.where(scores == max(scores))[0][0]
	return best_idx, best_idx_reg

	
def create_gff3(domain_type, ann_substring, unique_annotations, dom_start, dom_end, best_idx, annotation_best, strand, score, seq_id, db_seq, query_seq, domain_size, positions, OUTPUT_DOMAIN):
	with open(OUTPUT_DOMAIN, "a") as gff:	
			best_start = positions[best_idx][0]
			best_end = positions[best_idx][1]
			best_score = score[best_idx]
			# proportion of length of the best hit to the whole region length found by base
			length_proportion = int((best_end - best_start + 1)/(dom_end - dom_start + 1)*100)
			db_seq_best = db_seq[best_idx]
			query_seq_best = query_seq[best_idx]
			domain_size_best = domain_size[best_idx]
			[percent_ident, relat_align_len, relat_frameshifts] = filter_params(db_seq_best, query_seq_best, domain_size_best)
			ann_substring = "/".join(ann_substring.split("/")[1:])
			if ann_substring is '':
				ann_substring = "NONE(Annotations from different classes)"
			unique_annotations = ";".join(unique_annotations)
			if "/" in domain_type:
				gff.write("{}\t{}\t{}\t{}\t{}\t.\t{}\t{}\tName={},Classification=Ambiguous_domain,Region_Annotations={}\n".format(seq_id, configuration.SOURCE, configuration.DOMAINS_FEATURE, dom_start, dom_end, strand, configuration.PHASE, domain_type, unique_annotations))
			else:
				gff.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tName={},Classification={},Region_Annotations={},Best_HIT={}-{}[{}%],Best_annotation={},DB_Seq={},HIT_Seq={},Identity={},Relat_length={},Relat_frameshifts={}\n".format(seq_id, configuration.SOURCE, configuration.DOMAINS_FEATURE, dom_start, dom_end, best_score, strand, configuration.PHASE, domain_type, ann_substring, unique_annotations, best_start, best_end, length_proportion, annotation_best, query_seq_best, db_seq_best, percent_ident, relat_align_len, relat_frameshifts))
				

def filter_params(db, seq, protein_len):
	''' Calculate basic statistics of the quality of alignment '''
	num_ident = 0
	count_frm = 0
	alignment_len = 0
	for i,j in zip(db.upper(), seq):
		if i == j:
			num_ident += 1
		if j == "/" or j == "\\":
			count_frm += 1
	relat_align_len = round(len(db)/protein_len, 3) 
	align_identity = round(num_ident/len(seq), 2)
	relat_frameshifts = round(count_frm/math.ceil((len(seq)/100)),2)
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


def line_generator(tab_pipe, maf_pipe, start):
		if hasattr(line_generator, "a"):
			seq_id = line_generator.a.split("\t")[6]
			yield line_generator.a.encode("utf-8")
			del line_generator.a
		line_tab = ""
		for line_tab in tab_pipe:
			line_tab = line_tab.decode("utf-8")
			if not line_tab.startswith('#'):
					if start:
						seq_id = line_tab.split("\t")[6]
						start = False
					line_maf = [maf_pipe.readline() for line_count in range(4)]
					db_seq = line_maf[1].decode("utf-8").rstrip().split(" ")[-1]
					alignment_seq = line_maf[2].decode("utf-8").rstrip().split(" ")[-1]
					line = "{}\t{}\t{}".format(line_tab, db_seq, alignment_seq)
					line6 = line.split("\t")[6]
					if seq_id != line6:
						line_generator.a = line
						return 
					else:
						yield line.encode("utf-8")		
			else:
				maf_pipe.readline()	
		if line_tab == "":
			raise RuntimeError
		else:
			return

def domain_search(QUERY, LAST_DB, CLASSIFICATION, OUTPUT_DOMAIN):		
	''' Search for protein domains using our protein database and external tool LAST,
	stdout is parsed in real time and hits for a single sequence undergo further processing
	- tabular format(TAB) to get info about position, score, orientation 
	- MAF format to gain alignment and original sequence
	'''	
	
	#### samostatna funkcia
	tab = subprocess.Popen("lastal -F15 {} {} -L 10 -m 70 -f TAB".format(LAST_DB, QUERY), stdout=subprocess.PIPE, shell=True)
	maf = subprocess.Popen("lastal -F15 {} {} -L 10 -m 70 -f MAF".format(LAST_DB, QUERY), stdout=subprocess.PIPE, shell=True)

	tab_pipe = tab.stdout
	maf_pipe = maf.stdout
	maf_pipe.readline()


	seq_ids = []
	domain_reg = []
	xminimal_all = []
	xmaximal_all = []
	domains_all = []
	header_gff = "##gff-version 3"
	with open(OUTPUT_DOMAIN, "a") as gff:
		gff.write("{}\n".format(header_gff))
	
	start = True
	while True:	
		try:
			sequence_hits = np.genfromtxt(line_generator(tab_pipe, maf_pipe, start),
				names="score, name_db, start_db, al_size_db, strand_db, seq_size_db, name_q, start_q, al_size_q, strand_q, seq_size_q, block1, block2, block3, db_seq, q_seq",
				usecols="score, name_q, start_q, al_size_q, strand_q, seq_size_q, name_db, db_seq, q_seq, seq_size_db",
				dtype=None)
		
		except RuntimeError:
			break
		
		## PARSING LAST OUTPUT ##
		score = sequence_hits['score'].astype("int")
		seq_id = sequence_hits['name_q'][0].astype("str")
		start_hit = sequence_hits['start_q'].astype("int")
		end_hit = start_hit + sequence_hits['al_size_q'].astype("int")
		strand = sequence_hits['strand_q'].astype("str")
		seq_len = sequence_hits['seq_size_q'].astype("int")
		domain_db = sequence_hits['name_db'].astype("str")
		db_seq = sequence_hits['db_seq'].astype("str")
		query_seq = sequence_hits['q_seq'].astype("str") 
		domain_size = sequence_hits['seq_size_db'].astype("int")
		[reverse_strand_idx, positions_plus, positions_minus] = hits_processing(seq_len, start_hit, end_hit, strand)
		strand_gff = "+"
		[mins_plus, maxs_plus, data_plus, indices_plus] = overlapping_regions(positions_plus)
		[mins_minus, maxs_minus, data_minus, indices_minus] = overlapping_regions(positions_minus)
		positions = positions_plus + positions_minus
		indices_overal = indices_plus + [x + reverse_strand_idx for x in indices_minus]
		mins = mins_plus + mins_minus
		maxs = maxs_plus + maxs_minus
		# tuples dvojice min-max usekov ktore sa prekryvaju v ramci jedneho regionu=clustru
		data = data_plus + data_minus
		# process every region of overlapping hits sequentially
		count_region = 0
		for region in indices_overal:
			db_names = domain_db[np.array(region)]
			scores = score[np.array(region)]
			annotations = domain_annotation(db_names, CLASSIFICATION)
			[score_matrix, classes_dict] = score_table(mins[count_region], maxs[count_region], data[count_region], annotations, scores, CLASSIFICATION)	
			[ann_per_reg, max_scores_reg] = score_matrix_evaluation(score_matrix, classes_dict)
			[domain_type, ann_substring, unique_annotations] = group_annot_regs(ann_per_reg, mins[count_region], maxs[count_region])
			[best_idx, best_idx_reg] = best_score(scores, region)
			annotation_best = annotations[best_idx_reg]
			if count_region == len(indices_plus):
				strand_gff = "-"
			create_gff3(domain_type, ann_substring, unique_annotations, mins[count_region], maxs[count_region], best_idx, annotation_best, strand_gff, score, seq_id, db_seq, query_seq, domain_size, positions, OUTPUT_DOMAIN)
			count_region += 1
			domain_reg.append(domain_type)
		xminimal_all.append(mins)
		xmaximal_all.append(maxs)
		domains_all.append(domain_reg)
		domain_reg = []
		seq_ids.append(seq_id)
	return xminimal_all, xmaximal_all, domains_all, seq_ids
	
	
def main(args):
	
	t = time.time()
	
	QUERY = args.query
	LAST_DB = args.protein_database
	CLASSIFICATION = args.classification
	OUTPUT_DOMAIN = args.domain_gff
	NEW_PDB = args.new_pdb
	SUMMARY = args.summary_file
	OUTPUT_DIR = args.output_dir
	
	if NEW_PDB:
		subprocess.call("lastdb -p -cR01 {} {}".format(LAST_DB, LAST_DB), shell=True)
		
	if not os.path.exists(OUTPUT_DIR) and not os.path.exists(OUTPUT_DOMAIN):
		os.makedirs(OUTPUT_DIR)
		OUTPUT_DOMAIN = os.path.join(OUTPUT_DIR, os.path.basename(OUTPUT_DOMAIN))
		SUMMARY = os.path.join(OUTPUT_DIR, os.path.basename(SUMMARY))
	elif os.path.exists(OUTPUT_DIR) and not os.path.exists(OUTPUT_DOMAIN):
		OUTPUT_DOMAIN = os.path.join(OUTPUT_DIR, os.path.basename(OUTPUT_DOMAIN))
		SUMMARY = os.path.join(OUTPUT_DIR, os.path.basename(SUMMARY))
	
	
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
						help="reference sequence to find protein domains")
	parser.add_argument('-pdb', "--protein_database", type=str, default=LAST_DB, 
                        help='protein domains database file')
	parser.add_argument('-cs', '--classification', type=str, default=CLASSIFICATION, 
                        help='protein domains classification file')
	parser.add_argument("-oug", "--domain_gff",type=str, default=DOMAINS_GFF,
						help="output domains gff format")
	parser.add_argument("-npd","--new_pdb",type=str, default=False,
						help="create new protein database for lastal")
	parser.add_argument("-sum","--summary_file",type=str, default=DOM_SUMMARY,
						help=" output summary file containing overview of amount of domains in individual seqs")
	parser.add_argument("-dir","--output_dir",type=str, default=configuration.TMP,
						help="specify if you want to change the output directory")
	
	args = parser.parse_args()
	main(args)
