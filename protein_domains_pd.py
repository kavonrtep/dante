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
import sys

np.set_printoptions(threshold=np.nan)

def sc_hash():
	''' Create hash table for alignment similarity counting: for every 
	combination of aminoacids in alignment assign score from protein 
	scoring matrix defined in configuration file  '''
	sc_dict = {}
	with open(configuration.SC_MATRIX) as smatrix:
		count = 1
		for line in smatrix:
			if not line.startswith("#"):
				if count == 1:
					aa_all = line.rstrip().replace(" ","")
				else:
					count_aa = 1
					line = list(filter(None,line.rstrip().split(" ")))
					for aa in aa_all:
						sc_dict["{}{}".format(line[0],aa)] = line[count_aa]
						count_aa += 1
				count += 1
	return sc_dict
	
	
def characterize_fasta(QUERY, WIN_DOM):
	''' Find the sequences, their lengths, starts, ends and if 
	they exceed the window '''
	with open(QUERY) as query:
		headers = []
		fasta_lengths = []
		seq_starts = []
		seq_ends = []
		fasta_chunk_len = 0
		count_line = 1
		for line in query:
			line = line.rstrip()
			if line.startswith(">"):
				headers.append(line.rstrip())
				fasta_lengths.append(fasta_chunk_len)
				fasta_chunk_len = 0
				seq_starts.append(count_line + 1)
				seq_ends.append(count_line - 1)
			else:
				fasta_chunk_len += len(line)
			count_line += 1
		seq_ends.append(count_line)
		seq_ends = seq_ends[1:]
		fasta_lengths.append(fasta_chunk_len)
		fasta_lengths = fasta_lengths[1:]
		# control if there are correct (unique) names for individual seqs
		if len(headers) > len(set([header.split(" ")[0] for header in headers])):
			raise NameError('''Sequences in multifasta format are not named correctly:
							seq IDs(before the first space) are the same''')

	above_win = [idx for idx, value in enumerate(fasta_lengths) if value > WIN_DOM]
	below_win = [idx for idx, value in enumerate(fasta_lengths) if value <= WIN_DOM]
	lens_above_win = np.array(fasta_lengths)[above_win]
	return headers, above_win, below_win, lens_above_win, seq_starts, seq_ends


def split_fasta(QUERY, WIN_DOM, step, headers, above_win, below_win, lens_above_win, seq_starts, seq_ends):
	''' Create temporary file containing all sequences - the ones that exceed 
	the window are cut with a set overlap (greater than domain size with a reserve) '''
	with open(QUERY, "r") as query:
		count_fasta_divided = 0 
		count_fasta_not_divided = 0 
		ntf = NamedTemporaryFile(delete=False)
		divided = np.array(headers)[above_win]
		row_length = configuration.FASTA_LINE
		for line in query:
			line = line.rstrip()
			if line.startswith(">") and line in divided:
				stop_line = seq_ends[above_win[count_fasta_divided]] - seq_starts[above_win[count_fasta_divided]] + 1
				count_line = 0
				whole_seq = ""
				for line2 in query:
					whole_seq = "".join([whole_seq, line2.rstrip()])
					count_line += 1
					if count_line == stop_line:
						break
				## create list of starting positions for individual parts of a seq with a step given by a window and overlap
				windows_starts = list(range(0, lens_above_win[count_fasta_divided], step))
				## create list of ending positions (starting pos + window), the last element is the whole seq length
				windows_ends = [x + WIN_DOM if x + WIN_DOM < lens_above_win[count_fasta_divided] else  lens_above_win[count_fasta_divided] for x in windows_starts]
				count_part = 1
				for start_part, end_part in zip(windows_starts, windows_ends):
					seq_part = whole_seq[start_part:end_part]
					if count_part == len(windows_starts):
						ntf.write("{}_PART{}_LAST:{}-{}\n{}\n".format(line.split(" ")[0], count_part, start_part + 1, end_part, "\n".join([seq_part[i:i+row_length] for i in range(0, len(seq_part), row_length)])).encode("utf-8"))
					else:
						ntf.write("{}_PART{}:{}-{}\n{}\n".format(line.split(" ")[0], count_part, start_part + 1, end_part, "\n".join([seq_part[i:i+row_length] for i in range(0, len(seq_part), row_length)])).encode("utf-8"))
					count_part += 1
				count_fasta_divided += 1
			elif line.startswith(">") and line not in divided:
				length_seq = seq_ends[below_win[count_fasta_not_divided]] - seq_starts[below_win[count_fasta_not_divided]] + 1
				ntf.write("{}\n{}".format(line, "".join([query.readline() for x in range(length_seq)])).encode("utf-8"))
				count_fasta_not_divided += 1
		query_temp = ntf.name
		ntf.close()
	print(query_temp)
	return(query_temp)
				

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
	''' Join all overlapping intervals(hits) to clusters (potential domains),
	get list of start-end positions of individual hits within the interval, 
	list of minimus and maximums as well as the indices in the original 
	sequence_hits structure for the hits belonging to the same clusters '''
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
	''' Hash table where annotations of the hits within a clusters are the keys. 
	Each annotation has serial number assigned which indexes the row in the score_table '''
	classes_dict = {classes: idx for idx, classes in enumerate(set(annotations))}
	return classes_dict


def score_table(mins, maxs, data, annotations, scores, CLASSIFICATION):
	''' Score table is created based on the annotations occurance in the cluster.
	Matrix axis y corresponds to individual annotations (indexed according to classes_dict),
	axis x represents positions of analyzed seq in a given cluster.
	For every hit within cluster, array of scores on the corresponding position 
	is recorded to the table in case if the score on certain position is so far the highest 
	for the certain position and certain annotation '''
	classes_dict = annotations_dict(annotations)
	score_matrix = np.zeros((len(classes_dict),maxs-mins + 1), dtype=int)
	count = 0
	for item in annotations:
		saved_scores = score_matrix[classes_dict[item], data[count][0]-mins : data[count][1]-mins + 1]
		new_scores = [scores[count]]*len(saved_scores)
		score_matrix[classes_dict[item], data[count][0]-mins : data[count][1]-mins + 1] = [max(*pos_score) for pos_score in zip(saved_scores, new_scores)]	
		count += 1
	return score_matrix, classes_dict

	
def score_matrix_evaluation(score_matrix, classes_dict, THRESHOLD_SCORE):
	''' Score matrix is evaluated based on each position.
	For every position the list of annotations with a score which reaches 
	certain percentage of the overal best score of the cluster are stored '''
	ann_per_reg = []
	overal_best_score_reg = max((score_matrix.max(axis=1)))
	for position in score_matrix.T:
		## score threshold calculated as a percentage of the OVERALL best score in the cluster
		threshold = overal_best_score_reg * THRESHOLD_SCORE/100
		above_th = [idx for idx, score in enumerate(position) if position[idx] >= threshold]
		## select unique annotations in one position that are above threshold
		ann_per_pos = list(set([key for key, value in classes_dict.items() if value in above_th]))
		ann_per_reg.append(ann_per_pos)
	return ann_per_reg
	
	
def group_annot_regs(ann_per_reg):
	''' Get list of domains, annotations, longest common annotations and 
	counts of positions with certain annotation per regions '''
	## tranform list of lists (potential multiple annotations for every position ) to flat list of all annotations
	all_annotations = [item for sublist in ann_per_reg for item in sublist]
	unique_annotations = list(set(all_annotations))
	ann_pos_counts = [all_annotations.count(x) for x in unique_annotations]
	unique_annotations = list(set([item for sublist in ann_per_reg for item in sublist]))
	domain_type = list(set([annotation.split("/")[0] for annotation in unique_annotations]))
	classification_list = [os.path.join(annotation,"") for annotation in unique_annotations]
	ann_substring = os.path.commonprefix(classification_list).rpartition("/")[0]
	domain_type = "/".join(domain_type) 
	return domain_type, ann_substring, unique_annotations, ann_pos_counts
	
	
def best_score(scores, region):
	''' From overlapping intervals take the one with the highest score '''
	## if more hits have the same best score take only the first one
	best_idx = region[np.where(scores == max(scores))[0][0]]
	best_idx_reg = np.where(scores == max(scores))[0][0]
	return best_idx, best_idx_reg

	
def create_gff3(domain_type, ann_substring, unique_annotations, ann_pos_counts, dom_start, dom_end, step, best_idx, annotation_best, strand, score, seq_id, db_seq, query_seq, domain_size, positions, gff):
	''' Record obtained information about domain corresponding to individual cluster to common gff file '''
	best_start = positions[best_idx][0]
	best_end = positions[best_idx][1]
	best_score = score[best_idx]
	## proportion of length of the best hit to the whole region length found by base
	length_proportion = int((best_end - best_start + 1)/(dom_end - dom_start + 1)*100)
	db_seq_best = db_seq[best_idx]
	query_seq_best = query_seq[best_idx]
	domain_size_best = domain_size[best_idx]
	[percent_ident, align_similarity, relat_align_len, relat_frameshifts] = filter_params(db_seq_best, query_seq_best, domain_size_best)
	ann_substring = "/".join(ann_substring.split("/")[1:])
	if "PART" in seq_id:
		part = int(seq_id.split("PART")[1].split(":")[0].split("_")[0])
		dom_start = dom_start + (part-1)*step
		dom_end = dom_end + (part-1)*step
		best_start = best_start + (part-1)*step
		best_end = best_end + (part-1)*step
	if ann_substring is '':
		ann_substring = "NONE(Annotations from different classes)"
	if len(unique_annotations) > 1:
		unique_annotations = ";".join(["{}[{}bp]".format(ann, pos) for ann, pos in zip(unique_annotations, ann_pos_counts)])
	else:
		unique_annotations = unique_annotations[0]
	if "/" in domain_type:
		gff.write("{}\t{}\t{}\t{}\t{}\t.\t{}\t{}\tName={},Classification=Ambiguous_domain,Region_Annotations={}\n".format(seq_id, configuration.SOURCE, configuration.DOMAINS_FEATURE, dom_start, dom_end, strand, configuration.PHASE, domain_type, unique_annotations))
	else:
		gff.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tName={},Classification={},Region_Annotations={},Best_HIT={}:{}-{}[{}%],DB_Seq={},Query_Seq={},Identity={},Similarity={},Relat_length={},Relat_frameshifts={}\n".format(seq_id, configuration.SOURCE, configuration.DOMAINS_FEATURE, dom_start, dom_end, best_score, strand, configuration.PHASE, domain_type, ann_substring, unique_annotations, annotation_best, best_start, best_end, length_proportion, db_seq_best, query_seq_best, percent_ident, align_similarity, relat_align_len, relat_frameshifts))
			

def filter_params(db, query, protein_len):
	''' Calculate basic statistics of the quality of the alignment '''
	sc_dict = sc_hash()
	num_ident = 0
	count_frm = 0
	count_similarity = 0 
	alignment_len = 0
	for i,j in zip(db.upper(), query.upper()):
		if i == j and i != "X":
			num_ident += 1
		if j == "/" or j == "\\":
			count_frm += 1
		if (i.isalpha() or i == "*") and (j.isalpha() or j == "*"):
			if int(sc_dict["{}{}".format(i,j)]) > 0:
				count_similarity += 1
	## gapless alignment length proportional to the domain protein length
	relat_align_len = round((len(db) - db.count("-"))/protein_len, 3) 
	## proportional identical bases (except of X) to al.length
	align_identity = round(num_ident/len(db), 2)
	## proportional count of positive scores from scoring matrix to al. length 
	align_similarity = round(count_similarity/len(db),2)
	## number of frameshifts per 100 bp
	relat_frameshifts = round(count_frm/math.ceil((len(query)/100)),2)
	return align_identity, align_similarity, relat_align_len, relat_frameshifts, 	


def domains_stat(domains_all, seq_ids, SUMMARY):
	'''  Create a file containing amounts of individual domains types'''
	with open(SUMMARY, "w") as sumfile:
		count_seq = 0 
		if seq_ids == []:
			sumfile.write("NO DOMAINS")
		for seq_id in seq_ids:
			sumfile.write("{}\n".format(seq_id))
			dom_in_seq = Counter(domains_all[count_seq])
			[sumfile.write("\t{}:{}\n".format(k,v)) for k,v in dom_in_seq.items()]
			count_seq += 1


def line_generator(tab_pipe, maf_pipe, start):
	''' Yield individual lines of LASTAL stdout for single sequence '''
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
				line_id = line.split("\t")[6]
				if seq_id != line_id:
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


def domain_search(QUERY, LAST_DB, CLASSIFICATION, OUTPUT_DOMAIN, THRESHOLD_SCORE, WIN_DOM, OVERLAP_DOM):		
	''' Search for protein domains using our protein database and external tool LAST,
	stdout is parsed in real time and hits for a single sequence undergo further processing
	- tabular format(TAB) to get info about position, score, orientation 
	- MAF format to gain alignment and original sequence
	'''		
	
	step = WIN_DOM - OVERLAP_DOM
	[headers, above_win, below_win, lens_above_win, seq_starts, seq_ends] = characterize_fasta(QUERY, WIN_DOM)
	query_temp = split_fasta(QUERY, WIN_DOM, step, headers, above_win, below_win, lens_above_win, seq_starts, seq_ends)

	tab = subprocess.Popen("lastal -F15 {} {} -L 10 -m 70 -p BL80 -f TAB".format(LAST_DB, query_temp), stdout=subprocess.PIPE, shell=True)
	maf = subprocess.Popen("lastal -F15 {} {} -L 10 -m 70 -p BL80 -f MAF".format(LAST_DB, query_temp), stdout=subprocess.PIPE, shell=True)

	tab_pipe = tab.stdout
	maf_pipe = maf.stdout
	maf_pipe.readline()

	seq_ids = []
	domain_reg = []
	xminimal_all = []
	xmaximal_all = []
	domains_all = []
	header_gff = "##gff-version 3"
	with open(OUTPUT_DOMAIN, "w") as gff:
		gff.write("{}\n".format(header_gff))
	gff = open(OUTPUT_DOMAIN,"a")
	start = True
	while True:	
		try:
			sequence_hits = np.genfromtxt(line_generator(tab_pipe, maf_pipe, start),
				names="score, name_db, start_db, al_size_db, strand_db, seq_size_db, name_q, start_q, al_size_q, strand_q, seq_size_q, block1, block2, block3, db_seq, q_seq",
				usecols="score, name_q, start_q, al_size_q, strand_q, seq_size_q, name_db, db_seq, q_seq, seq_size_db",
				dtype=None)
		except RuntimeError:
			break
		## if there are no domains found
		if sequence_hits.size is 0:
			with open(OUTPUT_DOMAIN, "w") as gff:
				gff.write("NO DOMAINS")
			return [],[],[],[]
		
		############# PARSING LASTAL OUTPUT ############################
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
		data = data_plus + data_minus
		## process every region (cluster) of overlapping hits sequentially
		count_region = 0
		for region in indices_overal:
			db_names = domain_db[np.array(region)]
			scores = score[np.array(region)]
			annotations = domain_annotation(db_names, CLASSIFICATION)
			[score_matrix, classes_dict] = score_table(mins[count_region], maxs[count_region], data[count_region], annotations, scores, CLASSIFICATION)	
			ann_per_reg = score_matrix_evaluation(score_matrix, classes_dict, THRESHOLD_SCORE)
			[domain_type, ann_substring, unique_annotations, ann_pos_counts] = group_annot_regs(ann_per_reg)
			[best_idx, best_idx_reg] = best_score(scores, region)
			annotation_best = annotations[best_idx_reg]
			if count_region == len(indices_plus):
				strand_gff = "-"
			create_gff3(domain_type, ann_substring, unique_annotations, ann_pos_counts, mins[count_region], maxs[count_region], step, best_idx, annotation_best, strand_gff, score, seq_id, db_seq, query_seq, domain_size, positions, gff)
			count_region += 1
			domain_reg.append(domain_type)
		xminimal_all.append(mins)
		xmaximal_all.append(maxs)
		domains_all.append(domain_reg)
		domain_reg = []
		seq_ids.append(seq_id)
	os.unlink(query_temp)
	gff.close()
	if any("PART" in x for x in seq_ids):
		[xminimal_all, xmaximal_all, domains_all, seq_ids] = adjust_gff(OUTPUT_DOMAIN, WIN_DOM, OVERLAP_DOM, step)
	return xminimal_all, xmaximal_all, domains_all, seq_ids
	
	
def adjust_gff(OUTPUT_DOMAIN, WIN_DOM, OVERLAP_DOM, step):
	''' Original gff file is adjusted in case of containing cut parts 
	- for consecutive sequences overlap is divided to half with first half 
	of records(domains) belonging to the first sequence and second to the following one.
	Duplicate domains going through the middle of the overlap are removed.
	First and the last part (marked as LAST) of a certain sequence are 
	handled separately as the are overlapped from one side only '''
	xminimal_all = []
	xmaximal_all = []
	dom_all = []
	seq_id_all = []
	xminimal = []
	xmaximal = []
	dom  = []
	seen = set()
	adjusted_file = os.path.join(os.path.dirname(OUTPUT_DOMAIN), configuration.ADJUSTED_GFF)
	with open(adjusted_file, "w") as adjusted_gff:
		adjusted_gff.write("##gff-version 3\n")
		with open(OUTPUT_DOMAIN, "r") as primary_gff:
			next(primary_gff)
			start = True
			for line in primary_gff:
				split_line = line.split("\t")
				if start:
					seq_id_all.append(split_line[0].split("_PART")[0])
					start = False
				seq_id = split_line[0].split("_PART")[0]
				if "PART" in line:
					line_without_id = "\t".join(split_line[1:])
					part = int(split_line[0].split("_PART")[1].split(":")[0].split("_")[0])
					if seq_id != seq_id_all[-1]:
						seq_id_all.append(seq_id)
						xminimal_all.append(xminimal)
						xmaximal_all.append(xmaximal)
						dom_all.append(dom)
						xminimal = []
						xmaximal = []
						dom = []
						
					## first part of the sequence
					if part == 1:
						cut_end = WIN_DOM - OVERLAP_DOM/2 
						if int(split_line[3]) <= cut_end <= int(split_line[4]):
							if line_without_id not in seen:
								adjusted_gff.write("{}\t{}".format(seq_id, line_without_id))
								xminimal.append(split_line[3])
								xmaximal.append(split_line[4])
								dom.append(split_line[-1].split(",")[0].split("=")[1])
								seen.add(line_without_id)
						elif int(split_line[4]) < cut_end:
							adjusted_gff.write("{}\t{}".format(seq_id, line_without_id))
							xminimal.append(split_line[3])
							xmaximal.append(split_line[4])
							dom.append(split_line[-1].split(",")[0].split("=")[1])
								
					## last part of the sequence
					elif "LAST" in split_line[0]:
						cut_start = OVERLAP_DOM/2 + (part-1)*step 
						if int(split_line[3]) <= cut_start <= int(split_line[4]):
							if line_without_id not in seen:
								adjusted_gff.write("{}\t{}".format(seq_id, line_without_id))
								xminimal.append(split_line[3])
								xmaximal.append(split_line[4])
								dom.append(split_line[-1].split(",")[0].split("=")[1])
								seen.add(line_without_id)
						elif int(split_line[3]) > cut_start: 
							adjusted_gff.write("{}\t{}".format(seq_id, line_without_id))
							xminimal.append(split_line[3])
							xmaximal.append(split_line[4])
							dom.append(split_line[-1].split(",")[0].split("=")[1])
							
					## middle part of the sequence
					else:
						cut_start = OVERLAP_DOM/2 + (part-1)*step
						cut_end = WIN_DOM - OVERLAP_DOM/2 + (part-1)*step
						if int(split_line[3]) <= cut_start <= int(split_line[4]) or int(split_line[3]) <= cut_end <= int(split_line[4]):
							if line_without_id not in seen:
								adjusted_gff.write("{}\t{}".format(seq_id, line_without_id))
								xminimal.append(split_line[3])
								xmaximal.append(split_line[4])
								dom.append(split_line[-1].split(",")[0].split("=")[1])
								seen.add(line_without_id)
						elif int(split_line[3]) > cut_start and int(split_line[4]) < cut_end:
							adjusted_gff.write("{}\t{}".format(seq_id, line_without_id))
							xminimal.append(split_line[3])
							xmaximal.append(split_line[4])
							dom.append(split_line[-1].split(",")[0].split("=")[1])	
				
				## not divived
				else:
					if seq_id != seq_id_all[-1]:
						seq_id_all.append(seq_id)
						xminimal_all.append(xminimal)
						xmaximal_all.append(xmaximal)
						dom_all.append(dom)
						xminimal = []
						xmaximal = []
						dom = []
					adjusted_gff.write(line)
					xminimal.append(split_line[3])
					xmaximal.append(split_line[4])
					dom.append(split_line[-1].split(",")[0].split("=")[1])			
	xminimal_all.append(xminimal)
	xmaximal_all.append(xmaximal)
	dom_all.append(dom)
	os.remove(OUTPUT_DOMAIN)
	os.rename(adjusted_file, OUTPUT_DOMAIN)
	return xminimal_all, xmaximal_all, dom_all, seq_id_all 		
	
	
def main(args):
	
	t = time.time()
	
	QUERY = args.query
	LAST_DB = args.protein_database
	CLASSIFICATION = args.classification
	OUTPUT_DOMAIN = args.domain_gff
	NEW_LDB = args.new_ldb
	SUMMARY = args.summary_file
	OUTPUT_DIR = args.output_dir
	THRESHOLD_SCORE = args.threshold_score
	WIN_DOM = args.win_dom
	OVERLAP_DOM = args.overlap_dom
	
	if OUTPUT_DIR is None:
		OUTPUT_DIR = configuration.TMP
	if SUMMARY is None:
		SUMMARY = configuration.DOM_SUMMARY
	if OUTPUT_DOMAIN is None:
		OUTPUT_DOMAIN = configuration.DOMAINS_GFF
	
	
	if NEW_LDB:
		subprocess.call("lastdb -p -cR01 {} {}".format(LAST_DB, LAST_DB), shell=True)
		
	if not os.path.exists(OUTPUT_DIR) and not os.path.isabs(OUTPUT_DOMAIN):
		os.makedirs(OUTPUT_DIR)
		OUTPUT_DOMAIN = os.path.join(OUTPUT_DIR, os.path.basename(OUTPUT_DOMAIN))
		SUMMARY = os.path.join(OUTPUT_DIR, os.path.basename(SUMMARY))
	elif os.path.exists(OUTPUT_DIR) and not os.path.isabs(OUTPUT_DOMAIN):
		OUTPUT_DOMAIN = os.path.join(OUTPUT_DIR, os.path.basename(OUTPUT_DOMAIN))
		SUMMARY = os.path.join(OUTPUT_DIR, os.path.basename(SUMMARY))
	
	if not os.path.exists(LAST_DB):
		CLASSIFICATION = os.path.join(configuration.TOOL_DATA_DIR, CLASSIFICATION)
		LAST_DB = os.path.join(configuration.TOOL_DATA_DIR, LAST_DB) 
	

	
	[xminimal, xmaximal, domains_all, seq_ids] = domain_search(QUERY, LAST_DB, CLASSIFICATION, OUTPUT_DOMAIN, THRESHOLD_SCORE, WIN_DOM, OVERLAP_DOM)
	
	domains_stat(domains_all, seq_ids, SUMMARY)
	
	
	print("ELAPSED_TIME_DOMAINS = {} s".format(time.time() - t))
	
if __name__ == "__main__":
	import argparse
	from argparse import RawDescriptionHelpFormatter
	
	class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
		pass

	parser = argparse.ArgumentParser(
		description='''Script performs scanning of given DNA sequence(s) in (multi)fasta format in order to discover protein domains using our protein domains database. Domains are subsequently annotated and classified - in case certain domain has multiple annotations assigned, classifation is derived from the common classification level of all of them. Domains searching is accomplished engaging LASTAL alignment tool.
		
	DEPENDANCIES:
		- python 3.4 or higher
		- python3-numpy package
		- lastal 744 or higher [http://last.cbrc.jp/]
		- configuration.py module

	EXAMPLE OF USAGE:
		
		./protein_domains_pd.py -q PATH_TO_INPUT_SEQ -pdb PATH_TO_PROTEIN_DB -cs PATH_TO_CLASSIFICATION_FILE
		
	When running for the first time with a new database use -nld option allowing lastal to create indexed database files:
		
		-nld True
	
		''',
		epilog="""""",
		formatter_class=CustomFormatter)
	requiredNamed = parser.add_argument_group('required named arguments')
	requiredNamed.add_argument("-q","--query", type=str, required=True,
						help='input DNA sequence to search for protein domains in a fasta format. Multifasta format allowed.')
	requiredNamed.add_argument('-pdb', "--protein_database", type=str, required=True, 
                        help='protein domains database file')
	requiredNamed.add_argument('-cs', '--classification', type=str, required=True, 
                        help='protein domains classification file')
	parser.add_argument("-oug", "--domain_gff", type=str,
						help="output domains gff format")
	parser.add_argument("-nld","--new_ldb", type=str, default=False,
						help="create indexed database files for lastal in case of working with new protein db")
	parser.add_argument("-sum","--summary_file", type=str,
						help=" output summary file containing overview of amount of domains in individual seqs")
	parser.add_argument("-dir","--output_dir", type=str,
						help="specify if you want to change the output directory")
	parser.add_argument("-thsc","--threshold_score", type=int, default=80,
						help="percentage of the best score in the cluster to be tolerated when assigning annotations per base")
	parser.add_argument("-wd","--win_dom", type=int, default=10000000,
						help="window to process large input sequences sequentially")
	parser.add_argument("-od","--overlap_dom", type=int, default=10000,
						help="overlap of sequences in two consecutive windows")
	
	args = parser.parse_args()
	main(args)
