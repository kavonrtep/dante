#!/usr/bin/env python3

import numpy as np
import subprocess
import csv
import time
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import multiprocessing
import argparse
import os
from functools import partial
from multiprocessing import Pool
from tempfile import NamedTemporaryFile
from operator import itemgetter
from itertools import groupby
import gff
import protein_domains
import domains_filtering
import configuration
import visualization
import distutils
from distutils import dir_util
import tempfile
import re
from Bio import SeqIO
import sys
import pickle
import shutil
import warnings

t_profrep = time.time()
np.set_printoptions(threshold=np.nan)
warnings.filterwarnings("ignore", module="matplotlib")

def get_version(path):
	branch = subprocess.check_output("git rev-parse --abbrev-ref HEAD", shell=True, cwd=path).decode('ascii').strip()
	shorthash = subprocess.check_output("git log --pretty=format:'%h' -n 1  ", shell=True, cwd=path).decode('ascii').strip()
	revcount = len(subprocess.check_output("git log --oneline", shell=True,  cwd=path).decode('ascii').split())
	version_string = ("-------------------------------------"
		"-------------------------------------\n"
							  "PIPELINE VERSION         : "
		"{branch}-rv-{revcount}({shorthash})\n"
		"-------------------------------------"
		"-------------------------------------\n").format(
								  branch=branch,
								  shorthash=shorthash,
								  revcount=revcount,
                      )
	return(version_string)


def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected')


def check_fasta_id(QUERY):
	forbidden_ids = []
	headers = []
	for record in SeqIO.parse(QUERY, "fasta"):
		if any(x in record.id for x in configuration.FORBIDDEN_CHARS):
		 forbidden_ids.append(record.id)
		headers.append(record.id)
	if len(headers) > len(set([header.split(" ")[0] for header in headers])):
		raise NameError('''Sequences in multifasta format are not named correctly:
							seq IDs(before the first space) are the same''')
	return forbidden_ids, headers


def multifasta(QUERY):
	''' Create single fasta temporary files to be processed sequentially '''
	PATTERN = ">"
	fasta_list = []
	with open(QUERY, "r") as fasta:
		reader = fasta.read()
		splitter = reader.split(PATTERN)[1:]
		for fasta_num, part in enumerate(splitter):
			ntf = NamedTemporaryFile(delete=False)
			ntf.write("{}{}".format(PATTERN, part).encode("utf-8"))
			fasta_list.append(ntf.name)
			ntf.close()
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


def cluster_annotation(CL_ANNOTATION_TBL):
	''' Create dictionary of known annotations classes and related clusters '''
	cl_annotations = {} 			
	annot_table = np.genfromtxt(CL_ANNOTATION_TBL, dtype=str) 
	for line in annot_table:
		if line[1] in cl_annotations:
			cl_annotations[line[1]].append(line[0])
		else:
			cl_annotations[line[1]] = [line[0]]
	return list(cl_annotations.items()), list(cl_annotations.keys())


def read_annotation(CLS, cl_annotations_items):
	''' Dictionary of known repeat classes and related reads '''
	reads_annotations = {} 	
	with open(CLS, "r") as cls_file:
		count = 0
		for line in cls_file:
			line = line.rstrip()
			count += 1
			if count%2 == 0:
				reads = re.split("\s+", line) 
				for element in reads: 
					for key, value in cl_annotations_items:
						if clust in value:
							reads_annotations[element] = key 
			else:
				clust = re.split("\s+", line)[0].split(">CL")[1]
	return reads_annotations


def annot_profile(annotation_keys, part):
	''' Predefine dictionary of known annotations and partial sequence 
	repetitive profiles defined by parallel process '''
	subprofile = {} 				
	for key in annotation_keys:
		subprofile[key] = [np.zeros(part, dtype=int), np.zeros(part, dtype=int)]
	subprofile["ALL"] = [np.zeros(part, dtype=int), np.zeros(part, dtype=int)]
	return subprofile


def parallel_process(WINDOW, OVERLAP, seq_length, annotation_keys, reads_annotations, subfasta, BLAST_DB, E_VALUE, WORD_SIZE, BLAST_TASK, MAX_ALIGNMENTS, BITSCORE, DUST_FILTER, last_index, subsets_num, subset_index):
	''' Run parallel function to process the input sequence in windows
		Run blast for subsequence defined by the input index and window size
 		Create and increment subprofile vector based on reads aligned within window '''
	loc_start = subset_index + 1
	loc_end = subset_index + WINDOW
	if loc_end > seq_length:
		loc_end = seq_length
		subprofile = annot_profile(annotation_keys, seq_length - loc_start + 1)
	else:
		subprofile = annot_profile(annotation_keys, WINDOW + 1)
	
	# Find HSP records for every window defined by query location and parse the tabular stdout:
	# 1. query, 2. database read, 3. alignment start, 4. alignment end, 5. bitscore
	p = subprocess.Popen("blastn -query {} -query_loc {}-{} -db {} -evalue {} -word_size {} -dust {} -task {} -num_alignments {} -outfmt '6 qseqid sseqid qstart qend bitscore pident'".format(subfasta, loc_start, loc_end, BLAST_DB, E_VALUE, WORD_SIZE, DUST_FILTER, BLAST_TASK, MAX_ALIGNMENTS), stdout=subprocess.PIPE, shell=True)
	count_hits = 0
	for line in p.stdout:
		column = line.decode("utf-8").rstrip().split("\t")
		if float(column[4]) >= BITSCORE:
			count_hits += 1
			read = column[1]				# ID of individual aligned read
			if "reduce" in read:
				reads_representation = int(read.split("reduce")[-1])
			else:
				reads_representation = 1
			qstart = int(column[2])			# starting position of alignment
			qend = int(column[3])			# ending position of alignemnt
			if read in reads_annotations:
				annotation = reads_annotations[read]								
			else:
				annotation = "ALL"
			subprofile[annotation][0][qstart-subset_index-1 : qend-subset_index] = subprofile[annotation][0][qstart-subset_index-1 : qend-subset_index] + reads_representation
			subprofile[annotation][1][qstart-subset_index-1 : qend-subset_index] = subprofile[annotation][1][qstart-subset_index-1 : qend-subset_index] + float(column[5])*reads_representation
	subprofile["ALL"][0] = sum([item[0] for item in subprofile.values()])
	subprofile["ALL"][1] = sum([item[1] for item in subprofile.values()])
	for repeat in subprofile.keys():
		subprofile[repeat][1] = [int(round(quality/hits_num)) if hits_num !=0 else quality for hits_num, quality in zip(subprofile[repeat][0], subprofile[repeat][1])]
	print(subprofile[repeat][1])
	if subset_index == 0: 
		if subsets_num == 1:
			subprf_name = subprofile_single(subprofile, subset_index)
		else:
			subprf_name = subprofile_first(subprofile, subset_index, WINDOW, OVERLAP)
	elif subset_index == last_index:
		subprf_name = subprofile_last(subprofile, subset_index, OVERLAP)
	else:
		subprf_name = subprofiles_middle(subprofile, subset_index, WINDOW, OVERLAP)
	return subprf_name
	
	
def subprofile_single(subprofile, subset_index):
	subprofile['idx'] = list(range(1, len(subprofile["ALL"][0]) + 1))
	subprf_dict = NamedTemporaryFile(suffix='{}_.pickle'.format(subset_index),delete=False)
	with open(subprf_dict.name, 'wb') as handle:
		pickle.dump(subprofile, handle, protocol=pickle.HIGHEST_PROTOCOL)
	subprf_dict.close()
	return subprf_dict.name


def subprofile_first(subprofile, subset_index, WINDOW, OVERLAP):
	for key in subprofile.keys():
		subprofile[key][0] = subprofile[key][0][0 : -OVERLAP//2-1]
		subprofile[key][1] = subprofile[key][1][0 : -OVERLAP//2-1]
	subprofile['idx'] = list(range(subset_index + 1, subset_index + WINDOW-OVERLAP//2 + 1))
	subprf_dict = NamedTemporaryFile(suffix='{}_.pickle'.format(subset_index),delete=False)
	with open(subprf_dict.name, 'wb') as handle:
		pickle.dump(subprofile, handle, protocol=pickle.HIGHEST_PROTOCOL)
	subprf_dict.close()
	return subprf_dict.name
	
	
def subprofiles_middle(subprofile, subset_index, WINDOW, OVERLAP):
	for key in subprofile.keys():
		subprofile[key][0] = subprofile[key][0][OVERLAP//2 : -OVERLAP//2-1]
		subprofile[key][1] = subprofile[key][1][OVERLAP//2 : -OVERLAP//2-1]
	subprofile['idx'] = list(range(subset_index + OVERLAP//2 + 1, subset_index + WINDOW-OVERLAP//2 + 1))
	subprf_dict = NamedTemporaryFile(suffix='{}_.pickle'.format(subset_index),delete=False)
	with open(subprf_dict.name, 'wb') as handle:
		pickle.dump(subprofile, handle, protocol=pickle.HIGHEST_PROTOCOL)
	subprf_dict.close()
	return subprf_dict.name


def subprofile_last(subprofile, subset_index, OVERLAP):
	len_subprofile = len(subprofile['ALL'][0])
	for key in subprofile.keys():
		subprofile[key][0] = subprofile[key][0][OVERLAP//2:]
		subprofile[key][1] = subprofile[key][1][OVERLAP//2:]
	subprofile['idx'] = list(range(subset_index + OVERLAP//2 + 1, subset_index + len_subprofile +1))
	subprf_dict = NamedTemporaryFile(suffix='{}_.pickle'.format(subset_index),delete=False)
	with open(subprf_dict.name, 'wb') as handle:
		pickle.dump(subprofile, handle, protocol=pickle.HIGHEST_PROTOCOL)
	subprf_dict.close()
	return subprf_dict.name
	
	
def concatenate_prof(subprofiles_all, files_dict, seq_id, HTML_DATA):
	for subprofile in subprofiles_all:
		with open(subprofile, 'rb') as handle:
			individual_dict = pickle.load(handle)
			exclude = set(["idx"])
			for key in set(individual_dict.keys()).difference(exclude):
				if any(individual_dict[key][0]):
					indices = handle_zero_lines(individual_dict[key][0])
					if key not in files_dict.keys():
						prf_name = "{}/{}.wig".format(HTML_DATA, re.sub('[\/\|]','_',key))
						prf_qual_name = "{}/{}_qual.wig".format(HTML_DATA, re.sub('[\/\|]','_',key))
						with open(prf_name, "a") as prf_file,  open(prf_qual_name, "a") as prf_q_file:
							prf_file.write("{}{}\n".format(configuration.HEADER_WIG, seq_id))
							prf_q_file.write("{}{}\n".format(configuration.HEADER_WIG, seq_id))
							for i in indices:
								prf_file.write("{}\t{}\n".format(individual_dict['idx'][i], individual_dict[key][0][i]))
								prf_q_file.write("{}\t{}\n".format(individual_dict['idx'][i], int(individual_dict[key][1][i])))								
						files_dict[key] = [prf_name,[seq_id], prf_qual_name]
					else:
						prf_name = files_dict[key][0]
						prf_qual_name = files_dict[key][2]
						with open(prf_name, "a") as prf_file,  open(prf_qual_name, "a") as prf_q_file:
								if seq_id not in files_dict[key][1]:
									prf_file.write("{}{}\n".format(configuration.HEADER_WIG, seq_id))
									prf_q_file.write("{}{}\n".format(configuration.HEADER_WIG, seq_id))
									files_dict[key][1].append(seq_id)
								for i in indices:
									prf_file.write("{}\t{}\n".format(individual_dict['idx'][i], individual_dict[key][0][i]))
									prf_q_file.write("{}\t{}\n".format(individual_dict['idx'][i], int(individual_dict[key][1][i])))
	return files_dict


def concatenate_prof_CN(CV, subprofiles_all, files_dict, seq_id, HTML_DATA):
	for subprofile in subprofiles_all:
		with open(subprofile, 'rb') as handle:
			individual_dict = pickle.load(handle)
			exclude = set(["idx"])
			for key in set(individual_dict.keys()).difference(exclude):
				if any(individual_dict[key][0]):
					indices = handle_zero_lines(individual_dict[key][0])
					if key not in files_dict.keys():
						prf_name = "{}/{}.wig".format(HTML_DATA, re.sub('[\/\|]','_',key))
						prf_qual_name = "{}/{}_qual.wig".format(HTML_DATA, re.sub('[\/\|]','_',key))
						with open(prf_name, "a") as prf_file,  open(prf_qual_name, "a") as prf_q_file:
							prf_file.write("{}{}\n".format(configuration.HEADER_WIG, seq_id))
							prf_q_file.write("{}{}\n".format(configuration.HEADER_WIG, seq_id))
							for i in indices:
								prf_file.write("{}\t{}\n".format(individual_dict['idx'][i], int(individual_dict[key][0][i]/CV)))
								prf_q_file.write("{}\t{}\n".format(individual_dict['idx'][i], int(individual_dict[key][1][i])))
						files_dict[key] = [prf_name,[seq_id], prf_qual_name]
					else:
						prf_name = files_dict[key][0]
						prf_qual_name = files_dict[key][2]
						with open(prf_name, "a") as prf_file,  open(prf_qual_name, "a") as prf_q_file:
								if seq_id not in files_dict[key][1]:
									prf_file.write("{}{}\n".format(configuration.HEADER_WIG, seq_id))
									prf_q_file.write("{}{}\n".format(configuration.HEADER_WIG, seq_id))
									files_dict[key][1].append(seq_id)
								for i in indices:
									prf_file.write("{}\t{}\n".format(individual_dict['idx'][i], int(individual_dict[key][0][i]/CV)))
									prf_q_file.write("{}\t{}\n".format(individual_dict['idx'][i], int(individual_dict[key][1][i])))
	return files_dict
		

def handle_zero_lines(repeat_subhits):
	''' Clean lines which contains only zeros, i.e. positons which do not contain any hit. However border zero positions need to be preserved due to correct graphs plotting '''
	zero_idx = [idx for idx, val in enumerate(repeat_subhits) if val == 0]
	indices = [idx for idx, val in enumerate(repeat_subhits) if val != 0]
	zero_breakpoints = []
	for key, group in groupby(enumerate(zero_idx), lambda index_item: index_item[0] - index_item[1]):
		group = list(map(itemgetter(1),group))
		zero_breakpoints.append(group[0])
		zero_breakpoints.append(group[-1])
	if indices:
		indices.extend(zero_breakpoints)
		indices = sorted(set(indices), key=int)	
	else:
		indices = []
	return indices
	

def repeats_process_dom(OUTPUT_GFF, THRESHOLD, THRESHOLD_SEGMENT, HTML_DATA, xminimal, xmaximal, domains, seq_ids_dom, CN, seq_ids_all, seq_lengths_all, files_dict):
	''' Process the hits table separately for each fasta, create gff file and profile picture '''
	if files_dict:
		gff.create_gff(THRESHOLD, THRESHOLD_SEGMENT, OUTPUT_GFF, files_dict, seq_ids_all)
	else:
		with open(OUTPUT_GFF, "w") as gff_file:
			gff_file.write("{}\n".format(configuration.HEADER_GFF))
	seqs_all_part = seq_ids_all[0:configuration.MAX_PIC_NUM]
	graphs_dict = {}
	seqs_long = []
	if files_dict:	
		[graphs_dict, seqs_long] = visualization.vis_profrep(seq_ids_all, files_dict, seq_lengths_all, CN, HTML_DATA, seqs_all_part)
	count_seq = 0
	for seq in seqs_all_part:
		if seq in graphs_dict.keys():
			fig = graphs_dict[seq][0]
			ax = graphs_dict[seq][1]
			art = []
			lgd = ax.legend(bbox_to_anchor=(0.5,-0.1), loc=9, ncol=3)
			art.append(lgd)
			if seq in seq_ids_dom:
				dom_idx = seq_ids_dom.index(seq) 
				[fig, ax] = visualization.vis_domains(fig, ax, seq, xminimal[dom_idx], xmaximal[dom_idx], domains[dom_idx])
		elif seq in seqs_long:
			[fig, ax] = visualization.plot_figure(seq, seq_lengths_all[count_seq], CN)
			ax.text(0.3, 0.5, "Graphs are only displayed if sequence is not longer than {} bp".format(configuration.SEQ_LEN_VIZ),transform=ax.transAxes, fontsize=14, verticalalignment='center', color='blue')
		else:
			[fig, ax] = visualization.plot_figure(seq, seq_lengths_all[count_seq], CN)
			ax.hlines(0, 0, seq_lengths_all[count_seq], color="red", lw=4)
			if seq in seq_ids_dom:
				dom_idx = seq_ids_dom.index(seq) 
				[fig, ax] = visualization.vis_domains(fig, ax, seq, xminimal[dom_idx], xmaximal[dom_idx], domains[dom_idx])
		output_pic_png = "{}/{}.png".format(HTML_DATA, count_seq)
		fig.savefig(output_pic_png, bbox_inches="tight", format="png", dpi=configuration.IMAGE_RES)
		count_seq += 1	
	return None
	
	
def repeats_process(OUTPUT_GFF, THRESHOLD, THRESHOLD_SEGMENT, HTML_DATA, CN, seq_ids_all, seq_lengths_all, files_dict):
	''' Process the hits table separately for each fasta, create gff file and profile picture '''
	if files_dict:
		gff.create_gff(THRESHOLD, THRESHOLD_SEGMENT, OUTPUT_GFF, files_dict, seq_ids_all)
	else:
		with open(OUTPUT_GFF, "w") as gff_file:
			gff_file.write("{}\n".format(configuration.HEADER_GFF))
	seqs_all_part = seq_ids_all[0:configuration.MAX_PIC_NUM]
	graphs_dict = {}
	seqs_long = []
	if files_dict:	
		[graphs_dict, seqs_long] = visualization.vis_profrep(seq_ids_all, files_dict, seq_lengths_all, CN, HTML_DATA, seqs_all_part)
	count_seq = 0
	for seq in seqs_all_part:
		if seq in graphs_dict.keys():
			fig = graphs_dict[seq][0]
			ax = graphs_dict[seq][1]
			art = []
			lgd = ax.legend(bbox_to_anchor=(0.5,-0.1), loc=9, ncol=3)
			art.append(lgd)
		elif seq in seqs_long:
			[fig, ax] = visualization.plot_figure(seq, seq_lengths_all[count_seq], CN)
			ax.text(0.3, 0.5, "Graphs are only displayed if sequence is not longer than {} bp".format(configuration.SEQ_LEN_VIZ),transform=ax.transAxes, fontsize=14, verticalalignment='center', color='blue')
		else:
			[fig, ax] = visualization.plot_figure(seq, seq_lengths_all[count_seq], CN)
			ax.hlines(0, 0, seq_lengths_all[count_seq], color="red", lw=4)
		output_pic_png = "{}/{}.png".format(HTML_DATA, count_seq)
		fig.savefig(output_pic_png, bbox_inches="tight", format="png", dpi=configuration.IMAGE_RES)	
		plt.close()
		count_seq += 1	
	return None
	

def html_output(total_length, seq_lengths_all, seq_names, HTML, DB_NAME, REF, REF_LINK):
	''' Define html output with limited number of output pictures and link to JBrowse '''
	info = "\t\t".join(['<pre> {} [{} bp]</pre>'.format(seq_name, seq_length) for seq_name, seq_length in zip(seq_names, seq_lengths_all)])	
	if REF:
		ref_part_1 = REF.split("-")[0]
		ref_part_2 = "-".join(REF.split("-")[1:]).split(". ")[0]
		ref_part_3 = ". ".join("-".join(REF.split("-")[1:]).split(". ")[1:])
		ref_string = '''<h6> {} - <a href="{}" target="_blank" >{}</a>. {}'''.format(ref_part_1, REF_LINK, ref_part_2, ref_part_3)
	else:
		ref_string = "Custom Data"
	pictures = "\n\t\t".join(['<img src="{}.png" width=1800>'.format(pic)for pic in range(len(seq_names))[:configuration.MAX_PIC_NUM]])
	html_str = configuration.HTML_STR.format(info, total_length, DB_NAME, pictures, ref_string)
	with open(HTML,"w") as html_file:
		html_file.write(html_str)
		
				
def adjust_tracklist(jbrowse_data_path):
	starting_lines = []
	ending_lines = []
	end = False
	with open(os.path.join(jbrowse_data_path, "trackList.json"), "r") as track_list:
		for line in track_list:
			if "]" not in line and not end:
				starting_lines.append(line)
			else:
				end = True
				ending_lines.append(line)
	with open(os.path.join(jbrowse_data_path, "trackList.json"), "w") as track_list:
		for line in starting_lines:
			track_list.write(line)
	return ending_lines
	
				
def jbrowse_prep_dom(HTML_DATA, QUERY, OUT_DOMAIN_GFF, OUTPUT_GFF, N_GFF, total_length, JBROWSE_BIN, files_dict):
	''' Set up the paths, link and convert output data to be displayed as tracks in Jbrowse '''
	jbrowse_data_path = os.path.join(HTML_DATA, configuration.jbrowse_data_dir)
	with tempfile.TemporaryDirectory() as dirpath:
		subprocess.call(["{}/prepare-refseqs.pl".format(JBROWSE_BIN), "--fasta", QUERY, "--out", jbrowse_data_path])
		subprocess.call(["{}/flatfile-to-json.pl".format(JBROWSE_BIN), "--gff", OUT_DOMAIN_GFF, "--trackLabel", "GFF_domains", "--out",  jbrowse_data_path])
		subprocess.call(["{}/flatfile-to-json.pl".format(JBROWSE_BIN), "--gff", OUTPUT_GFF, "--trackLabel", "GFF_repeats", "--config", configuration.JSON_CONF_R, "--out",  jbrowse_data_path])
		subprocess.call(["{}/flatfile-to-json.pl".format(JBROWSE_BIN), "--gff", N_GFF, "--trackLabel", "N_regions", "--config", configuration.JSON_CONF_N, "--out",  jbrowse_data_path])		 
		count = 0
		# Control the total length processed, if above threshold, dont create wig image tracks 
		if files_dict:
			exclude = set(['ALL'])
			sorted_keys =  sorted(set(files_dict.keys()).difference(exclude))
			sorted_keys.insert(0, "ALL")
			ending_lines = adjust_tracklist(jbrowse_data_path)
			track_list = open(os.path.join(jbrowse_data_path, "trackList.json"), "a")
			for repeat_id in sorted_keys:
				color = configuration.COLORS_HEX[count]
				count += 1
				bw_name = "{}.bw".format(re.sub('[\/\|]','_',repeat_id))
				subprocess.call(["wigToBigWig", files_dict[repeat_id][0], os.path.join(HTML_DATA, configuration.CHROM_SIZES_FILE), os.path.join(jbrowse_data_path, bw_name)])
				track_list.write(configuration.TRACK_LIST.format("{", bw_name, repeat_id, repeat_id, "{", color, "}", "}"))
			for line in ending_lines:
				track_list.write(line)
		distutils.dir_util.copy_tree(dirpath,jbrowse_data_path)
	return None
	
	
def jbrowse_prep(HTML_DATA, QUERY, OUTPUT_GFF, N_GFF, total_length, JBROWSE_BIN, files_dict):
	''' Set up the paths, link and convert output data to be displayed as tracks in Jbrowse '''
	jbrowse_data_path = os.path.join(HTML_DATA, configuration.jbrowse_data_dir)
	with tempfile.TemporaryDirectory() as dirpath:
		subprocess.call(["{}/prepare-refseqs.pl".format(JBROWSE_BIN), "--fasta", QUERY, "--out", jbrowse_data_path])
		subprocess.call(["{}/flatfile-to-json.pl".format(JBROWSE_BIN), "--gff", OUTPUT_GFF, "--trackLabel", "GFF_repeats", "--config", configuration.JSON_CONF_R, "--out",  jbrowse_data_path])
		subprocess.call(["{}/flatfile-to-json.pl".format(JBROWSE_BIN), "--gff", N_GFF, "--trackLabel", "N_regions", "--config", configuration.JSON_CONF_N, "--out",  jbrowse_data_path])		 
		count = 0
		## Control the total length processed, if above threshold, dont create wig image tracks 
		if files_dict:
			exclude = set(['ALL'])
			sorted_keys =  sorted(set(files_dict.keys()).difference(exclude))
			sorted_keys.insert(0, "ALL")
			ending_lines = adjust_tracklist(jbrowse_data_path)
			track_list = open(os.path.join(jbrowse_data_path, "trackList.json"), "a")
			for repeat_id in sorted_keys:
				color = configuration.COLORS_HEX[count]
				count += 1
				bw_name = "{}.bw".format(re.sub('[\/\|]','_',repeat_id))
				subprocess.call(["wigToBigWig", files_dict[repeat_id][0], os.path.join(HTML_DATA, configuration.CHROM_SIZES_FILE), os.path.join(jbrowse_data_path, bw_name)])
				track_list.write(configuration.TRACK_LIST.format("{", bw_name, repeat_id, repeat_id, "{", color, "}", "}"))
			for line in ending_lines:
					track_list.write(line)
			track_list.close()
		distutils.dir_util.copy_tree(dirpath,jbrowse_data_path)
	return None
	
	
def genome2coverage(GS, BLAST_DB):
	''' Convert genome size to coverage '''
	num_of_reads = 0
	with open(BLAST_DB) as reads_all:
		first_line = reads_all.readline()
		if first_line.startswith(">"):
			num_of_reads += 1
			first_seq = reads_all.readline().rstrip()
		for line in reads_all:
			if line.startswith(">"):
				num_of_reads += 1
	len_of_read = len(first_seq)
	CV =  (num_of_reads *  len_of_read)/(GS*1000000) # GS in Mb
	return CV

	
def prepared_data(TBL, DB_ID, TOOL_DATA_DIR):
	''' Get prepared rep. annotation data from the table based on the selected species ID '''
	with open(TBL) as datasets:
		for line in datasets:
			if line.split("\t")[0] == DB_ID:
				DB_NAME = line.split("\t")[1]
				BLAST_DB = os.path.join(TOOL_DATA_DIR, line.split("\t")[2])
				print(BLAST_DB)
				CLS = os.path.join(TOOL_DATA_DIR, line.split("\t")[3])
				CL_ANNOTATION_TBL = os.path.join(TOOL_DATA_DIR, line.split("\t")[4])
				CV = float(line.split("\t")[5])
				REF = line.split("\t")[6]
				REF_LINK = line.split("\t")[7]
	return DB_NAME, BLAST_DB, CLS, CL_ANNOTATION_TBL, CV, REF, REF_LINK

def seq_sizes_file(seq_ids, seq_lengths_all, HTML_DATA):
	chrom_sizes = os.path.join(HTML_DATA, configuration.CHROM_SIZES_FILE)
	with open(chrom_sizes, "w") as chroms:
		for seq_id, seq_length in zip(seq_ids, seq_lengths_all):
			chroms.write("{}\t{}\n".format(seq_id, seq_length)) 
	
def main(args):
	## Command line arguments
	QUERY = args.query
	BLAST_DB = args.database
	CL_ANNOTATION_TBL = args.annotation_tbl 
	CLS = args.cls
	BITSCORE = args.bit_score
	E_VALUE = args.e_value
	WORD_SIZE = args.word_size
	WINDOW = args.window
	OVERLAP = args.overlap
	BLAST_TASK = args.task
	MAX_ALIGNMENTS = args.max_alignments
	NEW_DB = args.new_db
	THRESHOLD = args.threshold_repeat
	THRESHOLD_SEGMENT = args.threshold_segment
	OUTPUT_GFF = args.output_gff
	DOMAINS = args.protein_domains
	LAST_DB = args.protein_database
	CLASSIFICATION = args.classification
	OUT_DOMAIN_GFF = args.domain_gff	
	HTML = args.html_file
	HTML_DATA = args.html_path
	N_GFF = args.n_gff
	CN = args.copy_numbers
	GS = args.genome_size
	DB_ID = args.db_id
	TBL = args.datasets_tbl
	THRESHOLD_SCORE = args.threshold_score
	WIN_DOM = args.win_dom
	OVERLAP_DOM = args.overlap_dom
	TOOL_DATA_DIR = args.tool_dir
	JBROWSE_BIN = args.jbrowse_bin
	DUST_FILTER = args.dust_filter
	LOG_FILE = args.log_file
	
	
	if not JBROWSE_BIN: 
		try:
			JBROWSE_BIN = os.environ['JBROWSE_BIN']
		except KeyError:
			raise ValueError('There was no path to JBrowse bin found - set the enviroment variable JBROWSE_BIN or pass the argument explicitly')
	
	
	if CN and not DB_ID and not GS:
		raise ValueError("Genome size missing - if you want to convert hits to copy numbers please enter --genome_size parameter")


	## Check if there are forbidden characters in fasta IDs 
	[forbidden_ids, headers] = check_fasta_id(QUERY)
	if forbidden_ids:
		##################### USER ERROR ###############################
		raise UserWarning("The following IDs contain forbidden characters ('/' or '\\') - PLEASE REPLACE OR DELETE THEM:\n{}".format("\n".join(forbidden_ids)))
		
	if len(headers) > len(set([header.split(" ")[0] for header in headers])):
		raise NameError('''Sequences in multifasta format are not named correctly:
							seq IDs(before the first space) are the same''')

	
	## Create new blast database of reads
	if NEW_DB:
		subprocess.call("makeblastdb -in {} -dbtype nucl".format(BLAST_DB), shell=True)
	
	## Parse prepared annotation data table
	if TBL:
		[DB_NAME, BLAST_DB, CLS, CL_ANNOTATION_TBL, CV, REF, REF_LINK] = prepared_data(TBL, DB_ID, TOOL_DATA_DIR)
	else:
		REF = None
		REF_LINK = None
		DB_NAME = "CUSTOM"
		
	if TOOL_DATA_DIR:
		LAST_DB = os.path.join(LAST_DB, configuration.LAST_DB_FILE)
		CLASSIFICATION = os.path.join(CLASSIFICATION, configuration.CLASS_FILE)
	
	
	## Create dir to store outputs for html 
	if not os.path.exists(HTML_DATA):
		os.makedirs(HTML_DATA)
	
	if not os.path.isabs(HTML):
		HTML = os.path.join(HTML_DATA, HTML)
		
	if not os.path.isabs(OUT_DOMAIN_GFF):
		OUT_DOMAIN_GFF = os.path.join(HTML_DATA, OUT_DOMAIN_GFF)
	
	if not os.path.isabs(LOG_FILE):
		LOG_FILE = os.path.join(HTML_DATA, LOG_FILE)

	if not os.path.isabs(N_GFF):
		N_GFF = os.path.join(HTML_DATA, N_GFF)
		
	if not os.path.isabs(OUTPUT_GFF):
		OUTPUT_GFF = os.path.join(HTML_DATA, OUTPUT_GFF)

	
	path = os.path.dirname(os.path.realpath(__file__))
	version_string = get_version(path)
	
	log = os.open(LOG_FILE, os.O_RDWR|os.O_CREAT)
	
	os.write(log, version_string.encode("utf-8"))
	
	## Define parameters for parallel process
	STEP = WINDOW - OVERLAP		
	NUM_CORES = multiprocessing.cpu_count()	
	os.write(log, "NUM_OF_CORES = {}\n".format(NUM_CORES).encode("utf-8"))
	
	## Convert genome size to coverage
	if CN and GS:
		CV = genome2coverage(GS, BLAST_DB)
		os.write(log, "COVERAGE = {}\n".format(CV).encode("utf-8"))
	
	parallel_pool = Pool(NUM_CORES)

	## Assign clusters to repetitive classes
	[cl_annotations_items, annotation_keys] = cluster_annotation(CL_ANNOTATION_TBL)
	
	## Assign reads to repetitive classes
	reads_annotations = read_annotation(CLS, cl_annotations_items)
	
	## Detect all fasta sequences from input
	fasta_list = multifasta(QUERY)
	headers=[]
	files_dict = {}
	seq_count = 1
	start = 1
	total_length = 0
	seq_lengths_all = [] 
	Ngff = open(N_GFF,"w")
	Ngff.write("{}\n".format(configuration.HEADER_GFF))
	## Find hits for each fasta sequence separetely
	t_blast=time.time()	
	for subfasta in fasta_list:
		[header, sequence] = fasta_read(subfasta)
		os.write(log, "Sequence {} is being processed...\n".format(header).encode("utf-8"))
		os.fsync(log)
		indices_N = [indices + 1 for indices, n in enumerate(sequence) if n == "n" or n == "N"]
		if indices_N:
			gff.idx_ranges_N(indices_N, configuration.N_segment, header, Ngff, configuration.N_NAME, configuration.N_FEATURE)
		seq_length = len(sequence)
		headers.append(header)
		## Create parallel process																							
		subset_index = list(range(0, seq_length, STEP))
		## Situation when penultimal window is not complete but it is following by another one
		if len(subset_index) > 1 and subset_index[-2] + WINDOW >= seq_length:
			subset_index = subset_index[:-1]	
		last_index = subset_index[-1]
		index_range = range(len(subset_index))
		for chunk_index in index_range[0::configuration.MAX_FILES_SUBPROFILES]:
			multiple_param = partial(parallel_process, WINDOW, OVERLAP, seq_length, annotation_keys, reads_annotations, subfasta, BLAST_DB, E_VALUE, WORD_SIZE, BLAST_TASK, MAX_ALIGNMENTS, BITSCORE, DUST_FILTER, last_index, len(subset_index))
			subprofiles_all = parallel_pool.map(multiple_param, subset_index[chunk_index:chunk_index + configuration.MAX_FILES_SUBPROFILES])
			## Join partial profiles to the final profile of the sequence 
			if CN:							
				files_dict = concatenate_prof_CN(CV, subprofiles_all, files_dict, header, HTML_DATA)
			else:
				files_dict = concatenate_prof(subprofiles_all, files_dict, header, HTML_DATA)
			for subprofile in subprofiles_all:
				os.unlink(subprofile)
		total_length += seq_length 
		seq_lengths_all.append(seq_length)
	Ngff.close()
	os.write(log, "ELAPSED_TIME_BLAST = {} s\n".format(time.time() - t_blast).encode("utf-8"))
	os.write(log, "TOTAL_LENGHT_ANALYZED = {} bp\n".format(total_length).encode("utf-8"))
	
	## Create file containing size of sequences to convert wig to bigwig
	seq_sizes_file(headers, seq_lengths_all, HTML_DATA)
	
	## Protein domains module
	t_domains=time.time()
	if DOMAINS:
		os.write(log, "Domains module has started...\n". encode("utf-8"))
		os.fsync(log)
		domains_primary = NamedTemporaryFile(delete=False)
		protein_domains.domain_search(QUERY, LAST_DB, CLASSIFICATION, domains_primary.name, THRESHOLD_SCORE, WIN_DOM, OVERLAP_DOM)
		domains_primary.close()
		[xminimal, xmaximal, domains, seq_ids_dom] = domains_filtering.filter_qual_dom(domains_primary.name, OUT_DOMAIN_GFF, 0.35, 0.45, 0.8, 3, 'All', "")
		os.unlink(domains_primary.name)
		os.write(log, "ELAPSED_TIME_DOMAINS = {} s\n".format(time.time() - t_domains).encode("utf-8"))
		
		# Process individual sequences from the input file sequentially
		t_gff_vis = time.time() 
		repeats_process_dom(OUTPUT_GFF, THRESHOLD, THRESHOLD_SEGMENT, HTML_DATA, xminimal, xmaximal, domains, seq_ids_dom, CN, headers, seq_lengths_all, files_dict)
		os.write(log, "ELAPSED_TIME_GFF_VIS = {} s\n".format(time.time() - t_gff_vis).encode("utf-8"))
		
		# Prepare data for html output
		t_jbrowse=time.time()
		os.write(log, "JBrowse tracks are being prepared...\n".encode("utf-8"))
		os.fsync(log)
		jbrowse_prep_dom(HTML_DATA, QUERY, OUT_DOMAIN_GFF, OUTPUT_GFF, N_GFF, total_length, JBROWSE_BIN, files_dict)		
		os.write(log, "ELAPSED_TIME_JBROWSE_PREP = {} s\n".format(time.time() - t_jbrowse).encode("utf-8"))	
	else:
		# Process individual sequences from the input file sequentially
		t_gff_vis = time.time() 
		repeats_process(OUTPUT_GFF, THRESHOLD, THRESHOLD_SEGMENT, HTML_DATA, CN, headers, seq_lengths_all, files_dict)
		os.write(log, "ELAPSED_TIME_GFF_VIS = {} s\n".format(time.time() - t_gff_vis). encode("utf-8"))
		
		# Prepare data for html output
		t_jbrowse=time.time()
		jbrowse_prep(HTML_DATA, QUERY, OUTPUT_GFF, N_GFF, total_length, JBROWSE_BIN, files_dict)		
		os.write(log, "ELAPSED_TIME_JBROWSE_PREP = {} s\n".format(time.time() - t_jbrowse).encode("utf-8"))
	
	# Create HTML output
	t_html=time.time()
	os.write(log, "HTML output and JBrowse data structure are being prepared...\n".encode("utf-8"))
	os.fsync(log)
	html_output(total_length, seq_lengths_all, headers, HTML, DB_NAME, REF, REF_LINK)
	os.write(log, "ELAPSED_TIME_HTML = {} s\n".format(time.time() - t_html).encode("utf-8"))
	os.write(log, "ELAPSED_TIME_PROFREP = {} s\n".format(time.time() - t_profrep).encode("utf-8"))
	os.close(log)

	
	for subfasta in fasta_list:
		os.unlink(subfasta)
	
if __name__ == "__main__":
    import argparse  
    from argparse import RawDescriptionHelpFormatter  
    
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
	    pass
    
    # Default values(command line usage)
    HTML = configuration.HTML
    DOMAINS_GFF = configuration.DOMAINS_GFF
    REPEATS_GFF = configuration.REPEATS_GFF
    N_GFF = configuration.N_GFF
    LOG_FILE = configuration.LOG_FILE
    
    # Command line arguments
    parser = argparse.ArgumentParser(
    description='''
		
	DEPENDENCIES:
		- python 3.4 or higher with packages:
			- numpy
			- matplotlib
 		- [BLAST 2.2.28+](https://www.ncbi.nlm.nih.gov/books/NBK279690/) or higher
 		- ProfRep Modules:
			- gff.py
			- visualization.py
			- configuration.py 
			- protein_domains.py
			- domains_filtering.py

	EXAMPLE OF USAGE:
		
		./protein_domains_pd.py -q PATH_TO_INPUT_SEQ -pdb PATH_TO_PROTEIN_DB -cs PATH_TO_CLASSIFICATION_FILE
		''',
		epilog="""""", 
		formatter_class=CustomFormatter)
		
    Required = parser.add_argument_group('required arguments')
    altRequired = parser.add_argument_group('alternative required arguments - prepared datasets')
    blastOpt = parser.add_argument_group('optional arguments - BLAST Search')
    parallelOpt = parser.add_argument_group('optional arguments - Parallel Processing')
    protOpt = parser.add_argument_group('optional arguments - Protein Domains')
    outOpt = parser.add_argument_group('optional arguments - Output Paths')
    cnOpt = parser.add_argument_group('optional arguments - Copy Numbers/Hits ')
    galaxyOpt = parser.add_argument_group('optional arguments - Enviroment Variables')
    
    ################ INPUTS ############################################
    Required.add_argument('-q', '--query', type=str, required=True,
						help='input DNA sequence in (multi)fasta format')
    Required.add_argument('-d', '--database', type=str,
						help='blast database of all sequencing reads')
    Required.add_argument('-a', '--annotation_tbl', type=str,
						help='clusters annotation table, tab-separated number of cluster and its classification')
    Required.add_argument('-c', '--cls', type=str, 
						help='cls file containing reads assigned to clusters (hitsort.cls)')
    altRequired.add_argument('-tbl', '--datasets_tbl', type=str,
                        help='table with prepared annotation datasets') 
    altRequired.add_argument('-id', '--db_id', type=str,
                        help='annotation dataset ID (first column of datasets table)')  
                         					
	################ BLAST parameters ##################################
    blastOpt.add_argument('-bs', '--bit_score', type=float, default=50,
						help='bitscore threshold')
    blastOpt.add_argument('-m', '--max_alignments', type=int, default=10000000,
						help='blast filtering option: maximal number of alignments in the output')
    blastOpt.add_argument('-e', '--e_value', type=str, default=0.1,
						help='blast setting option: e-value')
    blastOpt.add_argument('-df', '--dust_filter', type=str, default="'20 64 1'",
						help='dust filters low-complexity regions during BLAST search')
    blastOpt.add_argument('-ws', '--word_size', type=int, default=11,
						help='blast search option: initial word size for alignment')
    blastOpt.add_argument('-t', '--task', type=str, default="blastn",
						help='type of blast to be triggered')
    blastOpt.add_argument('-n', '--new_db', type= str2bool, default=False,
						help='create a new blast database')	
						
	############### PARARELL PROCESSING ARGUMENTS ######################		
    parallelOpt.add_argument('-w', '--window', type=int, default=5000,
						help='sliding window size for parallel processing')
    parallelOpt.add_argument('-o', '--overlap', type=int, default=150,
						help='overlap for parallely processed regions, set greater than a read size')
						
	################ PROTEIN DOMAINS PARAMETERS ########################
    protOpt.add_argument('-pd', '--protein_domains', type=str2bool, default=True,
						help='use module for protein domains')
    protOpt.add_argument('-pdb', '--protein_database', type=str,
                        help='protein domains database')
    protOpt.add_argument('-cs', '--classification', type=str,
                        help='protein domains classification file')
    protOpt.add_argument('-wd', '--win_dom', type=int, default=10000000,
						help='protein domains module: sliding window to process large input sequences sequentially')
    protOpt.add_argument('-od', '--overlap_dom', type=int, default=10000,
						help='protein domains module: overlap of sequences in two consecutive windows')
    protOpt.add_argument('-thsc', '--threshold_score', type=int, default=80,
						help='protein domains module: percentage of the best score within the cluster to  significant domains')
		
	################ OUTPUTS ###########################################
    outOpt.add_argument('-lg', '--log_file', type=str, default=LOG_FILE,
                  		help='path to log file')
    outOpt.add_argument('-ouf', '--output_gff', type=str, default=REPEATS_GFF,
                        help='path to output gff of repetitive regions')
    outOpt.add_argument('-oug', '--domain_gff',type=str, default=DOMAINS_GFF,
						help='path to output gff of protein domains')
    outOpt.add_argument('-oun', '--n_gff',type=str, default=N_GFF,
						help='path to output gff of N regions')
    outOpt.add_argument('-hf', '--html_file', type=str, default=HTML,
                        help='path to output html file')
    outOpt.add_argument('-hp', '--html_path', type=str, default="output_dir",
                        help='path to html extra files')
								
	################ HITS/COPY NUMBERS ####################################
    cnOpt.add_argument('-cn', '--copy_numbers', type=str2bool, default=False,
                        help='convert hits to copy numbers')
    cnOpt.add_argument('-gs', '--genome_size', type=float,
                        help='genome size is required when converting hits to copy numbers and you use custom data')
    cnOpt.add_argument('-thr', '--threshold_repeat', type=int, default=5,
						help='threshold for hits/copy numbers per position to be considered repetitive')
    cnOpt.add_argument('-ths', '--threshold_segment', type=int, default=80,
                        help='threshold for the length of repetitive segment to be reported')                       
	
	################ GALAXY USAGE + JBrowse ##########################
    galaxyOpt.add_argument('-td', '--tool_dir', default=None,
                  		help='Galaxy tool data directory in galaxy')
    galaxyOpt.add_argument('-jb', '--jbrowse_bin', type=str,
                  		help='path to JBrowse bin directory')


    args = parser.parse_args()
    main(args)



