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
import pickle
import gff
import protein_domains_pd
import configuration
import visualization
import config_jbrowse
#####################
import distutils
from distutils import dir_util
import tempfile
#####################


t_profrep = time.time()
np.set_printoptions(threshold=np.nan)


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
		print(fasta_list)
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
			count += 1
			if count%2 == 0:
				reads = line.rstrip().split(" ") 
				for element in reads: 
					for key, value in cl_annotations_items:
						if clust in value:
							reads_annotations[element] = key 
			else:
				clust = line.split(" ")[0][3:]  
	return reads_annotations


def annot_profile(annotation_keys, part):
	''' Predefine dictionary of known annotations and partial sequence 
	repetitive profiles defined by parallel process '''
	subprofile = {} 				
	for key in annotation_keys:
		subprofile[key] = np.zeros(part, dtype=int)
	subprofile["ALL"] = np.zeros(part, dtype=int)
	return subprofile


def parallel_process(WINDOW, OVERLAP, seq_length, annotation_keys, reads_annotations, subfasta, BLAST_DB, E_VALUE, WORD_SIZE, BLAST_TASK, MAX_ALIGNMENTS, MIN_IDENTICAL, MIN_ALIGN_LENGTH, last_index, subsets_num, subset_index):
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
	# 1. query, 2. database read, 3. %identical, 4. alignment length, 5. alignment start, 6. alignment end
	p = subprocess.Popen("blastn -query {} -query_loc {}-{} -db {} -evalue {} -word_size {} -task {} -num_alignments {} -outfmt '6 qseqid sseqid pident length qstart qend'".format(subfasta, loc_start, loc_end, BLAST_DB, E_VALUE, WORD_SIZE, BLAST_TASK, MAX_ALIGNMENTS), stdout=subprocess.PIPE, shell=True)
	for line in p.stdout:
		column = line.decode("utf-8").rstrip().split("\t")
		if float(column[2]) >= MIN_IDENTICAL and int(column[3]) >= MIN_ALIGN_LENGTH:
			read = column[1]				# ID of individual aligned read
			qstart = int(column[4])			# starting position of alignment
			qend = int(column[5])			# ending position of alignemnt
			if read in reads_annotations:
				annotation = reads_annotations[read]								
			else:
				annotation = "ALL"
			subprofile[annotation][qstart-subset_index-1 : qend-subset_index] = subprofile[annotation][qstart-subset_index-1 : qend-subset_index] + 1
	#ntf = NamedTemporaryFile(suffix='_{}.pickle'.format(subset_index),delete=False)
	#with open(ntf.name, 'wb') as handle:
		#pickle.dump(subprofile, handle, protocol=pickle.HIGHEST_PROTOCOL)
	subprofile["ALL"] = sum(subprofile.values())
	if subset_index == 0: 
		if subsets_num == 1:
			[data, nonzero_len] = subprofile_single(subprofile, subset_index, annotation_keys, WINDOW, OVERLAP)
		else:
			[data, nonzero_len] = subprofile_first(subprofile, subset_index, annotation_keys, WINDOW, OVERLAP)
	elif subset_index == last_index:
		[data, nonzero_len] = subprofile_last(subprofile, subset_index, annotation_keys, WINDOW, OVERLAP)
	else:
		[data, nonzero_len] = subprofiles_middle(subprofile, subset_index, annotation_keys, WINDOW, OVERLAP)
	ntf = NamedTemporaryFile(suffix='_{}.tmp'.format(subset_index),delete=False)
	with open(ntf.name, 'wb') as tmp_subprof:
		np.savetxt(tmp_subprof, data, delimiter= "\t", fmt="%d")
	ntf.close()
	return ntf.name, nonzero_len
	
def subprofile_single(subprofile, subset_index, annotation_keys, WINDOW, OVERLAP):
	data = np.zeros((len(subprofile["ALL"]), len(annotation_keys) + 2), dtype=int)
	data[:,0] = np.array([list(range(1, len(subprofile["ALL"]) + 1))])
	exclude = set(["ALL"])
	#header = "{}\tALL\t{}".format(seq_id, "\t".join(sorted(set(annotation_keys).difference(exclude))))
	count = 2
	data[:,1] = subprofile["ALL"]
	for key in sorted(set(annotation_keys).difference(exclude)):
		data[:,count] = subprofile[key]
		count += 1
	# exclude rows containing only zeros
	data = data[~np.all(data[:,1:] == 0, axis=1)]
	nonzero_len = len(data)
	return data, nonzero_len
	
def subprofile_first(subprofile, subset_index, annotation_keys, WINDOW, OVERLAP):
	data = np.zeros((WINDOW-OVERLAP//2, len(annotation_keys) + 2), dtype=int)
	#if any(profile["ALL"]):
	data[:,0] = np.array([list(range(subset_index + 1, subset_index + WINDOW-OVERLAP//2 + 1))])
	#keys = set(nonzero_keys)
	exclude = set(["ALL"])
	#header = "{}\tALL\t{}".format(seq_id, "\t".join(sorted(set(annotation_keys).difference(exclude))))
	count = 2
	data[:,1] = subprofile["ALL"][0 : -OVERLAP//2-1]
	for key in sorted(set(annotation_keys).difference(exclude)):
		data[:,count] = subprofile[key][0 : -OVERLAP//2-1]
		count += 1
	# exclude rows containing only zeros
	data = data[~np.all(data[:,1:] == 0, axis=1)]
	nonzero_len = len(data)
	return data, nonzero_len
	
	
def subprofiles_middle(subprofile, subset_index, annotation_keys, WINDOW, OVERLAP):
	data = np.zeros((WINDOW-2*(OVERLAP//2), len(annotation_keys) + 2), dtype=int)
	#if any(profile["ALL"]):
	data[:,0] = np.array([list(range(subset_index + OVERLAP//2 + 1, subset_index + WINDOW-OVERLAP//2 + 1))])
	#keys = set(annotation_keys)
	exclude = set(["ALL"])
	#header = "{}\tALL\t{}".format(seq_id, "\t".join(sorted(set(annotation_keys).difference(exclude))))
	count = 2
	data[:,1] = subprofile["ALL"][OVERLAP//2 : -OVERLAP//2-1]
	for key in sorted(set(annotation_keys).difference(exclude)):
		data[:,count] = subprofile[key][OVERLAP//2 : -OVERLAP//2-1]
		count += 1
	# exclude rows containing only zeros
	data = data[~np.all(data[:,1:] == 0, axis=1)]
	nonzero_len = len(data)
	return data, nonzero_len
	

def subprofile_last(subprofile, subset_index, annotation_keys, WINDOW, OVERLAP):
	data = np.zeros((len(subprofile["ALL"])-OVERLAP//2, len(annotation_keys) + 2), dtype=int)
	#if any(profile["ALL"]):
	data[:,0] = np.array([list(range(subset_index + OVERLAP//2 + 1, subset_index + len(subprofile["ALL"])+1))])
	#keys = set(nonzero_keys)
	exclude = set(["ALL"])
	#header = "{}\tALL\t{}".format(seq_id, "\t".join(sorted(set(annotation_keys).difference(exclude))))
	count = 2
	data[:,1] = subprofile["ALL"][OVERLAP//2:]
	for key in sorted(set(annotation_keys).difference(exclude)):
		data[:,count] = subprofile[key][OVERLAP//2:]
		count += 1
	# exclude rows containing only zeros
	data = data[~np.all(data[:,1:] == 0, axis=1)]
	nonzero_len = len(data)
	return data, nonzero_len




#def concatenate_dict(profile_list, WINDOW, OVERLAP):
	#''' Concatenate sequential subprofiles to the final profile representing 
	#the whole sequence and deal with overlaping parts '''
	#profile = {}
	#with open(profile_list[0], 'rb') as handle_first:
		#first_dict = pickle.load(handle_first)
		#for key in first_dict.keys():
			#if len(profile_list) == 1:
				#profile[key] = first_dict[key]
			#else:
				##profile[key] = first_dict[key][0 : WINDOW-OVERLAP//2]
				#profile[key] = first_dict[key][0 : -OVERLAP//2-1]
				#for item in profile_list[1:-1]:
					#with open(item, 'rb') as handle:
						#individual_dict = pickle.load(handle)
						##profile[key] = np.append(profile[key], individual_dict[key][OVERLAP//2 : WINDOW-OVERLAP//2])
						#profile[key] = np.append(profile[key], individual_dict[key][OVERLAP//2 : -OVERLAP//2-1])
				#with open(profile_list[-1],'rb') as handle_last:
					#last_dict = pickle.load(handle_last)
					#profile[key] = np.append(profile[key], last_dict[key][OVERLAP//2:])
	#return profile

def concatenate_prof(OUTPUT, profile_list, CV):
	with open (OUTPUT, 'ab') as profile_tbl:
		for subprofile in profile_list:
			with open(subprofile, 'rb') as handle:
				if CV:
					for line in handle:
						#print((line.split("\t")[1:])/CV)
						counter = line.decode("utf-8").split("\t")[0]
						#cn_line = "\t".join([x for x in line.decode("utf-8").rstrip().split("\t")[1:]])
						cn_line = "\t".join([str(int((int(x)/CV))) for x in line.decode("utf-8").rstrip().split("\t")[1:]])
						profile_tbl.write("{}\t{}\n".format(counter, cn_line).encode("utf-8"))
				else:
					for line in handle:
						profile_tbl.write(line)
					


def hits_table(profile, OUTPUT, seq_id, seq_length, CV):
	''' Create table of blast hits of all fasta sequences '''
	nonzero_keys = []
	[nonzero_keys.append(k) for k,v in profile.items() if any(v)]
	data = np.zeros((seq_length, len(nonzero_keys) + 1), dtype=int)
	if any(profile["ALL"]):
		data[:,0] = np.array([list(range(1, seq_length + 1))])
		keys = set(nonzero_keys)
		exclude = set(["ALL"])
		header = "{}\tALL\t{}".format(seq_id, "\t".join(sorted(keys.difference(exclude))))
		count = 2
		if CV:
			data[:,1] = profile["ALL"]/CV
			for key in sorted(keys.difference(exclude)):
				data[:,count] = profile[key]/CV
				count += 1
		else:
			data[:,1] = profile["ALL"]
			for key in sorted(keys.difference(exclude)):
				data[:,count] = profile[key]
				count += 1
		# exclude rows containing only zeros
		data = data[~np.all(data[:,1:] == 0, axis=1)]
	if not np.any(data):
		data = [0]
		header = seq_id
	nonzero_len = len(data)
	################# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	with open(OUTPUT, "ab") as hits_tbl:
		np.savetxt(hits_tbl, data, delimiter= "\t", header=header, fmt="%d")
	return nonzero_len


def seq_process(OUTPUT, SEQ_INFO, OUTPUT_GFF, THRESHOLD, THRESHOLD_SEGMENT, GFF, HTML_DATA, OUT_DOMAIN_GFF, xminimal, xmaximal, domains, seq_ids_dom, CV):
	''' Process the hits table separately for each fasta, create gff file and profile picture '''
	with open(SEQ_INFO, "r") as s_info:
		next(s_info)
		repeats_all_seq = []
		#header_gff = "##gff-version 3"
		#with open(OUTPUT_GFF, "a") as gff_file:
			#gff_file.write("{}\n".format(header_gff))
		with open(OUTPUT_GFF, "w") as gff_file:
			gff_file.write("{}\n".format(configuration.HEADER_GFF))
		gff_repeats = open(OUTPUT_GFF, "a")	
		seq_count = 1
		for line in s_info:
			####################################
			present_repeats = []
			####################################
			line_parsed = line.strip().split("\t")
			fasta_start = int(line_parsed[3])
			fasta_end = int(line_parsed[4])
			seq_length = int(line_parsed[1])
			seq_repeats = np.genfromtxt(OUTPUT, names=True, dtype="int", skip_header=fasta_start-1, max_rows=fasta_end - fasta_start, delimiter="\t", deletechars="")
			seq_repeats = np.atleast_1d(seq_repeats)
			seq_id = seq_repeats.dtype.names[0]
			############################################################
			for repeat in seq_repeats.dtype.names[1:]:
				if not all(value == 0 for value in seq_repeats[repeat]):
					present_repeats.append(repeat)
					###
					if repeat not in repeats_all_seq:
						repeats_all_seq.append(repeat)
					### 
			############################################################
			if any(seq_repeats.shape):
				if GFF:
					gff.create_gff(seq_repeats, OUTPUT_GFF, THRESHOLD, THRESHOLD_SEGMENT, gff_repeats)
				#[repeats_all, max_wig] = create_wig(seq_repeats, seq_id, HTML_DATA, repeats_all)
				max_wig = create_wig(seq_id, present_repeats, seq_repeats, HTML_DATA)
			if seq_count <= configuration.MAX_PIC_NUM:
				[fig, ax] = visualization.vis_profrep(seq_id, present_repeats, seq_repeats, seq_length, CV)
				if seq_id in seq_ids_dom:
					dom_idx = seq_ids_dom.index(seq_id) 
					[fig, ax] = visualization.vis_domains(fig, ax, seq_id, xminimal[dom_idx], xmaximal[dom_idx], domains[dom_idx])
				output_pic_png = "{}/{}.png".format(HTML_DATA, seq_id)
				fig.savefig(output_pic_png, bbox_inches="tight", format="png", dpi=configuration.IMAGE_RES)	
			seq_count += 1
			plt.close()
		gff_repeats.close()
	print(repeats_all_seq)	
	return repeats_all_seq



def html_output(SEQ_INFO, total_length, seq_names, HTML, DB_NAME, REF, REF_LINK):
	''' Define html output with limited number of output pictures and link to JBrowse '''
	pictures = "\n\t\t".join(['<img src="{}.png" width=1800>'.format(pic)for pic in seq_names[:configuration.MAX_PIC_NUM] ])
	with open(SEQ_INFO, "r") as s_info:
		next(s_info)
		info = "\t\t".join(['<pre> {} [{} bp]</pre>'.format(line.split("\t")[0],line.split("\t")[1])for line in s_info])
	if REF:
		ref_part_1 = REF.split("-")[0]
		ref_part_2 = "-".join(REF.split("-")[1:]).split(". ")[0]
		ref_part_3 = ". ".join("-".join(REF.split("-")[1:]).split(". ")[1:])
		ref_string = '''<h6> {} - <a href="{}" target="_blank" >{}</a>. {}'''.format(ref_part_1, REF_LINK, ref_part_2, ref_part_3)
	else:
		ref_string = ""
	html_str = '''
	<!DOCTYPE html>
	<html>
	<body>
		<h2>PROFREP OUTPUT</h2>
		<h4> Sequences processed: </h4>
		{}
		<h4> Total length: </h4>
		<pre> {} bp </pre>
		<h4> Database: </h4>
		<pre> {} </pre>
		<h4> Interactive visualization: </h4>
		<hr>
		<h3> Repetitive profile(s)</h3> </br>
		{} <br/>
		<h4>References: </h4>
		{}
		</h6>
	</body>
	</html>
	'''.format(info, total_length, DB_NAME, pictures, ref_string)
	with open(HTML,"w") as html_file:
		html_file.write(html_str)


def jbrowse_prep(HTML_DATA, QUERY, OUT_DOMAIN_GFF, OUTPUT_GFF, repeats_all, N_GFF, total_length):
	''' Set up the paths, link and convert output data to be displayed as tracks in Jbrowse '''
	jbrowse_data_path = os.path.join(HTML_DATA, config_jbrowse.jbrowse_data_dir)
	#convert = "%2F"
	## Galaxy usage
	#if GALAXY:
		#results_path1 = "/".join(HTML.split("/")[HTML_DATA.split("/").index("database")+1:-1])
		#results_path2 = HTML_DATA.split("/")[-1]
		#link_part2 = os.path.join("data", config_jbrowse.link_to_files, results_path1, results_path2, configuration.jbrowse_data_dir).replace("/",convert)
	    ##JBROWSE_BIN = os.environ["JBROWSE_BIN"]
	    ##extra_data_path = "/".join(HTML_DATA.split("/")[-2:])
	    ##link_part2 = os.path.join(configuration.jbrowse_link_to_galaxy, extra_data_path, configuration.jbrowse_data_dir).replace("/",convert)
	    ##link = configuration.LINK_PART1 + link_part2

	## Local usage
	#else:            
		##JBROWSE_BIN = configuration.JBROWSE_BIN_PC
		##link_part2 = os.path.join(configuration.jbrowse_link_to_profrep, "jbrowse").replace("/", convert)
		##link = configuration.LINK_PART1_PC + link_part2
		#link_part2 = os.path.join("data", config_jbrowse.link_to_files, configuration.jbrowse_data_dir).replace("/", convert)
	#link_part1 = "http://{}/{}/index.html?data=".format(config_jbrowse.PC, config_jbrowse.JBROWSE)
	#link = link_part1 + link_part2
	jbrowse_bin = config_jbrowse.JBROWSE_BIN
	########################################################################################
	with tempfile.TemporaryDirectory() as dirpath:
	########################################################################################
		subprocess.call(["{}/prepare-refseqs.pl".format(jbrowse_bin), "--fasta", QUERY, "--out", jbrowse_data_path])
		subprocess.call(["{}/flatfile-to-json.pl".format(jbrowse_bin), "--gff", OUT_DOMAIN_GFF, "--trackLabel", "GFF_domains", "--out",  jbrowse_data_path])
		subprocess.call(["{}/flatfile-to-json.pl".format(jbrowse_bin), "--gff", OUTPUT_GFF, "--trackLabel", "GFF_repeats", "--config", configuration.JSON_CONF_R, "--out",  jbrowse_data_path])
		subprocess.call(["{}/flatfile-to-json.pl".format(jbrowse_bin), "--gff", N_GFF, "--trackLabel", "N_regions", "--config", configuration.JSON_CONF_N, "--out",  jbrowse_data_path])		 

		count = 0
		###############################################
		#control the total length processed, if above threshold, dont create wig image tracks 
		# default 1MB seq
		threshold_wig = 1000000
		if total_length <= threshold_wig:
		###############################################
			for repeat_id in repeats_all:
				color = configuration.COLORS_RGB[count]
				subprocess.call(["{}/wig-to-json.pl".format(jbrowse_bin), "--wig", "{}/{}.wig".format(HTML_DATA, repeat_id.split("/")[-1]), "--trackLabel", repeat_id, "--fgcolor", color, "--out",  jbrowse_data_path])
				count += 1
		distutils.dir_util.copy_tree(dirpath,jbrowse_data_path)
	return None


def create_wig(seq_id, present_repeats, seq_repeats, HTML_DATA):
	''' Create wiggle format to display individual profiles by JBrowse '''
	#header_repeats = seq_repeats.dtype.names
	#seq_id = header_repeats[0]
	max_wig = []
	#print(header_repeats)
	#for track in header_repeats[1:]:
	for track in present_repeats:
		header_wig = "variableStep\tchrom={}".format(seq_id)
		track_name = track.split("/")[-1]
		track_data = np.c_[seq_repeats[seq_id],seq_repeats[track]]
		track_nonzero = track_data[np.nonzero(track_data[:,1])]
		#if track not in repeats_all:
			#repeats_all.append(track)
		with open("{}/{}.wig".format(HTML_DATA, track_name), "ab") as track_file:
			np.savetxt(track_file, track_nonzero, fmt='%d', delimiter = "\t", header=header_wig, comments="")
		max_wig.append(max(seq_repeats[track]))
	return max_wig
	

def genome2coverage(GS, BLAST_DB):
	''' Convert genome size to coverage'''
	nr = subprocess.Popen('''cat {} | grep '>' | wc -l'''.format(BLAST_DB), stdout=subprocess.PIPE, shell=True)
	num_of_reads = int(nr.communicate()[0])
	########################################### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#########################################################
	lr = subprocess.Popen('''awk -v N=2 '{print}/>/&&--N<=0{exit}' ''' + BLAST_DB + '''| awk '$0 !~">"{print}' | awk '{sum+=length($0)}END{print sum}' ''', stdout=subprocess.PIPE, shell=True)
	len_of_read = int(lr.communicate()[0])
	CV = (num_of_reads*len_of_read)/(GS*1000000) # GS in Mbp
	print("COVERAGE = {}".format(CV))
	return CV

	
def prepared_data(TBL, DB_ID, TOOL_DATA_DIR):
	''' Get prepared rep. annotation data from the table based on the selected species ID '''
	with open(TBL) as datasets:
		for line in datasets:
			if DB_ID in line:
				DB_NAME = line.split("\t")[1]
				BLAST_DB = os.path.join(TOOL_DATA_DIR, line.split("\t")[2])
				CLS = os.path.join(TOOL_DATA_DIR, line.split("\t")[3])
				CL_ANNOTATION_TBL = os.path.join(TOOL_DATA_DIR, line.split("\t")[4])
				CV = float(line.split("\t")[5])
				REF = line.split("\t")[6]
				REF_LINK = line.split("\t")[7]
	return DB_NAME, BLAST_DB, CLS, CL_ANNOTATION_TBL, CV, REF, REF_LINK

	
def main(args):
	
	# Command line arguments
	QUERY = args.query
	BLAST_DB = args.database
	CL_ANNOTATION_TBL = args.annotation_tbl 
	CLS = args.cls
	MIN_IDENTICAL = args.identical
	MIN_ALIGN_LENGTH = args.align_length
	E_VALUE = args.e_value
	WORD_SIZE = args.word_size
	WINDOW = args.window
	OVERLAP = args.overlap
	BLAST_TASK = args.task
	MAX_ALIGNMENTS = args.max_alignments
	NEW_DB = args.new_db
	GFF = args.gff
	THRESHOLD = args.threshold
	THRESHOLD_SEGMENT = args.threshold_segment
	OUTPUT = args.output
	OUTPUT_GFF = args.output_gff
	DOMAINS = args.protein_domains
	LAST_DB = args.protein_database
	CLASSIFICATION = args.classification
	OUT_DOMAIN_GFF = args.domain_gff	
	HTML = args.html_file
	HTML_DATA = args.html_path
	N_GFF = args.n_gff
	CV = args.coverage
	CN = args.copy_numbers
	GS = args.genome_size
	DB_ID = args.db_id
	TBL = args.datasets_tbl
	DB_NAME = args.db_name
	THRESHOLD_SCORE = args.threshold_score
	WIN_DOM = args.win_dom
	OVERLAP_DOM = args.overlap_dom
	THRESHOLD_SCORE = args.threshold_score
	WIN_DOM = args.win_dom
	OVERLAP_DOM = args.overlap_dom
	GALAXY = args.galaxy_usage
	TOOL_DATA_DIR = args.tool_dir


	REF = None
	REF_LINK = None
	
	#if os.getenv("JBROWSE_BIN"):
		#GALAXY = True
		#CLASSIFICATION = os.path.join(configuration.TOOL_DATA_DIR, CLASSIFICATION)
		#LAST_DB = os.path.join(configuration.TOOL_DATA_DIR, LAST_DB)
		#TBL = os.path.join(configuration.TOOL_DATA_DIR, TBL)
	#else:
		#GALAXY = False
		
	
	
	# Parse prepared annotation data table
	if TBL:
		TBL = os.path.join(configuration.PROFREP_DATA, TBL)
		[DB_NAME, BLAST_DB, CLS, CL_ANNOTATION_TBL, CV, REF, REF_LINK] = prepared_data(TBL, DB_ID, TOOL_DATA_DIR)
	if GALAXY:
		LAST_DB = os.path.join(LAST_DB, configuration.LAST_DB_FILE)
		CLASSIFICATION = os.path.join(CLASSIFICATION, configuration.CLASS_FILE)

	
	# Calculate coverage 
	if not CN:
		CV = False
		
	# Create new blast database of reads
	if NEW_DB:
		subprocess.call("makeblastdb -in {} -dbtype nucl".format(BLAST_DB), shell=True)
	
	# Create dir to store outputs for html 
	if not os.path.exists(HTML_DATA):
		os.makedirs(HTML_DATA)
		
	if not os.path.isabs(OUT_DOMAIN_GFF):
		OUT_DOMAIN_GFF = os.path.join(HTML_DATA, OUT_DOMAIN_GFF)
	
	#if not os.path.isabs(N_GFF):
		#N_GFF = os.path.join(HTML_DATA, N_GFF)
	#############################################################################
	if not os.path.isabs(HTML):
		HTML = os.path.join(HTML_DATA, HTML)
	#############################################################################	
	# Define parameters for parallel process
	STEP = WINDOW - OVERLAP		
	NUM_CORES = multiprocessing.cpu_count()	
	print("NUM_OF_CORES = {}".format(NUM_CORES))
	parallel_pool = Pool(NUM_CORES)

	# Assign clusters to repetitive classes
	[cl_annotations_items, annotation_keys] = cluster_annotation(CL_ANNOTATION_TBL)
	
	# Assign reads to repetitive classes
	reads_annotations = read_annotation(CLS, cl_annotations_items)
	
	# Convert genome size to coverage
	if GS:
		CV = genome2coverage(GS, BLAST_DB)
	
	# Detect all fasta sequences from input
	fasta_list = multifasta(QUERY)
	headers=[]
	seq_count = 1
	start = 1
	total_length = 0
	with open(N_GFF, "w") as Ngff:
		Ngff.write("{}\n".format(configuration.HEADER_GFF))
	Ngff = open(N_GFF,"a")
	# Create file to record info about fasta sequences
	SEQ_INFO = "{}/{}".format(HTML_DATA, configuration.SEQ_INFO)
	with open(SEQ_INFO, "w") as s_info:
		s_info.write(configuration.s_info_header)
		# Find hits for each fasta sequence separetely
		t_blast=time.time()	
		for subfasta in fasta_list:
			[header, sequence] = fasta_read(subfasta)
			gff.N_gff(header, sequence, Ngff)
			seq_length = len(sequence)
			headers.append(header)
			with open(OUTPUT, 'a') as profile_tbl:
				profile_tbl.write( "{}\tALL\t{}\n".format(header, "\t".join(sorted(set(annotation_keys).difference(set(["ALL"]))))))
			# Create parallel process
			print(seq_length)																								
			subset_index = list(range(0, seq_length, STEP))
			# Situation when penultimal window is not complete but it is following by another one
			if len(subset_index) > 1 and subset_index[-2] + WINDOW >= seq_length:
				subset_index = subset_index[:-1]
			print(subset_index)	
			last_index = subset_index[-1]
			######## CONF #####################
			files_num = 1000
			###################################
			nonzero_total = 0
			index_range = range(len(subset_index))
			concatenated_prof = None
			for chunk_index in index_range[0::files_num]:
				multiple_param = partial(parallel_process, WINDOW, OVERLAP, seq_length, annotation_keys, reads_annotations, subfasta, BLAST_DB, E_VALUE, WORD_SIZE, BLAST_TASK, MAX_ALIGNMENTS, MIN_IDENTICAL, MIN_ALIGN_LENGTH, last_index, len(subset_index))
				parallel_data = parallel_pool.map(multiple_param, subset_index[chunk_index:chunk_index + files_num])
				profile_list, nonzero_len = zip(*parallel_data)
				nonzero_total = nonzero_total + sum(nonzero_len)
				# Join partial profiles to the final profile of the sequence 
				#if concatenated_prof:
					#profile_list.insert(0, concatenated_prof)	 							
				concatenate_prof(OUTPUT, profile_list, CV)
				#prof_temp = NamedTemporaryFile(suffix='.pickle',delete=False)
				#concatenated_prof = prof_temp.name
				#with open(concatenated_prof, 'wb') as handle:
					#pickle.dump(profile, handle, protocol=pickle.HIGHEST_PROTOCOL)
				#prof_temp.close()
				for subprofile in profile_list:
					os.unlink(subprofile)
			# Sum the profile counts to get "all" profile (including all repetitions and also hits not belonging anywhere)
			#profile["ALL"] = sum(profile.values())
			#nonzero_len = hits_table(profile, OUTPUT, header, seq_length, CV)
			if nonzero_total == 0:
				with open(OUTPUT, "a") as profile_tbl:
					no_repeats = ["0"] * (len(annotation_keys) + 2)
					profile_tbl.write("{}\n".format("\t".join(no_repeats)))
					nonzero_total = 1
			end = start + nonzero_total
			# Each line defines one sequence and its position in hits table 
			s_info.write("{}\t{}\t{}\t{}\t{}\n".format(header, seq_length, seq_count, start, end))
			start = end + 1
			seq_count += 1
			total_length += seq_length 
			#os.unlink(concatenated_prof)
	Ngff.close()
	print("ELAPSED_TIME_BLAST = {} s".format(time.time() - t_blast))
	print("TOTAL_LENGHT_ANALYZED = {} bp".format(total_length))
	
	# Protein domains module
	t_domains=time.time()
	if DOMAINS is True:
		[xminimal, xmaximal, domains, seq_ids_dom] = protein_domains_pd.domain_search(QUERY, LAST_DB, CLASSIFICATION, OUT_DOMAIN_GFF,  THRESHOLD_SCORE, WIN_DOM, OVERLAP_DOM)	
	print("ELAPSED_TIME_DOMAINS = {} s".format(time.time() - t_domains))
		
	t_gff_vis = time.time() 
	# Process individual sequences from the input file sequentially
	repeats_all = seq_process(OUTPUT, SEQ_INFO, OUTPUT_GFF, THRESHOLD, THRESHOLD_SEGMENT, GFF, HTML_DATA, OUT_DOMAIN_GFF, xminimal, xmaximal, domains, seq_ids_dom, CV)
	print("ELAPSED_TIME_GFF_VIS = {} s".format(time.time() - t_gff_vis))
	
	# Prepare data for html output
	t_jbrowse=time.time()
	jbrowse_prep(HTML_DATA, QUERY, OUT_DOMAIN_GFF, OUTPUT_GFF, repeats_all, N_GFF, total_length)		
	print("ELAPSED_TIME_JBROWSE_PREP = {} s".format(time.time() - t_jbrowse))
	t_html=time.time()
	html_output(SEQ_INFO, total_length, headers, HTML, DB_NAME, REF, REF_LINK)
	print("ELAPSED_TIME_HTML = {} s".format(time.time() - t_html))
	
	for subfasta in fasta_list:
		os.unlink(subfasta)
	
if __name__ == "__main__":
    
    # Default values(command line usage)
    HTML = configuration.HTML
    TMP = configuration.TMP
    DOMAINS_GFF = configuration.DOMAINS_GFF
    REPEATS_GFF = configuration.REPEATS_GFF
    N_REG = configuration.N_REG
    REPEATS_TABLE = configuration.REPEATS_TABLE
    #CLASSIFICATION = configuration.CLASSIFICATION
    #LAST_DB = configuration.LAST_DB

    
    # Command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-q', '--query', type=str, required=True,
						help='query sequence to be processed')
    parser.add_argument('-d', '--database', type=str,
						help='blast database of all reads')
    parser.add_argument('-a', '--annotation_tbl', type=str,
						help='clusters annotation table')
    parser.add_argument('-c', '--cls', type=str, 
						help='cls file containing reads assigned to clusters')
    parser.add_argument('-i', '--identical', type=float, default=95,
						help='blast filtering option: sequence indentity threshold between query and mapped read from db in %')
    parser.add_argument('-l', '--align_length', type=int, default=40,
						help='blast filtering option: minimal alignment length threshold in bp')
    parser.add_argument('-m', '--max_alignments', type=int, default=10000000,
						help='blast filtering option: maximal number of alignments in the output')
    parser.add_argument('-e', '--e_value', type=str, default=1e-15,
						help='blast setting option: e-value')
    parser.add_argument('-ws', '--word_size', type=int, default=11,
						help='blast setting option: initial word size for alignment')
    parser.add_argument('-t', '--task', type=str, default="blastn",
						help='type of blast to be triggered')
    parser.add_argument('-w', '--window', type=int, default=5000,
						help='window size for parallel processing')
    parser.add_argument('-o', '--overlap', type=int, default=150,
						help='overlap for parallely processed regions, set greater than read size')
    parser.add_argument('-n', '--new_db', default=False,
						help='create a new blast database')
    parser.add_argument('-g', '--gff', default=True,
						help='use module for gff')
    parser.add_argument('-th', '--threshold', type=int, default=5,
						help='threshold (number of hits) for report repetitive area in gff')
    parser.add_argument('-ths', '--threshold_segment', type=int, default=80,
                        help='threshold for a single segment length to be reported as repetitive reagion in gff')
    parser.add_argument('-pd', '--protein_domains', default=True,
						help='use module for protein domains')
    parser.add_argument('-pdb', '--protein_database', type=str,
                        help='protein domains database')
    parser.add_argument('-cs', '--classification', type=str,
                        help='protein domains classification file')
    parser.add_argument('-ou', '--output', type=str, default=REPEATS_TABLE,
						help='output profile table name')
    parser.add_argument('-ouf', '--output_gff', type=str, default=REPEATS_GFF,
                        help='output gff format')
    parser.add_argument("-oug", "--domain_gff",type=str, default=DOMAINS_GFF,
						help="output domains gff format")
    parser.add_argument("-oun", "--n_gff",type=str, default=N_REG,
						help="N regions gff format")
    parser.add_argument("-hf", "--html_file", type=str, default=HTML,
                        help="output html file name")
    parser.add_argument("-hp", "--html_path", type=str, default=TMP,
                        help="path to html extra files")
    parser.add_argument("-cv", "--coverage", type=float, 
                        help="coverage")
    parser.add_argument("-cn", "--copy_numbers", type=bool, 
                        help="convert hits to copy numbers")
    parser.add_argument("-gs", "--genome_size", type=float,
                        help="genome size")
    parser.add_argument("-id", "--db_id", type=str,
                        help="annotation database name")
    parser.add_argument("-tbl", "--datasets_tbl", type=str,
                        help="table with prepared anotation data")    
    parser.add_argument("-dbn", "--db_name", type=str,
                        help="custom database name")     
    parser.add_argument("-dir","--output_dir", type=str,
						help="specify if you want to change the output directory")
    parser.add_argument("-thsc","--threshold_score", type=int, default=80,
						help="percentage of the best score in the cluster to be tolerated when assigning annotations per base")
    parser.add_argument("-wd","--win_dom", type=int, default=10000000,
						help="window to process large input sequences sequentially")
    parser.add_argument("-od","--overlap_dom", type=int, default=10000,
						help="overlap of sequences in two consecutive windows")
    parser.add_argument("-gu", "--galaxy_usage", default=False,
                        help="option for galaxy usage only")
    parser.add_argument("-td", "--tool_dir", default=False,
                  		help="tool data directory in galaxy")



    args = parser.parse_args()
    main(args)

print("ELAPSED_TIME_PROFREP = {} s".format(time.time() - t_profrep))

