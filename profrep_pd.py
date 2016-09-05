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
import gff
import protein_domains_pd
import configuration
import visualization

t0 = time.time()
np.set_printoptions(threshold=np.nan)

# Create single fasta temporary files to be processed sequentially
def multifasta(QUERY):
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


# Read fasta, gain header and sequence without gaps
def fasta_read(subfasta):
	sequence_lines = []
	with open(subfasta, "r") as fasta:
		header = fasta.readline().strip().split(" ")[0]
		for line in fasta:
			clean_line = line.strip()			
			if clean_line:				
				sequence_lines.append(clean_line)
	sequence = "".join(sequence_lines)
	print(header)
	return header, sequence


# Create dictionary of known annotations classes and related clusters
def cluster_annotation(CL_ANNOTATION_TBL):
	cl_annotations = {} 			
	# Load annotation table as 2D array
	annot_table = np.genfromtxt(CL_ANNOTATION_TBL, dtype=str) 	
	for line in annot_table:
		if line[1] in cl_annotations:
			cl_annotations[line[1]].append(line[0])
		else:
			cl_annotations[line[1]] = [line[0]]
	return list(cl_annotations.items()), list(cl_annotations.keys())


# Create dictionary of known annotation classes and related reads
def read_annotation(CLS, cl_annotations_items):
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

# Predefine dictionary of known annotations and partial sequence repetitive profiles defined by parallel process
def annot_profile(annotation_keys, part):
	subprofile = {} 				
	for key in annotation_keys:
		subprofile[key] = np.zeros(part, dtype=int)
	subprofile["Unclassified_repetition"] = np.zeros(part, dtype=int)
	subprofile["all"] = np.zeros(part, dtype=int)
	return subprofile


# Run parallel function to process the input sequence in windows
# run blast for subsequence defined by the input index and window size
# create and increment subprofile vector based on reads aligned and position of alignment
def parallel_process(WINDOW, seq_length, annotation_keys, reads_annotations, subfasta, BLAST_DB, E_VALUE, WORD_SIZE, BLAST_TASK, MAX_ALIGNMENTS, MIN_IDENTICAL, MIN_ALIGN_LENGTH, subset_index):
	loc_start = subset_index + 1
	loc_end = subset_index + WINDOW
	if loc_end > seq_length:
		loc_end = seq_length
		subprofile = annot_profile(annotation_keys, seq_length - loc_start + 1)
	else:
		subprofile = annot_profile(annotation_keys, WINDOW + 1)
		
	# Find HSP records using blast for every window defined by query location and parse the tabular stdout -> 1. query, 2. database read, 3. %identical, 4. alignment length, 5. alignment start, 6. alignment end
	p = subprocess.Popen("blastn -query {} -query_loc {}-{} -db {} -evalue {} -word_size {} -task {} -num_alignments {} -outfmt '6 qseqid sseqid pident length qstart qend'".format(subfasta, loc_start, loc_end, BLAST_DB, E_VALUE, WORD_SIZE, BLAST_TASK, MAX_ALIGNMENTS), stdout=subprocess.PIPE, shell=True)
	for line in p.stdout:
		column = line.decode("utf-8").rstrip().split("\t")
		if float(column[2]) >= MIN_IDENTICAL and int(column[3]) >= MIN_ALIGN_LENGTH:
			read = column[1]				# ID of individual aligned read
			qstart = int(column[4])			# starting position of alignment
			qend = int(column[5])			# ending position of alignemnt
			# Assign repetition class to the read, if does not belong to any, assign to general "repeat" class
			if read in reads_annotations:
				annotation = reads_annotations[read]								
			else:
				annotation = "Unclassified_repetition"
			subprofile[annotation][qstart-subset_index-1 : qend-subset_index] = subprofile[annotation][qstart-subset_index-1 : qend-subset_index] + 1
	return subprofile


# Concatenate sequential subprofiles to the final profile representing the whole sequence and deal with overlaping parts
def concatenate_dict(profile_list, WINDOW, OVERLAP):
	profile = {}
	for key in profile_list[0].keys():
		if len(profile_list) == 1:
			profile[key] = profile_list[0][key]
		else:
			profile[key] = profile_list[0][key][0 : WINDOW-OVERLAP//2]
			for item in profile_list[1:-1]:
				profile[key] = np.append(profile[key], item[key][OVERLAP//2 : WINDOW-OVERLAP//2])
			profile[key] = np.append(profile[key], profile_list[-1][key][OVERLAP//2:])
	return profile


def hits_table(profile, OUTPUT, seq_id, seq_length):
	nonzero_keys = []
	[nonzero_keys.append(k) for k,v in profile.items() if any(v)]
	if any(profile["all"]):
		data = np.zeros((seq_length, len(nonzero_keys) + 1), dtype=int)
		data[:,0] = np.array([list(range(1, seq_length + 1))])
		data[:,1] = profile["all"] 
		keys = set(nonzero_keys)
		exclude = set(["all"])
		header = "{}\tall\t{}".format(seq_id, "\t".join(keys.difference(exclude)))
		count = 2
		for key in keys.difference(exclude):
			data[:,count] = profile[key]
			count += 1
	else:
		header = "{}\tno_hits".format(seq_id)    
		data = []
	with open(OUTPUT, "ab") as hits_tbl:
		np.savetxt(hits_tbl, data, delimiter= "\t", header=header, fmt="%d")
	
		
def seq_process(OUTPUT, SEQ_INFO, OUTPUT_GFF, THRESHOLD, THRESHOLD_SEGMENT, GFF, HTML_DATA, OUT_DOMAIN_GFF, xminimal, xmaximal, domains, seq_ids_dom):
	with open(SEQ_INFO, "r") as s_info:
		next(s_info)
		repeats_all = []
		header_gff = "##gff-version 3"
		with open(OUTPUT_GFF, "a") as gff_file:
			gff_file.write("{}\n".format(header_gff))
		for line in s_info:
			print(line)
			line_parsed = line.strip().split("\t")
			fasta_start = int(line_parsed[3])
			fasta_end = line_parsed[4]
			seq_length = int(line_parsed[1])
			seq_repeats = np.genfromtxt(OUTPUT, names=True, dtype="int", skip_header=fasta_start-1, max_rows=seq_length, delimiter="\t", deletechars="")
			seq_id = seq_repeats.dtype.names[0]
			if any(seq_repeats):
				if GFF:
					gff.create_gff(seq_repeats, OUTPUT_GFF, THRESHOLD, THRESHOLD_SEGMENT)
				repeats_all = create_wig(seq_repeats, seq_id, HTML_DATA, repeats_all)
			[fig, ax] = visualization.vis_profrep(seq_repeats, seq_length)
			if seq_id in seq_ids_dom:
				dom_idx = seq_ids_dom.index(seq_id) 
				[fig, ax] = visualization.vis_domains(fig, ax, seq_id, xminimal[dom_idx], xmaximal[dom_idx], domains[dom_idx])
			output_pic_png = "{}/{}.png".format(HTML_DATA, seq_id)
			fig.savefig(output_pic_png, bbox_inches="tight", format="png")	
	return repeats_all

	
def html_output(seq_names, link, HTML):
	pictures = "\n".join(["<img src={}.png>".format(pic)for pic in seq_names])
	html_str = """
	<!DOCTYPE html>
	<html>
	<body>
		<h1>Repetitive profile picture</h1>
		{}
		<a href={} target="_blank" >Link to JBrowser</a>	
	</body>
	</html>
	""".format(pictures, link)
	with open(HTML,"w") as html_file:
		html_file.write(html_str)


def jbrowse_prep(HTML_DATA, QUERY, OUT_DOMAIN_GFF, OUTPUT_GFF, repeats_all):
	print(repeats_all)
	jbrowse_data_path = os.path.join(HTML_DATA, configuration.jbrowse_data_dir)
	convert = "%2F"
	if os.getenv("JBROWSE_BIN"):
		JBROWSE_BIN = os.environ["JBROWSE_BIN"]
		extra_data_path = "/".join(HTML_DATA.split("/")[-2:])
		link_part2 = os.path.join(configuration.jbrowse_link_to_galaxy, extra_data_path, configuration.jbrowse_data_dir).replace("/",convert)
	else: 
		JBROWSE_BIN = configuration.JBROWSE_BIN_PC
		jbrowse_link_to_profrep = "data/profrep_data"
		link_part2 = os.path.join(jbrowse_link_to_profrep, "jbrowse").replace("/", convert)
		#QUERY = "/mnt/raid/users/ninah/profrep_git/tmp/test_seq_1"
	link = configuration.LINK_PART1_PC + link_part2
	subprocess.call("{}/prepare-refseqs.pl --fasta {} --out {}".format(JBROWSE_BIN, QUERY, jbrowse_data_path), shell=True)
	subprocess.call("{}/flatfile-to-json.pl --gff {} --trackLabel GFF_domains --out {}".format(JBROWSE_BIN, OUT_DOMAIN_GFF, jbrowse_data_path), shell=True)
	subprocess.call("{}/flatfile-to-json.pl --gff {} --trackLabel GFF_repetitions --out {}".format(JBROWSE_BIN, OUTPUT_GFF, jbrowse_data_path), shell=True)
	for repeat_id in repeats_all:
		subprocess.call("{}/wig-to-json.pl --trackLabel {} --wig {}/{}.wig --out {}".format(JBROWSE_BIN, repeat_id, HTML_DATA, repeat_id.split("/")[-1], jbrowse_data_path), shell=True)
	return link

	
def create_wig(seq_repeats, seq_id, HTML_DATA, repeats_all):
	header_repeats = seq_repeats.dtype.names
	seq_id = header_repeats[0]
	for track in header_repeats[1:]:
		header_wig = "variableStep\tchrom={}".format(seq_id)
		track_name = track.split("/")[-1]
		if track not in repeats_all:
			repeats_all.append(track)
		with open("{}/{}.wig".format(HTML_DATA, track_name), "ab") as track_file:
			np.savetxt(track_file, np.c_[seq_repeats[seq_id],seq_repeats[track]], fmt='%d', delimiter = "\t", header=header_wig, comments="")
	return repeats_all
	
	
def main(args):
	# Parse the command line arguments 
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
	SEQ_INFO = args.seq_info
	
	# Create new blast database of reads
	if NEW_DB:
		subprocess.call("makeblastdb -in {} -dbtype nucl".format(BLAST_DB), shell=True)
	
	if not os.path.exists(HTML_DATA):
		os.makedirs(HTML_DATA)
		
	# Define the parallel process
	STEP = WINDOW - OVERLAP		
	NUM_CORES = multiprocessing.cpu_count()	
	print("NUM_OF_CORES={}".format(NUM_CORES))
	parallel_pool = Pool(NUM_CORES)

	# Assign clusters to repetitive classes
	[cl_annotations_items, annotation_keys] = cluster_annotation(CL_ANNOTATION_TBL)
	
	# Assign reads to repetitive classes
	reads_annotations = read_annotation(CLS, cl_annotations_items)
	
	# Process every input fasta sequence sequentially
	fasta_list = multifasta(QUERY)
	headers=[]
	seq_count = 1
	start = 1
	with open(SEQ_INFO, "a") as s_info:
		s_info.write("seq_id\tseq_legth\tseq_count\tfile_start_pos\tfile_end_pos\n")
	for subfasta in fasta_list:
		[header, sequence] = fasta_read(subfasta)
		seq_length = len(sequence)
		end = start + seq_length
		seq_id = header[1:]
		##############################################################################################
		headers.append(header[1:])
		# Create parallel process																												
		subset_index = list(range(0, seq_length, STEP))	
		multiple_param = partial(parallel_process, WINDOW, seq_length, annotation_keys, reads_annotations, subfasta, BLAST_DB, E_VALUE, WORD_SIZE, BLAST_TASK, MAX_ALIGNMENTS, MIN_IDENTICAL, MIN_ALIGN_LENGTH)	
		profile_list = parallel_pool.map(multiple_param, subset_index)		 							

		# Join partial profiles to the final profile of the sequence 
		profile = concatenate_dict(profile_list, WINDOW, OVERLAP)

		# Sum the profile counts to get "all" repetitive profile (including all repetitions and also hits not belonging anywhere)
		profile["all"] = sum(profile.values())
		if not any(profile["all"]):
			end = start
		hits_table(profile, OUTPUT, seq_id, seq_length)
		with open(SEQ_INFO, "a") as s_info:
			s_info.write("{}\t{}\t{}\t{}\t{}\n".format(seq_id, seq_length, seq_count, start, end))
		start = end + 1
		seq_count += 1
		
	if DOMAINS:
		[xminimal, xmaximal, domains, seq_ids_dom] = protein_domains_pd.somain_search(QUERY, LAST_DB, CLASSIFICATION, OUT_DOMAIN_GFF, HTML_DATA)

	repeats_all = seq_process(OUTPUT, SEQ_INFO, OUTPUT_GFF, THRESHOLD, THRESHOLD_SEGMENT, GFF, HTML_DATA, OUT_DOMAIN_GFF, xminimal, xmaximal, domains, seq_ids_dom)
	link = jbrowse_prep(HTML_DATA, QUERY, OUT_DOMAIN_GFF, OUTPUT_GFF, repeats_all)		
	html_output(headers, link, HTML)
	
if __name__ == "__main__":
    
    HTML = configuration.HTML
    TMP = configuration.TMP
    DOMAINS_GFF = configuration.DOMAINS_GFF
    REPEATS_GFF = configuration.REPEATS_GFF
    REPEATS_TABLE = configuration.REPEATS_TABLE
    CLASSIFICATION = configuration.CLASSIFICATION
    LAST_DB = configuration.LAST_DB
    SEQ_INFO = configuration.SEQ_INFO
    
    
    #Define command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-q', '--query', type=str, required=True,
						help='query sequence to be processed')
    parser.add_argument('-d', '--database', type=str, required=True,
						help='blast database of all reads')
    parser.add_argument('-a', '--annotation_tbl', type=str, required=True,
						help='clusters annotation table')
    parser.add_argument('-c', '--cls', type=str, required=True,
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
    parser.add_argument('-th', '--threshold', type=int, default=50,
						help='threshold (number of hits) for report repetitive area in gff')
    parser.add_argument('-ths', '--threshold_segment', type=int, default=50,
                        help='threshold for a single segment length to be reported as repetitive reagion in gff')
    parser.add_argument('-pd', '--protein_domains', default=True,
						help='use module for protein domains')
    parser.add_argument('-pdb', '--protein_database', type=str, default=LAST_DB,
                        help='protein domains database')
    parser.add_argument('-cs', '--classification', type=str, default=CLASSIFICATION,
                        help='protein domains classification file')
    parser.add_argument('-ou', '--output', type=str, default=REPEATS_TABLE,
						help='output profile table name')
    parser.add_argument('-ouf', '--output_gff', type=str, default=REPEATS_GFF,
                        help='output gff format')
    parser.add_argument("-oug", "--domain_gff",type=str, default=DOMAINS_GFF,
						help="output domains gff format")
    parser.add_argument("-hf", "--html_file", type=str, default=HTML,
                        help="output html file name")
    parser.add_argument("-hp", "--html_path", type=str, default=TMP,
                        help="path to html extra files")
    parser.add_argument("-si", "--seq_info", type=str, default=SEQ_INFO,
                        help="file containg general info about sequence")

    args = parser.parse_args()
    main(args)

print("ELAPSED_TIME_PROFREP = {}".format(time.time() - t0))

