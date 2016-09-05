#!/usr/bin/env python3

import numpy as np
import subprocess
import csv
import time
import sys
import matplotlib.pyplot as plt
import multiprocessing
import argparse
import gff
from functools import partial
from multiprocessing import Pool
from tempfile import NamedTemporaryFile

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
		header = str(fasta.readline().strip())
		header = header.replace(",","_")	# prevent problems with csv parsing
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
	subprofile["repeat"] = np.zeros(part, dtype=int)
	subprofile["all"] = np.zeros(part, dtype=int)
	return subprofile


# Run parallel function to process the input sequence in windows
# run blast for subsequence defined by the input index and window size
# create and increment subprofile vector based on reads aligned and position of alignment
def parallel_process(WINDOW, query_length, annotation_keys, reads_annotations, subfasta, BLAST_DB, E_VALUE, WORD_SIZE, BLAST_TASK, MAX_ALIGNMENTS, MIN_IDENTICAL, MIN_ALIGN_LENGTH, subset_index):
	loc_start = subset_index + 1
	loc_end = subset_index + WINDOW
	if loc_end > query_length:
		loc_end = query_length
		subprofile = annot_profile(annotation_keys, query_length - loc_start + 1)
	else:
		subprofile = annot_profile(annotation_keys, WINDOW + 1)
		
	# Find HSP records using blast for every window defined by query location and parse the tabular stdout -> 1. query, 2. database read, 3. %identical, 4. alignment length, 5. alignment start, 6. alignment end
	p = subprocess.Popen("blastn -query {} -query_loc {}-{} -db ./reads -evalue {} -word_size {} -task {} -num_alignments {} -outfmt '6 qseqid sseqid pident length qstart qend'".format(subfasta, loc_start,loc_end, E_VALUE, WORD_SIZE, BLAST_TASK, MAX_ALIGNMENTS), stdout=subprocess.PIPE, shell=True)
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
				annotation = "repeat"
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


# Convert profile dictionary to output table and plot profile graphs
def profile_to_csv(profile, OUTPUT, OUTPUT_PIC, query_length, header, fig, ax, cm):
	data = np.zeros((len(profile), query_length), dtype=int)
	labels = []
	count = 0
	for key in list(profile.keys()):
		labels.append(key) 					
		data[count] = profile[key]		
		if np.any(data[count]):
			ax.plot(list(range(query_length)), data[count], label=labels[count])
		#ax.set_prop_cycle([cm(1.*count/len(profile))])
		if key == "all":
			all_position = count 				
		count += 1
	plt.legend(loc="upper left")
	OUPUT_PIC=fig.savefig("output.png")

	# Swap positions so that "all" record is the first in the table
	data[0], data[all_position] = data[all_position], data[0].copy()
	labels[0], labels[all_position] = labels[all_position], labels[0]
	
	# Write output to a csv table
	sequence_counter = np.array(list(range(1, query_length + 1)))				
	labels = np.append([header], labels)
	data = np.append([sequence_counter], data, 0)		
	output_table = open(OUTPUT, "a")
	writer = csv.writer(output_table, delimiter="\t")
	writer.writerow(labels)
	for values in np.transpose(data):		
		writer.writerow(values)
	return output_table


# Define profiles visualization
def set_visualization():
	fig = plt.figure()
	ax = fig.add_subplot(111)
	cm = plt.get_cmap("gist_rainbow")
	return fig, ax, cm
	

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
	NEW_DB=args.new_db
	GFF = args.gff
	THRESHOLD = args.threshold
	OUTPUT = args.output
	OUTPUT_GFF = args.output_gff
	OUTPUT_PIC = args.output_pic
	
	# Create new blast database of reads
	if NEW_DB:
		subprocess.call("makeblastdb -in {} -dbtype nucl -out ./reads".format(BLAST_DB), shell=True)
		
	# Define the parallel process
	STEP = WINDOW - OVERLAP		
	NUM_CORES = multiprocessing.cpu_count()	
	parallel_pool = Pool(NUM_CORES)

	# Assign clusters to repetitive classes
	[cl_annotations_items, annotation_keys] = cluster_annotation(CL_ANNOTATION_TBL)
	
	# Assign reads to repetitive classes
	reads_annotations = read_annotation(CLS, cl_annotations_items)
	
	# Process every input fasta sequence sequentially
	fasta_list = multifasta(QUERY)
	for subfasta in fasta_list:
		[header, sequence] = fasta_read(subfasta)
		query_length = len(sequence)
		
		# Create parallel process																												
		subset_index = list(range(0, query_length, STEP))	
		multiple_param = partial(parallel_process, WINDOW, query_length, annotation_keys, reads_annotations, subfasta, BLAST_DB, E_VALUE, WORD_SIZE, BLAST_TASK, MAX_ALIGNMENTS, MIN_IDENTICAL, MIN_ALIGN_LENGTH)	
		profile_list = parallel_pool.map(multiple_param, subset_index)		 							

		# Join partial profiles to the final profile of the sequence 
		profile = concatenate_dict(profile_list, WINDOW, OVERLAP)

		# Sum the profile counts to get "all" repetitive profile (including all repetitions and also hits not belonging anywhere)
		profile["all"] = sum(profile.values())

		# Set up the visualization
		[fig, ax, cm] = set_visualization()

		# Write output to a table and plot the profiles
		output_table = profile_to_csv(profile, OUTPUT, OUTPUT_PIC, query_length, header, fig, ax, cm)
	output_table.close()

	# Use gff module to convert profile table to GFF3 format
	if GFF:
		gff.main(OUTPUT, OUTPUT_GFF, THRESHOLD)		


if __name__ == "__main__":
    # Define command line arguments 
	parser = argparse.ArgumentParser()
	parser.add_argument('-q','--query', type=str, required=True,
						help='query sequence to be processed')
	parser.add_argument('-d','--database', type=str, required=True,
						help='blast database of all reads')
	parser.add_argument('-a','--annotation_tbl', type=str, required=True,
						help='clusters annotation table')
	parser.add_argument('-c','--cls', type=str, required=True,
						help='cls file containing reads assigned to clusters')
	parser.add_argument('-i','--identical', type=float, default=95,
						help='blast filtering option: sequence indentity threshold between query and mapped read from db in %')
	parser.add_argument('-l','--align_length', type=int, default=40,
						help='blast filtering option: minimal align length threshold in bp')
	parser.add_argument('-m','--max_alignments', type=int, default=10000000,
						help='blast filtering option: maximal number of alignments in the output')
	parser.add_argument('-e','--e_value', type=str, default=1e-15,
						help='blast setting option: e-value')
	parser.add_argument('-ws','--word_size', type=int, default=11,
						help='blast setting option: initial word size for alignment')
	parser.add_argument('-t','--task', type=str, default="blastn",
						help='type of blast to be triggered')
	parser.add_argument('-w','--window', type=int, default=5000,
						help='window size for parallel processing')
	parser.add_argument('-o','--overlap', type=int, default=150,
						help='overlap for parallely processed regions, set greater than read size')
	parser.add_argument('-n','--new_db', default=False,
						help='create a new blast database')
	parser.add_argument('-g','--gff', default=False,
						help='use module for gff')
	parser.add_argument('-th','--threshold',type=int, default=100,
						help='threshold (number of hits) for report repetitive area in gff')
	parser.add_argument('-ou','--output',type=str, default="output.csv",
						help='output profile table name')
	parser.add_argument('-ouf','--output_gff',type=str, default="output.gff",
                                                help='output gff format')
	parser.add_argument('-oup','--output_pic',type=str, default="output.png",
                                                help='output profile graph name')


	args = parser.parse_args()
	main(args)

print((time.time() - t0))
#plt.show()
