#!/usr/bin/env python3

import numpy as np
import subprocess
import csv
import time
import sys
import matplotlib.pyplot as plt
import multiprocessing
from multiprocessing import Pool
from functools import partial
import argparse
from tempfile import NamedTemporaryFile
#import GFF

t0 = time.time()

# setting to print out the whole range of arrays (for command line stdouts)
np.set_printoptions(threshold=np.nan)

# define arguments parsing
parser = argparse.ArgumentParser()
parser.add_argument('-q','--query', type=str, required=True,
                    help='query sequence to be processed')
parser.add_argument('-d','--database', type=str, required=True,
                    help='database for blast - all reads')
parser.add_argument('-a','--annotation_tbl', type=str, required=True,
                    help='annotation table')
parser.add_argument('-c','--cls', type=str, required=True,
                    help='cls file containing reads assigned to clusters')
parser.add_argument('-i','--identical', type=float, default=90,
                    help='blast filtering option: sequence indentity threshold between query and db in %')
parser.add_argument('-l','--align_length', type=int, default=40,
                    help='blast filtering option: minimal align length threshold in bp')
parser.add_argument('-e','--e_value', type=str, default=1e-15,
                    help='blast filtering option: e-value')
parser.add_argument('-ws','--word_size', type=int, default=11,
                    help='blast filtering option: initial word size')
parser.add_argument('-t','--task', type=str, default="blastn",
                    help='blast task to be triggered')
parser.add_argument('-m','--max_alignments', type=int, default=10000000,
                    help='blast filtering option: maximal number of alignments in the output')
parser.add_argument('-w','--window', type=int, default=5000,
                    help='window size for parallel processing')
parser.add_argument('-o','--overlap', type=int, default=150,
                    help='overlap for parallely processed regions, set greater than read size')
parser.add_argument('-n','--new_db', default=False,
                    help='create a new blast database -if running the program for the first time')
parser.add_argument('-g','--gff', default=False,
                    help='use module for gff')
parser.add_argument('-th','--threshold',type=int, default=100,
                    help='threshold (number of hits) for report repetitive area in gff')
parser.add_argument('-ou','--output',type=str, default="output2.csv",
                    help='output profile table name')


# command line arguments
args = parser.parse_args()

query = args.query
bl_database = args.database
cl_annotatation_table = args.annotation_tbl
CLS = args.cls
identical = args.identical
align_length = args.align_length
e_value = args.e_value
word_size = args.word_size
window = args.window
overlap = args.overlap
task = args.task
max_alignments = args.max_alignments
gff = args.gff
threshold = args.threshold
out = args.output

# create a blast database of all reads
if args.new_db:
	subprocess.call("makeblastdb -in {} -dbtype nucl".format(bl_database), shell=True)

# create subfasta files to be processed sequentially
def multifasta(query):
	pattern = ">"
	fasta_list = []
	handle_list = []
	with open(query, "r") as fasta:
		reader = fasta.read()
		splitter = reader.split(pattern)[1:]
		if len(splitter) > 1:
			for fasta_num, part in enumerate(splitter):
				ntf = NamedTemporaryFile(delete=False)
				print(ntf.name)
				ntf.write("{}{}".format(pattern, part).encode("utf-8"))
				handle_list.append(ntf)
				fasta_list.append(ntf.name)
				ntf.close()
			return fasta_list
		else:
			fasta_list.append(query)
			return fasta_list


# read fasta
def fasta_read(subfasta):
	print(subfasta)
	sequence_lines = []
	with open(subfasta, "r") as fasta:
		print(fasta)
		header = fasta.readline().strip()
		print(header)
		for line in fasta:
			clean_line = line.strip()		# strip empty spaces
			if clean_line:				# is not empty line
				sequence_lines.append(clean_line)
	sequence = "".join(sequence_lines)
	return header, sequence


# dictionary of cluster annotation
def cluster_annotation(cl_annotatation_table):
	cl_annotations = {} 		# input database of repetitions with coresponding clusters
	## loading tables as 2D arrays
	annot_table = np.genfromtxt(cl_annotatation_table, dtype=str) 	# table of cluster annotations
	## check the input
	for line in annot_table:
		if line[1] in cl_annotations:
			cl_annotations[line[1]].append(line[0])
		else:
			cl_annotations[line[1]] = [line[0]]
	return cl_annotations


# dictionary of read annotations
def read_annotation(CLS, cl_annotations):
	reads_annotations = {} 	# dictionary of reads and corresponding annotations
	with open(CLS, "r") as cls_file:
		count = 0
		for line in cls_file:
			count += 1
			if count%2 == 0:
				reads = line.rstrip().split(" ") # remove newline and devide on reads
				for element in reads: # for every read
					for key, value in list(cl_annotations.items()):
						if clust in value:
							reads_annotations[element] = key # assign annotation from anotation dictionary
			else:
				clust = line.split(" ")[0][3:] # take only the number of a cluster
	return reads_annotations


# dictionary of the sequence repetitive profile
def annot_profile(cl_annotations, part):
	profile = {} 				# output profile
	# fill profile dictionary with the repetition class names and corresponding vector of length of the query sequence for aligned reads mapping within each repetition class
	for key in list(cl_annotations.keys()):
		profile[key] = np.zeros(part, dtype=int)
	profile["repeat"] = np.zeros(part, dtype=int)
	profile["all"] = np.zeros(part, dtype=int)
	return profile


# parallel process
def processInput(subset):
	loc_start = subset + 1
	loc_end = subset + window
	if loc_end > query_len:
		loc_end = query_len
		print(query_len - loc_start + 1)
		profile = annot_profile(cl_annotations, query_len - loc_start + 1)
	else:
		profile = annot_profile(cl_annotations, window + 1)
	# find HSP records using blastn for every window defined by location and parse the tabular stdout -> 1. query, 2. database read, 3. %identical, 4. alignment length, 5. alignment start, 6. alignment end
	p = subprocess.Popen("blastn -query {} -query_loc {}-{} -db {} -evalue {} -word_size {} -task {} -num_alignments {} -outfmt '6 qseqid sseqid pident length qstart qend'".format(subfasta, loc_start,loc_end, bl_database, e_value,word_size, task, max_alignments), stdout=subprocess.PIPE, shell=True)
	for line in p.stdout:
		column = line.decode("utf8").rstrip().split("\t")
		if float(column[2]) >= identical and int(column[3]) >= align_length:
			read = column[1]				# ID of each individual aligned read
			qstart = int(column[4])		# starting position of alignment
			qend = int(column[5])			# ending position of alignemnt

			# assign repetition class to the cluster, if does not belong to any - assign to "all"
			if read in reads_annotations:
				annotation = reads_annotations[read]								# type of repetition
			else:
				annotation = "repeat"
			#print(annotation)
			profile[annotation][qstart-subset-1 : qend-subset-1] = profile[annotation][qstart-subset-1 : qend-subset-1] + 1
	return profile


# sum up values for the same keys on list of dictionaries
def conc_dict(profile_list):
	print(profile_list)
	profile_output = {}
	for k in profile_list[0].keys():
		if len(profile_list) == 1:
			profile_output[k] = profile_list[0][k]
		else:
			profile_output[k] = profile_list[0][k][0 : window-overlap//2]
			for x in profile_list[1:-1]:
				profile_output[k] = np.append(profile_output[k], x[k][overlap//2 : window-overlap//2])
				print("+++++++")
			profile_output[k] = np.append(profile_output[k], profile_list[-1][k][overlap//2:])
			print(len(profile_output[k]))
	return profile_output


# create matrix of values from the profile dictionary !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
def profile_to_csv(profile_output, ax, cm):
	data = np.zeros((len(profile_output), query_len), dtype=int)
	labels = []
	count = 0
	for key in list(profile_output.keys()):
		labels.append(key) 					# create list of labels from the dicitonary keys
		data[count] = profile_output[key]		# every line contains array of values for one repetition class
		if np.any(data[count]):
			ax.plot(list(range(query_len)), data[count], label=labels[count])
		#ax.set_color_cycle([cm(1.*count/len(profile))])
		if key == "all":
			all_position = count 				# save poisiton of "all" record
		count += 1
	plt.legend(loc="upper left")

	# swap positions so that "all" record is the first in the table
	data[0], data[all_position] = data[all_position], data[0].copy()
	labels[0], labels[all_position] = labels[all_position], labels[0]

	head_line = np.array(list(range(1, query_len+1)))
	labels = np.append([header], labels, 1)
	data = np.append([head_line], data, 0)

	# write output to a csv table
	output = open(out, "a")
	writer = csv.writer(output)
	writer.writerow(labels)
	for values in np.transpose(data):		# table transposed
		writer.writerow(values)
	return output


# create a profile graph
def set_visualization():
	fig = plt.figure()
	ax = fig.add_subplot(111)
	cm = plt.get_cmap("gist_rainbow")
	return ax, cm


fasta_list = multifasta(query)
#while fasta_counter<=len(fast_list):
for subfasta in fasta_list:
	[header,sequence] = fasta_read(subfasta)
	query_len = len(sequence)

	## cluster annotations
	cl_annotations = cluster_annotation(cl_annotatation_table)

	## dictionary of individual reads annotations to repetitive class
	reads_annotations = read_annotation(CLS, cl_annotations)

	## set up the parallel process
	step = window - overlap																# moving step
	num_cores = multiprocessing.cpu_count()												# number of processors
	par = Pool(num_cores)																# create pool
	subset = list(range(0, query_len, step))
	print(subset)																		# define vector of indeces for parallel functio
	#func=partial(processInput,len(subset))										# passing more input arguments to parallel function - empty profile and number of fragments
	profile_list = par.map(processInput, subset)			 										# create parallel process - output list of dictionaries for individual windows
	#print(profile_list)

	## merge output dictionaries of parallel process = sum up values fro the same keys
	profile_output = conc_dict(profile_list)

	## sum the profile counts to get "all" profile including all repetions and also hits not belonging anywhere
	profile_output["all"] = sum(profile_output.values())

	## visualization
	[ax, cm] = set_visualization()

	## output csv+plot
	output = profile_to_csv(profile_output, ax, cm)


# plot the profile graph
output.close()
if gff:
	gff_output = GFF.rep_regions(out, threshold)
	print(gff_output)
print((time.time() - t0))
plt.show()
