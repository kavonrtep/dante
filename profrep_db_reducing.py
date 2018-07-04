#!/usr/bin/env python3

import argparse
import subprocess
import re
from tempfile import NamedTemporaryFile
import os
import configuration


def group_reads(reads_files_list, IDENTITY_TH):
	''' Reduce the number of reads separately for each significant cluster based on similarities between them using cd-hit tool. 
		cd-hit produces reduced reads files containing only the representative reads.  
	'''
	reads_seqs_cl_reduced_list = []
	## Run cd-hit on each cluster separately
	for cluster_reads in reads_files_list:
		cl = cluster_reads.split("_")[-1]
		reads_seqs_cl_reduced = NamedTemporaryFile(suffix=".reduced_{}".format(cl), delete=False)
		subprocess.call("cd-hit-est -i {} -o {} -c {} -M {}".format(cluster_reads, reads_seqs_cl_reduced.name, IDENTITY_TH, configuration.MEM_LIM), shell=True)
		reads_seqs_cl_reduced_list.append(reads_seqs_cl_reduced.name)
		reads_seqs_cl_reduced.close()
	## Return the list of reduced reads files 
	return reads_seqs_cl_reduced_list
	
def representative_reads(READS_ALL, CLS_FILE, CL_SIZE_TH, CLS_REDUCED, IDENTITY_TH):
	''' Group the reads based on the sequences similarities. 
		Replace a group by only the one representative read preserving the quantitative info how much reads it represents
		1. Loop over the original cls file and find the significant clusters (min. number of reads) 
		2. Get the reads which are in individual significant clusters
		2. Get the reads sequences for individual clusters to run cd-hit which groups the reads for each cluster
		3. After getting all significant ones (clusters sorted by size) process the outputs from cd-hit and to get reads representations
		4. Create new cls file and write down significant clusters with the new reads IDs
		5. Continue reading unsignificant original cls file and copy the rest of clusters to the new cls unchanged 
	Find groups of similar reads and replace them with only one representative also preserving the number of reads it represents	
	'''
	reads_dict = {}
	cl = None
	line_cl_header = None
	modify_files = True
	reads_files_list = []
	cls_reduced_file = open(CLS_REDUCED, "w")
	## parse file of all clusters from RE 
	with open(CLS_FILE, "r") as cls_ori:
		for line in cls_ori:
			if line.startswith(">"):
				line_cl_header = line
				## cluster number
				cl = re.split('\t| ', line)[0].rstrip().lstrip(">")	
			else:
				reads_in_cl = line.rstrip().split("\t")
				## reduce only reads in the biggest clusters:
				## the threshold of cluster size is set as a minimum number of reads it has to contain
				if len(reads_in_cl) >= CL_SIZE_TH:
					## for significant cluster create a separate file to write reads sequences
					reads_seqs_cl_orig = NamedTemporaryFile(suffix="_{}".format(cl), delete=False)	
					reads_files_list.append(reads_seqs_cl_orig.name)	
					## for every read in the cluster create entry in reads_dict to which cluster it belongs and the file of the read sequence for this cluster
					for read in reads_in_cl:
						## Dictionary of reads from significant clusters -> KEY:read_id VALUE:[number of cluster, filename to reads sequences file]
						reads_dict[read] = [cl, reads_seqs_cl_orig.name]
						reads_seqs_cl_orig.close()
				## after getting all significant clusters to be reduced (original cls file sorted by size of clusters), process the reads reads in them and write to the modified reads and cls files
				elif modify_files:
					## get reads sequences for significant clusters from ALL reads
					get_read_sequences(READS_ALL, reads_dict, reads_files_list)
					## run cd-hit to reduce the reads for significant clusters
					reads_seqs_cl_reduced_list = group_reads(reads_files_list, IDENTITY_TH)
					reads_repre_dict = {}
					# for individual significant cluster, process the corresponding file of original and reduced reads 
					for reads_seqs_cl_orig, reads_seqs_cl_reduced in zip(reads_files_list, reads_seqs_cl_reduced_list):
						cl = reads_seqs_cl_reduced.split("_")[-1]
						## get reads quantitative represantion dictionary for individual cluster
						[reads_repre_dict, reads_in_cl_mod] = process_reads_groups(reads_seqs_cl_reduced, reads_repre_dict)
						## for each significant cluster write the new IDs of reduced reads to modified cls
						modify_cls(cls_reduced_file, reads_in_cl_mod, cl)
						os.unlink(reads_seqs_cl_orig)
						os.unlink(reads_seqs_cl_reduced)
					## write the last line that was chcecked but not reduced
					cls_reduced_file.write(line_cl_header)
					cls_reduced_file.write(line)
					modify_files = False
				## after reducing append the rest of clusters unchanged to the modified cls file
				else:
					cls_reduced_file.write(line_cl_header)
					cls_reduced_file.write(line)
	cls_reduced_file.close()
	return reads_repre_dict
		
		
def modify_cls(cls_reduced_file, reads_in_cl_mod, cl):
	''' For each significant cluster write down the new adjusted names of reads
	'''
	num_of_reads = len(reads_in_cl_mod)
	cls_reduced_file.write(">{}\t{}\n".format(cl, num_of_reads))
	cls_reduced_file.write("{}\n".format("\t".join(reads_in_cl_mod)))
	

def get_read_sequences(READS_ALL, reads_dict, reads_files_list):
	'''From file of ALL reads sequences take only the ones belonging to significant clusters.
	Distribute them to separate files based on the cluster they belong to
	'''
	with open(READS_ALL, "r") as reads_all:
		for line in reads_all:
			if line.startswith(">"):
				read = line.rstrip().lstrip(">")
				## check if the read belong to significant cluster
				if read in reads_dict.keys():	
					## then write it to file of reads for corresponding cluster 
					with open(reads_dict[read][1], "a") as reads_file:
						reads_file.write(line)
						reads_file.write(reads_all.readline())
	print(reads_files_list)
		

def process_reads_groups(reads_seqs_cl_reduced, reads_repre_dict):
	''' Process the .clstr output of cd-hit which contains groups of original reads
		Each group starts with > character, on separate lines are listed original reads IDs, the representative one is marked by *.
		Get the number of reads in every group to preserve the quantitative information
		Create dictionary of representative reads as keys and the amount of reads they represent as value:
		value = 0 : read is not representative and will not take place in the reduce database
		value > 0 : value indicates the number of reads it represents 
		Create list of new representative reads IDs encoding the original number of read in the group using 'reduce' tag:
				 e.g. 171freduce10 (10 original reads were reduced to one 171f representative)
	'''
	clstr_file = "{}.clstr".format(reads_seqs_cl_reduced)
	reads_in_cl_mod = []
	read_represent = ""
	with open(clstr_file, "r") as clstr:
		for line in clstr:
			count_reads = 0
			while not line.startswith(">"):
				if not line: 
					break
				read = line.split(" ")[1].split(".")[0]
				count_reads += 1
				if line.rstrip().endswith("*"):
					read_represent = read
				else:
					reads_repre_dict[read] = 0
				line = clstr.readline()
			if read_represent:	
				reads_repre_dict[read_represent] = count_reads
				reads_in_cl_mod.append("{}reduce{}".format(read_represent.rstrip().lstrip(">"), count_reads))
	return reads_repre_dict, reads_in_cl_mod
	
	
def reduce_reads(READS_ALL, READS_ALL_REDUCED, reads_repre_dict):
	''' Report a new file of reads sequences based on the original file of ALL reads using the reads representation dictionary.
		Loop over the reads in the original READS_ALL file 
		There are 3 options evaluated for the read:
			- the value in the dictionary equals to zero, read is not representative -> it will not take place in the new reads DB
			- the value is greater than zero, the read is representative -> in new read DB encode the number of representing reads using 'reduce' tag (<Original_read_ID>reduce<number_represented>)
			- the read is not in the dictionary -> add it unchanged from the original ALL reads database
	'''
	with open(READS_ALL_REDUCED, "w") as reads_all_red:
		with open(READS_ALL, "r") as reads_all_ori:
			for line in reads_all_ori:
				if line.startswith(">"):
					if line.rstrip() in reads_repre_dict:
						amount_represented = reads_repre_dict[line.rstrip()]
						if  amount_represented > 0:
							reads_all_red.write("{}reduce{}\n".format(line.rstrip(), amount_represented))
							reads_all_red.write(reads_all_ori.readline())
					else:
						reads_all_red.write(line)
						reads_all_red.write(reads_all_ori.readline())

def main(args):
	CLS_FILE = args.cls
	READS_ALL = args.reads_all
	CL_SIZE_TH = args.cluster_size
	IDENTITY_TH = args.identity_th
	CLS_REDUCED = args.cls_reduced
	READS_ALL_REDUCED = args.reads_reduced
	
	
	if not os.path.isabs(CLS_REDUCED):
		CLS_REDUCED = os.path.join(os.getcwd(), CLS_REDUCED)
		print(CLS_REDUCED)
		
	if not os.path.isabs(READS_ALL_REDUCED):
		READS_ALL_REDUCED = os.path.join(os.getcwd(), READS_ALL_REDUCED)
		print(READS_ALL_REDUCED)
	
	reads_repre_dict = representative_reads(READS_ALL, CLS_FILE, CL_SIZE_TH, CLS_REDUCED, IDENTITY_TH)
	reduce_reads(READS_ALL, READS_ALL_REDUCED, reads_repre_dict)


if __name__ == '__main__':

	
	# Command line arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('-r', '--reads_all', type=str, required=True,
						help='input file containing all reads sequences')
	parser.add_argument('-c', '--cls', type=str, required=True,
						help='input sorted cls file containing reads assigned to clusters')
	parser.add_argument('-rr', '--reads_reduced', type=str, default=configuration.READS_ALL_REDUCED,
						help='output file containing reduced number of reads')
	parser.add_argument('-cr', '--cls_reduced', type=str, default=configuration.CLS_REDUCED,
						help='output cls file containing adjusted clusters for the reduced reads database')					
	parser.add_argument('-i', '--identity_th', type=float, default=0.90,
						help='reads identity threshold for cdhit')
	parser.add_argument('-cs', '--cluster_size', type=int, default=1000,
						help='minimum cluster size to be included in reducing')
						
	args = parser.parse_args()
	main(args)
	
	
	
