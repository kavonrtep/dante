#!/usr/bin/env python3

import argparse
import subprocess
import re
from tempfile import NamedTemporaryFile
import os


def find_representative_reads(reads_files_list, IDENTITY_TH):
	reads_cl_reduced_files = []
	for cluster_reads in reads_files_list:
		cl = cluster_reads.split("_")[-1]
		reads_cl_reduced = NamedTemporaryFile(suffix=".reduced_{}".format(cl), delete=False)
		subprocess.call("cd-hit-est -i {} -o {} -c {} -M {}".format(cluster_reads, reads_cl_reduced.name, IDENTITY_TH, configuration.MEM_LIM), shell=True)
		reads_cl_reduced_files.append(reads_cl_reduced.name)
		reads_cl_reduced.close()
	return reads_cl_reduced_files
	
def reads_cl(READS_ALL, CLS_FILE, CL_SIZE_TH, CLS_REDUCED, IDENTITY_TH):
	cls_dict = {}
	cl = None
	line_cl_header = None
	modify_files = True
	reads_files_list = []
	# write or append mode?
	cls_reduced_file = open(CLS_REDUCED, "w")
	with open(CLS_FILE, "r") as cls_ori:
		for line in cls_ori:
			if line.startswith(">"):
				line_cl_header = line
				cl = re.split('\t| ', line)[0].rstrip().lstrip(">")	
			else:
				reads_cl = line.rstrip().split("\t")
				if len(reads_cl) >= CL_SIZE_TH:
					reads_cl_file = NamedTemporaryFile(suffix="_{}".format(cl), delete=False)	
					reads_files_list.append(reads_cl_file.name)		
					for read in reads_cl:
						cls_dict[read] = [cl, reads_cl_file.name]
						reads_cl_file.close()
					#[reads_repre_dict, reads_cl_mod] = process_similarity_clusters(reads_cl_reduced, reads_repre_dict)
					#modify_cls(cls_reduced_file, reads_cl_mod, cl)
				## while?
				elif modify_files:
					get_read_sequences(READS_ALL, cls_dict, reads_files_list)
					reads_cl_reduced_files = find_representative_reads(reads_files_list, IDENTITY_TH)
					reads_repre_dict = {}
					for reads_cl_file, reads_cl_reduced in zip(reads_files_list, reads_cl_reduced_files):
						cl = reads_cl_reduced.split("_")[-1]
						[reads_repre_dict, reads_cl_mod] = process_similarity_clusters(reads_cl_reduced, reads_repre_dict)
						modify_cls(cls_reduced_file, reads_cl_mod, cl)
						os.unlink(reads_cl_file)
						os.unlink(reads_cl_reduced)
					cls_reduced_file.write(line_cl_header)
					cls_reduced_file.write(line)
					modify_files = False
				else:
					cls_reduced_file.write(line_cl_header)
					cls_reduced_file.write(line)
	cls_reduced_file.close()
	return reads_repre_dict
		
		
def modify_cls(cls_reduced_file, reads_cl_mod, cl):
	num_of_reads = len(reads_cl_mod)
	cls_reduced_file.write(">{}\t{}\n".format(cl, num_of_reads))
	cls_reduced_file.write("{}\n".format("\t".join(reads_cl_mod)))
	

def get_read_sequences(READS_ALL, cls_dict, reads_files_list):
	with open(READS_ALL, "r") as reads_all:
		for line in reads_all:
			if line.startswith(">"):
				read = line.rstrip().lstrip(">")
				if read in cls_dict.keys():
					with open(cls_dict[read][1], "a") as reads_file:
						reads_file.write(line)
						reads_file.write(reads_all.readline())
					# delete record from dict?
	print(reads_files_list)
		

def process_similarity_clusters(reads_cl_reduced, reads_repre_dict):
	clstr_file = "{}.clstr".format(reads_cl_reduced)
	reads_cl_mod = []
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
				reads_cl_mod.append("{}_{}".format(read_represent.rstrip().lstrip(">"), count_reads))
	return reads_repre_dict, reads_cl_mod
	
	
def reduce_reads(READS_ALL, READS_ALL_REDUCED, reads_repre_dict):
	with open(READS_ALL_REDUCED, "w") as reads_all_red:
		with open(READS_ALL, "r") as reads_all_ori:
			for line in reads_all_ori:
				if line.startswith(">"):
					if line.rstrip() in reads_repre_dict:
						amount_represented = reads_repre_dict[line.rstrip()]
						if  amount_represented > 0:
							reads_all_red.write("{}_{}\n".format(line.rstrip(), amount_represented))
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
	
	reads_repre_dict = reads_cl(READS_ALL, CLS_FILE, CL_SIZE_TH, CLS_REDUCED, IDENTITY_TH)
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
	parser.add_argument('-cs', '--cluster_size', type=int, default=69000,
						help='minimum cluster size to be included in reducing')
						
	args = parser.parse_args()
	main(args)
	
	
	
