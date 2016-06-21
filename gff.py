#!/usr/bin/env python3
""" sequence repetitive profile conversion to GFF3 format """

import numpy as np
from operator import itemgetter
from itertools import groupby

def main(OUTPUT, THRESHOLD):
	# Predefine the standard columns of GFF format
	SOURCE = "Profrep"
	FEATURE = "Repeat"
	SCORE = "."
	STRAND = "."
	FRAME = "."
	# Find the position of headers in case of multifasta file
	PATTERN = ">"
	header = []
	with open(OUTPUT, "r") as input_profile:
		for num, line in enumerate(input_profile): 
			if PATTERN in line:
				header.append([num, line.rstrip().split(",")])
				print(num)
		csv_length = num + 1 
		print(csv_length)
	header.append([csv_length])
	annotation_output = np.genfromtxt(OUTPUT, dtype=int, delimiter=",")
	for fasta_counter in range(len(header) - 1):
		print(fasta_counter)
		fasta_start = header[fasta_counter][0]
		fasta_end = header[fasta_counter + 1][0]
		fasta_id = header[fasta_counter][1][0]
		# Choose the column, e.g. profile for certain repetitive class, ALL profile is omitted
		for column in range(2, annotation_output.shape[1]):
			repetitive_class = header[fasta_counter][1][column]
			nonzero_records = np.where(annotation_output[fasta_start + 1 : fasta_end, column] > THRESHOLD)[0] + 1
			if nonzero_records.any():						
				ranges = []
				# Take the vector of nonzero records and group it to vectors of consecutive sequences
				for key, group in groupby(enumerate(nonzero_records), lambda index_item: index_item[0] - index_item[1]):
					group = list(map(itemgetter(1), group))
					# Take boundaries of the group vectors
					ranges.append(group[0])
					ranges.append(group[-1])
				print(ranges)
				with open("{}.gff".format(OUTPUT.split(".")[0]), "a") as gff:	
					for i in range(len(ranges) - 1)[::2]:
						 gff.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tName={}\n".format(fasta_id, SOURCE, FEATURE, ranges[i], ranges[i + 1], SCORE, STRAND, FRAME, repetitive_class))
	
	#with open(OUTPUT,"r") as input_profile:
		#header=input_profile.readline().rstrip().rsplit(",")
		#annotation_output=np.loadtxt(OUTPUT, dtype=int, skiprows=1, delimiter=",")
		
		## Choose the column - profile for certain repetitive class, ALL profile is omitted
		#for column in range(2, annotation_output.shape[1]):		
			#nonzero_records=np.where(annotation_output[:, column] > THRESHOLD)[0]+1	
			#if nonzero_records.any():						
				#ranges = []
				## Take the vector of nonzero records and group it to the vectors of consecutive sequences
				#for key, group in groupby(enumerate(nonzero_records), lambda index_item: index_item[0] - index_item[1]):
					#group = list(map(itemgetter(1), group))
					## Take boundaries of the group vectors
					#ranges.append(group[0])
					#ranges.append(group[-1])
				#print(ranges)
				#with open("{}.gff".format(OUTPUT.split(".")[0]), "a") as gff:	
					#for i in range(len(ranges) - 1)[::2]:
						 #gff.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tName={}\n".format(header[0], SOURCE, FEATURE, ranges[i], ranges[i+1], SCORE, STRAND, FRAME, header[column]))


if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument("-ou","--output",type=str, required=True,
						help="output profile table name")
	parser.add_argument("-th","--threshold",type=int, default=100,
						help="threshold (number of hits) for report repetitive area in gff")
	OUTPUT = "output2.csv"
	THRESHOLD = 100
	main(OUTPUT, THRESHOLD)
