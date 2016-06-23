#!/usr/bin/env python3
""" sequence repetitive profile conversion to GFF3 format """

import numpy as np
from operator import itemgetter
from itertools import groupby
 
def main(OUTPUT, THRESHOLD):
	# Predefine the standard columns of GFF3 format
	SOURCE = "Profrep"
	FEATURE = "Repeat"
	SCORE = "."
	STRAND = "."
	FRAME = "."
	# Find the position of header(s) 
	PATTERN = ">"
	header = []
	with open(OUTPUT, "r") as input_profile:
		for num, line in enumerate(input_profile): 
			if PATTERN in line:
				header.append([num, line.rstrip().split(",")])
		csv_length = num + 1 
	header.append([csv_length])
	annotation_output = np.genfromtxt(OUTPUT, dtype=int, delimiter=",")
	for fasta_counter in range(len(header) - 1):
		fasta_start = header[fasta_counter][0]
		fasta_end = header[fasta_counter + 1][0]
		fasta_id = header[fasta_counter][1][0]
		# Choose the profile for certain repetitive class, "all" profile is omitted
		for column in range(2, annotation_output.shape[1]):
			repetitive_class = header[fasta_counter][1][column]
			nonzero_records = np.where(annotation_output[fasta_start+1 : fasta_end, column] > THRESHOLD)[0] + 1
			if nonzero_records.any():						
				ranges = []
				# Take the vector of nonzero records (above threshold) and group it to vectors of consecutive sequences, one group is one record in gff output
				for key, group in groupby(enumerate(nonzero_records), lambda index_item: index_item[0] - index_item[1]):
					group = list(map(itemgetter(1), group))
					# Take boundaries of the group vectors
					ranges.append(group[0])
					ranges.append(group[-1])
				with open("{}.gff".format(OUTPUT.split(".")[0]), "a") as gff:	
					for i in range(len(ranges) - 1)[::2]:
						 gff.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tName={}\n".format(fasta_id, SOURCE, FEATURE, ranges[i], ranges[i + 1], SCORE, STRAND, FRAME, repetitive_class))


if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument("-ou","--output",type=str, required=True,
						help="output profile table name")
	parser.add_argument("-th","--threshold",type=int, default=100,
						help="threshold (number of hits) for report repetitive area in gff")
	main(OUTPUT, THRESHOLD)
