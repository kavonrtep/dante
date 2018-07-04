#!/usr/bin/env python3
""" sequence repetitive profile conversion to GFF3 format """

import time
import configuration
from operator import itemgetter
from itertools import groupby
from tempfile import NamedTemporaryFile
import os


def create_gff(THRESHOLD, THRESHOLD_SEGMENT, OUTPUT_GFF, files_dict, headers):    
	seq_id = None
	exclude = set(['ALL'])
	gff_tmp_list = []
	for repeat in sorted(set(files_dict.keys()).difference(exclude)):
		gff_tmp = NamedTemporaryFile(delete=False)
		pos_seq_dict = {}
		with open(files_dict[repeat][0], "r") as repeat_f,  open(files_dict[repeat][2], "r") as quality_f:
			for line_r, line_q in zip(repeat_f, quality_f):
				if "chrom" in line_r:
					idx_ranges(THRESHOLD_SEGMENT, seq_id, gff_tmp, configuration.REPEATS_FEATURE, repeat, pos_seq_dict)
					seq_id = line_r.rstrip().split("chrom=")[1]
					pos_seq_dict = {}
				else:
					hits = int(line_r.rstrip().split("\t")[1])
					if hits >= THRESHOLD:
						position = int(line_r.rstrip().split("\t")[0])
						pos_seq_dict[position] = line_q.rstrip().split("\t")[1]
		idx_ranges(THRESHOLD_SEGMENT, seq_id, gff_tmp, configuration.REPEATS_FEATURE, repeat, pos_seq_dict)
		gff_tmp_list.append(gff_tmp.name)
		gff_tmp.close()
	sort_records(gff_tmp_list, headers, OUTPUT_GFF)
	for tmp in gff_tmp_list:
		os.unlink(tmp)
		
                                
def idx_ranges(THRESHOLD_SEGMENT, seq_id, gff_file, feature, repeat, pos_seq_dict):
	indices = sorted(pos_seq_dict.keys(), key=int)
	for key, group in groupby(enumerate(indices), lambda index_item: index_item[0] - index_item[1]):
		group = list(map(itemgetter(1), group))
		if len(group) > THRESHOLD_SEGMENT:
			sum_qual = 0
			for position in group:
				sum_qual += int(pos_seq_dict[position])
			qual_per_reg = sum_qual/len(group)
			# Take boundaries of the group vectors
			gff_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tName={};Average_PID={}\n".format(seq_id, configuration.SOURCE_PROFREP, feature, group[0], group[-1], configuration.GFF_EMPTY, configuration.GFF_EMPTY, configuration.GFF_EMPTY, repeat, round(qual_per_reg)).encode("utf-8"))
			
def idx_ranges_N(indices, THRESHOLD_SEGMENT, seq_id, gff_file, feature, att_name):
	for key, group in groupby(enumerate(indices), lambda index_item: index_item[0] - index_item[1]):
		group = list(map(itemgetter(1), group))
		if len(group) > THRESHOLD_SEGMENT:
			# Take boundaries of the group vectors
			gff_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tName={}\n".format(seq_id, configuration.SOURCE_PROFREP, feature, group[0], group[-1], configuration.GFF_EMPTY, configuration.GFF_EMPTY, configuration.GFF_EMPTY, att_name))


def sort_records(gff_tmp_list, headers, OUTPUT_GFF):
	opened_files = [open(i, "r") for i in gff_tmp_list]
	files_lines = dict((key,"") for key in opened_files)
	count_without_line = 0 
	count_seq = 0
	####################################################################
	present_seqs = headers
	####################################################################
	with open(OUTPUT_GFF, "w") as final_gff:
		final_gff.write("{}\n".format(configuration.HEADER_GFF))
		while True:
			for file_name in opened_files:
				if not files_lines[file_name]:
					line = file_name.readline()
					if line:
						files_lines[file_name] = line
					else:
						count_without_line += 1
			if count_without_line == len(opened_files):
				break
			count_without_line = 0
			count = 0
			lowest_pos = float("inf")
			for file_key in files_lines.keys():
				if files_lines[file_key].split("\t")[0] == present_seqs[count_seq]:
					count += 1
					start_pos = int(files_lines[file_key].split("\t")[3])
					if start_pos < lowest_pos:
						lowest_pos = start_pos
						record_to_write = files_lines[file_key]
						lowest_file_key = file_key
			if count == 0:
				count_seq += 1
			else:
				final_gff.write(record_to_write)
				files_lines[lowest_file_key] = ""

	
def main():
	pass

if __name__ == "__main__":
	main()
