#!/usr/bin/env python3

import argparse
import numpy as np
import os
from tempfile import NamedTemporaryFile
from collections import defaultdict
import configuration
import sys


def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected')


def check_dom_gff(DOM_GFF):
	with open(DOM_GFF) as domains_f:
		line = domains_f.readline()
		while line.startswith("#"):
			line = domains_f.readline()
		print(line)
		if len(line.split("\t")) == 9 and "Final_Classification" in line:
			pass
		else:
			raise IOError("There was detected an input GFF file that does not contain domains. Please check it and choose the domains GFF file")
			

def create_dom_dict(DOM_GFF):
	check_dom_gff(DOM_GFF)
	dict_domains = defaultdict(list)
	count_comment = check_file_start(DOM_GFF)
	with open(DOM_GFF, "r") as domains:
		for comment_idx in range(count_comment):
			next(domains)
		for line in domains:
			seq_id = line.split("\t")[0]
			ann_dom_lineage = line.split("\t")[8].split(";")[1].split("=")[-1]
			start_dom = int(line.split("\t")[3])      
			end_dom = int(line.split("\t")[4])
			strand_dom = line.split("\t")[6]
			dict_domains[seq_id].append((start_dom, end_dom, ann_dom_lineage, strand_dom))
	return dict_domains
	

def check_file_start(gff_file):
	count_comment = 0
	with open(gff_file, "r") as gff_all:
		line = gff_all.readline()
		while line.startswith("#"):
			line = gff_all.readline()
			count_comment += 1 
	return count_comment

	
def refining_intervals(OUT_REFINED, gff_removed, GAP_TH, DOM_NUM, dict_domains, domains_classes):
	with open(OUT_REFINED, "w") as joined_intervals:
		start_line = check_file_start(gff_removed)
		with open(gff_removed, "r") as repeats:
			for comment_idx in range(start_line):
				joined_intervals.write(repeats.readline())
			joined = False
			inicial = repeats.readline()
			print(inicial)
			if inicial is not "":
				seq_id_ini = inicial.rstrip().split("\t")[0] 
				start_ini = int(inicial.rstrip().split("\t")[3])
				end_ini = int(inicial.rstrip().split("\t")[4]) 
				ann_ini = inicial.rstrip().split("\t")[8].split(";")[0].split("=")[-1]
				starts_part = [start_ini]
				ends_part = [end_ini]
				for line in repeats:
					seq_id = line.rstrip().split("\t")[0]
					start = int(line.rstrip().split("\t")[3])
					end = int(line.rstrip().split("\t")[4])
					ann = line.rstrip().split("\t")[8].split(";")[0].split("=")[-1]
					if seq_id != seq_id_ini:
						if dict_domains and joined is True and configuration.WITH_DOMAINS in ann_ini:
							dict_domains = dom_refining(joined_intervals, dict_domains, seq_id_ini, start_ini, end_ini, ann_ini, starts_part, ends_part, DOM_NUM, domains_classes)
						else:
							joined_intervals.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tName={}\n".format(seq_id_ini, configuration.SOURCE_PROFREP, configuration.REPEATS_FEATURE, start_ini, end_ini, configuration.GFF_EMPTY, configuration.GFF_EMPTY, configuration.GFF_EMPTY, ann_ini))
						starts_part = [start]
						ends_part = [end]	
						joined = False
						seq_id_ini = seq_id
						start_ini = start
						end_ini = end
						ann_ini = ann
					elif start - end_ini <= GAP_TH and ann == ann_ini:
						starts_part.append(start)
						ends_part.append(end)
						end_ini = end
						joined = True
					else:
						if dict_domains and joined is True and configuration.WITH_DOMAINS in ann_ini:
							dict_domains = dom_refining(joined_intervals, dict_domains, seq_id, start_ini, end_ini, ann_ini, starts_part, ends_part, DOM_NUM, domains_classes)
						else:
							joined_intervals.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tName={}\n".format(seq_id, configuration.SOURCE_PROFREP, configuration.REPEATS_FEATURE, start_ini, end_ini, configuration.GFF_EMPTY, configuration.GFF_EMPTY, configuration.GFF_EMPTY, ann_ini))
						starts_part = [start] 
						ends_part = [end]
						joined = False
						start_ini = start
						end_ini = end
						ann_ini = ann	
				if dict_domains and joined is True and configuration.WITH_DOMAINS in ann_ini:
					dict_domains = dom_refining(joined_intervals, dict_domains, seq_id, start_ini, end_ini, ann_ini, starts_part, ends_part, DOM_NUM, domains_classes)
				else:
					joined_intervals.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tName={}\n".format(seq_id, configuration.SOURCE_PROFREP, configuration.REPEATS_FEATURE, start_ini, end_ini, configuration.GFF_EMPTY, configuration.GFF_EMPTY, configuration.GFF_EMPTY, ann_ini))
				del(starts_part)
				del(ends_part)


def create_clusters(REPEATS_GFF, BORDERS):
	start_line = check_file_start(REPEATS_GFF)
	cluster = []
	end_cluster = 0
	seq_id_ini = None
	gff_removed_parts = NamedTemporaryFile(delete=False)
	with open(REPEATS_GFF, "r") as repeats:
		for comment_idx in range(start_line):
			gff_removed_parts.write(repeats.readline().encode("utf-8"))
		for interval in repeats:
			seq_id = interval.split("\t")[0]
			start = int(interval.split("\t")[3])
			end = int(interval.split("\t")[4])
			ann = interval.split("\t")[8].split(";")[0].split("=")[-1]
			pid = int(interval.rstrip().split("\t")[8].split(";")[1].split("=")[-1].rstrip("%"))
			if seq_id != seq_id_ini:
				if len(cluster) > 1:
					remove_low_qual(cluster, BORDERS, gff_removed_parts)
				elif cluster:
					gff_removed_parts.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tName={};Average_PID={}%\n".format(cluster[0][0], configuration.SOURCE_PROFREP, configuration.REPEATS_FEATURE, cluster[0][1], cluster[0][2], configuration.GFF_EMPTY, configuration.GFF_EMPTY, configuration.GFF_EMPTY, cluster[0][3], cluster[0][4]).encode("utf-8"))
				seq_id_ini = seq_id
			if start <= end_cluster or end_cluster == 0:
				cluster.append([seq_id, start, end, ann, pid])
				if end > end_cluster:
					end_cluster = end
			else:
				if len(cluster) > 1:			
					remove_low_qual(cluster, BORDERS, gff_removed_parts)
				elif cluster:
					gff_removed_parts.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tName={};Average_PID={}%\n".format(cluster[0][0], configuration.SOURCE_PROFREP, configuration.REPEATS_FEATURE, cluster[0][1], cluster[0][2], configuration.GFF_EMPTY, configuration.GFF_EMPTY, configuration.GFF_EMPTY, cluster[0][3], cluster[0][4]).encode("utf-8"))
				end_cluster = end
				cluster = [[seq_id, start, end, ann, pid]]
	if cluster:
		remove_low_qual(cluster, BORDERS, gff_removed_parts)
	gff_removed_parts.close()
	return gff_removed_parts.name


def remove_low_qual(cluster, BORDERS, gff_removed_parts):
	score_sorted = sorted(cluster, key=lambda x: int(x[4]), reverse=True)
	for interval in score_sorted:
		start_toler = interval[1] - BORDERS
		end_toler =  interval[2] + BORDERS
		score_best = interval[4]
		## difference
		for item in score_sorted[:]:
			if item[1] >= start_toler and item[2] <= end_toler:
				if item[4] <= 0.95 * score_best: ### 5% off
					score_sorted.remove(item)
	position_sorted = sorted(score_sorted, key=lambda x: int(x[1]))
	for interval_refined in position_sorted:
		gff_removed_parts.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tName={};Average_PID={}%\n".format(interval_refined[0], configuration.SOURCE_PROFREP, configuration.REPEATS_FEATURE, interval_refined[1], interval_refined[2], configuration.GFF_EMPTY, configuration.GFF_EMPTY, configuration.GFF_EMPTY, interval_refined[3], interval_refined[4]).encode("utf-8"))

	

def dom_refining(joined_intervals, dict_domains, seq_id, start_ini, end_ini, ann_ini, starts_part, ends_part, DOM_NUM, domains_classes):
	count_dom = 0
	index_dom = 0
	strands = []
	indices_del = []
	for dom_attributes in dict_domains[seq_id]:
		ann_dom = dom_attributes[2]
		if dom_attributes[0] >= start_ini and dom_attributes[1] <= end_ini:
			repeat_class = "|".join(ann_ini.split("|")[1:])
			if ann_dom in domains_classes and ann_dom == repeat_class:
				strands.append(dom_attributes[3])
				count_dom += 1
			index_dom += 1 
	if len(set(strands)) <= 1 and count_dom >= DOM_NUM:
		joined_intervals.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tName={}\n".format(seq_id, configuration.SOURCE_PROFREP, configuration.REPEATS_FEATURE, start_ini, end_ini, configuration.GFF_EMPTY, configuration.GFF_EMPTY, configuration.GFF_EMPTY, ann_ini))
	else:
		for part_start, part_end in zip(starts_part, ends_part):
			joined_intervals.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tName={}\n".format(seq_id, configuration.SOURCE_PROFREP, configuration.REPEATS_FEATURE, part_start, part_end, configuration.GFF_EMPTY, configuration.GFF_EMPTY, configuration.GFF_EMPTY, ann_ini))
	del(dict_domains[seq_id][0:index_dom])
	return dict_domains
	
	
def get_unique_classes(CLASS_TBL):
	unique_classes = []
	with open(CLASS_TBL, "r") as class_tbl:
		for line in class_tbl:
			line_class = "|".join(line.rstrip().split("\t")[1:])
			if line_class not in unique_classes:
				unique_classes.append(line_class)
	return unique_classes

	
def main(args):
	# Command line arguments
	REPEATS_GFF = args.repeats_gff
	DOM_GFF = args.domains_gff
	GAP_TH = args.gap_threshold
	DOM_NUM = args.dom_number
	OUT_REFINED = args.out_refined 
	INCLUDE_DOM = args.include_dom
	CLASS_TBL = args.class_tbl
	BORDERS = args.borders
	
	#if CLASS_TBL and os.path.isdir(CLASS_TBL):
		#CLASS_TBL = os.path.join(CLASS_TBL, configuration.CLASS_FILE)

	gff_removed = create_clusters(REPEATS_GFF, BORDERS)

	if INCLUDE_DOM:
		unique_classes = get_unique_classes(CLASS_TBL)
		dict_domains = create_dom_dict(DOM_GFF)
		joined_intervals = refining_intervals(OUT_REFINED, gff_removed, GAP_TH, DOM_NUM, dict_domains, unique_classes)
	else:
		joined_intervals = refining_intervals(OUT_REFINED, gff_removed, GAP_TH, DOM_NUM, None, None)
	os.unlink(gff_removed)

		
if __name__ == "__main__":

    # Command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-rep_gff', '--repeats_gff', type=str, required=True,
						help='query sequence to be processed')
    parser.add_argument('-dom_gff', '--domains_gff', type=str,
						help='query sequence to be processed')
    parser.add_argument('-gth', '--gap_threshold', type=int, default=250,
						help='query sequence to be processed')
    parser.add_argument('-or', '--out_refined', type=str, default="output_refined.gff",
						help='query sequence to be processed')
    parser.add_argument('-dn', '--dom_number', type=int, default=2,
                        help='number of domains present to confirm one element type')
    parser.add_argument('-id', '--include_dom', type=str2bool, default=False,
						help='Include domains information to refine repeats regions output')
    parser.add_argument('-ct', '--class_tbl',
						help='Classification table to check the level of classification')
    parser.add_argument('-br', '--borders', type=int, default=10,
						help='Number of bp to be tolerated from one or the other side of two overlaping regions when comparing quality ')
    args = parser.parse_args()
    main(args)
