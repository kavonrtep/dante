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
	''' Check if the GFF file of protein domains was given '''
	with open(DOM_GFF) as domains_f:
		line = domains_f.readline()
		while line.startswith("#"):
			line = domains_f.readline()
		if len(line.split("\t")) == 9 and "Final_Classification" in line:
			pass
		else:
			raise IOError("There was detected an input GFF file that does not contain domains. Please check it and choose the domains GFF file")
			

def create_dom_dict(DOM_GFF):
	''' Create hash table of protein domains for individual sequences '''
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
	for seq_id in dict_domains.keys():
		dict_domains[seq_id] = sorted(dict_domains[seq_id], key=lambda x: x[0])
	return dict_domains
	

def check_file_start(gff_file):
	count_comment = 0
	with open(gff_file, "r") as gff_all:
		line = gff_all.readline()
		while line.startswith("#"):
			line = gff_all.readline()
			count_comment += 1 
	return count_comment

	
def interconnect_intervals(OUT_REFINED, gff_removed, GAP_TH, DOM_NUM, dict_domains, domains_classes):
	''' Second step of refining - INTERVALS INTERCONNECTING:
	Gradually checking regions of GFF which are sorted by starting point.
	Adding regions to be interconnected if the gap between the new and previous one is lower than threshold and theirclassification is the same as the one set by the first region.
	If the conditions are not fullfilled the adding is stopped and this whole expanded region is further evaluated:
	1. if domains are not included in joining (A. choosed as parameter or B. based on the classification the element does not belong to mobile elements):
		-> the fragments are joined reported as one region in refined gff
	2. if domains should be included:
		-> the fragments are joined if domains within the expanded region meets the criteria:
			1. they are at least certain number of domains 
			2. they have equal strand orientation 
			3. they are classified to the last classification level (checked from the class. table of the database) and this matches the classification of the expanded region
		-> otherwise the region is refragmented to previous parts, which are reported as they were in the original repeats gff
	'''
	with open(OUT_REFINED, "w") as joined_intervals:
		start_line = check_file_start(gff_removed)
		with open(gff_removed, "r") as repeats:
			for comment_idx in range(start_line):
				joined_intervals.write(repeats.readline())
			joined = False
			initial = repeats.readline()
			## if there are repeats in GFF, initialize
			if initial is not "":
				seq_id_ini = initial.rstrip().split("\t")[0] 
				start_ini = int(initial.rstrip().split("\t")[3])
				end_ini = int(initial.rstrip().split("\t")[4]) 
				ann_ini = initial.rstrip().split("\t")[8].split(";")[0].split("=")[-1]
				starts_part = [start_ini]
				ends_part = [end_ini]
				for line in repeats:
					seq_id = line.rstrip().split("\t")[0]
					start = int(line.rstrip().split("\t")[3])
					end = int(line.rstrip().split("\t")[4])
					ann = line.rstrip().split("\t")[8].split(";")[0].split("=")[-1]
					## if new sequence is detected the joining process is reset
					if seq_id != seq_id_ini:
						## handle the last expanded region from the previous sequence
						if dict_domains and joined is True and configuration.WITH_DOMAINS in ann_ini:
							# check the domains
							dict_domains = dom_refining(joined_intervals, dict_domains, seq_id_ini, start_ini, end_ini, ann_ini, starts_part, ends_part, DOM_NUM, domains_classes)
						else:
							# write the joined region if domains should not be considered
							joined_intervals.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tName={}\n".format(seq_id_ini, configuration.SOURCE_PROFREP, configuration.REPEATS_FEATURE, start_ini, end_ini, configuration.GFF_EMPTY, configuration.GFF_EMPTY, configuration.GFF_EMPTY, ann_ini))
						# initialize the search for expanded intervals
						starts_part = [start]
						ends_part = [end]	
						joined = False
						seq_id_ini = seq_id
						start_ini = start
						end_ini = end
						ann_ini = ann
					## prolonging the potential region to be expanded
					elif start - end_ini <= GAP_TH and ann == ann_ini:
						starts_part.append(start)
						ends_part.append(end)
						end_ini = end
						joined = True
					## end of expanding the interval
					else:
						if dict_domains and joined is True and configuration.WITH_DOMAINS in ann_ini:
							dict_domains = dom_refining(joined_intervals, dict_domains, seq_id, start_ini, end_ini, ann_ini, starts_part, ends_part, DOM_NUM, domains_classes)
						else:
							joined_intervals.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tName={}\n".format(seq_id, configuration.SOURCE_PROFREP, configuration.REPEATS_FEATURE, start_ini, end_ini, configuration.GFF_EMPTY, configuration.GFF_EMPTY, configuration.GFF_EMPTY, ann_ini))
						starts_part = [start] 
						ends_part = [end]
						# after handling the expanded interval,  set the "joined" indicator to be false, so the next region expanding and potential joining can begin  
						joined = False
						start_ini = start
						end_ini = end
						ann_ini = ann	
				## handle the last potentially expanded region
				if dict_domains and joined is True and configuration.WITH_DOMAINS in ann_ini:
					# check the domains
					dict_domains = dom_refining(joined_intervals, dict_domains, seq_id, start_ini, end_ini, ann_ini, starts_part, ends_part, DOM_NUM, domains_classes)
				else:
					# write the joined region
					joined_intervals.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tName={}\n".format(seq_id, configuration.SOURCE_PROFREP, configuration.REPEATS_FEATURE, start_ini, end_ini, configuration.GFF_EMPTY, configuration.GFF_EMPTY, configuration.GFF_EMPTY, ann_ini))
				del(starts_part)
				del(ends_part)


def cluster_regions_for_quality_check(REPEATS_GFF, BORDERS):
	''' First step of refining - REMOVING LOW CONFIDENCE REPEATS REGIONS
	Create clusters of overlapping regions from REPEATS_GFF. 
	These clusters subsequently serves for removing nested regions of different classification which have significantly lower quality comparing to the region they are nested in. 
	Overhang tolerance of couple of bases on the both sides is tolerated when evaluting the nested regions.	
	'''
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
			pid = int(interval.rstrip().split("\t")[8].split(";")[1].split("=")[-1])
			if seq_id != seq_id_ini:
				if len(cluster) > 1:
					remove_low_qual(cluster, BORDERS, gff_removed_parts)
				elif cluster:
					gff_removed_parts.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tName={};Average_PID={}\n".format(cluster[0][0], configuration.SOURCE_PROFREP, configuration.REPEATS_FEATURE, cluster[0][1], cluster[0][2], configuration.GFF_EMPTY, configuration.GFF_EMPTY, configuration.GFF_EMPTY, cluster[0][3], cluster[0][4]).encode("utf-8"))
				seq_id_ini = seq_id
				cluster = []
			## expanding the cluster
			if start <= end_cluster or end_cluster == 0:
				cluster.append([seq_id, start, end, ann, pid])
				if end > end_cluster:
					end_cluster = end
			## if the next region is not overlapping, evaluate the regions in cluster or write down the single region
			else:
				if len(cluster) > 1:			
					remove_low_qual(cluster, BORDERS, gff_removed_parts)
				elif cluster:
					gff_removed_parts.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tName={};Average_PID={}\n".format(cluster[0][0], configuration.SOURCE_PROFREP, configuration.REPEATS_FEATURE, cluster[0][1], cluster[0][2], configuration.GFF_EMPTY, configuration.GFF_EMPTY, configuration.GFF_EMPTY, cluster[0][3], cluster[0][4]).encode("utf-8"))
				end_cluster = end
				cluster = [[seq_id, start, end, ann, pid]]
	if cluster:
		remove_low_qual(cluster, BORDERS, gff_removed_parts)
	gff_removed_parts.close()
	return gff_removed_parts.name


def remove_low_qual(cluster, BORDERS, gff_removed_parts):
	''' Loop over regions in the cluster which are sorted based on descending quality (PID).
	Regions nested in the currently checked region (with some exceeding borders tolerance) are tested for quality. 
	If the quality difference between the currently highest quality and the nested region is more than a threshold nested one will be removed.
	Already checked region is removed from the list and procedure continues with the next best quality region in the cluster.
	'''
	score_sorted = sorted(cluster, key=lambda x: int(x[4]), reverse=True)
	for interval in score_sorted:
		start_toler = interval[1] - BORDERS
		end_toler =  interval[2] + BORDERS
		score_best = interval[4]
		## difference
		for item in score_sorted[:]:
			if item[1] >= start_toler and item[2] <= end_toler:
				if item[4] <= (1-configuration.QUALITY_DIFF_TO_REMOVE) * score_best: ### 5% difference tolerance
					score_sorted.remove(item)
	position_sorted = sorted(score_sorted, key=lambda x: int(x[1]))
	for interval_refined in position_sorted:
		gff_removed_parts.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tName={};Average_PID={}\n".format(interval_refined[0], configuration.SOURCE_PROFREP, configuration.REPEATS_FEATURE, interval_refined[1], interval_refined[2], configuration.GFF_EMPTY, configuration.GFF_EMPTY, configuration.GFF_EMPTY, interval_refined[3], interval_refined[4]).encode("utf-8"))

	

def dom_refining(joined_intervals, dict_domains, seq_id, start_ini, end_ini, ann_ini, starts_part, ends_part, DOM_NUM, domains_classes):
	''' Check the protein domains within the potential interval to be joined:
	1. appropriate domain has the same classification as interval of repeats to be joined and is classified to the last level 
	2. strands of appropriate domains are  uniform
	3. the count of such appropriate domains is more than a threshold
	'''
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
	## if criteria for domains are met, write the expanded interval
	if len(set(strands)) <= 1 and count_dom >= DOM_NUM:
		joined_intervals.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tName={}\n".format(seq_id, configuration.SOURCE_PROFREP, configuration.REPEATS_FEATURE, start_ini, end_ini, configuration.GFF_EMPTY, configuration.GFF_EMPTY, configuration.GFF_EMPTY, ann_ini))
	## if criteria are not met, refragment the expanded interval on original regions
	else:
		for part_start, part_end in zip(starts_part, ends_part):
			joined_intervals.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tName={}\n".format(seq_id, configuration.SOURCE_PROFREP, configuration.REPEATS_FEATURE, part_start, part_end, configuration.GFF_EMPTY, configuration.GFF_EMPTY, configuration.GFF_EMPTY, ann_ini))
	## delete already checked domains from the dict
	del(dict_domains[seq_id][0:index_dom])
	return dict_domains
	
	
def get_unique_classes(CLASS_TBL):
	''' Get all the lineages of current domains database classification table to subsequently check the protein domains if they are classified up to the last level.
	Only these domains will be considered as valid for interval joining. 
	If their classification is be finite (based on comparing to this list of unique classes) they will not be counted for minimum number of domains criterion within the segment to be joined
	'''
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
	
	# first step of refining - removing low confidence repeats regions
	gff_removed = cluster_regions_for_quality_check(REPEATS_GFF, BORDERS)
	
	# second step of refining - interconnecting repeats regions
	if INCLUDE_DOM:
		unique_classes = get_unique_classes(CLASS_TBL)
		dict_domains = create_dom_dict(DOM_GFF)
		joined_intervals = interconnect_intervals(OUT_REFINED, gff_removed, GAP_TH, DOM_NUM, dict_domains, unique_classes)
	else:
		joined_intervals = interconnect_intervals(OUT_REFINED, gff_removed, GAP_TH, DOM_NUM, None, None)
	os.unlink(gff_removed)

		
if __name__ == "__main__":

    # Command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-rep_gff', '--repeats_gff', type=str, required=True,
						help='original repeats regions GFF from PROFREP')
    parser.add_argument('-dom_gff', '--domains_gff', type=str,
						help='protein domains GFF if you want to support repeats joining by domains information')
    parser.add_argument('-gth', '--gap_threshold', type=int, default=250,
						help='gap tolerance between consecutive repeats to be interconnected')
    parser.add_argument('-our', '--out_refined', type=str, default="output_refined.gff",
						help='query sequence to be processed')
    parser.add_argument('-dn', '--dom_number', type=int, default=2,
                        help='min number of domains present to confirm the region joining')
    parser.add_argument('-id', '--include_dom', type=str2bool, default=False,
						help='include domains information to refine the repeats regions')
    parser.add_argument('-ct', '--class_tbl',
						help='classification table of protein domain database to check the level of classification')
    parser.add_argument('-br', '--borders', type=int, default=10,
						help='number of bp tolerated from one or the other side of two overlaping regions when evaluating quality')
    args = parser.parse_args()
    main(args)
