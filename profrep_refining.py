#!/usr/bin/env python3

import argparse
import numpy as np
import os
from tempfile import NamedTemporaryFile
from collections import defaultdict
import configuration


def create_dom_dict(DOM_GFF):
	dict_domains = defaultdict(list)
	with open(DOM_GFF, "r") as domains:
		next(domains)
		for line in domains:
			seq_id = line.split("\t")[0]
			ann_dom_lineage = line.split("\t")[8].split(";")[1].split("=")[-1].split("|")[-1]
			start_dom = int(line.split("\t")[3])      
			end_dom = int(line.split("\t")[4])
			strand_dom = line.split("\t")[6]
			dict_domains[seq_id].append((start_dom, end_dom, ann_dom_lineage, strand_dom))
	return dict_domains
	
		
def refining_intervals(OUT_REFINED, REPEATS_GFF, GAP_TH, DOM_NUM, INCLUDE_DOM, dict_domains):
	with open(OUT_REFINED, "w") as joined_intervals:
		joined_intervals.write("{}\n".format(configuration.HEADER_GFF))
		with open(REPEATS_GFF) as repeats:
			next(repeats)
			joined = False
			inicial = repeats.readline()
			if inicial is not "":
				seq_id_ini = inicial.rstrip().split("\t")[0] 
				start_ini = int(inicial.rstrip().split("\t")[3])
				end_ini = int(inicial.rstrip().split("\t")[4]) 
				ann_ini = inicial.rstrip().split("\t")[8].split("=")[-1].split("/")[-1]
				starts_part = [start_ini]
				ends_part = [end_ini]
				for line in repeats:
					seq_id = line.rstrip().split("\t")[0]
					start = int(line.rstrip().split("\t")[3])
					end = int(line.rstrip().split("\t")[4])
					ann = line.rstrip().split("\t")[8].split("=")[-1].split("/")[-1]
					if seq_id != seq_id_ini:
						if INCLUDE_DOM is True and joined is True and ann_ini not in configuration.WITHOUT_DOM:
							dict_domains = dom_refining(joined_intervals, dict_domains, seq_id_ini, start_ini, end_ini, ann_ini, starts_part, ends_part, DOM_NUM)
						else:
							joined_intervals.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tName={}\n".format(seq_id_ini, configuration.SOURCE, configuration.REPEATS_FEATURE, start_ini, end_ini, configuration.R_SCORE, configuration.R_STRAND, configuration.R_PHASE, ann_ini))
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
						if INCLUDE_DOM is True and joined is True and ann_ini not in configuration.WITHOUT_DOM:
							dict_domains = dom_refining(joined_intervals, dict_domains, seq_id, start_ini, end_ini, ann_ini, starts_part, ends_part, DOM_NUM)
						else:
							joined_intervals.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tName={}\n".format(seq_id, configuration.SOURCE, configuration.REPEATS_FEATURE, start_ini, end_ini, configuration.R_SCORE, configuration.R_STRAND, configuration.R_PHASE, ann_ini))
						starts_part = [start] 
						ends_part = [end]
						joined = False
						start_ini = start
						end_ini = end
						ann_ini = ann	
				if INCLUDE_DOM is True and joined is True and ann_ini not in configuration.WITHOUT_DOM:
					dict_domains = dom_refining(joined_intervals, dict_domains, seq_id, start_ini, end_ini, ann_ini, starts_part, ends_part, DOM_NUM)
				else:
					joined_intervals.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tName={}\n".format(seq_id, configuration.SOURCE, configuration.REPEATS_FEATURE, start_ini, end_ini, configuration.R_SCORE, configuration.R_STRAND, configuration.R_PHASE, ann_ini))
				del(starts_part)
				del(ends_part)


def dom_refining(joined_intervals, dict_domains, seq_id, start_ini, end_ini, ann_ini, starts_part, ends_part, DOM_NUM):
	count_dom = 0
	index_dom = 0
	strands = []
	indices_del = []
	for dom_attributes in dict_domains[seq_id]:
		ann_dom = dom_attributes[2]
		if dom_attributes[0] >= start_ini and dom_attributes[1] <= end_ini:
			if ann_dom == ann_ini:
				strands.append(dom_attributes[3])
				count_dom += 1
			index_dom += 1 
	if len(set(strands)) <= 1 and count_dom >= DOM_NUM:
		joined_intervals.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tName={}\n".format(seq_id, configuration.SOURCE, configuration.REPEATS_FEATURE, start_ini, end_ini, configuration.R_SCORE, configuration.R_STRAND, configuration.R_PHASE, ann_ini))
	else:
		for part_start, part_end in zip(starts_part, ends_part):
			joined_intervals.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tName={}\n".format(seq_id, configuration.SOURCE, configuration.REPEATS_FEATURE, part_start, part_end, configuration.R_SCORE, configuration.R_STRAND, configuration.R_PHASE, ann_ini))
	for idx in range(index_dom):
		del(dict_domains[seq_id][idx])
	return dict_domains
	
	
def main(args):
	# Command line arguments
	REPEATS_GFF = args.repeats_gff
	DOM_GFF = args.domains_gff
	GAP_TH = args.gap_threshold
	DOM_NUM = args.dom_number
	OUT_REFINED = args.out_refined 
	INCLUDE_DOM = args.include_dom
	
	if INCLUDE_DOM is True:
		dict_domains = create_dom_dict(DOM_GFF)
		joined_intervals = refining_intervals(OUT_REFINED, REPEATS_GFF, GAP_TH, DOM_NUM, INCLUDE_DOM, dict_domains)
	else:
		joined_intervals = refining_intervals(OUT_REFINED, REPEATS_GFF, GAP_TH, DOM_NUM, INCLUDE_DOM, None)

		
if __name__ == "__main__":

    # Command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-rep_gff', '--repeats_gff', type=str, required=True,
						help='query sequence to be processed')
    parser.add_argument('-dom_gff', '--domains_gff', type=str,
						help='query sequence to be processed')
    parser.add_argument('-gth', '--gap_threshold', type=int, default=500,
						help='query sequence to be processed')
    parser.add_argument('-or', '--out_refined', type=str, default="output_refined.gff",
						help='query sequence to be processed')
    parser.add_argument('-dn', '--dom_number', type=int, default=1,
                        help='number of domains present to confirm one element type')
    parser.add_argument('-id', '--include_dom', action='store_true', default=False,
						help='Include domains information to refine repeats regions output')
    args = parser.parse_args()
    main(args)
