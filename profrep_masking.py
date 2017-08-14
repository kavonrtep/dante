#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.Alphabet import generic_dna
import sys

def main(args):
	# Command line arguments
	QUERY = args.query
	MODE = args.mode
	REP_TBL = args.rep_table
	MASKED = args.output_masked
	R_TH = args.rp_threshold
		
	repeats_all = get_indices(REP_TBL, R_TH)

	if MODE == "lowercase":
		lower_mask(QUERY, repeats_all, MASKED)
	else:
		N_mask(QUERY, repeats_all, MASKED)

def get_indices(REP_TBL, R_TH):
	repeats_all = {}
	with open(REP_TBL, "r") as repeats_tbl:
		for line in repeats_tbl:
			index = line.split("\t")[0]
			value = line.split
			if line[0].isalpha():
				key = index
				repeats_all[key]=[]
			else:
				value = int(line.split("\t")[1])
				if value >= int(R_TH):
					repeats_all[key].append(int(index))
	return repeats_all
			

def lower_mask(QUERY, repeats_all, MASKED):
	allSeqs = list(SeqIO.parse(QUERY,'fasta'))
	for singleSeq in allSeqs:
		mutable = MutableSeq(str(singleSeq.seq),  generic_dna)
		for index in repeats_all[singleSeq.id]:
			mutable[index - 1] = mutable[index - 1].lower()
		singleSeq.seq = mutable
	with open(MASKED, "w") as handle:
		SeqIO.write(allSeqs, handle, 'fasta')

def N_mask(QUERY, repeats_all, MASKED):
	allSeqs = list(SeqIO.parse(QUERY,'fasta'))
	for singleSeq in allSeqs:
		mutable = MutableSeq(str(singleSeq.seq),  generic_dna)
		for index in repeats_all[singleSeq.id]:
			mutable[index - 1] = "N"
		singleSeq.seq = mutable
	with open(MASKED, "w") as handle:
		SeqIO.write(allSeqs, handle, 'fasta')


if __name__ == "__main__":

    # Command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-q', '--query', type=str, required=True,
						help='query sequence to be processed')
    parser.add_argument('-rt', '--rep_table', type=str, required=True,
						help='query sequence to be processed')
    parser.add_argument('-rth', '--rp_threshold', type=str, default=0,
						help='query sequence to be processed')
    parser.add_argument('-m', '--mode', default="lowercase", choices=['lowercase', 'N'],
						help='query sequence to be processed')
    parser.add_argument('-o', '--output_masked', type=str, default="output_masked",
						help='query sequence to be processed')
	
    args = parser.parse_args()
    main(args)
