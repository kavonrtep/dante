#!/usr/bin/env python3

import argparse


def main(args):
	# Command line arguments
	GFF_IN = args.gff_input
	GFF_OUT = args.gff_output
	INT_START = args.start
	INT_END = args.end
	NEW_SEQ_ID = args.new_seq_id
	SEQ_TO_CUT = args.seq_to_cut
	
	if GFF_OUT is None:
		GFF_OUT = "{}_cut{}:{}.gff3".format(GFF_IN, INT_START, INT_END)
		
	if not NEW_SEQ_ID:
		NEW_SEQ_ID = "{}_cut{}:{}".format(SEQ_TO_CUT, INT_START, INT_END)
	
	with open(GFF_OUT,"w") as gff_out:
		with open(GFF_IN, "r") as gff_in:
			gff_out.write(gff_in.readline())
			gff_out.write("##sequence region {} {} {}\n".format(SEQ_TO_CUT, INT_START, INT_END))
			for line in gff_in:
				if line.split("\t")[0] == SEQ_TO_CUT and int(line.split("\t")[3]) >= INT_START and int(line.split("\t")[4]) <= INT_END:
					new_start = int(line.split("\t")[3]) - INT_START + 1
					new_end = int(line.split("\t")[4]) - INT_START + 1
					gff_out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(NEW_SEQ_ID, line.split("\t")[1], line.split("\t")[2], new_start, new_end, line.split("\t")[5], line.split("\t")[6],line.split("\t")[7],line.split("\t")[8]))
		

		
if __name__ == "__main__":

    # Command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-gi', '--gff_input', type=str, required=True,
						help='choose gff file')
    parser.add_argument('-go', '--gff_output', type=str,
						help='choose gff file')
    parser.add_argument('-s', '--start', type=int, required=True,
						help='interval start')
    parser.add_argument('-e', '--end', type=int, required=True,
						help='interval end')
    parser.add_argument('-si', '--new_seq_id', type=str,
						help=' ')
    parser.add_argument('-sc', '--seq_to_cut', type=str, required=True,
						help=' ')
    args = parser.parse_args()
    main(args)
