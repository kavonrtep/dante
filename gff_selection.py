#!/usr/bin/env python3

import argparse


def main(args):
	# Command line arguments
	GFF_IN = args.gff_input
	GFF_OUT = args.gff_output
	REGION = args.region
	NEW_SEQ_ID = args.new_seq_id
	
	print(REGION)
	seq_to_cut = REGION.split(":")[0]
	int_start = int(REGION.split(":")[1].split("-")[0])
	int_end = int(REGION.split(":")[1].split("-")[1])
	
	if GFF_OUT is None:
		GFF_OUT = "{}_cut{}:{}.gff3".format(GFF_IN, int_start, int_end)
		
	if not NEW_SEQ_ID:
		NEW_SEQ_ID = "{}_cut{}:{}".format(seq_to_cut, int_start, int_end)
	
	with open(GFF_OUT,"w") as gff_out:
		with open(GFF_IN, "r") as gff_in:
			gff_out.write("##gff-version 3\n")
			gff_out.write("##sequence region {}\n".format(REGION))
			for line in gff_in:
				if not line.startswith("#") and line.split("\t")[0] == seq_to_cut and int(line.split("\t")[3]) >= int_start and int(line.split("\t")[4]) <= int_end:
					new_start = int(line.split("\t")[3]) - int_start + 1
					new_end = int(line.split("\t")[4]) - int_start + 1
					gff_out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(NEW_SEQ_ID, line.split("\t")[1], line.split("\t")[2], new_start, new_end, line.split("\t")[5], line.split("\t")[6],line.split("\t")[7],line.split("\t")[8]))
		

		
if __name__ == "__main__":

    # Command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-gi', '--gff_input', type=str, required=True,
						help='choose gff file')
    parser.add_argument('-go', '--gff_output', type=str,
						help='choose gff file')
    parser.add_argument('-si', '--new_seq_id', type=str,
						help=' ')
    parser.add_argument('-rg', '--region', type=str, required=True,
						help=' ')
    args = parser.parse_args()
    main(args)
