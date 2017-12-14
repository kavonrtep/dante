#!/usr/bin/env python3

import zipfile 
import tempfile
import argparse
from shutil import copyfile
import os
import configuration

def main(args):
	
	RE_ARCHIVE = args.re_archive
	OUTPUT_CLS = args.output_cls
	OUTPUT_READS_ALL = args.output_reads_all
	OUTPUT_ANNOTATION = args.output_annotation

	if not os.path.isabs(OUTPUT_CLS):
		OUTPUT_CLS = os.path.join(os.getcwd(), OUTPUT_CLS)
		
	if not os.path.isabs(OUTPUT_READS_ALL):
		OUTPUt_READS_ALL = os.path.join(os.getcwd(), OUTPUT_READS_ALL)
		
	if not os.path.isabs(OUTPUT_ANNOTATION):
		OUTPUT_ANNOTATION = os.path.join(os.getcwd(), OUTPUT_ANNOTATION)
		
	with tempfile.TemporaryDirectory() as dirpath:
		with zipfile.ZipFile(RE_ARCHIVE, 'r') as re_archive:
			re_archive.extractall(dirpath)
		copyfile(os.path.join(dirpath, configuration.HITSORT_CLS), OUTPUT_CLS)
		copyfile(os.path.join(dirpath, configuration.READS_ALL), OUTPUT_READS_ALL)
		copyfile(os.path.join(dirpath, configuration.ANNOTATION), OUTPUT_ANNOTATION)


if __name__ == '__main__':
	
	# Command line arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('-ar', '--re_archive', type=str, required=True,
						help='RepeatExplorer output data archive')		
	parser.add_argument('-oc', '--output_cls', type=str, default="output_hitsort_cls",
						help='Output cls file of all clusters')
	parser.add_argument('-or', '--output_reads_all', type=str, default="output_reads_all",
						help='Output file of all reads sequences')	
	parser.add_argument('-oa', '--output_annotation', type=str, default="output_annotation",
						help='Output file of clusters annotation')			
	args = parser.parse_args()
	main(args)
