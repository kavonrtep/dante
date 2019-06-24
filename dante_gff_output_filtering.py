#!/usr/bin/env python3

import time
import configuration
import os
import textwrap
import subprocess
from tempfile import NamedTemporaryFile
from collections import defaultdict


class Range():
    '''
    This class is used to check float range in argparse
    '''

    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __eq__(self, other):
        return self.start <= other <= self.end

    def __str__(self):
        return "float range {}..{}".format(self.start, self.end)

    def __repr__(self):
        return "float range {}..{}".format(self.start, self.end)


def check_file_start(gff_file):
    count_comment = 0
    with open(gff_file, "r") as gff_all:
        line = gff_all.readline()
        while line.startswith("#"):
            line = gff_all.readline()
            count_comment += 1
    return count_comment


def write_info(filt_dom_tmp, FILT_DOM_GFF, orig_class_dict, filt_class_dict,
               dom_dict, version_lines):
    '''
	Write domains statistics in beginning of filtered GFF
	'''
    with open(FILT_DOM_GFF, "w") as filt_gff:
        for line in version_lines:
            filt_gff.write(line)
        filt_gff.write("##CLASSIFICATION\tORIGINAL_COUNTS\tFILTERED_COUNTS\n")
        if not orig_class_dict:
            filt_gff.write("##NO DOMAINS CLASSIFICATIONS\n")
        for classification in sorted(orig_class_dict.keys()):
            if classification in filt_class_dict.keys():
                filt_gff.write("##{}\t{}\t{}\n".format(
                    classification, orig_class_dict[
                        classification], filt_class_dict[classification]))
            else:
                filt_gff.write("##{}\t{}\t{}\n".format(
                    classification, orig_class_dict[classification], 0))
        filt_gff.write("##-----------------------------------------------\n"
                       "##SEQ\tDOMAIN\tCOUNTS\n")
        if not dom_dict:
            filt_gff.write("##NO DOMAINS\n")
        for seq in sorted(dom_dict.keys()):
            for dom, count in sorted(dom_dict[seq].items()):
                filt_gff.write("##{}\t{}\t{}\n".format(seq, dom, count))
        filt_gff.write("##-----------------------------------------------\n")
        with open(filt_dom_tmp.name, "r") as filt_tmp:
            for line in filt_tmp:
                filt_gff.write(line)


def get_file_start(gff_file):
    count_comment = 0
    lines = []
    with open(gff_file, "r") as gff_all:
        line = gff_all.readline()
        while line.startswith("#"):
            lines.append(line)
            line = gff_all.readline()
            count_comment += 1
    return count_comment, lines


def filter_qual_dom(DOM_GFF, FILT_DOM_GFF, TH_IDENTITY, TH_SIMILARITY,
                    TH_LENGTH, TH_INTERRUPT, TH_LEN_RATIO, SELECTED_DOM,
                    ELEMENT):
    ''' Filter gff output based on domain and quality of alignment '''
    [count_comment, version_lines] = get_file_start(DOM_GFF)
    filt_dom_tmp = NamedTemporaryFile(delete=False)
    with open(DOM_GFF, "r") as gff_all, open(filt_dom_tmp.name,
                                             "w") as gff_filtered:
        for comment_idx in range(count_comment):
            next(gff_all)
        dom_dict = defaultdict(lambda: defaultdict(int))
        orig_class_dict = defaultdict(int)
        filt_class_dict = defaultdict(int)
        seq_ids_all = []
        xminimals = []
        xmaximals = []
        domains = []
        xminimals_all = []
        xmaximals_all = []
        domains_all = []
        start = True
        for line in gff_all:
            attributes = line.rstrip().split("\t")[-1]
            classification = attributes.split(";")[1].split("=")[1]
            orig_class_dict[classification] += 1
            ## ambiguous domains filtered out automatically
            if classification != configuration.AMBIGUOUS_TAG:
                al_identity = float(attributes.split(";")[-5].split("=")[1])
                al_similarity = float(attributes.split(";")[-4].split("=")[1])
                al_length = float(attributes.split(";")[-3].split("=")[1])
                relat_interrupt = float(attributes.split(";")[-2].split("=")[
                    1])
                db_len_proportion = float(attributes.split(";")[-1].split("=")[
                    1])
                dom_type = attributes.split(";")[0].split("=")[1]
                seq_id = line.split("\t")[0]
                xminimal = int(line.split("\t")[3])
                xmaximal = int(line.split("\t")[4])
                if al_identity >= TH_IDENTITY and al_similarity >= TH_SIMILARITY and al_length >= TH_LENGTH and relat_interrupt <= TH_INTERRUPT and db_len_proportion <= TH_LEN_RATIO and (
                        dom_type == SELECTED_DOM or
                        SELECTED_DOM == "All") and (ELEMENT in classification):
                    gff_filtered.writelines(line)
                    filt_class_dict[classification] += 1
                    dom_dict[seq_id][dom_type] += 1
                    if start:
                        seq_ids_all.append(line.split("\t")[0])
                        start = False
                    if seq_id != seq_ids_all[-1]:
                        seq_ids_all.append(seq_id)
                        xminimals_all.append(xminimals)
                        xmaximals_all.append(xmaximals)
                        domains_all.append(domains)
                        xminimals = []
                        xmaximals = []
                        domains = []
                    xminimals.append(xminimal)
                    xmaximals.append(xmaximal)
                    domains.append(dom_type)
    path = os.path.dirname(os.path.realpath(__file__))
    write_info(filt_dom_tmp, FILT_DOM_GFF, orig_class_dict, filt_class_dict,
               dom_dict, version_lines)
    os.unlink(filt_dom_tmp.name)
    xminimals_all.append(xminimals)
    xmaximals_all.append(xmaximals)
    domains_all.append(domains)
    return xminimals_all, xmaximals_all, domains_all, seq_ids_all


def get_domains_protseq(FILT_DOM_GFF, DOMAIN_PROT_SEQ):
    ''' Get the translated protein sequence of original DNA seq for all the filtered domains regions 
		The translated sequences are taken from alignment reported by LASTAL (Query_Seq attribute in GFF)	
	'''
    count_comment = check_file_start(FILT_DOM_GFF)
    with open(FILT_DOM_GFF, "r") as filt_gff:
        for comment_idx in range(count_comment):
            next(filt_gff)
        with open(DOMAIN_PROT_SEQ, "w") as dom_prot_file:
            for line in filt_gff:
                attributes = line.rstrip().split("\t")[8]
                positions = attributes.split(";")[3].split("=")[1].split(":")[
                    -1].split("[")[0]
                dom = attributes.split(";")[0].split("=")[1]
                dom_class = attributes.split(";")[1].split("=")[1]
                seq_id = line.rstrip().split("\t")[0]
                prot_seq_align = line.rstrip().split("\t")[8].split(";")[
                    6].split("=")[1]
                prot_seq = prot_seq_align.translate({ord(i): None
                                                     for i in '/\\-'})
                header_prot_seq = ">{}:{} {} {}".format(seq_id, positions, dom,
                                                        dom_class)
                dom_prot_file.write("{}\n{}\n".format(
                    header_prot_seq, textwrap.fill(prot_seq,
                                                   configuration.FASTA_LINE)))


def main(args):

    t = time.time()

    DOM_GFF = args.dom_gff
    DOMAIN_PROT_SEQ = args.domains_prot_seq
    TH_IDENTITY = args.th_identity
    TH_LENGTH = args.th_length
    TH_INTERRUPT = args.interruptions
    TH_SIMILARITY = args.th_similarity
    TH_LEN_RATIO = args.max_len_proportion
    FILT_DOM_GFF = args.domains_filtered
    SELECTED_DOM = args.selected_dom
    OUTPUT_DIR = args.output_dir
    # DELETE : ELEMENT = args.element_type.replace("_pipe_", "|")
    ELEMENT = args.element_type

    if DOMAIN_PROT_SEQ is None:
        DOMAIN_PROT_SEQ = configuration.DOM_PROT_SEQ
    if FILT_DOM_GFF is None:
        FILT_DOM_GFF = configuration.FILT_DOM_GFF

    if OUTPUT_DIR and not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    if not os.path.isabs(FILT_DOM_GFF):
        if OUTPUT_DIR is None:
            OUTPUT_DIR = os.path.dirname(os.path.abspath(DOM_GFF))
        FILT_DOM_GFF = os.path.join(OUTPUT_DIR, os.path.basename(FILT_DOM_GFF))
        DOMAIN_PROT_SEQ = os.path.join(OUTPUT_DIR,
                                       os.path.basename(DOMAIN_PROT_SEQ))

    [xminimals_all, xmaximals_all, domains_all, seq_ids_all] = filter_qual_dom(
        DOM_GFF, FILT_DOM_GFF, TH_IDENTITY, TH_SIMILARITY, TH_LENGTH,
        TH_INTERRUPT, TH_LEN_RATIO, SELECTED_DOM, ELEMENT)
    get_domains_protseq(FILT_DOM_GFF, DOMAIN_PROT_SEQ)

    print("ELAPSED_TIME_DOMAINS = {} s".format(time.time() - t))


if __name__ == "__main__":
    import argparse
    from argparse import RawDescriptionHelpFormatter

    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                          argparse.RawDescriptionHelpFormatter):
        pass

    parser = argparse.ArgumentParser(
        description=
        '''The script performs DANTE's output filtering for quality and/or extracting specific type of protein domain or mobile elements of origin. For the filtered domains it reports their translated protein sequence of original DNA.
		WHEN NO PARAMETERS GIVEN, IT PERFORMS QUALITY FILTERING USING THE DEFAULT PARAMETRES (optimized for Viridiplantae species)
		
		INPUTS:
			- GFF3 file produced by protein_domains.py OR already filtered GFF3
			
			FILTERING OPTIONS:
				> QUALITY: - Min relative length of alignemnt to the protein domain from DB (without gaps)
				   - Identity 
				   - Similarity (scoring matrix: BLOSUM82)
				   - Interruption in the reading frame (frameshifts + stop codons) per every starting 100 AA
				   - Max alignment proportion to the original length of database domain sequence 
				> DOMAIN TYPE: choose from choices ('Name' attribute in GFF)
				Records for ambiguous domain type (e.g. INT/RH) are filtered out automatically
				
				> MOBILE ELEMENT TYPE:
				arbitrary substring of the element classification ('Final_Classification' attribute in GFF)
				
		OUTPUTS:
			- filtered GFF3 file
			- fasta file of translated protein sequences (from original DNA) for the aligned domains that match the filtering criteria 
		
	DEPENDENCIES:
		- python 3.4 or higher
		> ProfRep modules:
			- configuration.py 

	EXAMPLE OF USAGE:
		Getting quality filtered integrase(INT) domains of all gypsy transposable elements:
		./domains_filtering.py -dom_gff PATH_TO_INPUT_GFF -pdb PATH_TO_PROTEIN_DB -cs PATH_TO_CLASSIFICATION_FILE --selected_dom INT --element_type Ty3/gypsy 

		''',
        epilog="""""",
        formatter_class=CustomFormatter)
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-dg",
                               "--dom_gff",
                               type=str,
                               required=True,
                               help="basic unfiltered gff file of all domains")
    parser.add_argument("-ouf",
                        "--domains_filtered",
                        type=str,
                        help="output filtered domains gff file")
    parser.add_argument("-dps",
                        "--domains_prot_seq",
                        type=str,
                        help="output file containg domains protein sequences")
    parser.add_argument("-thl",
                        "--th_length",
                        type=float,
                        choices=[Range(0.0, 1.0)],
                        default=0.8,
                        help="proportion of alignment length threshold")
    parser.add_argument("-thi",
                        "--th_identity",
                        type=float,
                        choices=[Range(0.0, 1.0)],
                        default=0.35,
                        help="proportion of alignment identity threshold")
    parser.add_argument("-ths",
                        "--th_similarity",
                        type=float,
                        choices=[Range(0.0, 1.0)],
                        default=0.45,
                        help="threshold for alignment proportional similarity")
    parser.add_argument(
        "-ir",
        "--interruptions",
        type=int,
        default=3,
        help=
        "interruptions (frameshifts + stop codons) tolerance threshold per 100 AA")
    parser.add_argument(
        "-mlen",
        "--max_len_proportion",
        type=float,
        default=1.2,
        help=
        "maximal proportion of alignment length to the original length of protein domain from database")
    parser.add_argument(
        "-sd",
        "--selected_dom",
        type=str,
        default="All",
        choices=[
            "All", "GAG", "INT", "PROT", "RH", "RT", "aRH", "CHDCR", "CHDII",
            "TPase", "YR", "HEL1", "HEL2", "ENDO"
        ],
        help="filter output domains based on the domain type")
    parser.add_argument(
        "-el",
        "--element_type",
        type=str,
        default="",
        help="filter output domains by typing substring from classification")
    parser.add_argument(
        "-dir",
        "--output_dir",
        type=str,
        default=None,
        help="specify if you want to change the output directory")
    args = parser.parse_args()
    main(args)
