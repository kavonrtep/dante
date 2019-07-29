#!/usr/bin/env python3

import argparse
import time
import os
import textwrap
from collections import defaultdict
from Bio import SeqIO
import configuration

t_nt_seqs_extraction = time.time()


def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected')


def check_file_start(gff_file):
    count_comment = 0
    with open(gff_file, "r") as gff_all:
        line = gff_all.readline()
        while line.startswith("#"):
            line = gff_all.readline()
            count_comment += 1
    return count_comment, line


def extract_nt_seqs(DNA_SEQ, DOM_GFF, OUT_DIR, CLASS_TBL, EXTENDED):
    ''' Extract nucleotide sequences of protein domains found by DANTE from input DNA seq.
		Sequences are saved in fasta files separately for each transposon lineage.
		Sequences extraction is based on position of Best_Hit alignment reported by LASTAL.
		The positions can be extended (optional) based on what part of database domain was aligned (Best_Hit_DB_Pos attribute).
		The strand orientation needs to be considered in extending and extracting the sequence itself
	'''
    [count_comment, first_line] = check_file_start(DOM_GFF)
    unique_classes = get_unique_classes(CLASS_TBL)
    files_dict = defaultdict(str)
    domains_counts_dict = defaultdict(int)
    allSeqs = SeqIO.to_dict(SeqIO.parse(DNA_SEQ, 'fasta'))
    with open(DOM_GFF, "r") as domains:
        for comment_idx in range(count_comment):
            next(domains)
        seq_id_stored = first_line.split("\t")[0]
        allSeqs = SeqIO.to_dict(SeqIO.parse(DNA_SEQ, 'fasta'))
        seq_nt = allSeqs[seq_id_stored]
        for line in domains:
            seq_id = line.split("\t")[0]
            dom_type = line.split("\t")[8].split(";")[0].split("=")[1]
            elem_type = line.split("\t")[8].split(";")[1].split("=")[1]
            strand = line.split("\t")[6]
            align_nt_start = int(line.split("\t")[8].split(";")[3].split(":")[
                -1].split("-")[0])
            align_nt_end = int(line.split("\t")[8].split(";")[3].split(":")[
                -1].split("-")[1].split("[")[0])
            if seq_id != seq_id_stored:
                seq_id_stored = seq_id
                seq_nt = allSeqs[seq_id_stored]
            if EXTENDED:
                ## which part of database sequence was aligned
                db_part = line.split("\t")[8].split(";")[4].split("=")[1]
                ## datatabse seq length
                dom_len = int(db_part.split("of")[1])
                ## start of alignment on database seq
                db_start = int(db_part.split("of")[0].split(":")[0])
                ## end of alignment on database seq
                db_end = int(db_part.split("of")[0].split(":")[1])
                ## number of nucleotides missing in the beginning
                dom_nt_prefix = (db_start - 1) * 3
                ## number of nucleotides missing in the end
                dom_nt_suffix = (dom_len - db_end) * 3
                if strand == "+":
                    dom_nt_start = align_nt_start - dom_nt_prefix
                    dom_nt_end = align_nt_end + dom_nt_suffix
                ## reverse extending for - strand
                else:
                    dom_nt_start = align_nt_start - dom_nt_suffix
                    dom_nt_end = align_nt_end + dom_nt_prefix
                ## correction for domain after extending having negative starting positon
                dom_nt_start = max(1, dom_nt_start)
            else:
                dom_nt_start = align_nt_start
                dom_nt_end = align_nt_end
            full_dom_nt = seq_nt.seq[dom_nt_start - 1:dom_nt_end]
            ## for - strand take reverse complement of the extracted sequence
            if strand == "-":
                full_dom_nt = full_dom_nt.reverse_complement()
            full_dom_nt = str(full_dom_nt)
            ## report when domain classified to the last level and no Ns in extracted seq
            if elem_type in unique_classes and "N" not in full_dom_nt:
                # lineages reported in separate fasta files
                if not elem_type in files_dict:
                    files_dict[elem_type] = os.path.join(
                        OUT_DIR, "{}.fasta".format(elem_type.split("|")[
                            -1].replace("/", "_")))
                with open(files_dict[elem_type], "a") as out_nt_seq:
                    out_nt_seq.write(">{}:{}-{}|{}[{}]\n{}\n".format(
                        seq_nt.id, dom_nt_start, dom_nt_end, dom_type,
                        elem_type, textwrap.fill(full_dom_nt,
                                                 configuration.FASTA_LINE)))
                domains_counts_dict[elem_type] += 1
    return domains_counts_dict


def get_unique_classes(CLASS_TBL):
    ''' Get all the lineages of current domains classification table to check if domains are classified to the last level.
		Only the sequences of unambiguous and completely classified domains will be extracted.
	'''
    unique_classes = []
    with open(CLASS_TBL, "r") as class_tbl:
        for line in class_tbl:
            line_class = "|".join(line.rstrip().split("\t")[1:])
            if line_class not in unique_classes:
                unique_classes.append(line_class)
    return unique_classes


def write_domains_stat(domains_counts_dict, OUT_DIR):
    ''' Report counts of domains for individual lineages
	'''
    total = 0
    with open(
            os.path.join(OUT_DIR,
                         configuration.EXTRACT_DOM_STAT), "w") as dom_stat:
        for domain, count in domains_counts_dict.items():
            dom_stat.write(";{}:{}\n".format(domain, count))
            total += count
        dom_stat.write(";TOTAL:{}\n".format(total))


def main(args):

    DNA_SEQ = args.input_dna
    DOM_GFF = args.domains_gff
    OUT_DIR = args.out_dir
    CLASS_TBL = args.classification
    EXTENDED = args.extended

    if not os.path.exists(OUT_DIR):
        os.makedirs(OUT_DIR)

    domains_counts_dict = extract_nt_seqs(DNA_SEQ, DOM_GFF, OUT_DIR, CLASS_TBL,
                                          EXTENDED)
    write_domains_stat(domains_counts_dict, OUT_DIR)

    print("ELAPSED_TIME_EXTRACTION = {} s\n".format(time.time() -
                                                    t_nt_seqs_extraction))


if __name__ == "__main__":

    # Command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',
                        '--input_dna',
                        type=str,
                        required=True,
                        help='path to input DNA sequence')
    parser.add_argument('-d',
                        '--domains_gff',
                        type=str,
                        required=True,
                        help='GFF file of protein domains')
    parser.add_argument('-cs',
                        '--classification',
                        type=str,
                        required=True,
                        help='protein domains classification file')
    parser.add_argument('-out',
                        '--out_dir',
                        type=str,
                        default=configuration.EXTRACT_OUT_DIR,
                        help='output directory')
    parser.add_argument(
        '-ex',
        '--extended',
        type=str2bool,
        default=True,
        help=
        'extend the domains edges if not the whole datatabase sequence was aligned')
    args = parser.parse_args()
    main(args)
