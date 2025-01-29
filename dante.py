#!/usr/bin/env python3

import numpy as np
import subprocess
import math
import time
from operator import itemgetter
from collections import Counter
from itertools import groupby
import os
import re
import configuration
from tempfile import NamedTemporaryFile
import sys
import warnings
import shutil
from collections import defaultdict

np.set_printoptions(threshold=sys.maxsize)

def alignment_scoring():
    ''' Create hash table for alignment similarity counting: for every 
	combination of aminoacids in alignment assign score from protein 
	scoring matrix defined in configuration file  '''
    score_dict = {}
    with open(configuration.SC_MATRIX) as smatrix:
        count = 1
        for line in smatrix:
            if not line.startswith("#"):
                if count == 1:
                    aa_all = line.rstrip().replace(" ", "")
                else:
                    count_aa = 1
                    line = list(filter(None, line.rstrip().split(" ")))
                    for aa in aa_all:
                        score_dict["{}{}".format(line[0], aa)] = line[count_aa]
                        count_aa += 1
                count += 1
    return score_dict


def characterize_fasta(QUERY, WIN_DOM):
    ''' Find the sequences, their lengths, starts, ends and if 
	they exceed the window '''
    with open(QUERY) as query:
        headers = []
        fasta_lengths = []
        seq_starts = []
        seq_ends = []
        fasta_chunk_len = 0
        count_line = 1
        for line in query:
            line = line.rstrip()
            if line.startswith(">"):
                headers.append(line.rstrip())
                fasta_lengths.append(fasta_chunk_len)
                fasta_chunk_len = 0
                seq_starts.append(count_line + 1)
                seq_ends.append(count_line - 1)
            else:
                fasta_chunk_len += len(line)
            count_line += 1
        seq_ends.append(count_line)
        seq_ends = seq_ends[1:]
        fasta_lengths.append(fasta_chunk_len)
        fasta_lengths = fasta_lengths[1:]
        # control if there are correct (unique) names for individual seqs:
        # LASTAL takes seqs IDs till the first space which can then create problems with ambiguous records
        if len(headers) > len(set([header.split(" ")[0] for header in headers
                                   ])):
            raise NameError(
                '''Sequences in multifasta format are not named correctly:
							seq IDs (before the first space) are the same''')

    above_win = [idx
                 for idx, value in enumerate(fasta_lengths) if value > WIN_DOM]
    below_win = [idx
                 for idx, value in enumerate(fasta_lengths)
                 if value <= WIN_DOM]
    lens_above_win = np.array(fasta_lengths)[above_win]
    return headers, above_win, below_win, lens_above_win, seq_starts, seq_ends


def split_fasta(QUERY, WIN_DOM, step, headers, above_win, below_win,
                lens_above_win, seq_starts, seq_ends):
    ''' Create temporary file containing all sequences - the ones that exceed 
	the window are cut with a set overlap (greater than domain size with a reserve) '''
    with open(QUERY, "r") as query:
        count_fasta_divided = 0
        count_fasta_not_divided = 0
        ntf = NamedTemporaryFile(delete=False)
        divided = np.array(headers)[above_win]
        row_length = configuration.FASTA_LINE
        for line in query:
            line = line.rstrip()
            if line.startswith(">") and line in divided:
                stop_line = seq_ends[above_win[
                    count_fasta_divided]] - seq_starts[above_win[
                        count_fasta_divided]] + 1
                count_line = 0
                whole_seq = []
                for line2 in query:
                    whole_seq.append(line2.rstrip())
                    count_line += 1
                    if count_line == stop_line:
                        break
                whole_seq = "".join(whole_seq)
                ## create list of starting positions for individual parts of a seq with a step given by a window and overlap
                windows_starts = list(range(0, lens_above_win[
                    count_fasta_divided], step))
                ## create list of ending positions (starting pos + window), the last element is the whole seq length
                windows_ends = [
                    x + WIN_DOM
                    if x + WIN_DOM < lens_above_win[count_fasta_divided] else
                    lens_above_win[count_fasta_divided] for x in windows_starts
                ]
                count_part = 1
                for start_part, end_part in zip(windows_starts, windows_ends):
                    seq_part = whole_seq[start_part:end_part]
                    if count_part == len(windows_starts):
                        ntf.write("{}_DANTE_PART{}_LAST:{}-{}\n{}\n".format(
                            line.split(" ")[0], count_part, start_part + 1,
                            end_part, "\n".join([seq_part[i:i + row_length]
                                                 for i in range(0, len(
                                                     seq_part), row_length)
                                                 ])).encode("utf-8"))
                    else:
                        ntf.write("{}_DANTE_PART{}:{}-{}\n{}\n".format(
                            line.split(" ")[0], count_part, start_part + 1,
                            end_part, "\n".join([seq_part[i:i + row_length]
                                                 for i in range(0, len(
                                                     seq_part), row_length)
                                                 ])).encode("utf-8"))
                    count_part += 1
                count_fasta_divided += 1
            elif line.startswith(">") and line not in divided:
                length_seq = seq_ends[below_win[
                    count_fasta_not_divided]] - seq_starts[below_win[
                        count_fasta_not_divided]] + 1
                ntf.write("{}\n{}".format(line, "".join([query.readline(
                ) for x in range(length_seq)])).encode("utf-8"))
                count_fasta_not_divided += 1
        query_temp = ntf.name
        ntf.close()
    return query_temp


def domain_annotation(elements, CLASSIFICATION):
    ''' Assign protein domain to each hit from protein database  '''
    domains = []
    annotations = []
    with open(CLASSIFICATION, "r") as cl_tbl:
        annotation = {}
        for line in cl_tbl:
            record = line.rstrip().split("\t")
            annotation[record[0]] = record[1:]
    for i in range(len(elements)):
        domains.append(elements[i].split("__")[0].split("-")[1])
        element_name = "__".join(elements[i].split("__")[1:])
        if element_name in annotation.keys():
            annotations.append("|".join([elements[i].split("__")[0].split("-")[
                1], ("|".join(annotation[element_name]))]))
        else:
            annotations.append("unknown|unknown")
    return annotations


def hits_processing(seq_len, start, end, strand):
    ''' Gain hits intervals separately for forward and reverse strand '''
    reverse_strand_idx = np.where(strand == "-")[0]
    if not reverse_strand_idx.any():
        start_pos_plus = start + 1
        end_pos_plus = end
        regions_plus = list(zip(start_pos_plus, end_pos_plus))
        regions_minus = []
    else:
        reverse_strand_idx = reverse_strand_idx[0]
        start_pos_plus = start[0:reverse_strand_idx] + 1
        end_pos_plus = end[0:reverse_strand_idx]
        start_pos_minus = seq_len[0] - end[reverse_strand_idx:] + 1
        end_pos_minus = seq_len[0] - start[reverse_strand_idx:]
        regions_plus = list(zip(start_pos_plus, end_pos_plus))
        regions_minus = list(zip(start_pos_minus, end_pos_minus))
    return reverse_strand_idx, regions_plus, regions_minus


def overlapping_regions(input_data):
    ''' Join all overlapping intervals(hits) to clusters (potential domains),
	get list of start-end positions of individual hits within the interval, 
	list of minimus and maximums as well as the indices in the original 
	sequence_hits structure for the hits belonging to the same clusters '''
    if input_data:
        sorted_idx, sorted_data = zip(*sorted(
            [(index, data) for index, data in enumerate(input_data)],
            key=itemgetter(1)))
        merged_ends = input_data[sorted_idx[0]][1]
        intervals = []
        data = []
        output_intervals = []
        output_data = []
        for i, j in zip(sorted_idx, sorted_data):
            if input_data[i][0] < merged_ends:
                merged_ends = max(input_data[i][1], merged_ends)
                intervals.append(i)
                data.append(j)
            else:
                output_intervals.append(intervals)
                output_data.append(data)
                intervals = []
                data = []
                intervals.append(i)
                data.append(j)
                merged_ends = input_data[i][1]
        output_intervals.append(intervals)
        output_data.append(data)
        mins = [x[0][0] for x in output_data]
        maxs = [max(x, key=itemgetter(1))[1] for x in output_data]
    else:
        mins = []
        maxs = []
        output_intervals = []
        output_data = []
    return mins, maxs, output_data, output_intervals


def annotations_dict(annotations):
    ''' Hash table where annotations of the hits within a clusters are the keys. 
	Each annotation has serial number assigned which indexes the row in the score_table '''
    classes_dict = {classes: idx
                    for idx, classes in enumerate(set(annotations))}
    return classes_dict


def score_table(mins, maxs, data, annotations, scores, CLASSIFICATION):
    ''' Score table is created based on the annotations occurance in the cluster.
	Matrix axis y corresponds to individual annotations (indexed according to classes_dict),
    axis x represents positions of analyzed seq in a given cluster.
    For every hit within cluster, array of scores on the corresponding position
    is recorded to the table in case if the score on certain position is so far the highest
	for the certain position and certain annotation '''
    classes_dict = annotations_dict(annotations)
    score_matrix = np.zeros((len(classes_dict), maxs - mins + 1), dtype=int)
    count = 0
    for item in annotations:
        saved_scores = score_matrix[classes_dict[item], data[count][0] - mins:
                                    data[count][1] - mins + 1]
        new_scores = [scores[count]] * len(saved_scores)
        score_matrix[classes_dict[item], data[count][0] - mins:data[count][
            1] - mins + 1] = [max(*pos_score)
                              for pos_score in zip(saved_scores, new_scores)]
        count += 1
    return score_matrix, classes_dict


def score_matrix_evaluation(score_matrix, classes_dict, THRESHOLD_SCORE):
    ''' Score matrix is evaluated based on each position.
	For every position the list of annotations with a score which reaches 
	certain percentage of the overal best score of the cluster are stored '''
    ann_per_reg = []
    overal_best_score_reg = max((score_matrix.max(axis=1)))
    for position in score_matrix.T:
        ## score threshold calculated as a percentage of the OVERALL best score in the cluster
        threshold = overal_best_score_reg * THRESHOLD_SCORE / 100
        above_th = [idx
                    for idx, score in enumerate(position)
                    if position[idx] >= threshold]
        ## select unique annotations in one position that are above threshold
        ann_per_pos = list(set(
            [key for key, value in classes_dict.items() if value in above_th]))
        ann_per_reg.append(ann_per_pos)
    return ann_per_reg


def group_annot_regs(ann_per_reg):
    ''' Get list of domains, annotations, longest common annotations and 
	counts of positions with certain annotation per regions '''
    ## tranform list of lists (potential multiple annotations for every position ) to flat list of all annotations
    all_annotations = [item for sublist in ann_per_reg for item in sublist]
    unique_annotations = list(set(all_annotations))
    ann_pos_counts = [all_annotations.count(x) for x in unique_annotations]
    unique_annotations = list(set(
        [item for sublist in ann_per_reg for item in sublist]))
    domain_type = list(set([annotation.split("|")[0]
                            for annotation in unique_annotations]))
    classification_list = [annotation.split("|")
                           for annotation in unique_annotations]
    ann_substring = "|".join(os.path.commonprefix(classification_list))
    domain_type = "/".join(domain_type)
    return domain_type, ann_substring, unique_annotations, ann_pos_counts


def best_score(scores, region):
    ''' From overlapping intervals take the one with the highest score '''
    ## if more hits have the same best score take only the first one
    best_idx = region[np.where(scores == max(scores))[0][0]]
    best_idx_reg = np.where(scores == max(scores))[0][0]
    return best_idx, best_idx_reg


def create_gff3(domain_type, ann_substring, unique_annotations, ann_pos_counts,
                dom_start, dom_end, step, best_idx, annotation_best,
                db_name_best, db_starts_best, db_ends_best, strand, score,
                seq_id, db_seq, query_seq, domain_size, positions, gff, consensus):
    ''' Record obtained information about domain corresponding to individual cluster to common gff file '''
    best_start = positions[best_idx][0]
    best_end = positions[best_idx][1]
    best_score = score[best_idx]
    ## proportion of length of the best hit to the whole region length found by base
    length_proportion = int((best_end - best_start + 1) /
                            (dom_end - dom_start + 1) * 100)
    db_seq_best = db_seq[best_idx]
    query_seq_best = query_seq[best_idx]
    domain_size_best = domain_size[best_idx]
    [percent_ident, align_similarity, relat_align_len, relat_interrupt,
     db_len_proportion
     ] = filter_params(db_seq_best, query_seq_best, domain_size_best)
    ann_substring = "|".join(ann_substring.split("|")[1:])
    annotation_best = "|".join([db_name_best] + annotation_best.split("|")[1:])
    if "DANTE_PART" in seq_id:
        part = int(seq_id.split("DANTE_PART")[1].split(":")[0].split("_")[0])
        dom_start = dom_start + (part - 1) * step
        dom_end = dom_end + (part - 1) * step
        best_start = best_start + (part - 1) * step
        best_end = best_end + (part - 1) * step
    if ann_substring == '':
        ann_substring = "NONE(Annotations from different classes)"
    if len(unique_annotations) > 1:
        unique_annotations = ",".join(["{}[{}bp]".format(
            ann, pos) for ann, pos in zip(unique_annotations, ann_pos_counts)])
    else:
        unique_annotations = unique_annotations[0]
    SOURCE = configuration.SOURCE_DANTE
    if "/" in domain_type:
        gff.write(
            "{}\t{}\t{}\t{}\t{}\t.\t{}\t{}\tName={};Final_Classification=Ambiguous_domain;Region_Hits_Classifications_={}\n".format(
                seq_id, SOURCE, configuration.DOMAINS_FEATURE, dom_start,
                dom_end, strand, configuration.PHASE, domain_type,
                unique_annotations))
    else:
        gff.write(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tName={};Final_Classification={};Region_Hits_Classifications={};Best_Hit={}:{}-{}[{}percent];Best_Hit_DB_Pos={}:{}of{};DB_Seq={};Region_Seq={};Query_Seq={};Identity={};Similarity={};Relat_Length={};Relat_Interruptions={};Hit_to_DB_Length={}\n".format(
                seq_id, SOURCE, configuration.DOMAINS_FEATURE, dom_start,
                dom_end, best_score, strand, configuration.PHASE, domain_type,
                ann_substring, unique_annotations, annotation_best, best_start,
                best_end, length_proportion, db_starts_best, db_ends_best,
                domain_size_best, db_seq_best, consensus, query_seq_best, percent_ident,
                align_similarity, relat_align_len, relat_interrupt,
                db_len_proportion))


def filter_params(db, query, protein_len):
    ''' Calculate basic statistics of the quality of the alignment '''
    score_dict = alignment_scoring()
    num_ident = 0
    count_interrupt = 0
    count_similarity = 0
    alignment_len = 0
    for i, j in zip(db.upper(), query.upper()):
        if i == j and i != "X":
            num_ident += 1
        if j == "/" or j == "\\" or j == "*":
            count_interrupt += 1
        if (i.isalpha() or i == "*") and (j.isalpha() or j == "*"):
            if int(score_dict["{}{}".format(i, j)]) > 0:
                count_similarity += 1
    ## gapless alignment length proportional to the domain protein length
    relat_align_len = round((len(db) - db.count("-")) / protein_len, 3)
    ## proportional identical bases (except of X) to al.length
    align_identity = round(num_ident / len(db), 2)
    ## proportional count of positive scores from scoring matrix to al. length
    align_similarity = round(count_similarity / len(db), 2)
    ## number of interruptions per 100 bp
    relat_interrupt = round(count_interrupt / math.ceil((len(query) / 100)), 2)
    ## Proportion of alignment to the original length of protein domain from database (indels included)
    db_len_proportion = round(len(db) / protein_len, 2)
    return align_identity, align_similarity, relat_align_len, relat_interrupt, db_len_proportion


def line_generator(tab_pipe, maf_pipe, start):
    ''' Yield individual lines of LASTAL stdout for single sequence '''
    if hasattr(line_generator, "dom"):
        seq_id = line_generator.dom.split("\t")[6]
        yield line_generator.dom.encode("utf-8")
        del line_generator.dom
    line_tab = ""
    for line_tab in tab_pipe:
        line_tab = line_tab.decode("utf-8")
        if not line_tab.startswith('#'):
            # is some versions of lastal, the record is split into multiple lines
            # first line contain only two columns, second line contain the rest
            # of the record. It must be merged into one line
            if len(line_tab.split("\t")) == 2:
                broken_line = True
                line_tab += tab_pipe.readline().decode("utf-8")
            else:
                broken_line = False
            if start:
                if not ('seq_id' in locals() and
                        seq_id != line_tab.split("\t")[6]):
                    seq_id = line_tab.split("\t")[6]
                    start = False
            Nlines = 5 if broken_line else 4
            line_maf = [maf_pipe.readline() for line_count in range(Nlines)]
            # if the record is broken, concatenate second and third line
            if broken_line:
                line_maf = [line_maf[0], line_maf[1] + line_maf[2], line_maf[3], line_maf[4]]
            db_seq = line_maf[1].decode("utf-8").rstrip().split(" ")[-1]
            alignment_seq = line_maf[2].decode("utf-8").rstrip().split(" ")[-1]
            line = "{}\t{}\t{}".format(line_tab, db_seq, alignment_seq)
            line_id = line.split("\t")[6]
            if seq_id != line_id:
                line_generator.dom = line
                return
            else:
                yield line.encode("utf-8")
        else:
            maf_pipe.readline()
    if line_tab == "":
        raise RuntimeError
    else:
        return


def get_version(path, LAST_DB):
    '''Return version is run from git repository '''
    version_string = (
        "##-----------------------------------------------\n"
        "##PROTEIN DATABASE VERSION : {PD}\n"
        "##-----------------------------------------------\n").format(
            PD=os.path.basename(LAST_DB)
        )
    if os.path.exists(".git"):
        try:
            branch = subprocess.check_output("git rev-parse --abbrev-ref HEAD",
                                             shell=True,
                                             cwd=path,
                                             stderr=subprocess.DEVNULL).decode(
                'ascii').strip()
            shorthash = subprocess.check_output("git log --pretty=format:'%h' -n 1  ",
                                                shell=True,
                                                cwd=path).decode('ascii').strip()
            revcount = len(subprocess.check_output("git log --oneline",
                                                   shell=True,
                                                   cwd=path).decode('ascii').split())
            version_string = (
                "##-----------------------------------------------\n"
                "##PIPELINE VERSION         : "
                "{branch}-rv-{revcount}({shorthash})\n"
                "##PROTEIN DATABASE VERSION : {PD}\n"
                "##-----------------------------------------------\n").format(
                    branch=branch,
                    shorthash=shorthash,
                    revcount=revcount,
                    PD=os.path.basename(LAST_DB))
        except:
            pass
    return version_string


def write_info(dom_gff_tmp, version_string):
    dom_gff_tmp.write("{}\n".format(configuration.HEADER_GFF))
    dom_gff_tmp.write(version_string)

def domain_search(QUERY, LAST_DB, CLASSIFICATION, OUTPUT_DOMAIN,
                  THRESHOLD_SCORE, WIN_DOM, OVERLAP_DOM, SCORING_MATRIX):
    ''' Search for protein domains using our protein database and external tool LAST,
	stdout is parsed in real time and hits for a single sequence undergo further processing
	- tabular format(TAB) to get info about position, score, orientation
	- MAF format to gain alignment and original sequence
	'''

    step = WIN_DOM - OVERLAP_DOM
    [headers, above_win, below_win, lens_above_win, seq_starts, seq_ends
     ] = characterize_fasta(QUERY, WIN_DOM)
    query_temp = split_fasta(QUERY, WIN_DOM, step, headers, above_win,
                             below_win, lens_above_win, seq_starts, seq_ends)

    ## TAB output contains all the alignment scores, positions, strands...
    lastal_columns = {
        "BL80" : ("score, name_db, start_db, al_size_db, strand_db,"
                   " seq_size_db, name_q, start_q, al_size_q, strand_q, seq_size_q,"
                   " block1, block2, block3, db_seq, q_seq"),
        "BL62" : ("score, name_db, start_db, al_size_db, strand_db,"
                   " seq_size_db, name_q, start_q, al_size_q, strand_q,"
                   " seq_size_q, block1, block2, block3, db_seq, q_seq"),
        "MIQS" : ("score, name_db, start_db, al_size_db, strand_db,"
                  " seq_size_db, name_q, start_q, al_size_q, strand_q,"
                  " seq_size_q, block1, db_seq, q_seq"),
    }
    tab = subprocess.Popen(
        "lastal -F15 {} {} -L 10 -m 70 -p {} -e 80 -f TAB".format(LAST_DB,
                                                                  query_temp,
                                                                  SCORING_MATRIX),
        stdout=subprocess.PIPE,
        shell=True)
    ## MAF output contains alignment sequences
    maf = subprocess.Popen(
        "lastal -F15 {} {} -L 10 -m 70 -p {}  -e 80 -f MAF".format(LAST_DB,
                                                                   query_temp,
                                                                   SCORING_MATRIX),
        stdout=subprocess.PIPE,
        shell=True)

    tab_pipe = tab.stdout
    maf_pipe = maf.stdout
    maf_pipe.readline()
    seq_ids = []
    dom_tmp = NamedTemporaryFile(delete=False)
    with open(dom_tmp.name, "w") as dom_gff_tmp:
        path = os.path.dirname(os.path.realpath(__file__))
        version_string = get_version(path, LAST_DB)
        write_info(dom_gff_tmp, version_string)
    gff = open(dom_tmp.name, "a")
    start = True
    while True:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                sequence_hits = np.genfromtxt(
                    line_generator(tab_pipe, maf_pipe, start),
                    names=lastal_columns[SCORING_MATRIX],
                    usecols=("score, name_q, start_q, al_size_q,"
                             " strand_q, seq_size_q, name_db, db_seq,"
                             " q_seq, seq_size_db, start_db, al_size_db"),
                    dtype=[('score', '<i8'), ('name_q', object), ('start_q', '<i8'),
                    ('al_size_q', '<i8'), ('strand_q', 'S1'), ('seq_size_q', '<i8'),
                    ('name_db', object), ('db_seq', object), ('q_seq', object),
                    ('seq_size_db', '<i8'), ('start_db', '<i8'), ('al_size_db', '<i8')],
                    comments=None)


        except RuntimeError:
            break
        ## if there are no domains found
        if sequence_hits.size == 0:
            gff.write("##No domains found\n")
            gff.close()
            shutil.copy2(dom_tmp.name, OUTPUT_DOMAIN)
            os.unlink(dom_tmp.name)
            return [], [], [], []

        ############# PARSING LASTAL OUTPUT ############################
        sequence_hits = np.atleast_1d(sequence_hits)
        score = sequence_hits['score'].astype("int")
        seq_id = sequence_hits['name_q'].astype("str")[0]
        start_hit = sequence_hits['start_q'].astype("int")
        end_hit = start_hit + sequence_hits['al_size_q'].astype("int")
        strand = sequence_hits['strand_q'].astype("str")
        seq_len = sequence_hits['seq_size_q'].astype("int")
        domain_db = sequence_hits['name_db'].astype("str")
        db_seq = sequence_hits['db_seq'].astype("str")
        query_seq = sequence_hits['q_seq'].astype("str")
        domain_size = sequence_hits['seq_size_db'].astype("int")
        db_start = sequence_hits['start_db'].astype("int") + 1
        db_end = sequence_hits['start_db'].astype("int") + sequence_hits[
            'al_size_db'].astype("int")

        [reverse_strand_idx, positions_plus, positions_minus
         ] = hits_processing(seq_len, start_hit, end_hit, strand)
        strand_gff = "+"
        [mins_plus, maxs_plus, data_plus, indices_plus
         ] = overlapping_regions(positions_plus)
        [mins_minus, maxs_minus, data_minus, indices_minus
         ] = overlapping_regions(positions_minus)
        positions = positions_plus + positions_minus
        indices_overal = indices_plus + [x + reverse_strand_idx
                                         for x in indices_minus]
        mins = mins_plus + mins_minus
        maxs = maxs_plus + maxs_minus
        data = data_plus + data_minus
        ## process every region (cluster) of overlapping hits sequentially
        count_region = 0
        for region in indices_overal:
            db_names = domain_db[np.array(region)]
            db_starts = db_start[np.array(region)]
            db_ends = db_end[np.array(region)]
            scores = score[np.array(region)]
            regions_above_threshold = [
                region[i]
                for i, _ in enumerate(region)
                if max(scores) / 100 * THRESHOLD_SCORE < scores[i]
            ]
            ## sort by score first:
            consensus = get_full_translation(
                translation_alignments(
                    query_seq=sortby(query_seq[regions_above_threshold], score[regions_above_threshold], True),
                    start_hit=sortby(start_hit[regions_above_threshold], score[regions_above_threshold], True),
                    end_hit=sortby(end_hit[regions_above_threshold], score[regions_above_threshold], True))
                )

            annotations = domain_annotation(db_names, CLASSIFICATION)
            [score_matrix, classes_dict] = score_table(
                mins[count_region], maxs[count_region], data[count_region],
                annotations, scores, CLASSIFICATION)
            ann_per_reg = score_matrix_evaluation(score_matrix, classes_dict,
                                                  THRESHOLD_SCORE)
            [domain_type, ann_substring, unique_annotations, ann_pos_counts
             ] = group_annot_regs(ann_per_reg)
            [best_idx, best_idx_reg] = best_score(scores, region)
            annotation_best = annotations[best_idx_reg]
            db_name_best = db_names[best_idx_reg]
            db_starts_best = db_starts[best_idx_reg]
            db_ends_best = db_ends[best_idx_reg]
            if count_region == len(indices_plus):
                strand_gff = "-"
            if strand_gff == "+":
                feature_start = min(start_hit[regions_above_threshold]) + 1
                feature_end = max(end_hit[regions_above_threshold])
            else:
                feature_end = seq_len[region][0] - min(start_hit[regions_above_threshold])
                feature_start = seq_len[region][0] - max(end_hit[regions_above_threshold]) + 1
            create_gff3(domain_type, ann_substring, unique_annotations,
                        ann_pos_counts, feature_start,feature_end,
                        step, best_idx, annotation_best, db_name_best,
                        db_starts_best, db_ends_best, strand_gff, score,
                        seq_id, db_seq, query_seq, domain_size, positions, gff, consensus)
            count_region += 1
        seq_ids.append(seq_id)
    os.unlink(query_temp)
    gff.close()
    dom_tmp.close()
    ## if any sequence from input data was split into windows, merge and adjust the data from individual windows
    if any("DANTE_PART" in x for x in seq_ids):
        adjust_gff(OUTPUT_DOMAIN, dom_tmp.name, WIN_DOM, OVERLAP_DOM, step)
    ## otherwise use the temporary output as the final domains gff
    else:
        shutil.copy2(dom_tmp.name, OUTPUT_DOMAIN)
    os.unlink(dom_tmp.name)

def domain_search(QUERY, LAST_DB, CLASSIFICATION, OUTPUT_DOMAIN,
                  THRESHOLD_SCORE, WIN_DOM, OVERLAP_DOM, SCORING_MATRIX):
    ''' Search for protein domains using our protein database and external tool LAST,
	stdout is parsed from temporary files instead of real-time streams
    '''

    step = WIN_DOM - OVERLAP_DOM
    [headers, above_win, below_win, lens_above_win, seq_starts, seq_ends
     ] = characterize_fasta(QUERY, WIN_DOM)
    query_temp = split_fasta(QUERY, WIN_DOM, step, headers, above_win,
                             below_win, lens_above_win, seq_starts, seq_ends)

    lastal_columns = {
        "BL80" : ("score, name_db, start_db, al_size_db, strand_db,"
                  " seq_size_db, name_q, start_q, al_size_q, strand_q, seq_size_q,"
                  " block1, block2, block3, db_seq, q_seq"),
        "BL62" : ("score, name_db, start_db, al_size_db, strand_db,"
                  " seq_size_db, name_q, start_q, al_size_q, strand_q,"
                  " seq_size_q, block1, block2, block3, db_seq, q_seq"),
        "MIQS" : ("score, name_db, start_db, al_size_db, strand_db,"
                  " seq_size_db, name_q, start_q, al_size_q, strand_q,"
                  " seq_size_q, block1, db_seq, q_seq"),
    }

    # Create temporary files for TAB and MAF outputs
    tab_tmp = NamedTemporaryFile(delete=False)
    tab_tmp_name = tab_tmp.name
    tab_tmp.close()

    maf_tmp = NamedTemporaryFile(delete=False)
    maf_tmp_name = maf_tmp.name
    maf_tmp.close()

    # Run LAST for TAB output, redirecting its stdout to the temp file
    subprocess.run(
        "lastal -F15 {} {} -L 10 -m 70 -p {} -e 80 -f TAB > {}".format(
            LAST_DB, query_temp, SCORING_MATRIX, tab_tmp_name
        ),
        shell=True,
        check=True
    )

    # Run LAST for MAF output, redirecting its stdout to the temp file
    subprocess.run(
        "lastal -F15 {} {} -L 10 -m 70 -p {} -e 80 -f MAF > {}".format(
            LAST_DB, query_temp, SCORING_MATRIX, maf_tmp_name
        ),
        shell=True,
        check=True
    )

    # Now open the resulting TAB and MAF files for parsing
    with open(tab_tmp_name, "rb") as tab_pipe, open(maf_tmp_name, "rb") as maf_pipe:

        # The original code had maf_pipe.readline() to skip the first line
        maf_pipe.readline()

        seq_ids = []
        dom_tmp = NamedTemporaryFile(delete=False)
        with open(dom_tmp.name, "w") as dom_gff_tmp:
            path = os.path.dirname(os.path.realpath(__file__))
            version_string = get_version(path, LAST_DB)
            write_info(dom_gff_tmp, version_string)
        gff = open(dom_tmp.name, "a")
        start = True
        while True:
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    sequence_hits = np.genfromtxt(
                        line_generator(tab_pipe, maf_pipe, start),
                        names=lastal_columns[SCORING_MATRIX],
                        usecols=("score, name_q, start_q, al_size_q,"
                                 " strand_q, seq_size_q, name_db, db_seq,"
                                 " q_seq, seq_size_db, start_db, al_size_db"),
                        dtype=[('score', '<i8'), ('name_q', object), ('start_q', '<i8'),
                               ('al_size_q', '<i8'), ('strand_q', 'S1'), ('seq_size_q', '<i8'),
                               ('name_db', object), ('db_seq', object), ('q_seq', object),
                               ('seq_size_db', '<i8'), ('start_db', '<i8'), ('al_size_db', '<i8')],
                        comments=None
                    )
            except RuntimeError:
                break
            if sequence_hits.size == 0:
                gff.write("##No domains found\n")
                gff.close()
                shutil.copy2(dom_tmp.name, OUTPUT_DOMAIN)
                os.unlink(dom_tmp.name)
                return [], [], [], []


            sequence_hits = np.atleast_1d(sequence_hits)
            score = sequence_hits['score'].astype("int")
            seq_id = sequence_hits['name_q'].astype("str")[0]
            start_hit = sequence_hits['start_q'].astype("int")
            end_hit = start_hit + sequence_hits['al_size_q'].astype("int")
            strand = sequence_hits['strand_q'].astype("str")
            seq_len = sequence_hits['seq_size_q'].astype("int")
            domain_db = sequence_hits['name_db'].astype("str")
            db_seq = sequence_hits['db_seq'].astype("str")
            query_seq = sequence_hits['q_seq'].astype("str")
            domain_size = sequence_hits['seq_size_db'].astype("int")
            db_start = sequence_hits['start_db'].astype("int") + 1
            db_end = sequence_hits['start_db'].astype("int") + sequence_hits[
                'al_size_db'].astype("int")

            [reverse_strand_idx, positions_plus, positions_minus
             ] = hits_processing(seq_len, start_hit, end_hit, strand)
            strand_gff = "+"
            [mins_plus, maxs_plus, data_plus, indices_plus
             ] = overlapping_regions(positions_plus)
            [mins_minus, maxs_minus, data_minus, indices_minus
             ] = overlapping_regions(positions_minus)
            positions = positions_plus + positions_minus
            indices_overal = indices_plus + [x + reverse_strand_idx
                                             for x in indices_minus]
            mins = mins_plus + mins_minus
            maxs = maxs_plus + maxs_minus
            data = data_plus + data_minus

            count_region = 0
            for region in indices_overal:
                db_names = domain_db[np.array(region)]
                db_starts = db_start[np.array(region)]
                db_ends = db_end[np.array(region)]
                scores = score[np.array(region)]
                regions_above_threshold = [
                    region[i]
                    for i, _ in enumerate(region)
                    if max(scores) / 100 * THRESHOLD_SCORE < scores[i]
                ]
                consensus = get_full_translation(
                    translation_alignments(
                        query_seq=sortby(query_seq[regions_above_threshold],
                                         score[regions_above_threshold], True),
                        start_hit=sortby(start_hit[regions_above_threshold],
                                         score[regions_above_threshold], True),
                        end_hit=sortby(end_hit[regions_above_threshold],
                                       score[regions_above_threshold], True))
                )

                annotations = domain_annotation(db_names, CLASSIFICATION)
                [score_matrix, classes_dict] = score_table(
                    mins[count_region], maxs[count_region], data[count_region],
                    annotations, scores, CLASSIFICATION)
                ann_per_reg = score_matrix_evaluation(score_matrix, classes_dict,
                                                      THRESHOLD_SCORE)
                [domain_type, ann_substring, unique_annotations, ann_pos_counts
                 ] = group_annot_regs(ann_per_reg)
                [best_idx, best_idx_reg] = best_score(scores, region)
                annotation_best = annotations[best_idx_reg]
                db_name_best = db_names[best_idx_reg]
                db_starts_best = db_starts[best_idx_reg]
                db_ends_best = db_ends[best_idx_reg]
                if count_region == len(indices_plus):
                    strand_gff = "-"
                if strand_gff == "+":
                    feature_start = min(start_hit[regions_above_threshold]) + 1
                    feature_end = max(end_hit[regions_above_threshold])
                else:
                    feature_end = seq_len[region][0] - min(start_hit[regions_above_threshold])
                    feature_start = (seq_len[region][0]
                                     - max(end_hit[regions_above_threshold]) + 1)
                create_gff3(domain_type, ann_substring, unique_annotations,
                            ann_pos_counts, feature_start,feature_end,
                            step, best_idx, annotation_best, db_name_best,
                            db_starts_best, db_ends_best, strand_gff, score,
                            seq_id, db_seq, query_seq, domain_size, positions, gff, consensus)
                count_region += 1
            seq_ids.append(seq_id)

        gff.close()
        dom_tmp.close()

    # Cleanup the temp TAB/MAF files
    # os.unlink(tab_tmp_name)
    # os.unlink(maf_tmp_name)
    # os.unlink(query_temp)

    # If any sequence from input data was split into windows, merge and adjust
    if any("DANTE_PART" in x for x in seq_ids):
        adjust_gff(OUTPUT_DOMAIN, dom_tmp.name, WIN_DOM, OVERLAP_DOM, step)
    else:
        shutil.copy2(dom_tmp.name, OUTPUT_DOMAIN)
    # os.unlink(dom_tmp.name)


def  sortby(a, by, reverse=False):
    ''' sort according values in the by list '''
    a_sorted = [i[0] for i in
                sorted(
                    zip(a, by),
                    key=lambda i: i[1],
                    reverse=reverse
                )]
    return a_sorted


def a2nnn(s):
    s1 = "".join([c if c in ['/', '\\'] else c + c + c
                  for c in s.replace("-", "")])
    # collapse frameshifts (/)
    s2 = re.sub("[a-zA-Z*]{2}//[a-zA-Z*]{2}", "//", s1)
    s3 = re.sub("[a-zA-Z*]/[a-zA-Z*]", "/", s2)
    return (s3)



def rle(s):
    '''run length encoding but max is set to 3 (codon)'''
    prev = ""
    count = 1
    char = []
    length = []
    L = 0
    for n in s:
        if n == prev and count < (3 - L):
            count += 1
        else:
            char.append(prev)
            length.append(count)
            L = 1 if prev == "/" else 0
            prev = n
            count = 1
    char.append(prev)
    length.append(count)
    sequence = char[1:]
    return sequence, length[1:]

def get_full_translation_old(translations):
    '''get one full length translation  from multiple partial
    aligned translation '''
    # find minimal set of alignments
    minimal_set = []
    not_filled_prev = len(translations[0])
    for s in translations:
        minimal_set.append(s)
        # check defined position - is there only '-' character?
        not_filled = sum([set(i) == {"-"} for i in  zip(*minimal_set)])
        if not_filled == 0:
            break
        if not_filled == not_filled_prev:
            # last added sequence did not improve coverage - remove it.
            minimal_set.pop()
        not_filled_prev = not_filled
    # merge translations
    final_translation = minimal_set[0]
    # record positions of joins to correct frameshifts reportings
    position_of_joins = set()
    position_of_joins_rle = set()
    if len(minimal_set) > 1:  # translation must be merged
        for s in minimal_set[1:]:
            s1 = re.search(r"-\w", final_translation)
            s2 = re.search(r"\w-", final_translation)
            if s1:
                position_of_joins.add(s1.start())
            if s2:
                position_of_joins.add((s2.end() - 1))
            final_translation = "".join(
                [b if a == "-" else a for a, b in zip(final_translation, s)])
    translation_rle = rle(final_translation)
    cumsumed_positions = np.cumsum(translation_rle[1])
    for p in position_of_joins:
        position_of_joins_rle.add(sum(cumsumed_positions <= p))
    # insert /\ when necessary
    for p in position_of_joins_rle:
        if translation_rle[0][p] not in ['/',"//","\\", "\\\\"]:
            if translation_rle[1][p] == 2:
                translation_rle[0][p] = translation_rle[0][p] + "/"
            if translation_rle[1][p] == 1:
                translation_rle[0][p] = "\\"
    consensus = "".join(translation_rle[0])
    return consensus


def get_full_translation(translations):
    """
    Create one full-length translation from multiple partial alignments.
    More efficient approach:
      1) pick only those partial sequences that add coverage,
      2) merge them in a single pass,
      3) do a single RLE pass to insert frameshift markers.
    """

    # 1) Track coverage to pick a minimal set of sequences
    M = len(translations[0])  # assume all have same length
    coverage = [False] * M
    chosen = []

    for s in translations:
        old_uncovered = coverage.count(False)
        # Update coverage
        for i, ch in enumerate(s):
            if ch != "-":
                coverage[i] = True
        new_uncovered = coverage.count(False)
        # Keep if coverage improved
        if new_uncovered < old_uncovered:
            chosen.append(s)
        # If fully covered, we can stop early
        if new_uncovered == 0:
            break

    # 2) Merge all chosen sequences in one pass
    # If none chosen, fallback to returning the first or empty
    if not chosen:
        return translations[0] if translations else ""

    merged = list(chosen[0])  # make it a list for mutability
    join_positions = set()

    for s in chosen[1:]:
        for i, (mch, sch) in enumerate(zip(merged, s)):
            # if we have a dash in merged and a real char in s => merge
            if mch == "-" and sch != "-":
                merged[i] = sch
                join_positions.add(i)

    final_str = "".join(merged)

    # 3) Run-length encode once
    val_list, run_list = rle(final_str)
    cumsum_positions = np.cumsum(run_list)

    # Convert join_positions into "RLE indices"
    join_rle_positions = set()
    for pos in join_positions:
        # find which RLE chunk pos belongs to
        chunk_index = np.searchsorted(cumsum_positions, pos, side='right')
        join_rle_positions.add(chunk_index)

    # 4) Insert slashes/backslashes
    for idx in join_rle_positions:
        if val_list[idx] not in ["/","//","\\","\\\\"]:
            if run_list[idx] == 2:
                val_list[idx] += "/"
            elif run_list[idx] == 1:
                val_list[idx] = "\\"

    # 5) Build final consensus
    return "".join(val_list)



# helper function for debugging
def translation_alignments(query_seq, start_hit, end_hit):
    pstart = min(start_hit)
    pend = max(end_hit)
    nnn = list()
    for s, start, end in zip(query_seq, start_hit, end_hit):
        nnn.append("-" * (start - pstart) + a2nnn(s) + "-" * (pend - end))
    return (nnn)



def adjust_gff(OUTPUT_DOMAIN, gff, WIN_DOM, OVERLAP_DOM, step):
    ''' Original gff file is adjusted in case of containing cut parts 
	- for consecutive sequences overlap is divided to half with first half 
	of records(domains) belonging to the first sequence and second to the following one.
	Duplicate domains going through the middle of the overlap are removed.
	First and the last part (marked as LAST) of a certain sequence are 
	handled separately as the are overlapped from one side only '''

    seq_id_all = []
    class_dict = defaultdict(int)
    seen = set()
    with open(OUTPUT_DOMAIN, "w") as adjusted_gff:
        with open(gff, "r") as primary_gff:
            start = True
            for line in primary_gff:
                if line.startswith("#"):
                    adjusted_gff.write(line)
                else:
                    split_line = line.split("\t")
                    classification = split_line[-1].split(";")[1].split("=")[1]
                    if start:
                        seq_id_all.append(split_line[0].split("_DANTE_PART")[
                            0])
                        start = False
                    seq_id = split_line[0].split("_DANTE_PART")[0]
                    if "DANTE_PART" in line:
                        line_without_id = "\t".join(split_line[1:])
                        part = int(split_line[0].split("_DANTE_PART")[1].split(
                            ":")[0].split("_")[0])
                        if seq_id != seq_id_all[-1]:
                            seq_id_all.append(seq_id)

                            ## first part of the sequence
                        if part == 1:
                            cut_end = WIN_DOM - OVERLAP_DOM / 2
                            if int(split_line[3]) <= cut_end <= int(split_line[
                                    4]):
                                if line_without_id not in seen:
                                    adjusted_gff.write("{}\t{}".format(
                                        seq_id, line_without_id))
                                    class_dict[classification] += 1
                                    seen.add(line_without_id)
                            elif int(split_line[4]) < cut_end:
                                adjusted_gff.write("{}\t{}".format(
                                    seq_id, line_without_id))
                                class_dict[classification] += 1

                                ## last part of the sequence
                        elif "LAST" in split_line[0]:
                            cut_start = OVERLAP_DOM / 2 + (part - 1) * step
                            if int(split_line[3]) <= cut_start <= int(
                                    split_line[4]):
                                if line_without_id not in seen:
                                    adjusted_gff.write("{}\t{}".format(
                                        seq_id, line_without_id))
                                    class_dict[classification] += 1
                                    seen.add(line_without_id)
                            elif int(split_line[3]) > cut_start:
                                adjusted_gff.write("{}\t{}".format(
                                    seq_id, line_without_id))
                                class_dict[classification] += 1

                        ## middle part of the sequence
                        else:
                            cut_start = OVERLAP_DOM / 2 + (part - 1) * step
                            cut_end = WIN_DOM - OVERLAP_DOM / 2 + (part -
                                                                   1) * step
                            if int(split_line[3]) <= cut_start <= int(
                                    split_line[4]) or int(split_line[
                                        3]) <= cut_end <= int(split_line[4]):
                                if line_without_id not in seen:
                                    adjusted_gff.write("{}\t{}".format(
                                        seq_id, line_without_id))
                                    class_dict[classification] += 1
                                    seen.add(line_without_id)
                            elif int(split_line[3]) > cut_start and int(
                                    split_line[4]) < cut_end:
                                adjusted_gff.write("{}\t{}".format(
                                    seq_id, line_without_id))
                                class_dict[classification] += 1
                    ## not divived
                    else:
                        if seq_id != seq_id_all[-1]:
                            seq_id_all.append(seq_id)
                        adjusted_gff.write(line)
                        class_dict[classification] += 1


def main(args):

    t = time.time()

    QUERY = args.query
    LAST_DB = args.protein_database
    CLASSIFICATION = args.classification
    OUTPUT_DOMAIN = args.domain_gff
    NEW_LDB = args.new_ldb
    OUTPUT_DIR = args.output_dir
    THRESHOLD_SCORE = args.threshold_score
    WIN_DOM = args.win_dom
    OVERLAP_DOM = args.overlap_dom
    SCORING_MATRIX = args.scoring_matrix
    configuration.SC_MATRIX = configuration.SC_MATRIX_SKELETON.format(SCORING_MATRIX)

    if OUTPUT_DOMAIN is None:
        OUTPUT_DOMAIN = configuration.DOMAINS_GFF
    if os.path.isdir(LAST_DB):
        LAST_DB = os.path.join(LAST_DB, configuration.LAST_DB_FILE)
    if os.path.isdir(CLASSIFICATION):
        CLASSIFICATION = os.path.join(CLASSIFICATION, configuration.CLASS_FILE)

    if NEW_LDB:
        subprocess.call("lastdb -p -cR01 {} {}".format(LAST_DB, LAST_DB),
                        shell=True)

    if OUTPUT_DIR and not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    if not os.path.isabs(OUTPUT_DOMAIN):
        if OUTPUT_DIR is None:
            OUTPUT_DIR = configuration.TMP
            if not os.path.exists(OUTPUT_DIR):
                os.makedirs(OUTPUT_DIR)
        OUTPUT_DOMAIN = os.path.join(OUTPUT_DIR,OUTPUT_DOMAIN)
    domain_search(QUERY, LAST_DB, CLASSIFICATION, OUTPUT_DOMAIN,
                  THRESHOLD_SCORE, WIN_DOM, OVERLAP_DOM, SCORING_MATRIX)

    print("ELAPSED_TIME_DOMAINS = {} s".format(time.time() - t))


if __name__ == "__main__":
    import argparse
    from argparse import RawDescriptionHelpFormatter

    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                          argparse.RawDescriptionHelpFormatter):
        pass

    parser = argparse.ArgumentParser(
        description=
        '''Script performs similarity search on given DNA sequence(s) in (multi)fasta against our protein domains database of all Transposable element for certain group of organisms (Viridiplantae or Metazoans). Domains are subsequently annotated and classified - in case certain domain has multiple annotations assigned, classifation is derived from the common classification level of all of them. Domains search is accomplished engaging LASTAL alignment tool.
		
	DEPENDENCIES:
		- python 3.4 or higher with packages:
			-numpy
		- lastal 744 or higher [http://last.cbrc.jp/]
		- configuration.py module

	EXAMPLE OF USAGE:
		
		./protein_domains_pd.py -q PATH_TO_INPUT_SEQ -pdb PATH_TO_PROTEIN_DB -cs PATH_TO_CLASSIFICATION_FILE
		
	When running for the first time with a new database use -nld option allowing lastal to create indexed database files:
		
		-nld True
	
		''',
        epilog="""""",
        formatter_class=CustomFormatter)

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        "-q",
        "--query",
        type=str,
        required=True,
        help=
        'input DNA sequence to search for protein domains in a fasta format. Multifasta format allowed.')
    requiredNamed.add_argument('-pdb',
                               "--protein_database",
                               type=str,
                               required=True,
                               help='protein domains database file')
    requiredNamed.add_argument('-cs',
                               '--classification',
                               type=str,
                               required=True,
                               help='protein domains classification file')
    parser.add_argument("-oug",
                        "--domain_gff",
                        type=str,
                        help="output domains gff format")
    parser.add_argument(
        "-nld",
        "--new_ldb",
        type=str,
        default=False,
        help=
        "create indexed database files for lastal in case of working with new protein db")
    parser.add_argument(
        "-dir",
        "--output_dir",
        type=str,
        help="specify if you want to change the output directory")
    parser.add_argument(
        "-M",
        "--scoring_matrix",
        type=str,
        default="BL80",
        choices=['BL80', 'BL62', 'MIQS'],
        help="specify scoring matrix to use for similarity search (BL80, BL62, MIQS)")
    parser.add_argument(
        "-thsc",
        "--threshold_score",
        type=int,
        default=80,
        help=
        "percentage of the best score in the cluster to be tolerated when assigning annotations per base")
    parser.add_argument(
        "-wd",
        "--win_dom",
        type=int,
        default=10000000,
        help="window to process large input sequences sequentially")
    parser.add_argument("-od",
                        "--overlap_dom",
                        type=int,
                        default=10000,
                        help="overlap of sequences in two consecutive windows")

    args = parser.parse_args()
    main(args)
