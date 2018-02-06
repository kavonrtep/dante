#!/usr/bin/env python3
""" visualization module """

import numpy as np
import configuration
import matplotlib.pyplot as plt
import matplotlib.lines as mlines


	
	
def vis_profrep(seq_ids_all, files_dict, seq_lengths_all, CN, HTML_DATA, seqs_all_part):
	''' visualization of repetitive profiles'''
	graphs_dict = {}	
	seq_id_repeats = []
	th_length = configuration.SEQ_LEN_VIZ
	print(th_length)
	exclude = set(['ALL'])
	sorted_keys =  sorted(set(files_dict.keys()).difference(exclude))
	sorted_keys.insert(0, "ALL")
	plot_num = 0 
	seqs_long = []
	seqs_count = 1
	seqs_max_limit = []
	for repeat in sorted_keys:
		with open(files_dict[repeat][0], "r") as repeat_f:
			positions_all = []
			hits_all = []
			include = True
			first_line = repeat_f.readline()
			seq_id_repeat = first_line.rstrip().split("chrom=")[1]
			seq_len_repeat = seq_lengths_all[seq_ids_all.index(seq_id_repeat)]
			if seq_id_repeat not in graphs_dict.keys():
				if seq_len_repeat > th_length:
					if seq_id_repeat not in seqs_long:
						seqs_long.append(seq_id_repeat)
					include = False
				else:
					[fig, ax] = plot_figure(seq_id_repeat, seq_len_repeat, CN)	
					graphs_dict[seq_id_repeat] = [fig, ax]
			seq_id_repeats.append(seq_id_repeat)
			for line in repeat_f:
				if "chrom" in line:
					seqs_count += 1
					if include:
						graphs_dict = plot_profile(graphs_dict, seq_id_repeats[-1], positions_all, hits_all, repeat, plot_num)
						positions_all = []
						hits_all = []
					seq_id_repeat = line.rstrip().split("chrom=")[1]
					seq_len_repeat = seq_lengths_all[seq_ids_all.index(seq_id_repeat)]
					if seq_id_repeat not in graphs_dict.keys():
						if seq_len_repeat > th_length:
							if seq_id_repeat not in seqs_long:
								seqs_long.append(seq_id_repeat)
							include = False
						else:	
							[fig, ax] = plot_figure(seq_id_repeat, seq_len_repeat, CN)	
							graphs_dict[seq_id_repeat] = [fig, ax]
					seq_id_repeats.append(seq_id_repeat)
					if seq_id_repeat not in seqs_all_part:
						break
				else:
					if include:
						positions_all.append(line.rstrip().split("\t")[0])
						hits_all.append(line.rstrip().split("\t")[1])
		if include:
			graphs_dict = plot_profile(graphs_dict, seq_id_repeats[-1], positions_all, hits_all, repeat, plot_num)
			seq_id_repeats.append(seq_id_repeat)
			positions_all = []
			hits_all = []	
		plot_num += 1	
	print(seqs_long)
	print(graphs_dict)	
	return graphs_dict, seqs_long

		#def vis_profrep(seq_ids_all, files_dict, seq_lengths_all, CN, HTML_DATA):
	#''' visualization of repetitive profiles'''
	##header = seq_repeats.dtype.names
	##seq_id = header[0]
	#graphs_dict = {}	
	#seq_id_repeats = []
	#th_length = configuration.SEQ_LEN_VIZ
	#exclude = set(['ALL'])
	#sorted_keys =  sorted(set(files_dict.keys()).difference(exclude))
	#sorted_keys.insert(0, "ALL")
	#plot_num = 0 
	#for repeat in sorted_keys:
		#with open(files_dict[repeat][0], "r") as repeat_f:
			#positions_all = []
			#hits_all = []
			#include = True
			###################################!!!!!!!!!!!!!!!!!!!!!!
			#first_line = repeat_f.readline()
			#seq_id_repeat = first_line.rstrip().split("chrom=")[1]
			#seq_len_repeat = seq_lengths_all[seq_ids_all.index(seq_id_repeat)]
			#if seq_id_repeat not in graphs_dict.keys():
				#if seq_len_repeat > th_length:
					#include = False
				#else:
					#[fig, ax] = plot_figure(seq_id_repeat, seq_len_repeat, CN)	
					#graphs_dict[seq_id_repeat] = [fig, ax]
			#seq_id_repeats.append(seq_id_repeat)
			#########################################################
			#for line in repeat_f:
				#if "chrom" in line:
					#if include:
						#graphs_dict = plot_profile(graphs_dict, seq_id_repeats[-1], positions_all, hits_all, repeat, plot_num)
						#positions_all = []
						#hits_all = []
					#####################################################
					#seq_id_repeat = line.rstrip().split("chrom=")[1]
					#seq_len_repeat = seq_lengths_all[seq_ids_all.index(seq_id_repeat)]
					#if seq_id_repeat not in graphs_dict.keys():
						#if seq_len_repeat > th_length:
							#include = False
						#else:	
							##if seq_id_repeat not in graphs_dict.keys():
							#[fig, ax] = plot_figure(seq_id_repeat, seq_len_repeat, CN)	
							#graphs_dict[seq_id_repeat] = [fig, ax]
					#seq_id_repeats.append(seq_id_repeat)
						##fig = graphs_dict[seq_id_repeats[-1]][0]
						##ax = graphs_dict[seq_id_repeats[-1]][1]
						##plot_num = graphs_dict[seq_id_repeats[-1]][2] 
						##print(seq_id_repeats)
						##if seq_id == seq_id_repeat:
						
						##plot_profile(fig, ax, positions_all, hits_all, plot_num)
						#####graphs_dict[seq_id_repeats[-1]][2] += 1	###################	
				#else:
					#if include:
						#positions_all.append(line.rstrip().split("\t")[0])
						#hits_all.append(line.rstrip().split("\t")[1])
		#if include:
			#graphs_dict = plot_profile(graphs_dict, seq_id_repeats[-1], positions_all, hits_all, repeat, plot_num)
			#seq_id_repeats.append(seq_id_repeat)
			##if seq_id == seq_id_repeat:
			#positions_all = []
			#hits_all = []
		##plot_profile(fig, ax, positions_all, hits_all, plot_num)
		##graphs_dict[seq_id_repeats[-1]][2] += 1		
		#plot_num += 1		
	#return graphs_dict
	
			
def plot_figure(seq_id, seq_length, CN):
	fig = plt.figure(figsize=(18, 8))
	ax = fig.add_subplot(111)
	ax.set_xlabel('sequence bp')
	if CN:
		ax.set_ylabel('copy numbers')
	else:
		ax.set_ylabel('hits')	
	ax.set_title(seq_id)
	plt.xlim([0, seq_length])	
	return fig, ax


def plot_profile(graphs_dict, seq_id_repeat, positions_all, hits_all, repeat, plot_num):
	if "|" in repeat:
		graphs_dict[seq_id_repeat][1].plot(positions_all, hits_all, label="|".join(repeat.split("|")[-2:]), color=configuration.COLORS_HEX[plot_num])
	else:
		graphs_dict[seq_id_repeat][1].plot(positions_all, hits_all, label=repeat, color=configuration.COLORS_HEX[plot_num])
	return graphs_dict

	
def vis_domains(fig, ax, seq_id, xminimal, xmaximal, domains):
	''' visualization of protein domains'''
	y_upper_lim = ax.get_ylim()[1]
	dom_uniq = list(set(domains)) 
	colors = [configuration.COLORS_HEX[dom_uniq.index(domain)] for domain in domains]
	colors_dom = [list(reversed(configuration.COLORS_HEX))[dom_uniq.index(domain)] for domain in domains]
	colors_legend = list(reversed(configuration.COLORS_HEX))[0:len(dom_uniq)]
	ax.hlines([y_upper_lim + y_upper_lim/10]*len(xminimal), xminimal, xmaximal, color=colors_dom, lw=2, label = dom_uniq)
	lines_legend = [] 
	ax2 = ax.twinx()	# add second axis for domains
	for count_uniq in list(range(len(dom_uniq))):
		lines_legend.append(mlines.Line2D([], [], color=colors_legend[count_uniq], markersize=15, label=dom_uniq[count_uniq]))
	ax2.legend(lines_legend, [line.get_label() for line in lines_legend], bbox_to_anchor=(1.05, 1),  loc='upper left', borderaxespad=0.)
	ax2.yaxis.set_visible(False)
	return fig, ax
	
def main():
	pass

if __name__ == "__main__":
	main()
	main()
	


