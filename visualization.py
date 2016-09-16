#!/usr/bin/env python3
""" visualization module """

import numpy as np
import configuration
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

def vis_profrep(seq_repeats, seq_length):
	header = seq_repeats.dtype.names
	seq_id = header[0]
	fig = plt.figure(figsize=(18, 8))
	ax = fig.add_subplot(111)
	ax.set_xlabel('sequence bp')
	ax.set_ylabel('copy numbers')
	ax.set_title(seq_id)
	plt.xlim([0, seq_length])	
	plot_num = 0
	if not any(seq_repeats.shape):
		ax.hlines(0, 0, seq_length, color="red", lw=4)
	else:
		ax.plot(seq_repeats[seq_id], seq_repeats[header[1]], label=header[1], color="#7F7F7F")
		for repeat in header[2:]: 
			ax.plot(seq_repeats[seq_id], seq_repeats[repeat], label=repeat, color=configuration.COLORS_HEX[plot_num])
			plot_num += 1
		art = []
		lgd = ax.legend(bbox_to_anchor=(0.5,-0.1 ), loc=9, ncol=3)
		art.append(lgd)
	return fig, ax
	
	
def vis_domains(fig, ax, seq_id, xminimal, xmaximal, domains):
	y_upper_lim = ax.get_ylim()[1]
	dom_uniq = list(set(domains))
	colors = [configuration.COLORS_HEX[dom_uniq.index(domain)] for domain in domains]
	colors_dom = [list(reversed(configuration.COLORS_HEX))[dom_uniq.index(domain)] for domain in domains]
	colors_legend = list(reversed(configuration.COLORS_HEX))[0:len(dom_uniq)]
	ax.hlines([y_upper_lim + y_upper_lim/10]*len(xminimal), xminimal, xmaximal, color=colors_dom, lw=2, label = dom_uniq)
	lines_legend = [] 
	ax2 = ax.twinx()
	for count_uniq in list(range(len(dom_uniq))):
		lines_legend.append(mlines.Line2D([], [], color=colors_legend[count_uniq], markersize=15, label=dom_uniq[count_uniq]))
	ax2.legend(lines_legend, [line.get_label() for line in lines_legend], bbox_to_anchor=(1.05, 1),  loc='upper left', borderaxespad=0.)
	ax2.yaxis.set_visible(False)
	return fig, ax
	


	


