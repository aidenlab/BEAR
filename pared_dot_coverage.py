import numpy as np
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import csv, math, argparse
import base64
from io import BytesIO
import pandas as pd
import itertools
 
bin_size = 300
width_of_bars=300
bar_color = "#3498DB"
line_palette = ['#FFEB3B', "#4CAF50", "#1B5E20"]
x_labels = "SARS-COV2 RefSeq Assembly"
y_labels = "De Novo SARS-COV2 Assembly"


def fill_blanks(f_txt, max_length):
	'''Fills in histogram blanks.

	Inputs: f_txt- .txt file string with position read counts
			max_length- int length of genome
	Outputs: filled- np array with 0s filled in
	'''
	raw = pd.read_csv(f_txt, sep="	", names=['first', 'second', 'third'], header=None)
	raw_np_array = np.zeros((len(raw),2))
	
	raw_np_array[:,0] = raw['second'].values
	raw_np_array[:,1] = raw['third'].values
	
	filled = np.zeros((max_length, 2))
	filled[:,0] = np.arange(1, max_length+1)
	filled[np.isin(filled[:,0], raw_np_array[:,0]),1] = raw_np_array[:,1]

	return filled


def make_hist(data, max_length):
	'''Makes a histogram from position read counts

	Inputs: data- np array of filled in posiiton read counts
			max_length- int length of genome
	Outputs: bin_pos- list of locations for histogram bars on x axis
			bin_counts- list of counts per bin
	'''
	bins = np.arange(1, max_length+bin_size, bin_size)
	digitized = np.digitize(data[:,0], bins)
	bin_counts = [data[digitized == i,1].sum() for i in range(1, len(bins))]
	bin_pos = [(bins[i]+bins[i+1]-1)/2.0 for i in range(len(bins)-1)]
	
	return bin_pos, bin_counts


def plot_alignment(ax_diagnostic, ax_align, f_csv):
	'''
	Plots horizontal bar graph of alignment percentages

	Inputs: ax_diagnostic- diagnostic symbol axis
			ax_align- alignment bar graph axis
			f_csv- string of .csv file for coronavirus alignment percentages
	Outputs: ax_diagnostic- diagnostic symbol axis
			ax_align- alignment bar graph axis
	'''
	align_data = pd.read_csv(f_csv)
	ax_align.barh(align_data['label'], align_data['percentage'], color="#D2B4DE")
	ax_align.set_xlim(0,100)
	ax_align.tick_params(axis="y",direction="in", left="true", pad=-5)
	plt.setp(ax_align.get_yticklabels(), ha="left")

	covid_val = align_data[align_data['label']==
		'wuhan_seafood_market_pneumonia_virus_isolate_Wuhan_Hu_1']['percentage'].values[0]
	
	ax_diagnostic.plot([.425, .575], [.5, .5], linewidth=20, color='black')
	if covid_val >= 75:
		ax_diagnostic.plot([.5, .5], [.1, .9], linewidth=20, color='black')
		ax_diagnostic.set_title("Positive",)
	else:
		ax_diagnostic.set_title("Negative")
	
	ax_diagnostic.set_xlim(0,1)
	ax_diagnostic.set_ylim(0,1)
	ax_diagnostic.axis('off')
	
	return ax_diagnostic, ax_align


def plot_dot_plot(ax_coverage, ax_dot, f_paf, max_length, bin_pos, bin_counts):
	'''Create dot plot and coverage track.

	Inputs: ax_coverage- coverage track axis
			ax_dot- dot plot axis
			f_paf- string of .paf file for dot plot
			max_length- int length of genome
			bin_pos- list of locations for histogram bars on x axis
			bin_counts- list of counts per bin
	Outputs: ax_coverage- coverage track axis
			ax_dot- dot plot axis
	'''
	col_names = [str(i) for i in range(18)]
	dot_data = pd.read_csv(f_paf, names=col_names, delimiter='	')
	
	dot_data = dot_data.sort_values('7')
	dot_data['cum_offset'] = dot_data['1'].cumsum()
	dot_data['upto_offset'] = dot_data['cum_offset'].values - dot_data['1'].values

	def dot_plot(x1, x2, y1, y2, segment_offset, num, denom):
		col_fract = float(num)/float(denom)
		col_idx = int(math.ceil(col_fract*len(line_palette)))-1
		col = line_palette[col_idx]
		ax_dot.plot([x1, x2], [y1+segment_offset, y2+segment_offset], linewidth=2, color=col)

	dot_data.apply(lambda row : dot_plot(row['7'], row['8'],
                     row['2'], row['3'], row['upto_offset'], row['9'], row['10']), axis = 1)

	ax_coverage.set_xlim(0, max_length)
	ax_dot.set_xlim((0, max_length))
	ax_dot.set_ylim((0, max_length))

	ax_dot.set_xlabel(x_labels)
	ax_dot.set_ylabel(y_labels)

	ax_coverage.bar(bin_pos, bin_counts, width=width_of_bars, color=bar_color)

	ax_coverage.yaxis.set_major_locator(plt.MaxNLocator(3))
	ax_coverage.xaxis.set_major_locator(plt.MaxNLocator(4))
	ax_dot.yaxis.set_major_locator(plt.MaxNLocator(4))
	ax_dot.xaxis.set_major_locator(plt.MaxNLocator(4))

	return ax_coverage, ax_dot


def plot(f_paf, f_csv, max_length, bin_pos, bin_counts, write_file):
	'''Plots coverage track on top of dot plot. Writes .pdf

	Inputs: f_paf- string of .paf file for dot plot
			f_csv- string of .csv file for alignment plot
			max_length- int length of genome
			bin_pos- list of locations for histogram bars on x axis
			bin_counts- list of counts per bin
			write_file- string name of file to write plots to
			
	Output: fig- matplotlib fig of coverage track, dot plot, alignment plot
	'''

	fig, axs = plt.subplots(2, 2, sharex=False, sharey=False, figsize=(14.35,10), 
		gridspec_kw={'height_ratios': [1, 9], 'width_ratios':[9, 5.35]})

	axs[0][0], axs[1][0] = plot_dot_plot(axs[0][0], axs[1][0], f_paf, max_length, 
		bin_pos, bin_counts)
	axs[0][1], axs[1][1] = plot_alignment(axs[0][1], axs[1][1], f_csv)

	fig.tight_layout()	
	plt.savefig(write_file+'.pdf',
		dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches=None, pad_inches=0.01,
        frameon=False, metadata=None)
	return fig


def main(f_txt, f_paf, f_csv, max_length, write_file):
	'''Reads .txt + .paf, plots coverage track + dot plot, writes plots to file

	Inputs: f_txt- .txt file string to read for histogram
			f_paf- .paf file string to read for dot plot
			f_csv- for alignment bar graph
			max_length- int length of viral genome
			write_file- string name of file to write plots to
	'''
	filled = fill_blanks(f_txt, max_length)
	bin_pos, bin_counts = make_hist(filled, max_length)
	fig = plot(f_paf, f_csv, max_length, bin_pos, bin_counts, write_file)


if __name__ == "__main__":
	'''
	Example usage:
	python pared_dot_coverage.py ".txt file for coverage" ".paf file for dot plot" ".csv for alignment" 
		"int genome length" "destination for plots"
	python pared_dot_coverage.py depth_per_base.txt contig_nCoV-2019.paf stats.csv 29903 covid_plots 
	'''
	parser = argparse.ArgumentParser(description='Plot histogram.')
	parser.add_argument('read_txt_file', metavar='rf_txt', type=str,
                     help='.txt file to read')
	parser.add_argument('read_paf_file', metavar='rf_paf', type=str,
                     help='.paf file to read')
	parser.add_argument('read_csv_file', metavar='rf_csv', type=str,
                     help='.csv file to read')

	parser.add_argument('max_length', metavar='m', type=int, 
                     help='Genome length')

	parser.add_argument('write_file', metavar='wf', type=str,
                     help='png/svg files to write')

	args = parser.parse_args()
	main(args.read_txt_file, args.read_paf_file, args.read_csv_file, 
		args.max_length, args.write_file)
	

