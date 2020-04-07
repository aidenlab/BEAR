import numpy as np
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import csv, math, argparse
import base64
from io import BytesIO
import pandas as pd
import itertools
 
bin_size = 5
width_of_bars = 4
bar_color = "#1565C0"
line_palette = ['#FFEB3B', "#4CAF50", "#1B5E20"]
x_labels = "SARS-COV2 RefSeq Assembly"
y_labels = "De Novo SARS-COV2 Assembly"


def fill_blanks(f_txt, x_length):
	'''Fills in histogram blanks.

	Inputs: f_txt- .txt file string with position read counts
			x_length- int length of x axis
	Outputs: filled- np array with 0s filled in
	'''
	raw = pd.read_csv(f_txt, sep="	", names=['first', 'second', 'third'], header=None)
	raw_np_array = np.zeros((len(raw),2))
	
	raw_np_array[:,0] = raw['second'].values
	raw_np_array[:,1] = raw['third'].values
	
	filled = np.zeros((x_length, 2))
	filled[:,0] = np.arange(1, x_length+1)
	filled[np.isin(filled[:,0], raw_np_array[:,0]),1] = raw_np_array[:,1]

	return filled


def make_hist(data, x_length):
	'''Makes a histogram from position read counts

	Inputs: data- np array of filled in posiiton read counts
			x_length- int length of genome
	Outputs: bin_pos- list of locations for histogram bars on x axis
			bin_counts- list of counts per bin
	'''
	bins = np.arange(1, x_length+bin_size, bin_size)
	digitized = np.digitize(data[:,0], bins)
	bin_counts = [data[digitized == i,1].sum() for i in range(1, len(bins))]
	bin_pos = [(bins[i]+bins[i+1]-1)/2.0 for i in range(len(bins)-1)]
	
	return bin_pos, bin_counts


def thicker_spines(ax, all_spines):
	'''Makes spines of axes thicker.

	Inputs: ax- matplot axis object
			all_spines- boolean, do we change all spines?
	Outputs: ax- matplot axis object
	'''
	spine_thickness = 2
	if all_spines:
		ax.spines['right'].set_linewidth(spine_thickness)
		ax.spines['top'].set_linewidth(spine_thickness)
	ax.spines['left'].set_linewidth(spine_thickness)
	ax.spines['bottom'].set_linewidth(spine_thickness)
	return ax


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
	
	ax_diagnostic.plot([.445, .555], [.5, .5], linewidth=20, color='black')
	if covid_val >= 75:
		ax_diagnostic.plot([.5, .5], [.1, .9], linewidth=20, color='black')
		ax_diagnostic.set_title("Positive",)
	else:
		ax_diagnostic.set_title("Negative")
	
	ax_diagnostic.set_xlim(0,1)
	ax_diagnostic.set_ylim(0,1)
	ax_diagnostic.axis('off')

	ax_align = thicker_spines(ax_align, True)
	
	return ax_diagnostic, ax_align


def plot_dot_plot(ax_coverage, ax_dot, filled, dot_data, x_length, y_length, bin_pos, bin_counts):
	'''Create dot plot and coverage track.

	Inputs: ax_coverage- coverage track axis
			ax_dot- dot plot axis
			filled- np array data for coverage track
			dot_data- np array data for dot plot
			x_length- int length of dot plot x axis
			y_length- int length of dot plot y axis
			bin_pos- list of locations for histogram bars on x axis
			bin_counts- list of counts per bin
	Outputs: ax_coverage- coverage track axis
			ax_dot- dot plot axis
	'''
	dot_data = dot_data.sort_values('7')
	dot_data['cum_offset'] = dot_data['1'].cumsum()
	dot_data['upto_offset'] = dot_data['cum_offset'].values - dot_data['1'].values

	def dot_plot(x1, x2, y1, y2, segment_offset, num, denom):
		col_fract = float(num)/float(denom)
		col_idx = int(math.ceil(col_fract*len(line_palette)))-1
		col = line_palette[col_idx]
		ax_dot.plot([x1, x2], [y1+segment_offset, y2+segment_offset], linewidth=3, color=col)

	dot_data.apply(lambda row : dot_plot(row['7'], row['8'],
                     row['2'], row['3'], row['upto_offset'], row['9'], row['10']), axis = 1)

	ax_coverage.set_xlim(0, x_length)
	ax_dot.set_xlim((0, x_length))
	ax_dot.set_ylim((0, y_length))

	ax_dot.set_xlabel(x_labels)
	ax_dot.set_ylabel(y_labels)

	ax_coverage.stackplot(filled[:,0], filled[:,1], color=bar_color, linewidth=.01)
	#ax_coverage.bar(bin_pos, bin_counts, width=width_of_bars, color=bar_color)

	(coverage_min, coverage_max) = ax_coverage.get_ylim()
	t = [int(0), int(math.ceil(coverage_max))]
	ax_coverage.set_yticks(t)
	ax_coverage.set_yticklabels(t)

	ax_coverage.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
	
	ax_dot.yaxis.set_major_locator(plt.MaxNLocator(4))
	ax_dot.xaxis.set_major_locator(plt.MaxNLocator(4))

	ax_coverage.spines['right'].set_visible(False)
	ax_coverage.spines['top'].set_visible(False)

	ax_coverage = thicker_spines(ax_coverage, False)
	ax_dot = thicker_spines(ax_dot, True)

	ax_coverage.spines['bottom'].set_position('zero')

	return ax_coverage, ax_dot


def plot(filled, dot_data, f_csv, bin_pos, bin_counts, x_length, y_length, write_file):
	'''Plots coverage track on top of dot plot. Writes .pdf

	Inputs: filled- np array data for coverage track
			dot_data- np array data for dot plot
			f_csv- string of .csv file for alignment plot
			bin_pos- list of locations for histogram bars on x axis
			bin_counts- list of counts per bin
			x_length- int length of dot plot x axis
			y_length- int length of dot plot y axis
			write_file- string name of file to write plots to
	Output: fig- matplotlib fig of coverage track, dot plot, alignment plot
	'''

	fig, axs = plt.subplots(2, 2, sharex=False, sharey=False, figsize=(15,10), 
		gridspec_kw={'height_ratios': [1, 9], 'width_ratios':[9, 6]})

	axs[0][0], axs[1][0] = plot_dot_plot(axs[0][0], axs[1][0], filled, dot_data, x_length, y_length,
		bin_pos, bin_counts)
	axs[0][1], axs[1][1] = plot_alignment(axs[0][1], axs[1][1], f_csv)

	fig.subplots_adjust(hspace=0.04, wspace = 0.05)	
	plt.savefig(write_file+'.pdf',
		dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches=None, pad_inches=0.01,
        frameon=False, metadata=None)
	return fig


def main(f_txt, dot_data, f_csv, x_length, y_length, write_file):
	'''Reads .txt + .paf, plots coverage track + dot plot, writes plots to file

	Inputs: f_txt- .txt file string to read for histogram
			dot_data- np array data for dot plot
			f_csv- for alignment bar graph
			x_length- int length of dot plot x axis
			y_length- int length of dot plot y axis
			write_file- string name of file to write plots to
	'''
	filled = fill_blanks(f_txt, x_length)
	bin_pos, bin_counts = make_hist(filled, x_length)
	fig = plot(filled, dot_data, f_csv, bin_pos, bin_counts, x_length, y_length, write_file)


if __name__ == "__main__":
	'''
	Example usage:
	python pared_dot_coverage.py ".txt file for coverage" ".paf file for dot plot" ".csv for alignment" 
		"y axis length" "destination for plots"
	python pared_dot_coverage.py depth_per_base.txt contig_nCoV-2019.paf stats.csv 29903 covid_plots 
	'''
	parser = argparse.ArgumentParser(description='Plot histogram.')
	parser.add_argument('read_txt_file', metavar='rf_txt', type=str,
                     help='.txt file to read')
	parser.add_argument('read_paf_file', metavar='rf_paf', type=str,
                     help='.paf file to read')
	parser.add_argument('read_csv_file', metavar='rf_csv', type=str,
                     help='.csv file to read')

	parser.add_argument('y_axis_length', metavar='y_length', type=int, 
                     help='y axis length')

	parser.add_argument('write_file', metavar='wf', type=str,
                     help='png/svg files to write')

	args = parser.parse_args()

	col_names = [str(i) for i in range(18)]
	dot_data = pd.read_csv(args.read_paf_file, names=col_names, delimiter='	') 
	x_length = np.sum(dot_data['6'].values[0])
	main(args.read_txt_file, dot_data, args.read_csv_file, 
		x_length, args.y_axis_length, args.write_file)
	

