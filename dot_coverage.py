import numpy as np
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import csv, math, argparse
import base64
from io import BytesIO
import pandas as pd
import itertools
from matplotlib import rcParams
import matplotlib.ticker as mticker
matplotlib.rcParams.update({'font.size': 10})

 
bin_size = 1
num_bins = 30000.0/bin_size
width_of_bars = 100000.0/num_bins
x_labels = "SARS-CoV-2 RefSeq Assembly"
y_labels = "de novo SARS-CoV-2 Assembly"
# col = "#D98880d" ###pink
col = "#80D1D9" ###blue
 

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
	#name_key = pd.read_csv('name_changes.csv')
	#align_data = align_data.merge(name_key, how="left")
	align_data['new'] = align_data['label'].str.replace('_',' ')
	align_data['new'] = align_data.new.str.title()


	def new_labels(new_name, perc):
		return "{}: {}%".format(new_name, perc)

	#align_data['new'] = align_data['new'] align_data['percentage']
	align_data['new'] = align_data.apply(lambda row : new_labels(row['new'], 
							row['percentage']), axis = 1)

	align_data['id'] = align_data['label'].str.lower()
	align_data = align_data.sort_values('percentage')

	align_data.plot(kind='barh', x='new', y='percentage', ax=ax_align, color=col, legend=False)
	
	#for i, v in enumerate(y):
    #	align_data.text(v + 3, i + .25, str(v))

	ax_align.set_xlabel('% of Reads that Align to Betacoronaviruses')
	ax_align.set_ylabel('')
	ax_align.set_xlim(0,100)
	ax_align.tick_params(axis="y",direction="in", left="true", pad=-5)
	plt.setp(ax_align.get_yticklabels(), ha="left")

	covid_val = align_data[align_data['id'].str.contains("wuhan")]['percentage'].values[0]
	
	ax_diagnostic.plot([.245, .355], [.5, .5], linewidth=20, color=col)
	if covid_val >= 75:
		ax_diagnostic.plot([.3, .3], [.1, .9], linewidth=20, color=col)
		ax_diagnostic.set_title("Test Result: Positive")
		ax_diagnostic.text(.4, .55, "SARS-CoV-2 was detected", color='black')
		ax_diagnostic.text(.4, .35, "in the sample", color='black')
	else:
		ax_diagnostic.set_title("Test Result: Negative")
		ax_diagnostic.text(.4, .55, "SARS-CoV-2 was not detected", color='black')
		ax_diagnostic.text(.4, .35, "in the sample", color='black')
	
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

	def dot_plot(x1, x2, y1, y2, segment_offset,):
		ax_dot.plot([x1, x2], [y1+segment_offset, y2+segment_offset], linewidth=3, color=col)

	dot_data.apply(lambda row : dot_plot(row['7'], row['8'],
                     row['2'], row['3'], row['upto_offset']), axis = 1)

	ax_coverage.set_xlim(0, x_length)
	ax_dot.set_xlim((0, x_length))
	ax_dot.set_ylim((0, y_length))

	ax_coverage.set_ylabel("Coverage")
	ax_dot.set_xlabel(x_labels)
	ax_dot.set_ylabel(y_labels)

	if bin_size == 1:
		ax_coverage.stackplot(filled[:,0], filled[:,1], color=col, linewidth=.01)
	else:
		ax_coverage.stackplot(filled[:,0], filled[:,1], color=col, linewidth=.01)

	(coverage_min, coverage_max) = ax_coverage.get_ylim()
	t_coverage = [int(0), int(math.ceil(coverage_max/ 100.0)) * 100]
	ax_coverage.set_yticks(t_coverage)
	#ax_coverage.set_yticklabels(t_coverage)
	ax_coverage.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

	def log_tick_formatter(val, pos=None):
		return "10$^{"+str(int(math.log10(val)))+"}$"

	t_loc = np.arange(0, max(x_length, y_length)+1, 10000)
	t_names = [0, 1, 2, 3, 4]
	ax_dot.set_xticks(t_loc)
	ax_dot.set_yticks(t_loc)

	ax_coverage.spines['right'].set_visible(False)
	ax_coverage.spines['top'].set_visible(False)
	ax_coverage.spines['left'].set_bounds(t_coverage[0], t_coverage[1])

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

	fig, axs = plt.subplots(2, 2, sharex=False, sharey=False, figsize=(17.75, 12), 
		gridspec_kw={'height_ratios': [1, 11], 'width_ratios':[11, 6.75]})

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
	

