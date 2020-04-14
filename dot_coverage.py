import numpy as np
import matplotlib, math, argparse
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import pandas as pd
matplotlib.rcParams.update({'font.size': 18})

 
bin_size = 1
num_bins = 30000.0/bin_size
width_of_bars = 100000.0/num_bins
x_labels = "SARS-CoV-2 RefSeq Assembly"
y_labels = "de novo SARS-CoV-2 Assembly"
col = "#80D1D9" ###blue


name_key = pd.DataFrame({ 
	'id': ['Bat_Coronavirus_HKU4_1',
	'Bat_Coronavirus_HKU5_1',
	'Bat_Coronavirus_HKU9_1', 
	'bat_Coronavirus_Parker', 
	'Bat_Hp_Betacoronavirus_Zhejiang2013',
	'Betacoronavirus_England', 
	'Betacoronavirus_Erinaceus', 
	'Betacoronavirus_HKU24_Strain_HKU24_R05005I',
	'Bovine_Coronavirus', 
	'Human_Coronavirus_HKU1', 
	'Human_Coronavirus_OC43_Strain_ATCC_VR_759',
	'Middle_East_Respiratory_Syndrome_Coronavirus', 
	'Mouse_Hepatitis_Virus_Strain_MHV_A59_C12_mutant',
	'Rabbit_Coronavirus_HKU14', 
	'Rousettus_Bat_Coronavirus_Isolate_GCCDC1_356', 
	'Severe_Acute_Respiratory_Syndrome_Coronavirus',
	'Wuhan_seafood_market_pneumonia_virus_isolate_Wuhan_Hu_1'], 
	'new': ['Bat coronavirus HKU4-1', 
	'Bat coronavirus HKU5-1', 
	'Bat coronavirus HKU9-1', 
	'Rat coronavirus Parker', 
	'Bat Hp-betacoronavirus', 
	'Betacoronavirus England 1', 
	'Betacoronavirus isolate Erinaceus CoV, 2012', 
	'Betacoronavirus HKU24 strain HKU24-R05005I', 
	'Bovine coronavirus', 
	'Human coronavirus HKU1', 
	'Human coronavirus OC43 strain ATCC VR-759', 
	'Middle East respiratory syndrome (MERS) coronavirus',
	'Mouse hepatitis virus strain MHV-A59 C12 mutant',
	'Rabbit coronavirus HKU14',
	'Rousettus bat coronavirus isolate GCCDC1 356',
	'Severe acute respiratory syndrome (SARS) coronavirus',
	'SARS-CoV-2']
	})
 

def fill_blanks(f_txt, x_length):
	'''Fills in histogram blanks.

	Inputs: f_txt- .txt file string with position read counts
			x_length- int length of x axis
	Outputs: filled- np array with 0s filled in
	'''
	if x_length == 0:
		return False
	raw = pd.read_csv(f_txt, sep="	", names=['first', 'second', 'third'], header=None)
	if len(raw.index) == 0:
		return False
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
	if data is False:
		return [], []
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


def strip_ticks(ax):
	ax.set_xticks([])
	ax.set_yticks([])
	ax = thicker_spines(ax, True)
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
	#Diagnostic axis style
	ax_diagnostic.axis('off')

	#Alignment axis style
	ax_align.set_xlabel('% of Reads that Align to Betacoronaviruses', labelpad=10)
	ax_align = thicker_spines(ax_align, True)

	#Read/format/sort alignment data
	try:
		align_data = pd.read_csv(f_csv)
	except:
		return strip_ticks(ax_diagnostic), strip_ticks(ax_align)
	
	name_key['id'] = name_key['id'].str.lower()
	align_data['id'] = align_data['label'].str.lower()
	align_data = align_data.merge(name_key, how="left")
	align_data = align_data.sort_values('percentage')
	align_data.loc[align_data['id'].str.contains('wuhan'), 'new'] = 'SARS-CoV-2'
	
	#Plot diagnostic symbol
	ax_diagnostic.plot([.27, .33], [.5, .5], linewidth=20, color=col)
	if align_data[align_data['new']=='SARS-CoV-2']['percentage'].values[0] >= 75:
		ax_diagnostic.plot([.3, .3], [.1, .9], linewidth=20, color=col)
		ax_diagnostic.text(.4, .83, "Test Result: Positive")
		ax_diagnostic.text(.4, .55, "SARS-CoV-2 was detected", color='black')
		ax_diagnostic.text(.4, .29, "in the sample", color='black')
	else:
		ax_diagnostic.text(.4, .83, "Test Result: Negative")
		ax_diagnostic.text(.4, .55, "SARS-COV2 was not detected", color='black')
		ax_diagnostic.text(.4, .29, "in the the sample", color='black')

	#Plot alignment bar chart
	def new_labels(new_name, perc):
		return "{}: {}%".format(new_name, perc)

	align_data['new'] = align_data.apply(lambda row : new_labels(row['new'], 
							row['percentage']), axis = 1)
	align_data.plot(kind='barh', x='new', y='percentage', ax=ax_align, color=col, legend=False)
	
	#Diagnostic axis limits
	ax_diagnostic.set_xlim(0,1)
	ax_diagnostic.set_ylim(0,1)
	
	#Alignment axis tick marks and limits
	ax_align.set_xlim(0,100)
	ax_align.tick_params(axis="y",direction="in", left="true", pad=-5)
	plt.setp(ax_align.get_yticklabels(), ha="left")
	ax_align.set_ylabel('')
	
	return ax_diagnostic, ax_align


def plot_coverage(ax_coverage, filled, bin_pos, bin_counts, x_length):
	'''Create dot plot and coverage track.

	Inputs: ax_coverage- coverage track axis
			filled- np array data for coverage track
			x_length- int length of dot plot x axis
			bin_pos- list of locations for histogram bars on x axis
			bin_counts- list of counts per bin
	Outputs: ax_coverage- coverage track axis
	'''
	#Coverage track axis style
	ax_coverage.set_yscale('log')
	ax_coverage.set_ylabel("Coverage", labelpad=10)
	ax_coverage.spines['right'].set_visible(False)
	ax_coverage.spines['top'].set_visible(False)
	ax_coverage.spines['bottom'].set_position(('axes', 0))#'zero')
	ax_coverage = thicker_spines(ax_coverage, False)
	
	if filled is False:
		return strip_ticks(ax_coverage)
	
	#Make coverage track
	if bin_size == 1:
		ax_coverage.stackplot(filled[:,0], filled[:,1], color=col, linewidth=.01)
	else:
		ax_coverage.stackplot(bin_pos, bin_counts, color=col, linewidth=.01)

	#Coverage track tick marks and limits
	ax_coverage.set_xlim(0, x_length)
	(coverage_min, coverage_max) = ax_coverage.get_ylim()
	t_coverage = [0, int(math.ceil(coverage_max/ 100.0)) * 100]
	t_log_coverage = [1, int(math.ceil(coverage_max/ 100.0)) * 100]
	#ax_coverage.set_yticks(t_coverage)
	ax_coverage.set_yticks(t_log_coverage)
	ax_coverage.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
	#ax_coverage.spines['left'].set_bounds(t_coverage[0], t_coverage[1])
	ax_coverage.spines['left'].set_bounds(t_log_coverage[0], t_log_coverage[1])
	

	return ax_coverage


def plot_dot_plot(ax_dot, dot_data, x_length, y_length):
	'''Create dot plot and coverage track.

	Inputs: ax_dot- dot plot axis
			dot_data- np array data for dot plot
			x_length- int length of dot plot x axis
			y_length- int length of dot plot y axis
	Outputs: ax_dot- dot plot axis
	'''
	#Dot plot axis style
	ax_dot = thicker_spines(ax_dot, True)
	ax_dot.set_xlabel(x_labels, labelpad=10)
	ax_dot.set_ylabel(y_labels, labelpad=10)

	if len(dot_data.index) == 0:
		return strip_ticks(ax_dot)
	
	#Make dot plot
	if y_length < 100000:
		dot_data = dot_data.sort_values('7')
	dot_data['cum_offset'] = dot_data['1'].cumsum()
	dot_data['upto_offset'] = dot_data['cum_offset'].values - dot_data['1'].values

	def dot_plot(x1, x2, y1, y2, segment_offset,):
		ax_dot.plot([x1, x2], [y1+segment_offset, y2+segment_offset], linewidth=3, color=col)

	dot_data.apply(lambda row : dot_plot(row['7'], row['8'],
                     row['2'], row['3'], row['upto_offset']), axis = 1)

	#Dot plot tick marks and limits
	ax_dot.set_xlim((0, x_length))
	ax_dot.set_ylim((0, y_length))
	
	def space_ticks(ax_length):
		tick_space = (ax_length+1)/3
		mult_10 = float(10**(math.floor(np.log10(tick_space))))
		return np.arange(0, ax_length+1, int(round(tick_space/ mult_10)) * mult_10)

	t_xx_loc = space_ticks(x_length)
	t_yy_loc = space_ticks(y_length)
	ax_dot.set_xticks(t_xx_loc)
	ax_dot.set_yticks(t_yy_loc)

	return ax_dot


def plot(filled, dot_data, f_csv, bin_pos, bin_counts, x_length, y_length, write_file):
	'''Plots coverage track, dot plot, alignment bar graph, diagnostic symbol. Writes .pdf

	Inputs: filled- np array data for coverage track
			dot_data- np array data for dot plot
			f_csv- string of .csv file for alignment plot
			bin_pos- list of locations for histogram bars on x axis
			bin_counts- list of counts per bin
			x_length- int length of dot plot x axis
			y_length- int length of dot plot y axis
			write_file- string name of file to write plots to
	'''
	fig, axs = plt.subplots(2, 2, sharex=False, sharey=False, figsize=(22, 12), 
		gridspec_kw={'height_ratios': [1, 11], 'width_ratios':[11, 11]})

	axs[0][0] = plot_coverage(axs[0][0], filled, bin_pos, bin_counts, x_length)
	axs[1][0] = plot_dot_plot(axs[1][0], dot_data, x_length, y_length)
	axs[0][1], axs[1][1] = plot_alignment(axs[0][1], axs[1][1], f_csv)

	fig.subplots_adjust(hspace=0.09, wspace = 0.05)	
	plt.savefig(write_file+'.pdf',
		dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches='tight', pad_inches=0.4,
        frameon=False, metadata=None)


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
	plot(filled, dot_data, f_csv, bin_pos, bin_counts, x_length, y_length, write_file)


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

	if len(dot_data.index) == 0:
		x_length = 0
	else:
		x_length = np.sum(dot_data['6'].values[0])
	
	y_length = args.y_axis_length

	main(args.read_txt_file, dot_data, args.read_csv_file, 
		x_length, y_length, args.write_file)
	

