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
	#raw = np.loadtxt(f_txt)
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


def make_new_tsv(f_paf, f_fa, new_tsv):
	'''Merges .paf and .fa and writes new .tsv

	Inputs: f_paf- string of .paf file for dot plot
			f_fa- string of .fa file for dot plot
			new_tsv- string name of file from merged .paf and .fa
	Outputs: int length of .paf row
	'''
	f1 = open(f_fa)
	f2 = open(f_paf)
	f3 = open(new_tsv, 'w+')

	rd_fa = csv.reader(f1, delimiter=" ", quotechar='"')
	rd_paf = csv.reader(f2, delimiter="\t", quotechar='"')
	new_file = csv.writer(f3, delimiter="\t", quotechar='"')

	def row_count(filename):
		with open(filename) as in_file:
			return sum(1 for _ in in_file)

	last_line_number = row_count(f_paf)
	row_paf = rd_paf.next()
	for row_fa in rd_fa:
		
		if row_fa[0][0] != '>':
			continue
		
		if row_fa[0][1:] == row_paf[0]:
			new_file.writerow(row_paf)
			if last_line_number != rd_paf.line_num:
				row_paf = rd_paf.next()
			
		else:
			fa_offset = int(row_fa[3][4:])
			new_row = [row_fa[0][1:], str(fa_offset)]
			new_file.writerow(new_row)

	f1.close()
	f2.close()
	f3.close()
	return len(row_paf)


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


def plot_dot_plot(ax_coverage, ax_dot, new_tsv, max_length, bin_pos, bin_counts, len_row):
	'''Create dot plot and coverage track.

	Inputs: ax_coverage- coverage track axis
			ax_dot- dot plot axis
			new_tsv- string name of file from merged .paf and .fa
			max_length- int length of genome
			bin_pos- list of locations for histogram bars on x axis
			bin_counts- list of counts per bin
			len_row- int length of .paf rows
	Outputs: ax_coverage- coverage track axis
			ax_dot- dot plot axis
	'''
	offset = 0
	col_names = [str(i) for i in range(len_row)]
	merged = pd.read_csv(new_tsv, names=col_names, delimiter='	')
	
	merged_sort = merged.sort_values('7')
	merged_sort['cum_offset'] = merged_sort['1'].cumsum()
	merged_sort['upto_offset'] = merged_sort['cum_offset'].values - merged_sort['1'].values

	plotted = merged_sort[merged_sort['2'].notnull()]
	not_plotted = merged_sort[merged_sort['2'].isnull()]
	
	def dot_plot(x1, x2, y1, y2, segment_offset, num, denom):
		col_fract = float(num)/float(denom)
		col_idx = int(math.ceil(col_fract*len(line_palette)))-1
		col = line_palette[col_idx]
		ax_dot.plot([x1, x2], [y1+segment_offset, y2+segment_offset], linewidth=2, color=col)

	plotted.apply(lambda row : dot_plot(row['7'], row['8'],
                     row['2'], row['3'], row['upto_offset'], row['9'], row['10']), axis = 1)

	offset += np.sum(merged_sort['1'].values)
	ax_coverage.set_xlim(0, max_length)
	ax_dot.set_xlim((0, max_length))
	ax_dot.set_ylim((0, offset))

	ax_dot.set_xlabel(x_labels)
	ax_dot.set_ylabel(y_labels)

	ax_coverage.bar(bin_pos, bin_counts, width=width_of_bars, color=bar_color)

	ax_coverage.yaxis.set_major_locator(plt.MaxNLocator(3))
	ax_coverage.xaxis.set_major_locator(plt.MaxNLocator(4))
	ax_dot.yaxis.set_major_locator(plt.MaxNLocator(4))
	ax_dot.xaxis.set_major_locator(plt.MaxNLocator(4))

	return ax_coverage, ax_dot


def plot(f_paf, f_fa, f_csv, max_length, bin_pos, bin_counts, write_file, new_tsv):
	'''Plots coverage track on top of dot plot. Writes .png, .svg, and .tsv

	Inputs: f_paf- string of .paf file for dot plot
			f_fa- string of .fa file for dot plot
			f_csv- string of .csv file for alignment plot
			max_length- int length of genome
			bin_pos- list of locations for histogram bars on x axis
			bin_counts- list of counts per bin
			write_file- string name of file to write plots to
			new_tsv- string name of file to write from .paf and .fa

	Output: fig- matplotlib fig of coverage track, dot plot, alignment plot
	'''

	fig, axs = plt.subplots(2, 2, sharex=False, sharey=False, figsize=(14.35,10), 
		gridspec_kw={'height_ratios': [1, 9], 'width_ratios':[9, 5.35]})

	len_row = make_new_tsv(f_paf, f_fa, new_tsv)
	axs[0][0], axs[1][0] = plot_dot_plot(axs[0][0], axs[1][0], new_tsv, max_length, 
		bin_pos, bin_counts, len_row)
	axs[0][1], axs[1][1] = plot_alignment(axs[0][1], axs[1][1], f_csv)

	fig.tight_layout()	
	plt.savefig(write_file+'.svg',
				dpi=None, facecolor='w', edgecolor='w',
		        orientation='portrait', papertype=None, format=None,
		        transparent=False, bbox_inches=None, pad_inches=0.01,
		        frameon=False, metadata=None)
	plt.savefig(write_file+'.png',
		dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches=None, pad_inches=0.01,
        frameon=False, metadata=None)
	plt.savefig(write_file+'.pdf',
		dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches=None, pad_inches=0.01,
        frameon=False, metadata=None)
	return fig


def append_html(fig, f_html):
	'''Opens html file, appends plots.

	Inputs: fig- matplotlib fig of coverage track over dot plot
			f_html- string name of html file to append plot to
	'''
	tmpfile = BytesIO()
	fig.savefig(tmpfile, format='png')
	encoded = base64.b64encode(tmpfile.getvalue()).decode('utf-8')
	html = '<img src=\'data:image/png;base64,{}\'>'.format(encoded)

	with open(f_html,'a') as f:
		f.write(html)


def html_to_pdf(f_html, new_pdf):
	'''Converts html to pdf.

	Inputs: f_html- string name of html file to append plot to
			new_pdf- strang name of pdf file to convert html to
	'''
	import pdfkit
	pdfkit.from_file(f_html, new_pdf)
	

def main(f_txt, f_paf, f_fa, f_csv, f_html, max_length, write_file, new_pdf, 
	new_tsv, convert_to_pdf,):
	'''Reads .txt + .paf +.fa file, plots coverage track + dot plot, writes plots to file
	Opens html file, appends plots, converts to pdf if convert_to_pdf is True

	Inputs: f_txt- .txt file string to read for histogram
			f_paf- .paf file string to read for dot plot
			f_fa- .fa file string to read for dot plot
			f_csv- for alignment bar graph
			f_html- string name of html file to append plot to
			max_length- int length of viral genome
			write_file- string name of file to write plots to
			new_pdf- string name of pdf file to convert html to
			new_tsv- string name of new tsv file from .paf and .fa
			convert_to_pdf- string, converts to pdf if == "convert"
	'''
	filled = fill_blanks(f_txt, max_length)
	bin_pos, bin_counts = make_hist(filled, max_length)
	make_new_tsv(f_paf, f_fa, new_tsv)
	fig = plot(f_paf, f_fa, f_csv, max_length, bin_pos, bin_counts, write_file, new_tsv)
	
	append_html(fig, f_html)

	if convert_to_pdf == "convert":
		html_to_pdf(f_html, new_pdf)



if __name__ == "__main__":
	'''
	Example usage:
	python dot_coverage.py ".txt file" ".paf file" ".fa file" "destination for plots" "bin size" 
		"genome length" "html stats" "new pdf" "pdf conversion?" "new tsv"
	python dot_coverage.py depth_per_base.txt contig_nCoV-2019.paf final.contigs.fa stats.csv stats.html 
		29903 covid_plots stats.pdf new.tsv False 
	'''
	parser = argparse.ArgumentParser(description='Plot histogram.')
	parser.add_argument('read_txt_file', metavar='rf_txt', type=str,
                     help='.txt file to read')
	parser.add_argument('read_paf_file', metavar='rf_paf', type=str,
                     help='.paf file to read')
	parser.add_argument('read_fa_file', metavar='rf_fa', type=str,
                     help='.fa file to read')
	parser.add_argument('read_csv_file', metavar='rf_csv', type=str,
                     help='.csv file to read')
	parser.add_argument('read_html_file', metavar='rf_html', type=str,
                     help='.html file to add plot to')

	parser.add_argument('max_length', metavar='m', type=int, 
                     help='Genome length')

	parser.add_argument('write_file', metavar='wf', type=str,
                     help='png/svg files to write')
	parser.add_argument('write_pdf_file', metavar='wf_pdf', type=str,
                     help='.pdf file to write html to')
	parser.add_argument('new_tsv', metavar='wf_tsv', type=str, 
                     help='new tsv from .paf and .fa')

	parser.add_argument('convert_to_pdf', metavar='pdf_bool', type=str, 
                     help='do we convert to pdf')

	
	args = parser.parse_args()
	main(args.read_txt_file, args.read_paf_file, args.read_fa_file, args.read_csv_file, 
		args.read_html_file, args.max_length, args.write_file, args.write_pdf_file, 
		args.new_tsv, args.convert_to_pdf)
	

