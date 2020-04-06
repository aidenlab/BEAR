import numpy as np
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import csv, math, argparse
import base64
from io import BytesIO
import pandas as pd
 

width_of_bars=500
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
	raw = pd.read_csv(f_txt, sep=" ", names=['first', 'second', 'third'], header=None)

	raw_np_array = np.zeros((len(raw),2))
	raw_np_array[:,0] = raw['second'].values
	raw_np_array[:,1] = raw['third'].values

	filled = np.zeros((max_length, 2))
	filled[:,0] = np.arange(1, max_length+1)
	filled[np.isin(filled[:,0], raw_np_array[:,0]),1] = raw_np_array[:,1]

	return filled


def make_hist(data, max_length, bin_size):
	'''Makes a histogram from position read counts

	Inputs: data- np array of filled in posiiton read counts
			max_length- int length of genome
			bin_size- int size of bin
	Outputs: bin_pos- list of locations for histogram bars on x axis
			bin_counts- list of counts per bin
	'''
	bins = np.arange(1, max_length+bin_size, bin_size)
	digitized = np.digitize(data[:,0], bins)
	bin_counts = [data[digitized == i,1].sum() for i in range(1, len(bins))]
	bin_pos = [(bins[i]+bins[i+1]-1)/2.0 for i in range(len(bins)-1)]
	
	return bin_pos, bin_counts


def plot(f_paf, max_length, bin_pos, bin_counts, write_file):
	'''Plots coverage track on top of dot plot. Writes a .png and .svg

	Inputs: f_paf- string of .paf file for dot plot
			max_length- int length of genome
			bin_pos- list of locations for histogram bars on x axis
			bin_counts- list of counts per bin
			write_file- string name of file to write plots to

	Output: fig- matplotlib fig of coverage track over dot plot
	'''

	fig, axs = plt.subplots(2, 1, sharex=True, sharey=False, figsize=(9,10), 
		gridspec_kw={'height_ratios': [1, 9]})

	offset = 0
	with open(f_paf) as f:
		rd = csv.reader(f, delimiter="\t", quotechar='"')
		for row in rd:
			col_fract = float(row[9])/float(row[10])
			col_idx = int(math.ceil(col_fract*len(line_palette)))-1
			col = line_palette[col_idx]
			axs[1].plot([float(row[7]), float(row[8])], [float(row[2])+offset, 
				float(row[3])+offset], linewidth=2, color=col)
			offset += float(row[1])
	
	axs[1].set_xlim((0, max_length))
	axs[1].set_ylim((0, offset))

	axs[1].set_xlabel(x_labels)
	axs[1].set_ylabel(y_labels)

	axs[0].bar(bin_pos, bin_counts, width=width_of_bars, color=bar_color)

	axs[0].yaxis.set_major_locator(plt.MaxNLocator(2))
	axs[1].yaxis.set_major_locator(plt.MaxNLocator(4))
	axs[1].xaxis.set_major_locator(plt.MaxNLocator(4))
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


def append_html(fig, html_file):
	'''Opens html file, appends plots.

	Inputs: fig- matplotlib fig of coverage track over dot plot
			html_file- string name of html file to append plot to
	'''
	tmpfile = BytesIO()
	fig.savefig(tmpfile, format='png')
	encoded = base64.b64encode(tmpfile.getvalue()).decode('utf-8')
	html = '<img src=\'data:image/png;base64,{}\'>'.format(encoded)

	with open(html_file,'a') as f:
		f.write(html)


def html_to_pdf(html_file, pdf_file):
	'''Converts html to pdf.

	Inputs: html_file- string name of html file to append plot to
			pdf_file- strang name of pdf file to convert html to
	'''
	import pdfkit
	pdfkit.from_file(html_file, pdf_file)
	

def main(f_txt, f_paf, max_length, bin_size, write_file, html_file, pdf_file, 
	convert_to_pdf):
	'''Reads .txt + .paf file, plots coverage track + dot plot, writes plots to file
	Opens html file, appends plots, converts to pdf if convert_to_pdf is True

	Inputs: f_txt- .txt file string to read for histogram
			f_paf- .paf file string to read for dot plot
			max_length- int length of viral genome
			bin_size- int bin size to use for histogram
			write_file- string name of file to write plots to
			html_file- string name of html file to append plot to
			pdf_file- strang name of pdf file to convert html to
	'''
	filled = fill_blanks(f_txt, max_length)
	bin_pos, bin_counts = make_hist(filled, max_length, bin_size)
	fig = plot(f_paf, max_length, bin_pos, bin_counts, write_file)
	
	append_html(fig, html_file)

	if convert_to_pdf == "convert":
		html_to_pdf(html_file, pdf_file)



if __name__ == "__main__":
	'''
	Example usage:
	python dot_coverage.py out.txt harder.paf /plots/covid_plots 500 29903 stats.html stats.pdf convert
	'''
	parser = argparse.ArgumentParser(description='Plot histogram.')
	parser.add_argument('read_txt_file', metavar='rf_txt', type=str,
                     help='.txt file to read')
	parser.add_argument('read_paf_file', metavar='rf_paf', type=str,
                     help='.paf file to read')
	
	parser.add_argument('write_file', metavar='wf', type=str,
                     help='png/svg files to write')
	parser.add_argument('bin_size', metavar='b', type=int,
                     help='Histogram bin size')
	parser.add_argument('max_length', metavar='m', type=int, 
                     help='Genome length')

	parser.add_argument('read_html_file', metavar='rf_html', type=str,
                     help='.html file to add plot to')
	parser.add_argument('write_pdf_file', metavar='wf_pdf', type=str,
                     help='.pdf file to write html to')

	parser.add_argument('convert_to_pdf', metavar='pdf_bool', type=str, 
                     help='do we convert to pdf')
    
	args = parser.parse_args()
	main(args.read_txt_file, args.read_paf_file, args.max_length, args.bin_size,
		args.write_file, args.read_html_file, args.write_pdf_file, args.convert_to_pdf)
	

