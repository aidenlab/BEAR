import numpy as np
import matplotlib.pyplot as plt
import argparse


def fill_blanks(max_length):
	max_length = int(np.max(raw[:,0]))

	filled = np.zeros((max_length, 2))
	filled[:,0] = np.arange(1, max_length+1)
	filled[np.isin(filled[:,0], raw[:,0]),1] = raw[:,1]

	return filled


def make_hist(data, max_length, bin_size, write_file):
	max_length = int(np.max(raw[:,0]))

	bins = np.arange(1, max_length+1, bin_size)
	digitized = np.digitize(data[:,0], bins)
	bin_counts = [data[digitized == i,1].sum() for i in range(1, len(bins))]
	bin_pos = [(bins[i]+bins[i+1]-1)/2.0 for i in range(len(bins)-1)]
	
	plt.bar(bin_pos, bin_counts, width=width_of_bars, color=bar_color)
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
	plt.show()
	

def main(raw, max_length, bin_size, write_file):
	filled = fill_blanks(max_length)
	make_hist(filled, max_length, bin_size, write_file)


if __name__ == "__main__":
	#Example usage:
	#python coverage_track.py /Users/ABlackburn/Downloads/out.txt /Users/ABlackburn/Documents/covid_hist 500 30000
	width_of_bars=400
	bar_color = "#3498DB"
	
	parser = argparse.ArgumentParser(description='Plot histogram.')
	parser.add_argument('read_file', metavar='rf', type=str, 
                     help='.txt file to read')
	parser.add_argument('write_file', metavar='wf', type=str, 
                     help='png/svg files to write')
	parser.add_argument('bin_size', metavar='b', type=int, 
                     help='Histogram bin size')
	parser.add_argument('max_length', metavar='m', type=int,
                     help='Genome length')
    
	args = parser.parse_args()
	raw = np.loadtxt(args.read_file)

	max_length = int(np.max(raw[:,0]))
	main(raw, args.max_length, args.bin_size, args.write_file)
	

