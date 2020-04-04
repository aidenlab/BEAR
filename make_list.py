import csv
import argparse
import numpy as np

def main(read_file, write_file):
	xx = []
	yy = []
	with open(read_file) as f:
		rd = csv.reader(f, delimiter="\t", quotechar='"')
		for row in rd:
			xx.append(row[2])
			xx.append(row[7])
			yy.append(row[3])
			yy.append(row[8])
	plt_mat = np.zeros((len(xx), 2))
	plt_mat[:,0] = xx
	plt_mat[:,1] = yy
	np.savetxt(write_file, plt_mat)


if __name__ == "__main__":
	#Example usage
	#python make_list.py /Users/ABlackburn/Downloads/contig_nCoV-2019.paf /Users/ABlackburn/Documents/list_destination.txt
	parser = argparse.ArgumentParser(description='Read .paf file, make list of 3rd and 8th fields')
	parser.add_argument('read_file', metavar='rf', type=str, 
                     help='.paf file to read')
	parser.add_argument('write_file', metavar='wf', type=str, 
                     help='file to dump to')

	args = parser.parse_args()
	main(args.read_file, args.write_file)
    