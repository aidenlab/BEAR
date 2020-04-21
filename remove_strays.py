import numpy as np
import pandas as pd
import argparse


max_primer_size = 50 #largest primer size included in viral enrichment PCR step

def remove(cov_pd, number_list, wf):
	new_list = np.zeros(number_list.shape)
	zero_idx = np.where(number_list == 0)[0]
	hole_length = np.ediff1d(zero_idx)
	hole_loc = np.where(hole_length > max_primer_size)[0] 


	data_start = np.add(zero_idx[hole_loc], 1)
	data_end = np.add(data_start, hole_length[hole_loc])
	
	for i in range(data_start.shape[0]):
		new_list[data_start[i]:data_end[i]] = number_list[data_start[i]:data_end[i]]

	cov_pd['third'] = new_list
	cov_pd['third'] = cov_pd['third'].astype(int)
	np.savetxt(wf, cov_pd.values, fmt='%s', delimiter='	')


if __name__ == "__main__":
	#python remove_strays.py depth_per_base.txt new.txt
	parser = argparse.ArgumentParser(description=
		'Remove non zero values with fewer than 3 consecutive occurences')
	
	parser.add_argument('read_txt_file', type=str,
                     help='.txt file to read')
	parser.add_argument('write_txt_file', type=str,
                     help='.txt file to read')
	args = parser.parse_args()

	cov_pd = pd.read_csv(args.read_txt_file, sep="	", names=['first', 'second', 'third'], header=None)
	number_list = cov_pd['third'].values

	#number_list = np.asarray([2, 2, 0, 0, 0, 1, 2, 3, 0, 0, 1, 3 ,4, 0, 0, 5])
	remove(cov_pd, number_list, args.write_txt_file)