import pandas as pd
import argparse
import os 

def get_stats(save_file):
	file_name = os.path.basename(os.getcwd())

	data	=	str(lib_name), \
				str(breadth_of_coverage.iloc[0,1]),	\
				str(samtools_coverage_mapq30.iloc[0,3]), 	\
				str(samtools_coverage_mapq30.iloc[0,6]), 	\
				str(samtools_coverage_mapq30.iloc[1,3]), 	\
				str(samtools_coverage_mapq30.iloc[1,6]), 	\
				str(samtools_coverage_mapq30.iloc[0,3]), 	\
				str(samtools_coverage_mapq30.iloc[0,6]), 	\
				str(samtools_coverage_mapq30.iloc[1,3]), 	\
				str(samtools_coverage_mapq30.iloc[1,6]) 

	with open(save_file, 'w') as file:
		for num in data:
			file.write(num + ',')
	file.close()

	

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description=
		'Compile POLAR-BEAR Results')
	parser.add_argument('top_dir_name', type=str,
		help='File name needed')
	parser.add_argument('path_to_map30_tsv_file', type=str,
		help='PATH to csv files needed')
	parser.add_argument('path_to_map0_tsv_file', type=str,
		help='PATH to csv files needed')
	parser.add_argument('path_to_csv_file', type=str,
		help='PATH to tsv files needed')
	parser.add_argument('path_to_save_file', type=str,
		help='PATH to save output')
	args = parser.parse_args()

	lib_name = args.top_dir_name
	samtools_coverage_mapq30 = pd.read_csv(args.path_to_map30_tsv_file, sep='\t')
	samtools_coverage_mapq0 = pd.read_csv(args.path_to_map0_tsv_file, sep='\t')
	breadth_of_coverage = pd.read_csv(args.path_to_csv_file)
	get_stats(args.path_to_save_file)