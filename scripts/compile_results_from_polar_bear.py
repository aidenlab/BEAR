import pandas as pd
import argparse


def compile_results_into_file(library_name, top_directory):
    all_align_stats_file = str(top_directory) + "/aligned/all_alignment_stats.txt"
    viral_align_stats_file = str(top_directory) + "/aligned/viral_alignment_stats.txt"
    viral_depth_per_base = str(top_directory) + "/aligned/viral_depth_per_base.txt"
    path_to_qc_result = str(top_directory) + "/aligned/qc_stats.txt"
    path_to_clinical_result = str(top_directory) + "/final/result.csv"

    with open(viral_align_stats_file) as viral_stat_file:
        for line in viral_stat_file:
            if line.startswith('SN'):
                if 'reads mapped:' in line:
                    mapped = int(line.split()[-1])
                if 'reads duplicated:' in line:
                    dups = int(line.split()[-8])
                if 'insert size average:' in line:
                    insert_var = float(line.split()[-1])
                if 'insert size standard deviation:' in line:
                    insert_dev_var = float(line.split()[-1])

    with open(all_align_stats_file) as all_stat_file:
        for line in all_stat_file:
            if 'paired in sequencing' in line:
                reads = int(line.split()[0])

    total_reads_var = f"{int((reads / 2)):,}"
    mapped_reads_var = str(round((mapped * 100 / reads), 1)) + '% (' + str(f"{int(mapped):,}") + ')'
    dup_reads_var = str(round((dups * 100 / mapped), 1)) + '% (' + str(f"{int(dups):,}") + ')'

    depth_per_base = pd.read_csv(viral_depth_per_base, sep='\t', names=["RF", "B", "CV"],
                                 dtype={'RF': 'object', 'B': 'int64', 'CV': 'int64'})
    total_num_of_bases = depth_per_base.loc[depth_per_base.CV >= 0, 'CV'].count()
    base_with_at_lest_5x_cov = depth_per_base.loc[depth_per_base.CV >= 5, 'CV'].count()
    breadth_of_coverage_var = str(round((100 * base_with_at_lest_5x_cov / total_num_of_bases), 1))

    if float(breadth_of_coverage_var) >= 5:
        clinical_result_var = "Positive"
    else:
        clinical_result_var = "Negative"

    qc_stats_to_write = '{0}\n{1}\n{2}\n{3}\n{4}\n{5}\n{6}'.format(
        "Library name ," + str(library_name),
        "Total reads ," + str(total_reads_var),
        "Aligned reads ," + str(mapped_reads_var),
        "BoC ," + str(breadth_of_coverage_var),
        "Duplicates ," + str(dup_reads_var),
        "Insert ," + str(insert_var),
        "Insert deviation ," + str(insert_dev_var)
    )

    clinical_result_to_write = str(library_name) + "," + str(clinical_result_var)

    with open(path_to_qc_result, 'w') as qc_file:
        qc_file.write(qc_stats_to_write)
    qc_file.close()

    with open(path_to_clinical_result, 'w') as clinical_file:
        clinical_file.write(clinical_result_to_write)
    clinical_file.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compile POLAR-BEAR results and QC data')
    parser.add_argument('library_name', type=str)
    parser.add_argument('top_directory', type=str)
    arg = parser.parse_args()

    compile_results_into_file(arg.library_name, arg.top_directory)
