import pandas as pd
import argparse
import os.path
import numpy

polar_control_index_ids = ['A:H6'
                           'A:H12'
                           'B:H6'
                           'B:H12'
                           'C:H6'
                           'C:H12'
                           'D:H6'
                           'D:H12']


def compile_results_into_file(libname_var, all_align_stats_file, viral_align_stats_file,
                              viral_depth_per_base, path_to_clinical_result, path_to_qc_result,
                              paf_file):
    if os.path.isfile(paf_file):
        lengths = []
        with open(paf_file, 'r') as data:
            for line in data:
                ls = line.split()
                lengths.append(int(ls[1]))

        sorted_lengths = sorted(lengths, reverse=True)
        csum = numpy.cumsum(sorted_lengths)
        n2 = int(sum(lengths) / 2)
        csumn2 = min(csum[csum >= n2])
        ind = numpy.where(csum == csumn2)

        n50_var = sorted_lengths[int(ind[0])]
        min_length_var = sorted_lengths[0]
        max_length_var = sorted_lengths[1]
        sum_length_var = sum(sorted_lengths)
        number_of_contigs_var = len(sorted_lengths)
    else:
        n50_var = 'NA'
        min_length_var = 'NA'
        max_length_var = 'NA'
        sum_length_var = 'NA'
        number_of_contigs_var = 'NA'

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
            if 'Index ID:' in line:
                index_id = str(line.split()[2])

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

    if index_id in polar_control_index_ids:
        sample_type_var = "Negative Control"
    else:
        sample_type_var = "Clinical Sample"

    qc_stats_to_write = str(libname_var) + '\n' + \
                        str(total_reads_var) + '\n' + \
                        str(mapped_reads_var) + '\n' + \
                        str(breadth_of_coverage_var) + '\n' + \
                        str(dup_reads_var) + '\n' + \
                        str(insert_var) + '\n' + \
                        str(insert_dev_var) + '\n' + \
                        str(number_of_contigs_var) + '\n' + \
                        str(sum_length_var) + '\n' + \
                        str(min_length_var) + '\n' + \
                        str(max_length_var) + '\n' + \
                        str(n50_var)

    clinical_result_to_write = str(libname_var) + "," + str(clinical_result_var) + "," + str(sample_type_var)

    with open(path_to_qc_result, 'w') as qc_file:
        qc_file.write(qc_stats_to_write)
    qc_file.close()

    with open(path_to_clinical_result, 'w') as clinical_file:
        clinical_file.write(clinical_result_to_write)
    clinical_file.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compile POLAR-BEAR results and QC data to file')
    parser.add_argument('libname_var', type=str)
    parser.add_argument('all_align_stats_file', type=str)
    parser.add_argument('viral_align_stats_file', type=str)
    parser.add_argument('viral_depth_per_base', type=str)
    parser.add_argument('path_to_clinical_result', type=str)
    parser.add_argument('path_to_qc_result', type=str)
    parser.add_argument('paf_file', type=str)
    arg = parser.parse_args()

    compile_results_into_file(arg.libname_var, arg.all_align_stats_file, arg.viral_align_stats_file,
                              arg.viral_depth_per_base, arg.path_to_clinical_result, arg.path_to_qc_result,
                              arg.paf_file)
