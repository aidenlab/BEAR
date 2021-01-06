import pandas as pd
import argparse
import os.path

polar_control_index_ids = [ 'A:H6'
                            'A:H12'
                            'B:H6'
                            'B:H12'
                            'C:H6'
                            'C:H12'
                            'D:H6'
                            'D:H12' ]


def stats(all_align_stats_file, viral_align_stats_file, viral_depth_per_base, path_to_qc_stats, path_to_result, megahit_log_file, libname):

    with open(viral_align_stats_file) as viral_stat_file:
        for line in viral_stat_file:
            if line.startswith('SN'):
                if 'reads mapped:' in line:
                    mapped = int(line.split()[-1])
                if 'reads duplicated:' in line:
                    dups = int(line.split()[-8])
                if 'insert size average:' in line:
                    insert = float(line.split()[-1])
                if 'insert size standard deviation:' in line:
                    insert_dev = float(line.split()[-1])

    with open(all_align_stats_file) as all_stat_file:
        for line in all_stat_file:
            if 'paired in sequencing' in line:
                reads = int(line.split()[0])
            if 'Index ID:' in line:
                index_id = str(line.split()[2])

    total_reads = f"{int((reads/2)):,}"
    mapped_reads = str(round((mapped*100/reads), 1)) + '% (' + str(f"{int(mapped):,}") + ')'
    dup_reads = str(round((dups*100/mapped), 1)) + '% (' + str(f"{int(dups):,}") + ')'

    depth_per_base = pd.read_csv(viral_depth_per_base, sep='\t', names=["RF", "B", "CV"], dtype={'RF': 'object', 'B': 'int64', 'CV': 'int64'})
    total_num_of_bases=depth_per_base.loc[depth_per_base.CV >= 0, 'CV'].count()
    base_with_at_lest_5x_cov = depth_per_base.loc[depth_per_base.CV >= 5, 'CV'].count()
    breadth_of_coverage = str(round((100* base_with_at_lest_5x_cov / total_num_of_bases), 1))
    arq_mapped = 'NA'

    if os.path.isfile(megahit_log_file):
        with open(megahit_log_file) as asm_log_file:
            second_to_last_line = asm_log_file.readlines()[-2]
            asm_stats = second_to_last_line.split('-')[3]
            contigs = asm_stats.split()[0]
            asm_size = asm_stats.split()[3]
            max = asm_stats.split()[9]
            n50 = asm_stats.split()[15]
    else:
        contigs = 'NA'
        asm_size = 'NA'
        max = 'NA'
        n50 = 'NA'

    QC_STAT_DATA_TO_WRITE =  str(libname) + '\n' + \
                             str(total_reads) + '\n' + \
                             str(mapped_reads) + '\n' + \
                             str(arq_mapped) + '\n' + \
                             str(breadth_of_coverage) + '\n' + \
                             str(dup_reads) + '\n' + \
                             str(insert) + '\n' + \
                             str(insert_dev) + '\n' + \
                             str(contigs) + '\n' + \
                             str(asm_size) + '\n' + \
                             str(max) + '\n' + \
                             str(n50)



    with open(path_to_qc_stats, 'w') as file:
        file.write(QC_STAT_DATA_TO_WRITE)
    file.close()

    if float(breadth_of_coverage) >= 5:
        result = "Positive"
    else:
        result = "Negative"

    if index_id in polar_control_index_ids:
        sample_type = "Negative Control"
    else:
        sample_type = "Clinical Sample"

    RESULT_DATA_TO_WRITE = str(libname) + "," + str(result) + ","  + str(sample_type) 

    with open(path_to_result, 'w') as file:
        file.write(RESULT_DATA_TO_WRITE)
    file.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=
        'Compile POLAR-BEAR results and QC data to file')
    parser.add_argument('all_align_stats_file', type=str,
        help='Path to Samtools stats file for all data')
    parser.add_argument('viral_align_stats_file', type=str,
        help='Path to Samtools stats file for viral data')
    parser.add_argument('viral_depth_per_base', type=str,
        help='Path to Samtools depth per base file')
    parser.add_argument('path_to_qc_stats', type=str,
        help='Path to result file')
    parser.add_argument('path_to_result', type=str,
        help='Path to result file')
    parser.add_argument('megahit_log_file', type=str,
        help='Path to MEGAHIT log file')
    parser.add_argument('libID', type=str,
        help='Library ID')
    args = parser.parse_args()

stats(args.all_align_stats_file, args.viral_align_stats_file, args.viral_depth_per_base, args.path_to_qc_stats, args.path_to_result, args.megahit_log_file, args.libID)


