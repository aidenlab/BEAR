import pandas as pd
import argparse

def format_number_w_commads(number):
    return f"{int((number)):,}"

def return_breadth_of_coverage(viral_depth_per_base):
    depth_per_base = pd.read_csv(viral_depth_per_base,
                                 sep='\t', names=["RF", "B", "CV"],
                                 dtype={'RF': 'object', 'B': 'int64', 'CV': 'int64'})
    total_num_of_bases = depth_per_base.loc[depth_per_base.CV >= 0, 'CV'].count()
    base_with_at_lest_5x_cov = depth_per_base.loc[depth_per_base.CV >= 5, 'CV'].count()
    breadth_of_coverage = str(round((100 * base_with_at_lest_5x_cov / total_num_of_bases), 1))
    return breadth_of_coverage

def return_number_of_mapped_reads(number_of_mapped_reads, number_of_reads):
    mapped_reads_percentage = round( ((number_of_mapped_reads/2)/number_of_reads * 100), 2)
    mapped_reads_var = str(format_number_w_commads(number_of_mapped_reads)) + str(' (') + str(mapped_reads_percentage) + str('%)')
    return mapped_reads_var

def return_viral_read_dup(number_of_viral_duplicates, number_of_viral_reads):
    viral_dup_percentage = round((int(number_of_viral_duplicates) * 100 / int(number_of_viral_reads)), 1)
    viral_read_dup = str(format_number_w_commads(number_of_viral_duplicates)) + str(' (') + str(viral_dup_percentage) + str('%)')
    return viral_read_dup

def determine_diagnostic_result(breadth_of_coverage, control_count, read_length, reads):
    if float(breadth_of_coverage) >= 5:
        clinical_result_var = "positive"
    elif reads < 150000 or read_length < 90:
        clinical_result_var = "Inconclusive (read depth/length of sample not sufficient for analysis)"
    elif int(control_count) < 20000:
        clinical_result_var = "Inconclusive (control not adequately present in sample)"
    else:
        clinical_result_var = "negative"
    return clinical_result_var

def write_stats_file(path_to_qc_result, library_name, total_reads_var, mapped_reads_var, breadth_of_coverage_var, dup_reads_var, insert_var, insert_dev_var):
    qc_stats_to_write = '{0}\n{1}\n{2}\n{3}\n{4}\n{5}\n{6}'.format(
        "Library name: " + str(library_name),
        "Total reads: " + str(total_reads_var),
        "Aligned reads: " + str(mapped_reads_var),
        "BoC: " + str(breadth_of_coverage_var),
        "Duplicates: " + str(dup_reads_var),
        "Insert: " + str(insert_var),
        "Insert deviation: " + str(insert_dev_var)
    )

    with open(path_to_qc_result, 'w') as qc_file:
        qc_file.write(qc_stats_to_write)
    qc_file.close()

def write_result_file(path_to_clinical_result, library_name, clinical_result_var):
    result_header = "sample,result\n"
    clinical_result_to_write = str(library_name) + "," + str(clinical_result_var)

    with open(path_to_clinical_result, 'w') as clinical_file:
        clinical_file.write(result_header)
        clinical_file.write(clinical_result_to_write)
    clinical_file.close()

def return_data_from_viral_align_stats(viral_align_stats_file_path):
    with open(viral_align_stats_file_path) as viral_stat_file:
        for line in viral_stat_file:
            if line.startswith('SN'):
                if 'reads mapped and paired:' in line:
                    viral_reads = float(line.split(':')[1].split('#')[0])
                if 'reads duplicated:' in line:
                    viral_duplicates = float(line.split(':')[1].split('#')[0])
                if 'insert size average:' in line:
                    insert_var = float(line.split(':')[1])
                if 'insert size standard deviation:' in line:
                    insert_dev_var = float(line.split(':')[1])

    return viral_reads, viral_duplicates, insert_var, insert_dev_var


def return_data_from_all_align_stats(all_align_stats_file_path):
    reads = []
    mapped_reads_count = []
    with open(all_align_stats_file_path) as all_stat_file:
        for line in all_stat_file:
            if not "SN" in line:
                if 'mapped (' in line:
                    mapped_reads_count = float(line.split()[0])
                if 'paired in sequencing' in line:
                    reads = int(line.split()[0])/2

    return mapped_reads_count, reads

def return_data_from_amplicon_counts(amplicon_control_counts_path):
    with open(amplicon_control_counts_path) as amplicon_control_counts_file:
        lines = amplicon_control_counts_file.readlines()
        control_count = (lines[1].split()[1])
    return control_count

def compile_results_into_file(library_name, top_directory):
    all_align_stats_file = str(top_directory) + "/aligned/all_alignment_stats.txt"
    viral_align_stats_file = str(top_directory) + "/aligned/viral_alignment_stats.txt"
    viral_depth_per_base_file = str(top_directory) + "/aligned/viral_depth_per_base.txt"
    path_to_qc_result = str(top_directory) + "/aligned/qc_stats.txt"
    path_to_clinical_result = str(top_directory) + "/final/result.csv"
    amplicon_control_counts = str(top_directory) + "/aligned/ampliconCoverage.txt"

    mapped_reads_count, reads = return_data_from_all_align_stats(all_align_stats_file)
    viral_reads, viral_duplicates, insert_var, insert_dev_var = return_data_from_viral_align_stats(viral_align_stats_file)
    control_count = return_data_from_amplicon_counts(amplicon_control_counts)

    total_reads_var = format_number_w_commads(reads)
    mapped_reads_var = return_number_of_mapped_reads(mapped_reads_count, reads)
    dup_reads_var = return_viral_read_dup(viral_duplicates, viral_reads)
    breadth_of_coverage_var = return_breadth_of_coverage(viral_depth_per_base_file)
    diagnostic_result_var = determine_diagnostic_result(breadth_of_coverage_var, control_count, reads)

    write_stats_file(path_to_qc_result,
                     library_name,
                     total_reads_var,
                     mapped_reads_var,
                     breadth_of_coverage_var,
                     dup_reads_var,
                     insert_var,
                     insert_dev_var)

    write_result_file(path_to_clinical_result,
                      library_name,
                      diagnostic_result_var)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compile POLAR-BEAR results and QC data')
    parser.add_argument('library_name', type=str)
    parser.add_argument('top_directory', type=str)
    arg = parser.parse_args()

    compile_results_into_file(arg.library_name, arg.top_directory)

