#!/usr/local/bin/bash
## Polar BEAR FDA EUA pipeline

printHelpAndExit() {
    cat <<PRINTHELPANDEXIT
Program: POLAR-BEAR (POLAR Bioinformatics Evaluation of Assembly and Resequencing)
Version: 2.0
Contact: Neva Durand <neva@broadinstitute.org> & Per Adastra <adastra.aspera.per@gmail.com>
Usage:
        $0 [options]
Options:

   -d  Top level directory which must contain a subdirectory (fastq/) with fastq files
   -t  Number of threads for BWA alignment (Default: $THREADS)
   -h  Print this help and exit

PRINTHELPANDEXIT
exit
}

while getopts "d:t:sh" opt;
do
    case $opt in
        d) TOP_DIR=$OPTARG ;;
        t) THREADS=$OPTARG ;;
        h) printHelpAndExit ;;
    esac
done

### Determine if running in BaseSpace
if [ -f POLAR-BEAR/native.app.txt ];
then
  APP_MODE=1
  basespace_output_path_for_sample=$(awk 'NR==1 {print; exit}' POLAR-BEAR/native.app.txt)
  PYTHON=python3
else
  PYTHON=python
fi

### Threads
if [ -z $THREADS ];
then
    THREADS=16
fi

### PATHS
if [ -z $TOP_DIR ];
then
    TOP_DIR=$(pwd)
fi

# Make sure top directory does not include a '/'
END_OF_TOP_DIR=$(echo $TOP_DIR | rev | head -c 1)
if  [ "$END_OF_TOP_DIR" = "/" ]
then
  TOP_DIR=$(echo $TOP_DIR | rev | cut -c2- | rev)
fi

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd )"
LIB_NAME=$(echo $TOP_DIR | awk -F "/" '{print $NF}')

### Scripts
REMRECOMBO="${PIPELINE_DIR}/scripts/accugenomics/remRecombo.sh"
NT_TO_IS="${PIPELINE_DIR}/scripts/accugenomics/NT_IS_LOOKUP_TABLE_v0.4.2_seperate.txt"
AMPLICONS="${PIPELINE_DIR}/scripts/accugenomics/VarDict-amplicon.v2.1.bed"
NON_CROSS_REACT_REGIONS="${PIPELINE_DIR}/references/non_sars_cross_reactive_sars_cov_2_regions.bed"
COMPILE_RESULT="${PIPELINE_DIR}/scripts/compile_results_from_polar_bear.py"

### SSQC Settings
QCUTOFF=0
GOODBASECHANGE=1
PE=0
SEQSPLIT=1

### Misc vars
PATHOGEN_NAME="Sars-CoV-2"
REFERENCE="${PIPELINE_DIR}/references/sars_cov_2_accukit_ISv0.4.1/sars_cov_2_accukit_ISv0.4.1.fasta"

# Check for required installed software
echo "ʕ·ᴥ·ʔ : Checking dependencies..."


command -v bwa >/dev/null 2>&1 || { echo >&2 "ʕ·ᴥ·ʔ : BWA required but it's not installed!"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo >&2 "ʕ·ᴥ·ʔ : Samtools required but it's not installed!"; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo >&2 "ʕ·ᴥ·ʔ : Bedtools required but it's not installed!"; exit 1; }
command -v "${PYTHON}" >/dev/null 2>&1 || { echo >&2 "ʕ·ᴥ·ʔ : Python required but it's not installed!"; exit 1; }

# Check for data (FASTQ) files
# We assume the files exist in a "fastq" directory
FASTQ_DIR=${TOP_DIR}"/fastq/*_R*.fastq*"
READ1_STR="R1"
READ2_STR="R2"

if ls $FASTQ_DIR 1> /dev/null 2>&1;
then
    echo "ʕ·ᴥ·ʔ : Looking for fastq files...fastq files exist"
    testname=$(ls -l ${FASTQ_DIR} | awk 'NR==1{print $9}')
    if [ "${testname: -3}" == ".gz" ]
    then
        read1=${TOP_DIR}"/fastq/*${READ1_STR}*.fastq.gz"
    else
        read1=${TOP_DIR}"/fastq/*${READ1_STR}*.fastq"
    fi
else
    echo "ʕ·ᴥ·ʔ : Failed to find any files matching ${FASTQ_DIR}"
    exit
fi

if  [ "$APP_MODE" = 1 ]
then
  export WORK_DIR="data/scratch/polar-bear-fda-eua"
else
  export WORK_DIR="${TOP_DIR}/polar-bear-fda-eua"
fi

# Check to make sure output folders do not already exist

if ! mkdir "${WORK_DIR}" >/dev/null 2>&1; then echo "ʕ·ᴥ·ʔ : Unable to create ${WORK_DIR}! Exiting!"; exit 1; fi
if ! mkdir "${WORK_DIR}/aligned">/dev/null 2>&1; then echo "ʕ·ᴥ·ʔ : Unable to create ${WORK_DIR}/aligned! Exiting!"; exit 1; fi
if ! mkdir "${WORK_DIR}/debug">/dev/null 2>&1; then echo "ʕ·ᴥ·ʔ : Unable to create ${WORK_DIR}/debug! Exiting!"; exit 1; fi
if ! mkdir "${WORK_DIR}/final">/dev/null 2>&1; then echo "ʕ·ᴥ·ʔ : Unable to create ${WORK_DIR}/final! Exiting!"; exit 1; fi

# Create an array comprised of a FASTQ files
declare -a read1files=()
declare -a read2files=()

for i in ${read1}
do
    ext=${i#*$READ1_STR}
    name=${i%$READ1_STR*}
    name1=${name}${READ1_STR}
    name2=${name}${READ2_STR}
    read1files+=($name1$ext)
    read2files+=($name2$ext)
done

####### First block of work: Alignment of reads to reference
echo "ʕ·ᴥ·ʔ : Aligning files matching $FASTQ_DIR to $PATHOGEN_NAME reference assembly"

for ((i = 0; i < ${#read1files[@]}; ++i)); do
    file1=${read1files[$i]}
    file2=${read2files[$i]}

    FILE=$(basename ${file1%$read1str})
    ALIGNED_FILE=${WORK_DIR}/aligned/${FILE}"_mapped"

    # Align reads to viral reference
    bwa mem -t $THREADS $REFERENCE $file1 $file2 > $ALIGNED_FILE".sam" 2> ${WORK_DIR}/debug/align.out

    mv
    # Samtools fixmate fills in mate coordinates and insert size fields for deduping
    # Samtools fixmate is also converting SAM to BAM
    samtools fixmate -m $ALIGNED_FILE".sam" $ALIGNED_FILE"_matefixd.sam"

    # Sort reads based on position for deduping
    samtools sort -o $ALIGNED_FILE"_matefixd_sorted.sam" $ALIGNED_FILE"_matefixd.sam" 2> ${WORK_DIR}/debug/merge.out

done

# Merge BAMs if multiple BAMs were generated
samtools merge ${WORK_DIR}/aligned/sorted_merged.sam ${WORK_DIR}/aligned/*_matefixd_sorted.sam 2> data/logs/marge_out.txt

ls "${WORK_DIR}/aligned/" > data/logs/sanity_check_04.txt

######## Second block of work: Seperate viral data from control data
#echo "ʕ·ᴥ·ʔ : Removing Recombinants..."
#
#"${REMRECOMBO}" "${NT_TO_IS}" "${WORK_DIR}/aligned/sorted_merged.sam" "${QCUTOFF}" "${GOODBASECHANGE}" "${PE}" "${SEQSPLIT}" 2> ${WORK_DIR}/debug/recombo.out
#
#for SAM in ${WORK_DIR}/aligned/*sam;
#do
#    samtools view -hb $SAM > ${SAM%.sam}".bam"
#    rm $SAM
#done
#
## Get coverage of viral reference from read catagories
#echo "ʕ·ᴥ·ʔ : Analyzing Coverage..."
#
#samtools index "${WORK_DIR}/aligned/sorted_merged-good.bam"
#samtools index "${WORK_DIR}/aligned/sorted_merged-IS.bam"
#samtools index "${WORK_DIR}/aligned/sorted_merged-bad.bam"
#echo $'virus\taccukit\tchimeras'  > ${WORK_DIR}/aligned/ampliconCoverage.txt
#samtools bedcov -Q 4 "$AMPLICONS" "${WORK_DIR}/aligned/sorted_merged-good.bam" | awk '$1=="MN908947.3" { ar=int($9/($3-$2)); nt+=ar}END{printf ("%i\t",  nt)}' >> ${WORK_DIR}/aligned/ampliconCoverage.txt
#samtools bedcov -Q 4 "$AMPLICONS" "${WORK_DIR}/aligned/sorted_merged-IS.bam" | awk '$1 ~ /-SNAQ$/ { ar=int($9/($3-$2)); nt+=ar }END{printf ("%i\t",  nt)}' >> ${WORK_DIR}/aligned/ampliconCoverage.txt
#samtools bedcov -Q 4 "$AMPLICONS" "${WORK_DIR}/aligned/sorted_merged-bad.bam" | awk '{ ar=int($9/($3-$2)); nt+=ar}END{printf ("%i\n",  nt)}' >> ${WORK_DIR}/aligned/ampliconCoverage.txt
#
## Mark dups
#samtools markdup "${WORK_DIR}/aligned/sorted_merged-good.bam" "${WORK_DIR}/aligned/sorted_merged_dups_marked_viral.bam" 2> ${WORK_DIR}/debug/good_dedup.out
#samtools markdup "${WORK_DIR}/aligned/sorted_merged-IS.bam" "${WORK_DIR}/aligned/sorted_merged_dups_marked_IS.bam" 2> ${WORK_DIR}/debug/IS_dedup.out
#
## Get BoC
## To avoid cross-reaction with SARS a few regions are excluded from analysis
#samtools depth -a -b $NON_CROSS_REACT_REGIONS -Q 4 "${WORK_DIR}/aligned/sorted_merged_dups_marked_viral.bam" | awk '$1=="MN908947.3"' > ${WORK_DIR}/aligned/viral_depth_per_base.txt 2> ${WORK_DIR}/debug/viral_depth.out
#
## Gather alignment qc statistics
#echo "ʕ·ᴥ·ʔ :samtools flagstat result" > ${WORK_DIR}/aligned/all_alignment_stats.txt
#samtools flagstat "${WORK_DIR}/aligned/sorted_merged.bam"  >> ${WORK_DIR}/aligned/all_alignment_stats.txt
#
#echo "ʕ·ᴥ·ʔ :samtools flagstat result" > ${WORK_DIR}/aligned/viral_alignment_stats.txt
#samtools flagstat "${WORK_DIR}/aligned/sorted_merged_dups_marked_viral.bam"  >> ${WORK_DIR}/aligned/viral_alignment_stats.txt
#
#echo "ʕ·ᴥ·ʔ : samtools stats result " > ${WORK_DIR}/aligned/viral_alignment_stats.txt
#samtools stats "${WORK_DIR}/aligned/sorted_merged_dups_marked_viral.bam" >> ${WORK_DIR}/aligned/viral_alignment_stats.txt
#
## Write results to a file
#echo "ʕ·ᴥ·ʔ : Compiling results"
#"${PYTHON}" $COMPILE_RESULT $LIB_NAME $WORK_DIR
#echo "ʕ·ᴥ·ʔ : Pipeline completed, check ${WORK_DIR}/final for diagnositc result"
#
#if  [ "$APP_MODE" = 1 ]
#then
#  mkdir "${basespace_output_path_for_sample}/alignments/"
#  mv "${WORK_DIR}/aligned/sorted_merged-good.bam" "${basespace_output_path_for_sample}/alignments/"
#  mv "${WORK_DIR}/aligned/sorted_merged-IS.bam" "${basespace_output_path_for_sample}/alignments/"
#  mv "${WORK_DIR}/aligned/sorted_merged-bad.bam" "${basespace_output_path_for_sample}/alignments/"
#
#  mkdir "${basespace_output_path_for_sample}/results/"
#  mv "${WORK_DIR}/final/*" "${basespace_output_path_for_sample}/results/"

#  mv ${WORK_DIR}"/debug/*" data/logs/

#fi