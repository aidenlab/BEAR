#!/bin/bash
### Polar BEAR FDA EUA pipeline
### NEED TO DO

### REQUIRED SOFTWARE 
### NEED TO DO


## Threads
threads=32

## Target pathogen
VIRUS="SARS-CoV-2"

### VARIABLES: Automatically set
TOP_DIR=$(pwd)
LIB_NAME=$(pwd | awk -F "/" '{print $NF}')
PIPELINE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
VIRAL_REFERENCE="${PIPELINE_DIR}/reference_files/SARS-CoV_ISv0.4.1_seperate.fasta.gz"
NT_TO_IS="${PIPELINE_DIR}/reference_files/NT_IS_LOOKUP_TABLE_v0.4.2_seperate.txt"
AMPLICONS="${PIPELINE_DIR}/reference_files/VarDict-amplicon.v2.1.bed"
REMRECOMBO="${PIPELINE_DIR}/accuGenomics/remRecombo"
COMPILE_RESULT="${PIPELINE_DIR}/compile_results.py"
NON_CROSS_REACT_REGIONS="${PIPELINE_DIR}/reference_files/noncross_reactive_sars_cov_2_regions.bed"

# Usage and commands
usageHelp="Usage: ${0##*/} [-d TOP_DIR] [-t THREADS] -jkrh"
dirHelp="* [TOP_DIR] is the top level directory (default \"$TOP_DIR\")\n\
  [TOP_DIR]/fastq must contain the fastq files"
threadsHelp="* [THREADS] is number of threads for BWA alignment"
indexHelp="* -j produce index file for aligned files"
stageHelp="* -k start pipeline after alignment"
helpHelp="* -h: print this help and exit"

printHelpAndExit() {
    echo -e "$usageHelp"
    echo -e "$dirHelp"
    echo -e "$threadsHelp"
    echo -e "$indexHelp"
    echo -e "$reducedHelp"
    echo -e "$stageHelp"
    echo "$helpHelp"
    exit "$1"
}

while getopts "d:t:hjkn" opt; do
    case $opt in
    h) printHelpAndExit 0;;
        d) TOP_DIR=$OPTARG ;;
    j) produceIndex=1 ;;
    t) threads=$OPTARG ;;
    k) afteralignment=1 ;;
    n) afteralignment=compileresultonly=1 ;;
    [?]) printHelpAndExit 1;;
    esac
done

# We assume the files exist in a fastq directory
FASTQ_DIR=${TOP_DIR}"/fastq/*_R*.fastq*"
READ1_STR="_R1"
READ2_STR="_R2"

# Check for installed software
command -v bwa >/dev/null 2>&1 || { echo >&2 "!*** BWA required but it's not installed."; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo >&2 "!*** Samtools required but it's not installed."; exit 1; }
command -v python >/dev/null 2>&1 || { echo >&2 "!*** Python required but it's not installed."; exit 1; }

# Check for fastq files
if [ ! -d "$TOP_DIR/fastq" ]; then
    echo "Directory \"$TOP_DIR/fastq\" does not exist."
    echo "Create \"$TOP_DIR/fastq\" and put fastq files to be aligned there."
    printHelpAndExit 1
else
    if stat -t ${FASTQ_DIR} >/dev/null 2>&1
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
        echo "***! Failed to find any files matching ${FASTQ_DIR}"
    printHelpAndExit 1
    fi
fi

declare -a read1files=()
declare -a read2files=()
for i in ${read1}
do
    ext=${i#*$READ1_STR}
    name=${i%$READ1_STR*}
    # these names have to be right or it'll break                                                                            
    name1=${name}${READ1_STR}
    name2=${name}${READ2_STR}
    read1files+=($name1$ext)
    read2files+=($name2$ext)
done

export WORK_DIR=${TOP_DIR}/polar-bear-fda-eua

######################################################################
########## Align viral reference
######################################################################
if [[ "$afteralignment" -ne 1 ]]
then
    if ! mkdir "${WORK_DIR}"; then echo "***! Unable to create ${WORK_DIR}! Exiting"; exit 1; fi

    echo -e "ʕ·ᴥ·ʔ : Aligning files matching $FASTQ_DIR\n to $VIRUS reference assembly and Accukit control sequence"

    if ! mkdir "${WORK_DIR}/virus"; then echo "***! Unable to create ${WORK_DIR}/virus! Exiting"; exit 1; fi
    if ! mkdir "${WORK_DIR}/virus/aligned"; then echo "***! Unable to create ${WORK_DIR}/virus/aligned! Exiting"; exit 1; fi
    if ! mkdir "${WORK_DIR}/virus/debug"; then echo "***! Unable to create ${WORK_DIR}/virus/debug! Exiting"; exit 1; fi
    
    for ((i = 0; i < ${#read1files[@]}; ++i)); do
    file1=${read1files[$i]}
    file2=${read2files[$i]}
    
    FILE=$(basename ${file1%$read1str})
    ALIGNED_FILE=${WORK_DIR}/virus/aligned/${FILE}"_mapped"
    
        # Align reads to virus
    bwa mem -t $threads $VIRAL_REFERENCE $file1 $file2 > $ALIGNED_FILE".sam" 2> ${WORK_DIR}/virus/debug/align.out
    
        # Samtools fixmate and sort, output as BAM
    samtools fixmate -m $ALIGNED_FILE".sam" $ALIGNED_FILE".bam"
    samtools sort -@ $threads -o $ALIGNED_FILE"_matefixd_sorted.bam" $ALIGNED_FILE".bam"  2> ${WORK_DIR}/virus/debug/sort.out

    done

        # Merge sorted BAMs
    if samtools merge ${WORK_DIR}/virus/aligned/sorted_merged.bam ${WORK_DIR}/virus/aligned/*_matefixd_sorted.bam 2> ${WORK_DIR}/virus/debug/merge.out
    then
    rm ${WORK_DIR}/virus/aligned/*_sorted.bam 
    rm ${WORK_DIR}/virus/aligned/*.sam  
    rm ${WORK_DIR}/virus/aligned/*"_mapped"*bam  
    fi

fi

if [[ "$compileresultonly" -ne 1 ]]
then
    echo "ʕ·ᴥ·ʔ : Done with alignment" 

    echo "ʕ·ᴥ·ʔ : Removing Recombinants..."

    "${REMRECOMBO}" "${NT_TO_IS}" ${WORK_DIR}/virus/aligned/sorted_merged.bam 0  2> ${WORK_DIR}/virus/debug/recombo.out

    echo "ʕ·ᴥ·ʔ : Analyzing Coverage..."

    echo $'virus\taccukit\tchimeras'  > ${WORK_DIR}/virus/aligned/ampliconCoverage.txt 
    samtools bedcov -Q 4 "$AMPLICONS" "${WORK_DIR}/virus/aligned/sorted_merged-good.bam" | awk '$1=="MN908947.3" { ar=int($9/($3-$2)); nt+=ar}END{printf ("%i\t",  nt)}' > ${WORK_DIR}/virus/aligned/ampliconCoverage.txt 
    samtools bedcov -Q 4 "$AMPLICONS" "${WORK_DIR}/virus/aligned/sorted_merged-IS.bam" | awk '$1 ~ /-SNAQ$/ { ar=int($9/($3-$2)); nt+=ar }END{printf ("%i\t",  nt)}' >> ${WORK_DIR}/virus/aligned/ampliconCoverage.txt 
    samtools bedcov -Q 4 "$AMPLICONS" "${WORK_DIR}/virus/aligned/sorted_merged-bad.bam" | awk '{ ar=int($9/($3-$2)); nt+=ar}END{printf ("%i\n",  nt)}' >> ${WORK_DIR}/virus/aligned/ampliconCoverage.txt

    samtools markdup "${WORK_DIR}/virus/aligned/sorted_merged-good.bam" "${WORK_DIR}/virus/aligned/sorted_merged_dups_marked_good.bam" 2> ${WORK_DIR}/virus/debug/good_dedup.out
    samtools markdup "${WORK_DIR}/virus/aligned/sorted_merged-IS.bam" "${WORK_DIR}/virus/aligned/sorted_merged_dups_marked_IS.bam" 2> ${WORK_DIR}/virus/debug/IS_dedup.out

    # Gather alignment statistics

    samtools depth -a -b $NON_CROSS_REACT_REGIONS -Q 4 "${WORK_DIR}/virus/aligned/sorted_merged_dups_marked_good.bam" | awk '$1=="MN908947.3"' > ${WORK_DIR}/virus/aligned/viral_depth_per_base.txt 2> ${WORK_DIR}/virus/debug/viral_depth.out
    samtools depth -a -Q 4 "${WORK_DIR}/virus/aligned/sorted_merged_dups_marked_IS.bam" | awk '$1!="MN908947.3"' > ${WORK_DIR}/virus/aligned/control_depth_per_base.txt 2> ${WORK_DIR}/virus/debug/control_depth.out

    # In case you want to visualize the bams, index them. 
    if [ -n "$produceIndex" ]
    then
        samtools index ${WORK_DIR}/virus/aligned/sorted_merged_dups_marked-good.bam 2> ${WORK_DIR}/virus/debug/index.out
    fi

    # Statistics 
    echo "ʕ·ᴥ·ʔ :samtools flagstat result" > ${WORK_DIR}/virus/aligned/alignment_stats.txt
    samtools flagstat "${WORK_DIR}/virus/aligned/sorted_merged_dups_marked_good.bam"  >> ${WORK_DIR}/virus/aligned/alignment_stats.txt

    echo "ʕ·ᴥ·ʔ : samtools stats result " >> ${WORK_DIR}/virus/aligned/alignment_stats.txt
    samtools stats "${WORK_DIR}/virus/aligned/sorted_merged_dups_marked_good.bam" >> ${WORK_DIR}/virus/aligned/alignment_stats.txt

fi

echo "ʕ·ᴥ·ʔ : Compiling results" 

LIB_NAME=$(echo $TOP_DIR | awk -F "/" '{print $NF}')
python $COMPILE_RESULT ${WORK_DIR}/virus/aligned/viral_depth_per_base.txt ${WORK_DIR}/virus/result.txt $LIB_NAME

FIN=$(less ${WORK_DIR}/virus/result.txt)
echo "ʕ·ᴥ·ʔ " $FIN

# echo "ʕ·ᴥ·ʔ : Pipeline completed, check ${WORK_DIR} for diagnositc result"
