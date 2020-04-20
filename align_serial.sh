#!/bin/bash
##### 
### Polar pipeline script for viral diagnostic
### Serial version
### Given paired end sequencing of putative viral data
###   -> Aligns to a selection of potential viruses, sorts and merges
###   -> Assembles the data into contigs
###   -> After alignment, runs samtools depth to create stats.csv
###   -> After contigging, creates pairwise alignment
###   -> Final report created as PDF with dotplot from paf and stats
#####

### REQUIRED SOFTWARE 
## You must have the following software installed
## and available in your PATH
## BWA; Samtools; Minimap2; Megahit; Python

## Threads
threads=1

### VARIABLES: Automatically set
TOP_DIR=$(pwd)
PIPELINE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BETACORONA_REF_DIR="${PIPELINE_DIR}/betacoronaviruses/*/*/*.fasta"
BETACORONA_SMALL="${PIPELINE_DIR}/betacoronaviruses/close/*/*.fasta \
		  ${PIPELINE_DIR}/betacoronaviruses/match/*/*.fasta" 

# Usage and commands
usageHelp="Usage: ${0##*/} [-d TOP_DIR] [-t THREADS] -jrh"
dirHelp="* [TOP_DIR] is the top level directory (default \"$TOP_DIR\")\n\
  [TOP_DIR]/fastq must contain the fastq files"
threadsHelp="* [THREADS] is number of threads for BWA alignment"
indexHelp="* -j produce index file for aligned files"
reducedHelp="* -r reduced set for alignment"
helpHelp="* -h: print this help and exit"

printHelpAndExit() {
    echo -e "$usageHelp"
    echo -e "$dirHelp"
    echo -e "$threadsHelp"
    echo -e "$indexHelp"
    echo -e "$reducedHelp"
    echo "$helpHelp"
    exit "$1"
}

while getopts "d:t:hrj" opt; do
    case $opt in
	h) printHelpAndExit 0;;
        d) TOP_DIR=$OPTARG ;;
	j) produceIndex=1 ;;
	r) reducedSet=1 ;;
	t) threads=$OPTARG ;;
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
command -v minimap2 >/dev/null 2>&1 || { echo >&2 "!*** Minimap2 required but it's not installed."; exit 1; }
command -v megahit >/dev/null 2>&1 || { echo >&2 "!*** Megahit required but it's not installed."; exit 1; }
command -v python >/dev/null 2>&1 || { echo >&2 "!*** Python required but it's not installed."; exit 1; }

# Check for fastq files
if [ ! -d "$TOP_DIR/fastq" ]; then
    echo "Directory \"$TOP_DIR/fastq\" does not exist."
    echo "Create \"$TOP_DIR/fastq\" and put fastq files to be aligned there."
    printHelpAndExit 1
else
    if stat -t ${FASTQ_DIR} >/dev/null 2>&1
    then
        echo "(-: Looking for fastq files...fastq files exist"
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

read1files=()
read2files=()
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

# replace spaces with commas for megahit
read1filescomma=$(echo "${read1files[*]}" | sed 's/ /,/g;s/,$//')
read2filescomma=$(echo "${read2files[*]}" | sed 's/ /,/g;s/,$//')

if [[ "$reducedSet" -eq 1 ]]
then
    REFERENCES=$BETACORONA_SMALL
else
    REFERENCES=$BETACORONA_REF_DIR
fi

WORK_DIR=${TOP_DIR}/work
LOG_DIR=${TOP_DIR}/log
FINAL_DIR=${TOP_DIR}/final

if ! mkdir "${WORK_DIR}"; then echo "***! Unable to create ${WORK_DIR}! Exiting"; exit 1; fi
if ! mkdir "${LOG_DIR}"; then echo "***! Unable to create ${LOG_DIR}! Exiting"; exit 1; fi
if ! mkdir "${FINAL_DIR}"; then echo "***! Unable to create ${FINAL_DIR}! Exiting"; exit 1; fi

for REFERENCE in $REFERENCES
do
    ######################################################################
    ########## Align 
    ######################################################################
    REFERENCE_NAME=$(echo $REFERENCE | sed 's:.*/::' | rev | cut -c7- | rev )
    if [[ "$REFERENCE" == *match* ]]
    then
	MATCH_REF=${REFERENCE}
	MATCH_NAME=${REFERENCE_NAME}
    fi

    echo -e "(-: Aligning files matching $FASTQ_DIR\n to genome $REFERENCE_NAME"

    if ! mkdir "${WORK_DIR}/${REFERENCE_NAME}"; then echo "***! Unable to create ${WORK_DIR}/${REFERENCE_NAME}! Exiting"; exit 1; fi
    if ! mkdir "${WORK_DIR}/${REFERENCE_NAME}/aligned"; then echo "***! Unable to create ${WORK_DIR}/${REFERENCE_NAME}/aligned! Exiting"; exit 1; fi
    if ! mkdir "${WORK_DIR}/${REFERENCE_NAME}/debug"; then echo "***! Unable to create ${WORK_DIR}/${REFERENCE_NAME}/debug! Exiting"; exit 1; fi

    for ((i = 0; i < ${#read1files[@]}; ++i)); do
        file1=${read1files[$i]}
        file2=${read2files[$i]}

	FILE=$(basename ${file1%$read1str})
	ALIGNED_FILE=${WORK_DIR}/${REFERENCE_NAME}/aligned/${FILE}"_mapped"

        # Align reads
	bwa mem -t $threads $REFERENCE $file1 $file2 > $ALIGNED_FILE".sam" 2> ${WORK_DIR}/${REFERENCE_NAME}/debug/align.out

	# Samtools fixmate and sort, output as BAM
	samtools fixmate -m $ALIGNED_FILE".sam" $ALIGNED_FILE".bam"
	samtools sort -@ $threads -o $ALIGNED_FILE"_matefixd_sorted.bam" $ALIGNED_FILE".bam"  2> ${WORK_DIR}/${REFERENCE_NAME}/debug/sort.out
    done

    # Merge sorted BAMs
    if samtools merge ${WORK_DIR}/${REFERENCE_NAME}/aligned/sorted_merged.bam ${WORK_DIR}/${REFERENCE_NAME}/aligned/*_matefixd_sorted.bam 2> ${WORK_DIR}/${REFERENCE_NAME}/debug/merge.out
    then
	rm ${WORK_DIR}/${REFERENCE_NAME}/aligned/*_sorted.bam 
	rm ${WORK_DIR}/${REFERENCE_NAME}/aligned/*.sam  
    fi
    
    if samtools markdup ${WORK_DIR}/${REFERENCE_NAME}/aligned/sorted_merged.bam ${WORK_DIR}/${REFERENCE_NAME}/aligned/sorted_merged_dups_marked.bam 2> ${WORK_DIR}/${REFERENCE_NAME}/debug/dedup.out
    then
	rm ${WORK_DIR}/${REFERENCE_NAME}/aligned/sorted_merged.bam
    fi

    samtools depth -a ${WORK_DIR}/${REFERENCE_NAME}/aligned/sorted_merged_dups_marked.bam > ${WORK_DIR}/${REFERENCE_NAME}/aligned/depth_per_base.txt 2> ${WORK_DIR}/${REFERENCE_NAME}/debug/coverage.out

    # In case you want to visualize the bams, index them. 
    if [ -n "$produceIndex" ]
    then
	samtools index ${WORK_DIR}/${REFERENCE_NAME}/aligned/sorted_merged_dups_marked.bam 2> ${WORK_DIR}/${REFERENCE_NAME}/debug/index.out
    fi

    # Statistics 
    samtools stats ${WORK_DIR}/${REFERENCE_NAME}/aligned/sorted_merged_dups_marked.bam | grep ^SN | cut -f 2- > ${WORK_DIR}/${REFERENCE_NAME}/aligned/stats.txt
done

echo "(-: Done with alignment" 
# Gather alignment statistics (coverage %)
echo "label,percentage" > ${WORK_DIR}/stats.csv
for f in ${WORK_DIR}/*/aligned/depth_per_base.txt
do
    awk -v fname=$(basename ${f%%/aligned*}) '$3>0{count++}END{if (NR==0){NR=1} printf("%s,%0.02f\n", fname, count*100/NR)}' $f >> ${WORK_DIR}/stats.csv
done

# Produce contigs - this can happen concurrently with alignment
megahit -1 $read1filescomma -2 $read2filescomma -o ${WORK_DIR}/contigs -m 750 --min-contig-len 100  &> ${LOG_DIR}/contig.out
mv ${WORK_DIR}/contigs/final.contigs.fa ${FINAL_DIR}/.
CONTIG_LENGTH=$(tail -n2 ${WORK_DIR}/contigs/log |grep -o 'total.*' | awk '{print $2}')
echo "(-: Done with contigs" 
# Dot plot - after contig and alignment
minimap2 -x asm5 $MATCH_REF ${FINAL_DIR}/final.contigs.fa > ${WORK_DIR}/contig.paf 2> ${LOG_DIR}/minimap.out
if [ ! -s "${WORK_DIR}/contig.paf" ]
then
    echo "!*** Pairwise alignment by minimap2 failed."
    touch ${WORK_DIR}/contig.paf
fi
echo "(-: Done with pairwise comparison" 
python ${PIPELINE_DIR}/dot_coverage.py ${WORK_DIR}/${MATCH_NAME}/aligned/depth_per_base.txt ${WORK_DIR}/contig.paf ${WORK_DIR}/stats.csv $CONTIG_LENGTH ${FINAL_DIR}/report --crop_y &> ${LOG_DIR}/dotplot.out

echo "(-: Done with dotplot"
echo "(-: Pipeline completed, check ${FINAL_DIR} for the report"
