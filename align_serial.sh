#!/bin/bash
##### 
### Polar pipeline script for viral diagnostic
### Serial version
### Given paired end sequencing of putative viral data
###   -> Aligns to a selection of potential viruses, sorts and merges
###   -> Assembles the data into contigs
###   -> After alignment, runs flagstat aligned data to create stats.html
###   -> After contigging, creates dotplot to add to stats.html
###   -> Final report created as PDF from stats.html
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

# We assume the files exist in a fastq directory
FASTQ_DIR=${TOP_DIR}"/fastq/*_R*.fastq*"
READ1_STR="_R1"
READ2_STR="_R2"

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

while getopts "d:t:hr" opt; do
    case $opt in
	h) printHelpAndExit 0;;
        d) TOP_DIR=$OPTARG ;;
	j) produceIndex=1 ;;
	r) reducedSet=1 ;;
	t) threads=$OPTARG ;;
	[?]) printHelpAndExit 1;;
    esac
done

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
    echo -e "(-: Aligning files matching $FASTQ_DIR\n to genome $REFERENCE_NAME"

    if ! mkdir "${WORK_DIR}/${REFERENCE_NAME}"; then echo "***! Unable to create ${WORK_DIR}/${REFERENCE_NAME}! Exiting"; exit 1; fi
    if ! mkdir "${WORK_DIR}/${REFERENCE_NAME}/aligned"; then echo "***! Unable to create ${WORK_DIR}/${REFERENCE_NAME}/aligned! Exiting"; exit 1; fi
    if ! mkdir "${WORK_DIR}/${REFERENCE_NAME}/debug"; then echo "***! Unable to create ${WORK_DIR}/${REFERENCE_NAME}/debug! Exiting"; exit 1; fi

    for ((i = 0; i < ${#read1files[@]}; ++i)); do
        usegzip=0
        file1=${read1files[$i]}
        file2=${read2files[$i]}

	FILE=$(basename ${file1%$read1str})
	ALIGNED_FILE=${WORK_DIR}/${REFERENCE_NAME}/aligned/${FILE}"_mapped.sam"

        # Align reads
	bwa mem -t $threads $REFERENCE $file1 $file2 > $ALIGNED_FILE

	# Sort SAM and convert to BAM
	samtools sort -@ $threads $ALIGNED_FILE -o ${ALIGNED_FILE}"_sorted.bam"
    done

    # Merge sorted BAMs
    if samtools merge ${WORK_DIR}/${REFERENCE_NAME}/aligned/sorted_merged.bam ${WORK_DIR}/${REFERENCE_NAME}/aligned/*_sorted.bam
    then
	rm ${WORK_DIR}/${REFERENCE_NAME}/aligned/*_sorted.bam 
	rm ${WORK_DIR}/${REFERENCE_NAME}/aligned/*.sam  
    fi

    if [[ "$REFERENCE" == *match* ]]
    then
	matchname=${REFERENCE_NAME}
	matchref=${REFERENCE}
    fi

    # In case you want to visualize the bams, index them. 
    if [ -n "$produceIndex" ]
    then
	samtools index ${WORK_DIR}/${REFERENCE_NAME}/aligned/sorted_merged.bam
    fi

    # Statistics (in particular percentage of mapped reads)
    samtools flagstat ${WORK_DIR}/${REFERENCE_NAME}/aligned/sorted_merged.bam > ${WORK_DIR}/${REFERENCE_NAME}/aligned/stats.txt 
done

# Gather alignment statistics (mapping %)
echo "label,percentage" > ${WORK_DIR}/stats.csv
for f in ${WORK_DIR}/*/aligned/stats.txt
do
    awk -v fname=$(basename ${f%%/aligned*}) 'BEGIN{OFS=","}$4=="mapped"{split($5, a, "("); split(a[2],b, "%"); print fname, b[1]}' $f >> ${WORK_DIR}/stats.csv
done

# Produce contigs - this can happen concurrently with alignment
megahit -1 $read1filescomma -2 $read2filescomma -m 750 -o ${WORK_DIR}/contigs
mv ${WORK_DIR}/contigs/final.contigs.fa ${FINAL_DIR}/.
CONTIG_LENGTH=$(tail -n2 ${WORK_DIR}/contigs/log |grep -o 'total.*' | awk '{print $2}')

# Dot plot - after contig and alignment
minimap2 -x asm5 $matchref ${FINAL_DIR}/final.contigs.fa > ${WORK_DIR}/contig_${matchname}.paf
if [ ! -s "${WORK_DIR}/contig_${matchname}.paf" ]
then
    echo "!*** Pairwise alignment by minimap2 failed."
    exit 1
fi
samtools rmdup ${WORK_DIR}/${matchname}/aligned/sorted_merged.bam ${WORK_DIR}/${matchname}/aligned/sorted_merged_dedup.bam 
samtools depth ${WORK_DIR}/${matchname}/aligned/sorted_merged_dedup.bam > ${WORK_DIR}/${matchname}/aligned/depth_per_base.txt	
python ${PIPELINE_DIR}/dot_coverage.py ${WORK_DIR}/${matchname}/aligned/depth_per_base.txt ${WORK_DIR}/contig_${matchname}.paf ${WORK_DIR}/stats.csv $CONTIG_LENGTH ${FINAL_DIR}/report
