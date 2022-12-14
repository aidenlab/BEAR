#!/bin/bash
: '
██████  ███████  █████  ██████
██   ██ ██      ██   ██ ██   ██
██████  █████   ███████ ██████
██   ██ ██      ██   ██ ██   ██
██████  ███████ ██   ██ ██   ██

Bioinformatics Evaluation of Assembly and Resequencing (BEAR) Pipeline
Serial version

By default, the BEAR pipeline will...
> Align data against a selection of potential viruses, sorts, merges and marks duplicates.
> Calculate alignment statistics.
> Assemble the data into contigs.
> Perform a pairwise alignment between the assembled contigs and the intended viral target.
> Create a final report with alignment statistics, a genome by genome dotplot, and the diagnostic result.
'

TOP_DIR=$(pwd)
threads=4
stage="all"

# Usage and commands
usageHelp="Usage: ${0##*/} [-d TOP_DIR] [-t THREADS] -jkrh"
dirHelp="* [TOP_DIR] is the top level directory (default: \"$TOP_DIR\")\n\
        Note: [TOP_DIR]/fastq must contain the fastq files"
threadsHelp="* [THREADS] is number of threads for BWA alignment (default: \"$threads\")"
reducedHelp="* -r Only align to SARS-CoV-2"
stageHelp="* -k Run specific stage of pipeline"
forceHelp="* -f Overwrite files if they already exist"
helpHelp="* -h Print this help and exit"

printHelpAndExit() {
    echo -e "$usageHelp"
    echo -e "$dirHelp"
    echo -e "$threadsHelp"
    echo -e "$reducedHelp"
    echo -e "$stageHelp"
    echo "$helpHelp"
    exit "$1"
}


while getopts "d:t:k:rfh" opt; do
    case $opt in
        d) TOP_DIR=$OPTARG ;;
        t) threads=$OPTARG ;;
        k) stage=$OPTARG ;;
        r) reducedSet=1 ;;
        f) forceOverwrite=1 ;;
        h) printHelpAndExit 0;;
       [?]) printHelpAndExit 1;;
    esac
done


### Set directories
PIPELINE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BETACORONA_REF_DIR="${PIPELINE_DIR}/betacoronaviruses/*/*/*.fasta"
BETACORONA_SMALL="${PIPELINE_DIR}/betacoronaviruses/close/*/*.fasta \
          ${PIPELINE_DIR}/betacoronaviruses/match/*/*.fasta"
MATCH_REF="${PIPELINE_DIR}/betacoronaviruses/match/*/*.fasta"
MATCH_NAME=$(echo $MATCH_REF | sed 's:.*/::' | rev | cut -c7- | rev )

export WORK_DIR=${TOP_DIR}/work
export LOG_DIR=${TOP_DIR}/log
export FINAL_DIR=${TOP_DIR}/final


Align_Reference ()
{
    REFERENCE=$1

    # get name of FASTA file and remove last ".fasta" to use as the reference name
    REFERENCE_NAME=$(echo $REFERENCE | sed 's:.*/::' | rev | cut -c7- | rev )
    REFERENCE_DIR=${WORK_DIR}/${REFERENCE_NAME}
    echo -e "ʕ•ᴥ•ʔ < Aligning files matching ${FASTQ_DIR}\n to the \"${REFERENCE_NAME}\" genome.)"

    # create a reference specific directory and subdirectories for log file and alignments
    if [[ "$forceOverwrite" -ne 1 ]]; then
      if ! mkdir "${REFERENCE_DIR}"; then echo "ʕ•ᴥ•ʔ < Error! Unable to create \"${REFERENCE_DIR}!\" Exiting.)"; exit 1; fi
      if ! mkdir "${REFERENCE_DIR}/aligned"; then echo "ʕ•ᴥ•ʔ < Error! Unable to create \"${REFERENCE_DIR}/aligned!\" Exiting.)"; exit 1; fi
      if ! mkdir "${REFERENCE_DIR}/debug"; then echo "ʕ•ᴥ•ʔ < Error! Unable to create \"${REFERENCE_DIR}/debug!\" Exiting.)"; exit 1; fi
    fi

    # loop through FASTQs and align against reference
    for ((i = 0; i < ${#read1files[@]}; ++i)); do
            file1=${read1files[$i]}
            file2=${read2files[$i]}

        FILE=$(basename ${file1%$read1str})
        ALIGNED_FILE="${REFERENCE_DIR}/aligned/${FILE}_mapped"

        # align using BWA mem w/ default paramaters
        bwa mem -t $threads $REFERENCE $file1 $file2 > "${ALIGNED_FILE}_temp.sam" 2> "${REFERENCE_DIR}/debug/align.out"

        # fill in mate coordinates and insert size fields, sort and compress alignment using Samtools
        samtools fixmate -m "${ALIGNED_FILE}_temp.sam" "${ALIGNED_FILE}_temp.bam"
        samtools sort -@ $threads -o "${ALIGNED_FILE}_matefixd_temp.bam" "${ALIGNED_FILE}_temp.bam"  2> "${REFERENCE_DIR}/debug/sort.out"
    done

    # if more than one, merge alignment files using Samtools
    if samtools merge "${REFERENCE_DIR}/aligned/sorted_merged.bam" "${REFERENCE_DIR}/aligned/"*"matefixd"* 2> "${REFERENCE_DIR}/debug/merge.out"
    then
        rm "${REFERENCE_DIR}/aligned/"*"temp"*
    fi

    # mark duplicates
    if samtools markdup "${REFERENCE_DIR}/aligned/sorted_merged.bam" "${REFERENCE_DIR}/aligned/sorted_dedupd.bam" 2> "${REFERENCE_DIR}/debug/dedup.out"
    then
        rm "${REFERENCE_DIR}/aligned/sorted_merged.bam"
    fi

    # calculate depth per base
    samtools depth -a "${REFERENCE_DIR}/aligned/sorted_dedupd.bam" > "${REFERENCE_DIR}/aligned/depth_per_base.txt" 2> "${REFERENCE_DIR}/debug/coverage.out"

    # calculate alignment statistics
    samtools stats "${REFERENCE_DIR}/aligned/sorted_dedupd.bam" | grep ^SN | cut -f 2- > "${REFERENCE_DIR}/aligned/stats.txt"

}

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

# Step 1: check for required software.

if [[ "$stage" != "assembly" ]]; then
    command -v bwa >/dev/null 2>&1 || { echo >&2 "ʕ•ᴥ•ʔ < Error! BWA required but it's not installed.)"; exit 1; }
    command -v samtools >/dev/null 2>&1 || { echo >&2 "ʕ•ᴥ•ʔ < Error! Samtools required but it's not installed.)"; exit 1; }
fi
command -v minimap2 >/dev/null 2>&1 || { echo >&2 "ʕ•ᴥ•ʔ < Error! Minimap2 required but it's not installed.)"; exit 1; }
command -v megahit >/dev/null 2>&1 || { echo >&2 "ʕ•ᴥ•ʔ < Error! Megahit required but it's not installed.)"; exit 1; }
command -v python >/dev/null 2>&1 || { echo >&2 "ʕ•ᴥ•ʔ < Error! Python required but it's not installed.)"; exit 1; }

# Step 2: make sure folders created by pipeline do not already exist.
if [[ "$forceOverwrite" -ne 1 && "$stage" != "assembly" && "$stage" != "plot" ]]; then
  if ! mkdir "${WORK_DIR}"; then echo "ʕ•ᴥ•ʔ < Error! Unable to create ${WORK_DIR}! Exiting.)"; exit 1; fi
  if ! mkdir "${LOG_DIR}"; then echo "ʕ•ᴥ•ʔ < Error! Unable to create ${LOG_DIR}! Exiting.)"; exit 1; fi
  if ! mkdir "${FINAL_DIR}"; then echo "ʕ•ᴥ•ʔ < Error! Unable to create ${FINAL_DIR}! Exiting.)"; exit 1; fi
fi 

# Step 3: check that FASTQs exist.

if [[ "$stage" != "plot" ]]; then
    FASTQ_DIR="${TOP_DIR}/fastq/*_R*.fastq*"
    READ1_STR="_R1"
    READ2_STR="_R2"

    if [[ ! -d "$TOP_DIR/fastq" ]]; then
        echo "ʕ•ᴥ•ʔ < Error! Directory \"${TOP_DIR}/fastq\" does not exist.)"
        echo "ʕ•ᴥ•ʔ < Create \"${TOP_DIR}/fastq\" and put fastq files to be aligned there. Exiting.)"
        printHelpAndExit 1
    else
        if stat -t ${FASTQ_DIR} >/dev/null 2>&1
        then
            echo "ʕ•ᴥ•ʔ < Looking for fastq files...fastq files exist.)"
            testname=$(ls -l ${FASTQ_DIR} | awk 'NR==1{print $9}')
            if [ "${testname: -3}" == ".gz" ]
            then
                read1="${TOP_DIR}/fastq/*${READ1_STR}*.fastq.gz"
            else
                read1="${TOP_DIR}/fastq/*${READ1_STR}*.fastq"
            fi
        else
            echo "ʕ•ᴥ•ʔ < Error! Failed to find any files matching \"${FASTQ_DIR}.\" Exiting.)"
            printHelpAndExit 1
        fi
    fi
fi

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


# Step 4: define what reference(s) to align against.
if [[ "$reducedSet" -eq 1 ]]
then
    REFERENCES=$MATCH_REF
else
    REFERENCES=$BETACORONA_REF_DIR
fi

# Step 5: align reads against references.
if [[ "$stage" == "align" || "$stage" == "all" ]]; then

    export -f Align_Reference
    for REFERENCE in $REFERENCES
    do
        Align_Reference $REFERENCE &
    done
    wait

    echo "ʕ•ᴥ•ʔ < Done with alignment.)"
fi

# Step 6: compile breadth of coverage values
echo "label,percentage" > "${WORK_DIR}/stats.csv"
for DEPTH_PER_BASE in "${WORK_DIR}/"*"/aligned/depth_per_base.txt"; do
    awk -v fname=$(basename ${DEPTH_PER_BASE%%/aligned*}) 'BEGIN{count=0; onisland=0}$3>0{if (!onisland){onisland=1; island_start=$2}}$3<=0{if (onisland){island_end=$2; if (island_end-island_start>=20){count=count+island_end-island_start}} onisland=0}END{if (onisland){island_end=$2; if (island_end-island_start>=50){count=count+island_end-island_start}} if (NR==0){NR=1} printf("%s,%0.02f\n", fname, count*100/NR) }' $DEPTH_PER_BASE >> ${WORK_DIR}/stats.csv
done


if [[ "$stage" == "assembly" || "$stage" == "all" ]]; then

    read_count=$(less "${WORK_DIR}/${MATCH_NAME}/aligned/stats.txt" | grep "raw total sequences:" | awk '{print $4}')

    if [ "$read_count" -gt 500 ]; then
        # Step 7: produce assembly
        if [[ "$noReport" -ne 1 ]]; then
          read1filescomma=$(echo "${read1files[*]}" | sed 's/ /,/g;s/,$//')
          read2filescomma=$(echo "${read2files[*]}" | sed 's/ /,/g;s/,$//')

          echo "ʕ•ᴥ•ʔ < Assembling contigs.)"
          megahit -m 750 -1 $read1filescomma -2 $read2filescomma -o "${WORK_DIR}/contigs"   &> ${LOG_DIR}/contig.out
          if [[ -f "${WORK_DIR}/contigs/done" ]]; then

            mv "${WORK_DIR}/contigs/final.contigs.fa" ${FINAL_DIR}/.
            echo "ʕ•ᴥ•ʔ < Performing pairwise comparison)"
            minimap2 -x asm5 $MATCH_REF "${FINAL_DIR}/final.contigs.fa" > "${WORK_DIR}/contigs/contig.paf" 2> "${LOG_DIR}/minimap.out"
            if [[ ! -s "${WORK_DIR}/contigs/contig.paf" ]]; then
              echo "ʕ•ᴥ•ʔ < Error! Pairwise alignment by minimap2 failed.)"
            fi
          else
            echo "ʕ•ᴥ•ʔ < Error! Assembly failed.)"
          fi
        fi
    else
        echo "ʕ•ᴥ•ʔ < Error! Insufficient reads for assembly!)"
    fi
fi

if [[ "$stage" == "plot" || "$stage" == "all" ]]; then
    CONTIG_LENGTH=$(tail -n2 ${WORK_DIR}/contigs/log |grep -o 'total.*' | awk '{print $2}')
    python "${PIPELINE_DIR}/remove_strays.py" "${WORK_DIR}/${MATCH_NAME}/aligned/depth_per_base.txt" "${WORK_DIR}/${MATCH_NAME}/aligned/depth_per_base_without_primers_islands.txt"
    python "${PIPELINE_DIR}/dot_coverage.py" "${WORK_DIR}/${MATCH_NAME}/aligned/depth_per_base_without_primers_islands.txt" "${WORK_DIR}/contigs/contig.paf" "${WORK_DIR}/stats.csv" $CONTIG_LENGTH "${FINAL_DIR}/report" &> "${LOG_DIR}/dotplot.out"
fi

echo "ʕ•ᴥ•ʔ < Pipeline completed, check \"${FINAL_DIR}\" for the report.)"
