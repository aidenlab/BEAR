#!/bin/bash
##### 
### Polar pipeline script for viral diagnostic
### Serial version
### Given paired end sequencing of putative viral data
###   -> Aligns to a selection of potential viruses, sorts and merges
###   -> Assembles the data into contigs
###   -> After alignment, runs samtools depth to create result.csv
###   -> After contigging, creates pairwise alignment
###   -> Final report created as PDF with dotplot from paf and stats
#####

### REQUIRED SOFTWARE 
## You must have the following software installed
## and available in your PATH
## BWA; Samtools

## Threads
threads=16
VIRUS="SARS-CoV-2"

### VARIABLES: Automatically set
TOP_DIR=$(pwd)
PIPELINE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REFERENCE="${PIPELINE_DIR}/reference_files/SARS-CoV_ISv0.4.1_seperate.fasta.gz"
NT_TO_IS="${PIPELINE_DIR}/reference_files/"
REMRECOMBO="${PIPELINE_DIR}/accuGenomics/remRecombo"
COVERAGE="${PIPELINE_DIR}/accuGenomics/coverage"
HUMAN_REF="${PIPELINE_DIR}/reference_files/chrM.fa.gz"

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

while getopts "d:t:hjk" opt; do
    case $opt in
	h) printHelpAndExit 0;;
        d) TOP_DIR=$OPTARG ;;
	j) produceIndex=1 ;;
	t) threads=$OPTARG ;;
	k) afteralignment=1 ;;
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

Align_Reference ()
{
    ######################################################################
    ########## Align 
    ######################################################################
    REFERENCE_NAME=$(echo $1 | sed 's:.*/::' | rev | cut -c7- | rev )

    echo -e "ʕ·ᴥ·ʔ : Aligning files matching $FASTQ_DIR\n to $VIRUS reference assembly and Accukit control sequence"

    if ! mkdir "${WORK_DIR}/${REFERENCE_NAME}"; then echo "***! Unable to create ${WORK_DIR}/${REFERENCE_NAME}! Exiting"; exit 1; fi
    if ! mkdir "${WORK_DIR}/${REFERENCE_NAME}/aligned"; then echo "***! Unable to create ${WORK_DIR}/${REFERENCE_NAME}/aligned! Exiting"; exit 1; fi
    if ! mkdir "${WORK_DIR}/${REFERENCE_NAME}/debug"; then echo "***! Unable to create ${WORK_DIR}/${REFERENCE_NAME}/debug! Exiting"; exit 1; fi
    
    for ((i = 0; i < ${#read1files[@]}; ++i)); do
        file1=${read1files[$i]}
        file2=${read2files[$i]}
        
        FILE=$(basename ${file1%$read1str})
        ALIGNED_FILE=${WORK_DIR}/${REFERENCE_NAME}/aligned/${FILE}"_mapped"
        
        # Align reads
        bwa mem -k 32 -t $threads $1 $file1 $file2 > $ALIGNED_FILE".sam" 2> ${WORK_DIR}/${REFERENCE_NAME}/debug/align.out
        # Samtools fixmate and sort, output as BAM
        samtools fixmate -m $ALIGNED_FILE".sam" $ALIGNED_FILE".bam"
        samtools sort -@ $threads -o $ALIGNED_FILE"_matefixd_sorted.bam" $ALIGNED_FILE".bam"  2> ${WORK_DIR}/${REFERENCE_NAME}/debug/sort.out
    done

        # Merge sorted BAMs
    if samtools merge ${WORK_DIR}/${REFERENCE_NAME}/aligned/sorted_merged.bam ${WORK_DIR}/${REFERENCE_NAME}/aligned/*_matefixd_sorted.bam 2> ${WORK_DIR}/${REFERENCE_NAME}/debug/merge.out
    then
        rm ${WORK_DIR}/${REFERENCE_NAME}/aligned/*_sorted.bam 
        rm ${WORK_DIR}/${REFERENCE_NAME}/aligned/*.sam  
        rm ${WORK_DIR}/${REFERENCE_NAME}/aligned/*"_mapped"*bam  
    fi

}

if [[ "$afteralignment" -ne 1 ]]
then
    if ! mkdir "${WORK_DIR}"; then echo "***! Unable to create ${WORK_DIR}! Exiting"; exit 1; fi
    
    export -f Align_Reference
    for REF in $REFERENCE $HUMAN_REF
    do
        Align_Reference $REF &
    done
    wait
fi
echo "ʕ·ᴥ·ʔ : Done with alignment" 

echo "ʕ·ᴥ·ʔ : Removing Recombinants..."

REFERENCE_NAME=$(echo $1 | sed 's:.*/::' | rev | cut -c7- | rev )
"${REMRECOMBO}" "${NT_TO_IS}" ${WORK_DIR}/${REFERENCE_NAME}/aligned/sorted_merged.bam 0 

#samtools index "${pi}/${fn}.bam"
#rm -f "${pi}/${fn}-rsort.temp"

echo "ʕ·ᴥ·ʔ : Calculating Viral Load..."

#average reads per amplicon = total bases / expected fragment length 
# print coverage, median value
#echo -e "SAMPLE\tNT_READS\tIS_READS\tREC_READS\tNT_AVG\tIS_AVG\tREC_AVG" > "$coverage"
samtools bedcov -Q 4 "$NT_TO_IS" "${stem}-good.bam" | awk '$1=="MN908947.3" { ar=int($9/($3-$2)); nt+=ar}END{printf ("%i\t",  nt)}' > ${WORK_DIR}/aligned/ampliconCoverage.txt 
samtools bedcov -Q 4 "$NT_TO_IS" "${stem}-IS.bam" | awk '$1 ~ /-SNAQ$/) { ar=int($9/($3-$2)); nt+=ar }END{printf ("%i\t",  nt)}' >> ${WORK_DIR}/aligned/ampliconCoverage.txt 
samtools bedcov -Q 4 "$NT_TO_IS" "${stem}-bad.bam" | awk '{ ar=int($9/($3-$2)); nt+=ar}END{printf ("%i\n",  nt)}' >> ${WORK_DIR}/aligned/ampliconCoverage.txt 

awk -v threshold=$viralThreshold -v controlThreshold=$controlThreshold -v humanThreshold=$humanThreshold '{
 recomb = $3/2; 
 viral_load = 100*($1 + recomb)/($2 + recomb);
 if (viral_load > threshold) { 
   print "Positive";
 }
 else {
   if (100*$1 > controlThreshold) {
     if (humanThreshold) {
       print "Negative";
     }
     else {
       print "Undetermined, no human detected";
     }
   }
   else {
     print "Undetermined, no control detected";
   }
}'  ${WORK_DIR}/aligned/ampliconCoverage.txt


## NEVA HERE PAUSED
##
# Gather alignment statistics
echo "ʕ·ᴥ·ʔ : Compiling results" 
    
if samtools markdup ${WORK_DIR}/${REFERENCE_NAME}/aligned/sorted_merged.bam ${WORK_DIR}/${REFERENCE_NAME}/aligned/sorted_merged_dups_marked.bam 2> ${WORK_DIR}/${REFERENCE_NAME}/debug/dedup.out
then
    rm ${WORK_DIR}/${REFERENCE_NAME}/aligned/sorted_merged.bam
fi

samtools depth -a -Q 1 ${WORK_DIR}/${REFERENCE_NAME}/aligned/sorted_merged_dups_marked.bam | awk '$1=="MN908947.3"' > ${WORK_DIR}/${REFERENCE_NAME}/aligned/depth_per_base.txt 2> ${WORK_DIR}/${REFERENCE_NAME}/debug/depth.out

#    samtools coverage -q 1 ${WORK_DIR}/${REFERENCE_NAME}/aligned/sorted_merged_dups_marked.bam > ${WORK_DIR}/${REFERENCE_NAME}/aligned/coverage_1.tsv 2> ${WORK_DIR}/${REFERENCE_NAME}/debug/coverage_1.out
#    samtools coverage -q 0 ${WORK_DIR}/${REFERENCE_NAME}/aligned/sorted_merged_dups_marked.bam > ${WORK_DIR}/${REFERENCE_NAME}/aligned/coverage_0.tsv 2> ${WORK_DIR}/${REFERENCE_NAME}/debug/coverage_0.out

# In case you want to visualize the bams, index them. 
if [ -n "$produceIndex" ]
then
    samtools index ${WORK_DIR}/${REFERENCE_NAME}/aligned/sorted_merged_dups_marked.bam 2> ${WORK_DIR}/${REFERENCE_NAME}/debug/index.out
fi

# Statistics 
echo "ʕ·ᴥ·ʔ :samtools flagstat result" > ${WORK_DIR}/${REFERENCE_NAME}/aligned/alignment_stats.txt
samtools flagstat ${WORK_DIR}/${REFERENCE_NAME}/aligned/sorted_merged_dups_marked.bam  >> ${WORK_DIR}/${REFERENCE_NAME}/aligned/alignment_stats.txt

echo "ʕ·ᴥ·ʔ : samtools stats result " >> ${WORK_DIR}/${REFERENCE_NAME}/aligned/alignment_stats.txt
samtools stats ${WORK_DIR}/${REFERENCE_NAME}/aligned/sorted_merged_dups_marked.bam >> ${WORK_DIR}/${REFERENCE_NAME}/aligned/alignment_stats.txt


REFERENCE_NAME='Wuhan_seafood_market_pneumonia_virus_isolate_Wuhan_Hu_1'

echo "label,breadth_of_coverage" > ${WORK_DIR}/${REFERENCE_NAME}/aligned/stats.csv
for f in ${WORK_DIR}/*/aligned/depth_per_base.txt
do
    awk -v fname=$(basename ${f%%/aligned*}) 'BEGIN{count=0; onisland=0}$3>3{if (!onisland){onisland=1; island_start=$2}}$3<=3{if (onisland){island_end=$2; if (island_end-island_start>=50){count=count+island_end-island_start}} onisland=0}END{if (onisland){island_end=$2; if (island_end-island_start>=50){count=count+island_end-island_start}} if (NR==0){NR=1} printf("%s,%0.02f\n", fname, count*100/NR) }' $f >> ${WORK_DIR}/${REFERENCE_NAME}/aligned/stats.csv
done

python ${PIPELINE_DIR}/compile_results.py ${TOP_DIR} ${WORK_DIR}/${REFERENCE_NAME}/aligned/coverage_1.tsv ${WORK_DIR}/${REFERENCE_NAME}/aligned/coverage_0.tsv ${WORK_DIR}/${REFERENCE_NAME}/aligned/stats.csv ${WORK_DIR}/result.csv 


echo "ʕ·ᴥ·ʔ : Pipeline completed, check ${WORK_DIR} for diagnositc result"
