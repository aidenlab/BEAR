#!/bin/bash
##### 
### Polar pipeline script for viral diagnostic
### Given paired end sequencing of putative viral data
###   -> Aligns to a selection of potential viruses, sorts and merges
###   -> In parallel, assembles the data into contigs
###   -> After alignment, runs flagstat aligned data to create stats.html
###   -> After contigging, creates dotplot to add to stats.html
###   -> Final report created as PDF from stats.html
#####

### VARIABLES: SET THESE FOR YOUR SYSTEM

## Software commands; add load modules as necessary
# BWA for alignment to different genomes
LOAD_BWA="spack load bwa@0.7.17 arch=\`spack arch\`"
BWA_CMD="bwa"
# Samtools for file manipulation
LOAD_SAMTOOLS=""
SAMTOOLS_CMD="samtools"
# Megahit for contig creation
LOAD_MEGAHIT=""
MEGAHIT_CMD="/gpfs0/work/brian/scripts/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit"
# Minimap2 for creating dotplot
LOAD_MINIMAP2=""
MINIMAP2_CMD="minimap2"
# Python for creating dotplot
LOAD_PYTHON=""
PYTHON_CMD="/gpfs0/apps/x86/anaconda3/bin/python"

## Queues
QUEUE="commons"
QUEUE_X86="x86"

## Threads
threads=8
threadstring="-t \$SLURM_JOB_CPUS_PER_NODE"

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
usageHelp="Usage: ${0##*/} [-d TOP_DIR] -jrh"
dirHelp="* [TOP_DIR] is the top level directory (default \"$TOP_DIR\")\n\
  [TOP_DIR]/fastq must contain the fastq files"
indexHelp="* -j produce index file for aligned files"
reducedSet="* -r reduced set for alignment"
helpHelp="* -h: print this help and exit"

printHelpAndExit() {
    echo -e "$usageHelp"
    echo -e "$dirHelp"
    echo -e "$indexHelp"
    echo -e "$reducedSet"
    echo "$helpHelp"
    exit "$1"
}

while getopts "d:hr" opt; do
    case $opt in
	h) printHelpAndExit 0;;
        d) TOP_DIR=$OPTARG ;;
	j) produceIndex=1 ;;
	r) reducedSet=1 ;;
	[?]) printHelpAndExit 1;;
    esac
done

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

	dependsort="afterok"

        # Align reads
	jid=`sbatch <<- ALGNR | egrep -o -e "\b[0-9]+$"
		#!/bin/bash -l
		#SBATCH -p $QUEUE
		#SBATCH -o ${WORK_DIR}/${REFERENCE_NAME}/debug/align-%j.out
		#SBATCH -e ${WORK_DIR}/${REFERENCE_NAME}/debug/align-%j.err
		#SBATCH -t 2880
		#SBATCH -n 1
		#SBATCH -c $threads
		#SBATCH --mem-per-cpu=4G
                #SBATCH -J "align_${FILE}"
		#SBATCH --threads-per-core=1

		$LOAD_BWA
		$BWA_CMD 2>&1 | awk '\\\$1=="Version:"{printf(" BWA %s; ", \\\$2)}'
		echo "Running command $BWA_CMD mem $threadstring $REFERENCE $file1 $file2 > $ALIGNED_FILE"
		srun --ntasks=1 $BWA_CMD mem $threadstring $REFERENCE $file1 $file2 > $ALIGNED_FILE
		if [ \$? -ne 0 ]                      
		then  
			exit 1                                         
		else  
			echo "(-: Mem align of $name$ext.sam done successfully"             
		fi
ALGNR`
	dependalign="afterok:$jid"

        ######################################################################
        ##########Sort SAMs, convert to BAM         
        ######################################################################
	jid=`sbatch <<- SORTSAM | egrep -o -e "\b[0-9]+$"
		#!/bin/bash -l
		#SBATCH -p $QUEUE
		#SBATCH -o ${WORK_DIR}/${REFERENCE_NAME}/debug/sortsam-%j.out
		#SBATCH -e ${WORK_DIR}/${REFERENCE_NAME}/debug/sortsam-%j.err
		#SBATCH -t 2880 
		#SBATCH -n 1
		#SBATCH -c $threads
		#SBATCH --mem-per-cpu=4G
		#SBATCH --threads-per-core=1
		#SBATCH -d $dependalign

		$LOAD_SAMTOOLS 
		$SAMTOOLS_CMD sort -m 4G -@ $threads $ALIGNED_FILE -o ${ALIGNED_FILE}"_sorted.bam"
SORTSAM`
	dependsort="${dependsort}:$jid"
    done


    ######################################################################
    ##########Merge sorted BAMs, get stats, and index if flag set
    ######################################################################
    jid=`sbatch <<- MERGESAM | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $QUEUE
	#SBATCH -o ${WORK_DIR}/${REFERENCE_NAME}/debug/mergesam-%j.out
	#SBATCH -e ${WORK_DIR}/${REFERENCE_NAME}/debug/mergesam-%j.err
	#SBATCH -t 2880 
	#SBATCH -n 1 
	#SBATCH -c 1
	#SBATCH --mem-per-cpu=4G
	#SBATCH --threads-per-core=1 
	#SBATCH -d $dependsort

	$LOAD_SAMTOOLS
	if $SAMTOOLS_CMD merge ${WORK_DIR}/${REFERENCE_NAME}/aligned/sorted_merged.bam ${WORK_DIR}/${REFERENCE_NAME}/aligned/*_sorted.bam
	then
		rm ${WORK_DIR}/${REFERENCE_NAME}/aligned/*_sorted.bam
		rm ${WORK_DIR}/${REFERENCE_NAME}/aligned/*.sam
	fi

MERGESAM`

    dependmerge="afterok:$jid"

    if [[ "$REFERENCE" == *match* ]]
    then
	dependmatchdone="afterok:$jid"
	matchname=${REFERENCE_NAME}
	matchref=${REFERENCE}
    fi

    if [ -n "$produceIndex" ]
    then
	jid=`sbatch <<- INDEXSAM | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $QUEUE
	#SBATCH -o ${WORK_DIR}/${REFERENCE_NAME}/debug/indexsam-%j.out
	#SBATCH -e ${WORK_DIR}/${REFERENCE_NAME}/debug/indexsam-%j.err
	#SBATCH -t 2880 
	#SBATCH -n 1 
	#SBATCH -c 1
	#SBATCH --mem-per-cpu=10G
	#SBATCH --threads-per-core=1 
	#SBATCH -d $dependmerge

	$LOAD_SAMTOOLS
	$SAMTOOLS_CMD index ${WORK_DIR}/${REFERENCE_NAME}/aligned/sorted_merged.bam

INDEXSAM`
    fi

    dependstats="afterok"
    jid=`sbatch <<- SAMSTATS | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $QUEUE
	#SBATCH -o ${WORK_DIR}/${REFERENCE_NAME}/debug/samstats-%j.out
	#SBATCH -e ${WORK_DIR}/${REFERENCE_NAME}/debug/samstats-%j.err
	#SBATCH -t 2880 
	#SBATCH -n 1 
	#SBATCH -c 1
	#SBATCH --mem-per-cpu=10G
	#SBATCH --threads-per-core=1 
	#SBATCH -d $dependmerge

	$LOAD_SAMTOOLS
	$SAMTOOLS_CMD flagstat ${WORK_DIR}/${REFERENCE_NAME}/aligned/sorted_merged.bam > ${WORK_DIR}/${REFERENCE_NAME}/aligned/stats.txt

SAMSTATS`

    dependstats="${dependstats}:$jid"
done

echo "#!/bin/bash -l" > "$WORK_DIR"/collect_stats.sh
echo "#SBATCH -p $QUEUE" >> "$WORK_DIR"/collect_stats.sh
echo "#SBATCH -o ${LOG_DIR}/collectstats-%j.out"  >> "$WORK_DIR"/collect_stats.sh
echo "#SBATCH -e ${LOG_DIR}/collectstats-%j.err" >> "$WORK_DIR"/collect_stats.sh
echo "#SBATCH -t 30" >> "$WORK_DIR"/collect_stats.sh
echo "#SBATCH -n 1 " >> "$WORK_DIR"/collect_stats.sh
echo "#SBATCH -c 1" >> "$WORK_DIR"/collect_stats.sh
echo "#SBATCH --mem=200" >> "$WORK_DIR"/collect_stats.sh
echo "#SBATCH --threads-per-core=1 " >> "$WORK_DIR"/collect_stats.sh
echo "#SBATCH -d $dependstats"  >> "$WORK_DIR"/collect_stats.sh 
echo "echo \"label,percentage\" > $WORK_DIR/stats.csv " >> "$WORK_DIR"/collect_stats.sh
echo "for f in $WORK_DIR/*/aligned/stats.txt; do"  >> "$WORK_DIR"/collect_stats.sh
echo  "awk -v fname=\$(basename \${f%%/aligned*}) 'BEGIN{OFS=\",\"}\$4==\"mapped\"{split(\$5,a,\"(\"); print fname, a[2]}' \$f >> ${WORK_DIR}/stats.csv"  >> $WORK_DIR/collect_stats.sh 
echo "	done "  >> "$WORK_DIR"/collect_stats.sh

sbatch < "$WORK_DIR"/collect_stats.sh

######################################################################
######################################################################
########## Produce contigs - can happen in parallel
######################################################################
######################################################################

jid=`sbatch <<- CONTIG | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $QUEUE_X86
	#SBATCH -o ${LOG_DIR}/contig-%j.out
	#SBATCH -e ${LOG_DIR}/contig-%j.err
	#SBATCH -t 600
	#SBATCH -n 1 
	#SBATCH -c 1
	#SBATCH --mem-per-cpu=10G
	#SBATCH --threads-per-core=1
	#SBATCH -d $dependmerge

	$LOAD_MEGAHIT
	echo "Running $MEGAHIT_CMD -1 $read1filescomma -2 $read2filescomma -o ${WORK_DIR}/contigs"
	$MEGAHIT_CMD -1 $read1filescomma -2 $read2filescomma -o ${WORK_DIR}/contigs
	mv ${WORK_DIR}/contigs/final.contigs.fa ${FINAL_DIR}/.
CONTIG`

# need to wait for match alignment
dependcontig="${dependmatchdone}:$jid"

######################################################################
######################################################################
########## Dot plot - after contig and alignment
######################################################################
######################################################################

jid=`sbatch <<- DOTPLOT | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH --partition=$QUEUE_X86
	#SBATCH -o ${TOP_DIR}/dotplot-%j.out
	#SBATCH -e ${TOP_DIR}/dotplot-%j.err
	#SBATCH -t 600
	#SBATCH -n 1 
	#SBATCH -c 2
	#SBATCH --mem-per-cpu=2G
	#SBATCH --threads-per-core=1 
	#SBATCH -d $dependcontig
	
	$LOAD_PYTHON
	$LOAD_MINIMAP2
	$LOAD_SAMTOOLS
	$MINIMAP2 -x asm5 $matchref ${FINAL_DIR}/final.contigs.fa > ${FINAL_DIR}/contig_${matchname}.paf
	if [ -s ${FINAL_DIR}/contig_${matchname}.paf ]
	then
    		echo "!*** Pairwise alignment by minimap2 failed."
		exit 1
	fi
	$SAMTOOLS_CMD rmdup ${WORK_DIR}/${matchname}/aligned/sorted_merged.bam ${WORK_DIR}/${matchname}/aligned/sorted_merged_dedup.bam 
    	$SAMTOOLS_CMD depth ${WORK_DIR}/${matchname}/aligned/sorted_merged_dedup.bam > ${WORK_DIR}/${matchname}/aligned/depth_per_base.txt	

	CONTIG_LENGTH=$(tail -n2 ${WORK_DIR}/contigs/log |grep -o 'total.*' | awk '{print \$2}')

	$PYTHON_CMD ${PIPELINE_DIR}/dot_coverage.py ${WORK_DIR}/${matchname}/aligned/depth_per_base.txt ${FINAL_DIR}/contig_${matchname}.paf ${WORK_DIR}/stats.csv $CONTIG_LENGTH ${FINAL_DIR}/report.pdf
DOTPLOT`

echo "(-: Finished adding all jobs... Now is a good time to get that cup of coffee... Last job id $jid"
