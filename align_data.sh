#!/bin/bash
### Script to align sample to viral genomes
### Eventually will do entire pipline

# Variables: do we want to be able to set betacorona ref dir?
TOP_DIR=$(pwd)
BETACORONA_REF_DIR="/gpfs0/work/brian/references/betacoronaviruses/*/*.fasta"
BETACORONA_SMALL="/gpfs0/work/brian/references/betacoronaviruses/wuhan*/*.fasta /gpfs0/work/brian/references/betacoronaviruses/sars*/*.fasta /gpfs0/work/brian/references/betacoronaviruses/severe_acute_respiratory_syndrome_coronavirus_2_strain_USA_WA1/*.fasta"
FASTQ_DIR=${TOP_DIR}"/fastq/*_R*.fastq*"
READ1_STR="_R1"
READ2_STR="_R2"

usageHelp="Usage: ${0##*/} [-d TOP_DIR] -jrh"
dirHelp="* [TOP_DIR] is the top level directory (default\n  \"$TOP_DIR\")\n     [TOP_DIR]/fastq must contain the fastq files"
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
    read1filescomma+=$name1$ext","
    read2filescomma+=$name2$ext","
done

# replace commas with spaces for iteration, put in array
read1files=($(echo $read1filescomma | tr ',' ' '))
read2files=($(echo $read2filescomma | tr ',' ' '))

threads=8
threadstring="-t \$SLURM_JOB_CPUS_PER_NODE"

if [ -n "$reducedSet" ]
then
    REFERENCES=$BETACORONA_SMALL
else
    REFERENCES=$BETACORONA_REF_DIR
fi

for REFERENCE in $REFERENCES
do

    ######################################################################
    ######################################################################
    ##########Step #1: Align 
    ######################################################################
    ######################################################################

    REFERENCE_NAME=$(echo $REFERENCE | sed 's:.*/::' | rev | cut -c7- | rev )
    echo -e "(-: Aligning files matching $FASTQ_DIR\n to genome $REFERENCE_NAME"

    if ! mkdir ${TOP_DIR}/${REFERENCE_NAME}; then echo "***! Unable to create ${TOP_DIR}/${REFERENCE_NAME}! Exiting"; exit 1; fi
    if ! mkdir ${TOP_DIR}/${REFERENCE_NAME}/aligned; then echo "***! Unable to create ${TOP_DIR}/${REFERENCE_NAME}_aligned! Exiting"; exit 1; fi
    if ! mkdir ${TOP_DIR}/${REFERENCE_NAME}/debug; then echo "***! Unable to create ${TOP_DIR}/${REFERENCE_NAME}/debug! Exiting"; exit 1; fi
    errorfile=${REFERENCE_NAME}/debug/alignerror

    for ((i = 0; i < ${#read1files[@]}; ++i)); do
        usegzip=0
        file1=${read1files[$i]}
        file2=${read2files[$i]}

	FILE=$(basename ${file1%$read1str})
	ALIGNED_FILE=${TOP_DIR}/${REFERENCE_NAME}/aligned/${FILE}"_mapped.sam"

	dependsort="afterok"

        # Align reads
	jid=`sbatch <<- ALGNR | egrep -o -e "\b[0-9]+$"
		#!/bin/bash -l
		#SBATCH -p commons
		#SBATCH -o ${TOP_DIR}/${REFERENCE_NAME}/debug/align-%j.out
		#SBATCH -e ${TOP_DIR}/${REFERENCE_NAME}/debug/align-%j.err
		#SBATCH -t 2880
		#SBATCH -n 1
		#SBATCH -c $threads
		#SBATCH --mem-per-cpu=10G
                #SBATCH -J "align_${FILE}"
		#SBATCH --threads-per-core=1

		spack load bwa@0.7.17 arch=\`spack arch\`
		bwa 2>&1 | awk '\\\$1=="Version:"{printf(" BWA %s; ", \\\$2)}' 
		echo "Running command bwa mem $threadstring $REFERENCE $file1 $file2 > $ALIGNED_FILE"
		srun --ntasks=1 bwa mem $threadstring $REFERENCE $file1 $file2 > $ALIGNED_FILE
		if [ \$? -ne 0 ]                      
		then  
			touch $errorfile            
			exit 1                                         
		else  
			echo "(-: Mem align of $name$ext.sam done successfully"             
		fi
ALGNR`
	dependalign="afterok:$jid"

        ######################################################################
        ######################################################################
        ##########Step #2: Sort SAMs                                       
        ######################################################################
        ######################################################################
	jid=`sbatch <<- SORTSAM | egrep -o -e "\b[0-9]+$"
		#!/bin/bash -l
		#SBATCH -p commons
		#SBATCH -o ${TOP_DIR}/${REFERENCE_NAME}/debug/sortsam-%j.out
		#SBATCH -e ${TOP_DIR}/${REFERENCE_NAME}/debug/sortsam-%j.err
		#SBATCH -t 2880 
		#SBATCH -n 1
		#SBATCH -c 8
		#SBATCH --mem-per-cpu=10G
		#SBATCH --threads-per-core=1
		#SBATCH -d $dependalign 
		samtools sort -m 4G -@ 8 $ALIGNED_FILE -o ${ALIGNED_FILE}"_sorted.bam"
SORTSAM`
	dependsort="${dependsort}:$jid"
    done


    ######################################################################
    ######################################################################
    ##########Step #3: Merge sorted SAMs into a BAM, get stats, and index
    ######################################################################
    ######################################################################
    jid=`sbatch <<- MERGESAM | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p commons
	#SBATCH -o ${TOP_DIR}/${REFERENCE_NAME}/debug/mergesam-%j.out
	#SBATCH -e ${TOP_DIR}/${REFERENCE_NAME}/debug/mergesam-%j.err
	#SBATCH -t 2880 
	#SBATCH -n 1 
	#SBATCH -c 1
	#SBATCH --mem-per-cpu=10G
	#SBATCH --threads-per-core=1 
	#SBATCH -d $dependsort

	if samtools merge ${TOP_DIR}/${REFERENCE_NAME}/aligned/sorted_merged.bam ${TOP_DIR}/${REFERENCE_NAME}/aligned/*_sorted.bam
	then
		rm ${TOP_DIR}/${REFERENCE_NAME}/aligned/*_sorted.bam
		rm ${TOP_DIR}/${REFERENCE_NAME}/aligned/*.sam
	fi

	if [ ${REFERENCE_NAME} = "severe_acute_respiratory_syndrome_coronavirus_2_strain_USA_WA1" ]; then
    	samtools depth ${TOP_DIR}/${REFERENCE_NAME}/aligned/sorted_merged.bam > ${TOP_DIR}/${REFERENCE_NAME}/aligned/depth_per_base.txt
	fi
	
MERGESAM`

    dependmerge="afterok:$jid"

    if [ -n "$produceIndex" ]
    then
	jid=`sbatch <<- INDEXSAM | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p commons
	#SBATCH -o ${TOP_DIR}/${REFERENCE_NAME}/debug/indexsam-%j.out
	#SBATCH -e ${TOP_DIR}/${REFERENCE_NAME}/debug/indexsam-%j.err
	#SBATCH -t 2880 
	#SBATCH -n 1 
	#SBATCH -c 1
	#SBATCH --mem-per-cpu=10G
	#SBATCH --threads-per-core=1 
	#SBATCH -d $dependmerge

	samtools index ${TOP_DIR}/${REFERENCE_NAME}/aligned/sorted_merged.bam

INDEXSAM`
    fi

    dependstats="afterok"
    jid=`sbatch <<- SAMSTATS | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p commons
	#SBATCH -o ${TOP_DIR}/${REFERENCE_NAME}/debug/samstats-%j.out
	#SBATCH -e ${TOP_DIR}/${REFERENCE_NAME}/debug/samstats-%j.err
	#SBATCH -t 2880 
	#SBATCH -n 1 
	#SBATCH -c 1
	#SBATCH --mem-per-cpu=10G
	#SBATCH --threads-per-core=1 
	#SBATCH -d $dependmerge

	samtools flagstat ${TOP_DIR}/${REFERENCE_NAME}/aligned/sorted_merged.bam > ${TOP_DIR}/${REFERENCE_NAME}/aligned/stats.txt

SAMSTATS`

    dependstats="${dependstats}:$jid"
done

echo "#!/bin/bash -l" > $TOP_DIR/collect_stats.sh
echo "#SBATCH -p commons" >> $TOP_DIR/collect_stats.sh
echo "#SBATCH -o ${TOP_DIR}/collectstats-%j.out"  >> $TOP_DIR/collect_stats.sh
echo "#SBATCH -e ${TOP_DIR}/collectstats-%j.err" >> $TOP_DIR/collect_stats.sh
echo "#SBATCH -t 30" >> $TOP_DIR/collect_stats.sh
echo "#SBATCH -n 1 " >> $TOP_DIR/collect_stats.sh
echo "#SBATCH -c 1" >> $TOP_DIR/collect_stats.sh
echo "#SBATCH --mem=200" >> $TOP_DIR/collect_stats.sh
echo "#SBATCH --threads-per-core=1 " >> $TOP_DIR/collect_stats.sh
echo "#SBATCH -d $dependstats"  >> $TOP_DIR/collect_stats.sh 
echo "echo \"<table>\" > $TOP_DIR/stats.html " >> $TOP_DIR/collect_stats.sh
echo "for f in $TOP_DIR/*/aligned/stats.txt; do"  >> $TOP_DIR/collect_stats.sh
echo  "awk -v fname=\${f%%aligned*} '\$4==\"mapped\"{split(\$5,a,\"(\"); print \"<tr><td> \"fname\" </td>\", \"<td> \"a[2]\" </td></tr>\"}' \$f >> ${TOP_DIR}/stats.html"  >> $TOP_DIR/collect_stats.sh 
echo "	done "  >> $TOP_DIR/collect_stats.sh
echo "echo \"</table>\" >> $TOP_DIR/stats.html " >> $TOP_DIR/collect_stats.sh

sbatch < $TOP_DIR/collect_stats.sh

######################################################################
######################################################################
##########Step #4: Produce contigs - can happen in parallel
######################################################################
######################################################################

jid=`sbatch <<- CONTIG | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p x86
	#SBATCH -o ${TOP_DIR}/contig-%j.out
	#SBATCH -e ${TOP_DIR}/contig-%j.err
	#SBATCH -t 600
	#SBATCH -n 1 
	#SBATCH -c 1
	#SBATCH --mem-per-cpu=10G
	#SBATCH --threads-per-core=1
	#SBATCH -d $dependmerge

	echo "Running /gpfs0/work/brian/scripts/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 $read1filescomma -2 $read2filescomma -o ${TOP_DIR}/contigs"
	/gpfs0/work/brian/scripts/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 $read1filescomma -2 $read2filescomma -o ${TOP_DIR}/contigs

	mv contig-* contigs/
CONTIG`

# need to wait for alignment, though this is waiting on all alignments
# and on sort instead of waiting on just the alignments we need
dependcontig="$dependsort:$jid"

######################################################################
######################################################################
##########Step #5: Dot plot - after contig and alignment
######################################################################
######################################################################

jid=`sbatch <<- DOTPLOT | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH --partition=x86
	#SBATCH -o ${TOP_DIR}/dotplot-%j.out
	#SBATCH -e ${TOP_DIR}/dotplot-%j.err
	#SBATCH -t 600
	#SBATCH -n 1 
	#SBATCH -c 2
	#SBATCH --mem-per-cpu=2G
	#SBATCH --threads-per-core=1 
	#SBATCH -d $dependcontig
	
	minimap2 -x asm5 /gpfs0/work/brian/references/betacoronaviruses/severe_acute_respiratory_syndrome_coronavirus_2_strain_USA_WA1/*.fasta ${TOP_DIR}/contigs/final.contigs.fa > ${TOP_DIR}/contig_nCoV-2019.paf
	
	/gpfs0/apps/x86/anaconda3/bin/python ../COVID19/dot_coverage.py severe_acute_respiratory_syndrome_coronavirus_2_strain_USA_WA1/aligned/depth_per_base.txt ${TOP_DIR}/contig_nCoV-2019.paf dotplot 500 29867 stats.html stats.pdf False

   
DOTPLOT`

echo "(-: Finished adding all jobs... Now is a good time to get that cup of coffee... Last job id $jid"