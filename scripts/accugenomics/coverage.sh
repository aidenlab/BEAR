#!/bin/bash -e
##coverage </path/rootFileName> <path/to/run/coverage/output/file> <full path to bed file describing amplicon start/stop> <seqSplit>
## <seqSplit> indicates if the sequneces were split into part, false, just count SNAQ and NM
##Program flow: use bamstats to count each amplicon in NT, IS, recombinant and CC bams ($1 contains prefix to bam files).
##	awk script median amplicon reads + sums all amplicon counts for each bam, makes table indicated by $2
## may want to add mapq and length filters for the reads.

fp=$(dirname "$1")
fn=$(basename "$1")
if [[ $4 -ne 0 ]];then
	if [ ! -f "${1}-NT.coverage" ];then
	       samtools bedcov -Q 4 "$3" "${1}-good.bam" > "${1}-NT.coverage"
	       fi
	       
	if [ ! -f "${1}-IS.coverage" ];then
	       samtools bedcov -Q 4 "$3" "${1}-IS.bam" > "${1}-IS.coverage"
	       fi
	
	if [ ! -f "${1}-bad.coverage" ];then
	       samtools bedcov -Q 4 "$3" "${1}-bad.bam" > "${1}-bad.coverage"
	fi
	if [ ! -f "${1}-cc.coverage" ];then
		samtools bedcov -Q 4 "$3" "${1}-cc.bam" > "${1}-cc.coverage"
	fi
	if [ ! -f "${1}-ukn.coverage" ];then
	        samtools bedcov -Q 4 "$3" "${1}-ukn.bam" > "${1}-ukn.coverage"
	fi

	awk -v filestr=$fn -v fo="$2" '
	FILENAME ~ /-NT.coverage$/ {
		if ($1=="MN908947.3") {
			ar = int($9 / ( $3 - $2 ))#average reads per amplicon = total bases / expected fragment lenthg
			na[nc++]=ar;nt+=ar
		} else {next}
	}	
	FILENAME ~ /-IS.coverage$/ {
		if ($1 ~ /-SNAQ$/) {
			ar = int($9/($3 -$2))
			 ia[ic++]=ar;it+=ar
		} else {next}
	}
	FILENAME ~ /-ukn.coverage$/ {
	        if ($1 ~ /-SNAQ$/) {
	                ar = int($9/($3 -$2))
	                 ia[ic++]=ar;it+=ar
	        } else {next}
	}
	FILENAME ~ /-bad.coverage$/ {
		ar = int($9/($3 -$2))
		ba[bc++]=ar;bt+=ar
	}END {printf ("%s\t%i\t%i\t%i\t", filestr, nt, it, bt) >> fo
		asort(na);asort (ia);asort (ib)
		if ((nc % 2) == 1){nm=na[ int(nc/2) ]} else {nm=( na[nc/2] + na[nc/2-1])/2};
		if ((ic % 2) == 1){im=ia[ int(ic/2) ]} else {im=( ia[ic/2] + ia[ic/2-1])/2};
		if ((bc % 2) == 1){bm=ba[ int(bc/2) ]} else {bm=( ba[bc/2] + ba[bc/2-1])/2};
		 printf ("%i\t%i\t%i\n", nm, im, bm ) >> fo
	}' "${1}-NT.coverage" "${1}-IS.coverage" "${1}-bad.coverage" "${1}-ukn.coverage"
	
	#clean up
	#rm -f "${1}-good.bam" "${1}-IS.bam" "${1}-bad.bam" "${1}-cc.bam"
	#rm -f "${1}-NT.coverage" "${1}-IS.coverage" "${1}-bad.coverage" "${1}-cc.coverage"
else
	if [ ! -f "${1}.coverage" ];then
		samtools bedcov -Q 4 "$3" "${1}.bam" > "${1}.coverage"
	fi
	awk -v filestr=$fn -v fo="$2" '{
        if ($1=="MN908947.3") {
	        ar = int($9 / ( $3 - $2 ))#average reads per amplicon = total bases / expected fragment lenthg
                na[nc++]=ar;nt+=ar
        } else if ($1 ~ /-SNAQ$/) {
        	ar = int($9/($3 -$2))
                ia[ic++]=ar;it+=ar
        }} END {
	printf ("%s\t%i\t%i\t%i\t", filestr, nt, it, 0) >> fo
        asort(na);asort (ia);asort (ib)
        if ((nc % 2) == 1){nm=na[ int(nc/2) ]} else {nm=( na[nc/2] + na[nc/2-1])/2};
        if ((ic % 2) == 1){im=ia[ int(ic/2) ]} else {im=( ia[ic/2] + ia[ic/2-1])/2};
        printf ("%i\t%i\t%i\n", nm, im, 0 ) >> fo
	}
	' "${1}.coverage"
fi
