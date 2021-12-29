#!/bin/bash -e

##usage: remRecombo <lookupTable> <sam> <qScore> <gbc> <PE> <seqSplit>
## <lookupTable> file path to lookup table (see comment below)
## <sam> file path to sam file (not bam because script uses two cycles through file)
## <qScore> minimum qscore a base change position must have to be accepted (zero will skip test and script will run faster)
## <gbc> number of base change positions required per insert
## <PE> 0 SE, 2 PE sequence
## For SNAQ analysis requiring reads containing both control and native sequences.
## two passes through SAM, first to ID bins for reads, second pass to sort,
## awk Import NT IS lookup table, create recombinant hash table RefName : POS : base, where base is switched between IS and NT
##    Use cigar for each read and look at each position to see if found in hash, if base indicates recombo between NT and IS,
##    place qname:chrom:pos++ in badQname (recombination), or ISQname or NTQname if base reference (a qScore test determines if added base accepterd)
##    After first pass sort into bins based on Qname hashes
##Complexity Control extraction.  If RefName contains CC, hard coded position extracted into a file: SampleID matchQuality/PoolID Sequence.
##		a or b indicates bad CC from pool1 or pool2, other wise 1 or 2 indicates CC sequence passed QC and indicates it source.
##	CC QC uses the bases flanking the Ns to ensure alignment was correct.
##	NT and IS CC-regions counted to allow for future relative yield calculations.
##Clean up: covert sam to bam, index, create fastq, tabulate CC results
## Version:
##	0.1 First release, uses concatenated IS
##	0.2 Add removal of pairs, instead of either R1 or R2
##	0.3 Add single read support, count CC region in NT, IS, and CC, for
##		reads indicate if bad match to pool 1 =a and pool 2 =b, move CC to their own bam
##	0.4 A good IS or NT read requires n base change positions with qscore > x; removed splitSeq boolean.
## 	0.5 Create new ukn bin to handle inserts w/o bc position.  Attempt to address suplemental alignments, 
##		PE chimeras, require good bc position SE or PE seq, SAM file input due to two pass,   
##For Support: tmorrison@accuGenomics.com

pi=$(dirname "${2}")
fp=$(basename "$2")
fe=${fp##*.}
fn=${fp%.*}
if [[ $fe -ne "bam" ]];then
	echo "Input file not bam"
	exit 1
fi
complexityFile="${pi}/${fn}-CC.txt"
goodFile="${pi}/${fn}-good.sam"
isFile="${pi}/${fn}-IS.sam"
badFile="${pi}/${fn}-bad.sam"
ccFile="${pi}/${fn}-cc.sam"
unkFile="${pi}/${fn}-ukn.sam"
hashTallies="${pi}/${fn}-tallies.txt"
tempFile="${pi}/tempFile.sam"

if [[ $6 -eq 1 ]];then
	samtools view -H "${pi}/${fn}.${fe}" > $goodFile
	cp $goodFile $badFile
	cp $goodFile $ccFile
	cp $goodFile $isFile
	cp $goodFile $unkFile
fi
##the awk script uses sam file when two passes needed when splitting sequence
##If no splitting, one can change the script to take a bam file (i.e., don't need to make a sam file) by samtools view -h <bam> | awk '....' $1 -
##Further, if CC counting only, the FNR==1 & FNUM==1 & FNUM==3 can be removed and file input reduced to just awk '....' -
awk -v badFile="$badFile" -v ccFile="$ccFile" -v unkFile="$unkFile" -v goodFile="$goodFile" -v hashTallies="$hashTallies" -v complexityFile="$complexityFile" -v sampleName="$fn" -v isFile="$isFile" -v qCutoff="$3" -v gbc="$4" -v seqSplit="$6" ' 
BEGIN {
	xs="!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
}
FNR==1{
	FNUM++
}
FNUM==1{ #process lookup table
	if(NR==1 && ($1 != "NT_CHROM" || $2 != "NT_POS" || $3 != "NT_REF" || $4 != "IS_CHROM" || $5 != "IS_POS" || $6 != "IS_REF" )) {
		print "ERROR: NT-IS lookup file header incorrect" > "/dev/stderr"
		print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 > "/dev/stderr"
		exit 1
	}
	if($6 ~/[gatc]/) {
		ISs[$1 ":" $2] = toupper($3) #CHROM:POS = good base
		ISs[$4 ":" $5] = toupper($6)
		NTs[$1 ":" $2] = toupper($6) #CHROM:POS = bad base
		NTs[$4 ":" $5] = toupper($3)
	}
	next
}
FNUM==3 && !/^@/ && seqSplit {
	## Capture CC info and sort read into NT, IS, CC, REComination, UNKnown bin
	if($1 ":" $3 ":" $4 in ccQname){
		print $0 >> ccFile #Complexity controls
	} else if($1 ":" $3 ":" $4 in badQname ){ #IS or NT read with a base change indicating recombination
		print $0 >> badFile #this will ensure R2, which could be recombo, does not end up in good bam
	} else if ($7 != "="  && $7 != $3){
		print $0 >> badFile #this read pair is a chimera is a recobinant assumes NT is single reference genome contig
	} else if ( (ISQname[$1 ":" $3 ":" $4] + ISQname[$1 ":" $7 ":" $8]) >= gbc) { #sum of base change positions in read pairs
		print $0 >> isFile #IS sequence with at least gbc number of base change postions
	} else if ( (NTQname[$1 ":" $3 ":" $4] + NTQname[$1 ":" $7 ":" $8]) >= gbc) {
		print $0 >> goodFile #NT sequences with at least gbc number of base change positions
	} else {
		print $0 >> unkFile #read does not have base change position with sufficient qscore
	}
	next
}
!/^@/ {
	if (seqSplit || $3 ~ /-CC$/ ) {
		tmp=$10
		tmpq=$11
		cmp=$6
		n=patsplit(cmp,a,/[0-9]*[MSIDH]/) 
		pnt=0
		for(i=1;i<=n;i++){
			cd=substr(a[i],length(a[i]),1) 
			nm=substr(a[i],1,length(a[i])-1) 
			if(cd=="M" || cd=="=" || cd=="X"){
				pnt=pnt+nm
			}
			else if (cd=="D" || cd=="N"){
				tmp=substr(tmp,1,pnt) substr(xs,1,nm) substr(tmp,pnt+1,length(tmp)) 
				if (qCutoff) {tmpq=substr(tmpq,1,pnt) substr(xs,1,nm) substr(tmpq,pnt+1,length(tmpq))}
				pnt=pnt+nm 
			}
			else if (cd=="I"){
		                tmp=substr(tmp,1,pnt) substr(tmp,pnt+nm+1,length(tmp)-nm-pnt )
				if (qCutoff) {tmpq=substr(tmpq,1,pnt) substr(tmpq,pnt+nm+1,length(tmpq)-nm-pnt )}
			}	
	       	        else if (cd=="S"){
		                tmp=substr(tmp,1,pnt) substr(tmp,pnt+nm+1,length(tmp)-nm-pnt )
				if (qCutoff) {tmpq=substr(tmpq,1,pnt) substr(tmpq,pnt+nm+1,length(tmpq)-nm-pnt )} 
	       	        }
			else if (cd=="H" || cd=="P"){
			}
			else{
				print $1 " cigar Term " cd " not recognized" > "/dev/stderr"
				exit 1
			}
		}
	}
	if ($3 ~ /-CC$/ || $3 ~ /-CC.[0-9]*$/ ) { #put && ! ntf if eliminating recombinant CC
		if(!($1 ":" $3 ":" $4 in ccQname)) {
			t2 = 655 - $4 + 1
			if (t2 >0 && t2 + 9 < length(tmp) && !($1 in CCc2)) { #ensures no double counting
				CCc2[$1]=0
				ccQname[$1 ":" $3 ":" $4]=0
				if ($7 ~ !/[\*=]/) {ccQname[$1 ":" $7 ":" $8]=0} #captured paired read too
				t1=substr(tmp,t2, 10)
				if(t1 ~/^T[AGTC]{8}G/) {
					printf("%s\t%s\t%s\n", sampleName, "2", t1) > complexityFile
				} else {
					printf("%s\t%s\t%s\n", sampleName, "b", t1) > complexityFile
				}
			} else {
				t2=433 - $4 +1 
				if (t2 >0 && t2 + 9 < length(tmp) && !($1 in CCc1)) {
					CCc1[$1]=0
					ccQname[$1 ":" $3 ":" $4]=0
					if ($7 ~ !/[\*=]/) {ccQname[$1 ":" $7 ":" $8]=0} #capture paired read too
					t1=substr(tmp,t2,10)
					if (t1 ~/G[AGTC]{8}T/) {
						printf("%s\t%s\t%s\n", sampleName, "1", t1) > complexityFile 
                        		} else {
						printf("%s\t%s\t%s\n", sampleName, "a", t1) > complexityFile
					}
				}
			}
		}
	} else {
		if ($3 == "nCov-19-5110-7127_20ABYFCP.0-SNAQ" && $4 <= 443 && ($4 + length($10)) >= 452 && !($1 in ISc1)) {
                        ISc1[$1]=0 #cc region
                } else if ($3 == "nCov-19-5110-7127_20ABYFCP.0-SNAQ" && $4 <= 665 && ($4 + length($10)) >= 674 && !($1 in ISc2)) {
                        ISc2[$1]=0
                } else if ($3=="MN908947.3" && $4 <= 5532 && ($4 + length($10)) >= 5541 && !($1 in NTc1)) {
                        NTc1[$1]=0
                } else if ($3=="MN908947.3" && $4 <= 5754 && ($4 + length($10)) >= 5763 && !($1 in NTc2)) {
                        NTc2[$1]=0
		} 
		if (seqSplit) {
			ntf=0
			qPass=1 #if qCutoff=0, this will bipass qscore test
			n=split(tmp,a,"")
			if (qCutoff) {split(tmpq,aq,"")} #only needed if qScore test is implemented
			for (i=1;i<=n;i++){
				if (qCutoff) { 
					if (aq[i] >= qCutoff) {
						qPass=1
					} else {
						qPass=0
					} 
				}
				if ($3 ":" $4+i-1 in NTs && qPass){ #Basechange position with passing qscore otherwise not interested in POS
					if (NTs[$3 ":" $4+i-1] == a[i]){ #base match indicates recombination
						badQname[$1 ":" $3 ":" $4]++ #QNAME : CHROM : POS -- each aligned read is tracked
						recPos[$3 ":" $4+i-1 ":" NTs[$3 ":" $4+i-1]]++ #track recombinant position (possible future use)
						if ($7 ~ !/[\*=]/) {badQname[$1 ":" $7 ":" $8]++} #paired end marked as possible recombinant too
					} else if ($3 ~ /-SNAQ$/ && ISs[$3 ":" $4+i-1] == a[i]) { #count expected ref base
						ISQname[$1 ":" $3 ":" $4]++
					} else if (ISs[$3 ":" $4+i-1] == a[i]){
						NTQname[$1 ":" $3 ":" $4]++
					}
				}
        		        
        		}
		}
	}
}
END {
print length(NTc1) "\t" length(ISc1) "\t" length(CCc1) "\t" length(NTc2) "\t" length(ISc2) "\t" length(CCc2) > complexityFile 
	for (key in recPos){
		print key "\t" recPos[key] > hashTallies #position and number of events where recombination detected...not used 
	}
}' "$1" "${pi}/${fn}.${fe}" "${pi}/${fn}.${fe}" 

if [ $6 != 0 ];then #were split files created?
	echo "Creating BAM..."
	samtools view -b1h  "${goodFile}" -o "${pi}/${fn}-good.bam"
	samtools index "${pi}/${fn}-good.bam"
	rm -f "${goodFile}"
	samtools view -b1h "${pi}/${fn}-IS.sam" -o "${pi}/${fn}-IS.bam"
	samtools index "${pi}/${fn}-IS.bam"
	rm -f "${pi}/${fn}-IS.sam"
	samtools view -b1h "${pi}/${fn}-bad.sam" -o "${pi}/${fn}-bad.bam"
	samtools index "${pi}/${fn}-bad.bam"
	rm -f "${pi}/${fn}-bad.sam"
	samtools view -b1h "${pi}/${fn}-cc.sam" -o "${pi}/${fn}-cc.bam"
	samtools index "${pi}/${fn}-cc.bam"
	rm -f "${pi}/${fn}-cc.sam"
	samtools  view -b1h "${pi}/${fn}-ukn.sam" -o "${pi}/${fn}-ukn.bam"
	samtools index "${pi}/${fn}-ukn.bam"
	rm -f "${pi}/${fn}-ukn.sam"
	
	echo "Creating fastq..."
	if [[ $5 == 1 ]];then
		if [ -f "${pi}/${fn}-good.bam" ];then	
			samtools sort -n "${pi}/${fn}-good.bam" -o "${pi}/${fn}-rsort.temp" -O bam -@ 32 -m 1G >/dev/null 2>&1
			samtools fastq --verbosity 3 -1 "${pi}/${fn}-NT-R1.fastq.gz" -2 "${pi}/${fn}-NT-R2.fastq.gz" "${pi}/${fn}-rsort.temp" &>/dev/null
		fi
	        if [ -f "${pi}/${fn}-IS.bam" ];then
			samtools sort -n "${pi}/${fn}-IS.bam" -o "${pi}/${fn}-rsort.temp" -O bam -@ 32 -m 1G >/dev/null 2>&1
			samtools fastq --verbosity 3 -1 "${pi}/${fn}-IS-R1.fastq.gz" -2 "${pi}/${fn}-IS-R2.fastq.gz" "${pi}/${fn}-rsort.temp" &>/dev/null
		fi
	        if [ -f "${pi}/${fn}-bad.bam" ];then
			samtools sort -n "${pi}/${fn}-bad.bam" -o "${pi}/${fn}-rsort.temp" -O bam -@ 32 -m 1G >/dev/null 2>&1
			samtools fastq --verbosity 3 -1 "${pi}/${fn}-bad-R1.fastq.gz" -2 "${pi}/${fn}-bad-R2.fastq.gz" "${pi}/${fn}-rsort.temp" &>/dev/null
		fi
		if [ -f "${pi}/${fn}-CC.bam" ];then
			samtools sort -n "${pi}/${fn}-CC.bam" -o "${pi}/${fn}-rsort.temp" -O bam -@ 32 -m 1G >/dev/null 2>&1
			samtools fastq --verbosity 3 -1 "${pi}/${fn}-CC-R1.fastq.gz" -2 "${pi}/${fn}-CC-R2.fastq.gz" "${pi}/${fn}-rsort.temp" &>/dev/null
		fi
	else
	     if [ -f "${pi}/${fn}-good.bam" ];then
		 samtools fastq -0 "${pi}/${fn}-NT-R1.fastq.gz" "${pi}/${fn}-good.bam"
	     fi
	     if [ -f "${pi}/${fn}-IS.bam" ];then
		 samtools fastq -0 "${pi}/${fn}-IS-R1.fastq.gz" "${pi}/${fn}-IS.bam"
	     fi
	     if [ -f "${pi}/${fn}-bad.bam" ];then
		 samtools fastq -0 "${pi}/${fn}-bad-R1.fastq.gz" "${pi}/${fn}-bad.bam"
		 fi
		 if [ -f "${pi}/${fn}-CC.bam" ];then
		 samtools fastq -0 "${pi}/${fn}-CC-R1.fastq.gz" "${pi}/${fn}-CC.bam"
		 fi
	fi
	
	echo "Removing scratch files..."
	#rm -f "${pi}/${fn}-IS.sam"
	#rm -f "${pi}/${fn}-good.sam"
	#rm -f "${pi}/${fn}-bad.sam"
	#rm -f "${pi}/${fn}-cc.sam"
	rm -f "${pi}/tempFile.sam"
fi


#if [[$(stat -c%s $complexityFile) > 0 ]];then

## Read through <prefix>-CC.txt
##   Create count of each unique <ID><type><seq>
##   Last row of <prefix>-CC.txt has NT1 IS1 NT2 IS2 totCC_reads
awk -v outFile1="${pi}/${fn}-CC-counts.txt" -v outFile2="${2}_CC-read-counts.txt" 'FNR==NR{
	last++;next
}
FNR != last {
	uSEQ[$1 "\t" $2 "\t" $3]++
}
FNR == last {
	coverage=$1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6	
}
END {
	print "SampleID\tTarget\tSequence\tCount\tNT_P1\tIS_P1\tCC_P1\tNT_P2\tIS_P2\tCC_P2" > outFile1
	print "SampleID\tTarget\tUniqueCount\tNT_P1\tIS_P1\tCC_P1\tNT_P2\tIS_P2\tCC_P2" > outFile2
	for (key in uSEQ) {
		print key "\t" uSEQ[key] "\t" coverage > outFile1
		split(key,arr,"\t")
		cUniq[arr[1]"\t"arr[2]]++
	}
	
	for (key in cUniq) {
		print key "\t" cUniq[key] "\t" coverage > outFile2
	}
}' "$complexityFile" "$complexityFile"

#rm -f "$complexityFile"

