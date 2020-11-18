
touch compiled_results.csv

echo "lib,breadth,viral_reads_mapq1,viral_coverage_mapq1,control_reads_mapq1,control_coverage_mapq1,viral_reads_mapq0,viral_coverage_mapq0,control_reads_mapq0,control_coverage_mapq0" > compiled_results.csv

for file in CTD*; do 
	NEW_ENTRY=$(head -1 $file/polar-bear-fda-eua/result.csv)
	echo $NEW_ENTRY >> compiled_results.csv
done
