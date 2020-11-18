dir=$1

SARSCOV2_ALGN_STATS=$dir"/polar-bear-fda-eua/Wuhan_seafood_market_pneumonia_virus_isolate_Wuhan_Hu_1/aligned/alignment_stats.txt"
BREADTH_STAT=$dir"/polar-bear-fda-eua/Wuhan_seafood_market_pneumonia_virus_isolate_Wuhan_Hu_1/aligned/stats.csv"


grep "reads paired" $SARSCOV2_ALGN_STATS | awk '{print $4/2}'
grep "mapped (" $SARSCOV2_ALGN_STATS | head -1 | awk -F ":" '{print $1}' | awk -F "(" '{print $2}' | awk '{print $1}' 
echo "NA"
grep "Wuhan_seafood_market_pneumonia_virus_isolate_Wuhan_Hu_1" $BREADTH_STAT | awk -F "," '{print $2}' 
echo "NA"
grep "average quality" $SARSCOV2_ALGN_STATS | awk '{print $4}'
grep "reads duplicated" $SARSCOV2_ALGN_STATS | awk '{print $4}'
grep "insert size average" $SARSCOV2_ALGN_STATS | awk '{print $5}'
grep "insert size standard deviation" $SARSCOV2_ALGN_STATS | awk '{print $6}'
