for i in "BRNACC" "BRNAMY" "BRNCDT" "BRNCTXB24" "BRNCTXBA9" "BRNHIP" "BRNHYP" "BRNPUT" "BRNSNA"; do
#for i in "BRNAMY" "BRNCTXB24"; do
sbatch filter_results.sh $i;
done 
