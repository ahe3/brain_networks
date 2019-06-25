for i in "BRNACC" "BRNAMY" "BRNCDT" "BRNCTXB24" "BRNCTXBA9" "BRNHIP" "BRNHYP" "BRNPUT" "BRNSNA"; do
sbatch filter_fdr.sh $i;
done 
