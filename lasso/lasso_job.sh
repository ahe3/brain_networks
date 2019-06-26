#for i in "BRNACC" "BRNAMY" "BRNCDT" "BRNCTXB24" "BRNHYP" "BRNCTXBA9" "BRNHIP" "BRNPUT" "BRNSNA"; do
for i in "BRNAMY"; do
sbatch lasso.sh $i;
done
