#for i in "BRNACC" "BRNCDT" "BRNCTXB24" "BRNHYP" "BRNCTXBA9" "BRNHIP" "BRNPUT" "BRNSNA"; do
for i in "BRNCTXBA9" "BRNHIP" "BRNPUT" "BRNSNA"; do
#for i in "BRNAMY"; do
sbatch quic_STARS.sh $i;
done
