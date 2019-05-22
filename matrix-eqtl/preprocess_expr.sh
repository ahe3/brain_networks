if [[ ! -d data ]]; then
	mkdir data
	for i in "BRNACC" "BRNAMY" "BRNCDT" "BRNCTXB24" "BRNCTXBA9" "BRNHIP" "BRNHYP" "BRNPUT" "BRNSNA"; do
		mkdir data/$i;
	done	
fi

#python filter_expr_samples.py
#Rscript qnormalize_expr.R
#python split_expr.py
