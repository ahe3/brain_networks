import os

expr_folder = "../wgcna/filtered_expression_by_tissue/"
samples_folder = "../wgcna/samples/"
output_folder = "data/"

for file in os.listdir(expr_folder):
	tissue = file.split('_')[0]
	print("Preprocessing " + tissue)
	samples = []

	samples_file = open(samples_folder + tissue + ".samples.txt", "r")
	for line in samples_file:
		samples.append(line.rstrip())
	samples_file.close()

	expr_file = open(expr_folder + tissue + "_expression.txt", "r")
	myfile = open(output_folder + tissue + "/GE.txt", "w")
	linenum = 1
	indexes = []
	for line in expr_file:
		line = line.rstrip()
		line = line.split('\t')
		if linenum == 1:
			myfile.write("id")

			for i in range(len(line)):
				if line[i] in samples:
					s = line[i].split('.')
					samp = s[0] + "-" + s[1]
					indexes.append(i)
					myfile.write("\t" + samp)
			myfile.write("\n")
			linenum = 2

		else:
			myfile.write(line[0])
			expr = line[1:]
			for i in range(len(expr)):
				if i in indexes:
					myfile.write("\t" + expr[i])
			myfile.write("\n")
	expr_file.close()
	myfile.close()
