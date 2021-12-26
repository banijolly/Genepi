import sys,os

find_k="NC_045512.2"

inDir = sys.argv[1]


listOfFiles = os.listdir(inDir)
for l in listOfFiles:
	ID = l.split('.')[0]
	rep = ID
	outFile = ID + ".fa"
	fasta = open(l,"rt")
	fasta_data = fasta.read()
	fasta_data = fasta_data.replace(find_k,rep)
	fout = open(outFile, "wt")
	fout.write(fasta_data)
	fout.close()
