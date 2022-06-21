import pandas as pd
import sys
#Read the RPKM data
rpkm = pd.read_csv(sys.argv[1], sep = '\t', skiprows = [0,1], index_col = 0)
#Remove description row and transpose
rpkm = rpkm.drop('Description', axis=1).transpose()
#Make an empty output DF
GNrpkm = pd.DataFrame(index = rpkm.index)
#Normalize RPKM data to the sample-wide sum
c = 0
for gene in rpkm.columns:
	if sum(list(rpkm[gene])) > 0:
		GNrpkm[gene] = rpkm[gene] * (1000000 / (sum(list(rpkm[gene]))))
	else:
		GNrpkm[gene] = rpkm[gene]
	c += 1
	print (str(c) + ' out of ' + str(len(rpkm.columns)))
#Save normalized expression data
GNrpkm.transpose().to_csv(sys.argv[2], sep = '\t', compression = 'gzip')
