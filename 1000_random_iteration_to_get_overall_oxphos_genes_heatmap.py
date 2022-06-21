import matplotlib
matplotlib.use('SVG')
import matplotlib.pyplot as plt
import pandas as pd
import sys
import seaborn as sns
import numpy as np
import scipy.stats  as stats
sns.set(font_scale=1.8)
#read residuals data
data = pd.read_csv(sys.argv[1], sep = '\t', index_col = 0)
#Leave only stable gene IDs
data.columns = [x[:15] for x in list(data.columns)]
#Reas sample Info
sample_info = pd.read_csv(sys.argv[2])[['SAMPID', 'SMTSD']]
sample_info.columns = ['sample', 'tissue']
#Read list of OXPHOS-related genes
genes = pd.read_csv(sys.argv[3])
#Focus only on structural OXPHOS genes
oxphosGenes = genes.loc[genes['Type'] != 'Assembly']
#Leave only OXPHOS genes data
data = data[list(oxphosGenes['ENSGnumber'])]
#Assign tissue identity to each sample
data = data.reindex(list(sample_info['sample']))
data['tissue'] = list(sample_info['tissue'])
#Remove NA values
data = data.dropna()
#Leave only data for tissues with 50 samples or more
tissues = [x for x in list(set(list(data['tissue']))) if list(data['tissue']).count(x) >= 50]
data = data.loc[data['tissue'].isin(tissues)]

#Make dictionaries to hold correlation and p-values
corralation_dict = {}
for geneA in data.columns[:-1]:
	corralation_dict[geneA] = {}
	for geneB in data.columns[:-1]:
		corralation_dict[geneA][geneB] = []

pval_dict = {}
for geneA in data.columns[:-1]:
	pval_dict[geneA] = {}
	for geneB in data.columns[:-1]:
		pval_dict[geneA][geneB] = []

#Run 1000 iteration of sampling 50 sample per tissue and calculate gene-gene expression correlation for each iteration
c = 0
while c < 1000:
	df = data.groupby('tissue').sample(n=50)
	for geneA in data.columns[:-1]:
		for geneB in data.columns[:-1]:
			corralation_dict[geneA][geneB].append(stats.spearmanr(df[geneA], df[geneB])[0])
			pval_dict[geneA][geneB].append(stats.spearmanr(df[geneA], df[geneB])[-1])
	c += 1
	print (str(c) + ' random inerations done out of 1000')

#Generate dicts of the median corelations and p-values for all 1000 iterations
median_corralation_dict = {}
median_pval_dict = {}
for geneA in pval_dict:
	median_corralation_dict[geneA] = {}
	median_pval_dict[geneA] = {}
	for geneB in pval_dict:
		median_pval_dict[geneA][geneB] = np.percentile(pval_dict[geneA][geneB], 95)
		if median_pval_dict[geneA][geneB] <= 0.05:
			median_corralation_dict[geneA][geneB] = np.median(corralation_dict[geneA][geneB])
		else:
			median_corralation_dict[geneA][geneB] = 0
#Convert the dicts to dataframes
pval = pd.DataFrame.from_dict(median_pval_dict).fillna(0)[list(oxphosGenes['ENSGnumber'])].reindex(oxphosGenes['ENSGnumber'])
pval.columns = list(oxphosGenes['GeneName'])
pval.index = list(oxphosGenes['GeneName'])

rho = pd.DataFrame.from_dict(median_corralation_dict).fillna(0)[list(oxphosGenes['ENSGnumber'])].reindex(oxphosGenes['ENSGnumber'])
rho.columns = list(oxphosGenes['GeneName'])
rho.index = list(oxphosGenes['GeneName'])
#Save the resulting dataframes
pval.to_csv(sys.argv[4])
rho.to_csv(sys.argv[5])
#color mtDNA and nDNA genes differently
row_colors = []
for i in genes['Genome']:
	if i == 'mtDNA':
		row_colors.append('orange')
	else:
		row_colors.append('green')

#Plot heatmaps
n=sns.clustermap(rho, figsize=(45, 45), cmap = 'coolwarm', center = 0, vmin = -1, vmax = 1)
l=sns.clustermap(pval, figsize=(45, 45))
plt.setp(n.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
plt.setp(n.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
plt.setp(l.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
plt.setp(l.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
figu=l.savefig(sys.argv[6])
fig=n.savefig(sys.argv[7])

