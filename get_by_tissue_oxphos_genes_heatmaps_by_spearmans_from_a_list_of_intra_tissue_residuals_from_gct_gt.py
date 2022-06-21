import matplotlib
matplotlib.use('SVG')
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
import seaborn as sns
import numpy as np
import scipy.stats  as stats

genes = pd.read_csv(sys.argv[2])
outputDir = str(sys.argv[3])
oxphosGenes = genes.loc[genes['Type'] != 'Assembly']

for m in os.listdir(sys.argv[1]):
	print(m.split('.')[0])
	data = pd.read_csv(sys.argv[1] + m, sep = '\t', index_col = 0)
	data.columns = [str(x)[:15] for x in list(data.columns)]
	data = data[list(oxphosGenes['ENSGnumber'])].dropna()
	if len(data.index) >=50:
		row_colors = []
		sns.set(font_scale=1.8)
		for i in genes['Genome']:
			if i == 'mtDNA':
				row_colors.append('orange')
			else:
				row_colors.append('green')

		correlation_dict = {}
		for geneA in data.columns:
			correlation_dict[geneA] = {}
			for geneB in data.columns:
				correlation_dict[geneA][geneB] = stats.spearmanr(data[geneA], data[geneB])[0]

		rho = pd.DataFrame.from_dict(correlation_dict)[list(oxphosGenes['ENSGnumber'])].reindex(list(oxphosGenes['ENSGnumber']))
		rho.columns = list(oxphosGenes['GeneName'])
		rho.index = list(oxphosGenes['GeneName'])
		rho.to_csv(outputDir + str(m.split('.')[0] + '_oxphos_genes_correlation_matrix.csv'))
		n=sns.clustermap(rho, figsize=(45, 50), cmap = 'coolwarm', center = 0, vmin = -1, vmax = 1, row_colors = row_colors)
		plt.setp(n.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
		plt.setp(n.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
		n.fig.suptitle(str(m.split('.')[0]), size = 40)
		figu=n.savefig(outputDir + str(m.split('.')[0] + '_oxphos_genes_correlation_clustermap.svg'))



