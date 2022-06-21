import pandas as pd
import sys
#Read residuals data
redisals = pd.read_csv(sys.argv[1], sep = '\t', index_col = 0)
#Reas sample Info
sample_info = pd.read_csv(sys.argv[2])[['SAMPID', 'SMTSD']]
sample_info.columns = ['sample', 'tissue']
#Define ouput directory
output_dir = str(sys.argv[3])
#Leave on tissues with 50 or over samples
tissues = list(set([x for x in list(sample_info['tissue']) if list(sample_info['tissue']).count(x) >= 50]))
sample_info = sample_info.loc[sample_info['tissue'].isin(tissues)]
redisals = redisals.loc[redisals.index.isin(list(sample_info['sample']))].reindex(list(sample_info['sample']))
redisals['tissue'] = list(sample_info['tissue'])
#Generate per tissue residuals list in output directory
for tissue in list(set(redisals['tissue'])):
	redisals.loc[redisals['tissue'] == tissue][list(redisals.columns[:-1])].to_csv(output_dir + str(tissue) + '.txt.gz', sep = '\t', compression = 'gzip')


