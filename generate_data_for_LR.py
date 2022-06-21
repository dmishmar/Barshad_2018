import pandas as pd
import sys
#Read Normalized data
data = pd.read_csv(sys.argv[1], sep = '\t', index_col = 0)
#Read donor and sample info
donorInfo = pd.read_csv(sys.argv[2])[['SUBJID', 'GENDER', 'AGE', 'DTHHRDY']]
sampleInfo = pd.read_csv(sys.argv[3])[['SAMPID', 'SMTSD']]

#Filter out sample withoud info assignment
samples_with_tissue_assingment = [x for x in list(data.columns) if x in list(sampleInfo['SAMPID'])]
samples_to_keep = [x for x in samples_with_tissue_assingment if x.split('-')[0] + '-' + x.split('-')[1] in list(donorInfo['SUBJID'])]
lmData = data[samples_to_keep].transpose()

#Assign age, gender, COD and tissu data to each sample
age = []
gender = []
cod = []
tissue = []
c = 0
for sample in samples_to_keep:
	donorInfo_sample = donorInfo.loc[donorInfo['SUBJID'] == sample.split('-')[0] + '-' + sample.split('-')[1]]
	sampleInfo_sample = sampleInfo.loc[sampleInfo['SAMPID'] == sample]
	age.append(list(donorInfo_sample['AGE'])[0])
	gender.append(list(donorInfo_sample['GENDER'])[0])
	cod.append(list(donorInfo_sample['DTHHRDY'])[0])
	tissue.append(list(sampleInfo_sample['SMTSD'])[0])
	c += 1
	print (str(c) + ' out of ' + str(len(samples_to_keep)))
lmData['age'] = age
lmData['gender'] = gender
lmData['cod'] = cod
lmData['tissue'] = tissue
#Save the data
lmData.to_csv(sys.argv[4], sep = '\t', compression = 'gzip')

