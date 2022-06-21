import numpy as np
import pandas as pd
import sys
from sklearn.linear_model import LinearRegression
import warnings
warnings.filterwarnings("ignore")
#Read data
df = pd.read_csv(sys.argv[1], sep = '\t', index_col = 0)
#Generate gender dummies
df['gender_dummy'] = df.gender.map({1:0, 2:1})
#Generate COD dummies
cod_dummies = pd.get_dummies(df.cod, prefix='cod_dummy').iloc[:, 1:]
df = pd.concat([df, cod_dummies], axis=1)
#Calculate linear regression model residuals
variables = ['age','gender_dummy','cod_dummy_1.0', 'cod_dummy_2.0', 'cod_dummy_3.0', 'cod_dummy_4.0']
lm = LinearRegression()
residuals = pd.DataFrame(index = df.index)
x = df[variables]
c = 0
for gene in df.columns[:-9]:
	y = df[gene]
	lm.fit(x,y)
	residuals[str(gene)] = df[gene] - lm.predict(x)
	c += 1
	print ('gene ' + str(c) + ' out of ' + str(len(df.columns[:-9])))
#Save the residuals data
residuals.to_csv(sys.argv[2], sep = '\t', compression = 'gzip')
