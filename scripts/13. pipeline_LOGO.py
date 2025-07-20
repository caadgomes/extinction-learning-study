
######### LOGO - python code

import pandas as pd, numpy as np
import os
from sklearn import linear_model
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import LeaveOneGroupOut
from scipy.stats import pearsonr

base_dir = '/home/data/'
learning_dir = os.path.join(base_dir, 'group','learning')

# specify task
task = 'acquisition'

if task == "acquisition":
    cols = ['lCEB_rPFC', 'lHIP_lACC', 'lHIP_lPFC', 'rACC_rPFC', 'rAMY_rACC',
            'rCEB_lHIP', 'rHIP_rACC']
elif task == "extinction":
    cols = ['lAMY_lACC','lAMY_lHIP','lCEB_rACC','lHIP_lPFC','rAMY_rACC','rAMY_rPFC',
            'rHIP_rACC']
elif task == "renewal":
    cols = ['lAMY_lHIP','lPFC_lACC','rPFC_rACC','rPFC_rHIP']

# EC
df_EC = pd.read_csv(os.path.join(learning_dir,'task-AER_desc-EC_df.tsv'), sep='\t')
df_EC2 = df_EC.loc[(df_EC.task==task),:]

x = df_EC2[cols].reset_index(drop=True)
y = df_EC2.loc[:,'learning'].reset_index(drop=True)

groups = list(df_EC2.AG)

lm = linear_model.LinearRegression(fit_intercept=True)
m = lm.fit(x, y)

feature_names = x.columns
model_coefficients = m.coef_

coefficients_df = pd.DataFrame(data = model_coefficients, 
                              index = feature_names, 
                              columns = ['beta'])

logo=LeaveOneGroupOut()
scaler = StandardScaler()
yhat = np.zeros(len(y))
for train, test in logo.split(x, y, groups):
    x_std = scaler.fit_transform(x.loc[train,:])
    y_std = scaler.fit_transform(np.array(y[train]).reshape(-1, 1))
    m = lm.fit(X=x_std, y=y_std.flatten())
    yhat[test] = m.predict(x.values[test, :])
    

df_all = pd.DataFrame({'y':y,'yhat':yhat,'modality':'EC', 'task':task})
 
# FC #
df_FC = pd.read_csv(os.path.join(learning_dir,'task-AER_desc-FC_df.tsv'), sep='\t')
df_FC2 = df_FC.loc[(df_FC.task==task),:]

# no directionality so check if roi pair is not reversed
colsr = ['_'.join(s.split('_')[::-1]) for s in cols]
cols = list(set.intersection(set(cols+colsr), set(df_FC2.columns)))

x = df_FC2[cols].reset_index(drop=True)
y = df_FC2.loc[:,'learning'].reset_index(drop=True)

groups = list(df_FC2.AG)

lm = linear_model.LinearRegression(fit_intercept=True)
m = lm.fit(x, y)

feature_names = x.columns
model_coefficients = m.coef_

coefficients_df = pd.DataFrame(data = model_coefficients, 
                              index = feature_names, 
                              columns = ['beta'])

logo=LeaveOneGroupOut()
scaler = StandardScaler()
yhat = np.zeros(len(y))
for train, test in logo.split(x, y, groups):
    # standardize x[train]
    x_std = scaler.fit_transform(x.loc[train,:])
    y_std = scaler.fit_transform(np.array(y[train]).reshape(-1, 1))
    m = lm.fit(X=x_std, y=y_std.flatten())
    yhat[test] = m.predict(x.values[test, :])
    
df_all = pd.concat([df_all,
                    pd.DataFrame({'y':y,'yhat':yhat,'modality':'FC', 'task':task})])


# SC #
df_SC = pd.read_csv(os.path.join(learning_dir,'task-AER_desc-SC_df.tsv'), sep='\t')
df_SC2 = df_SC.loc[(df_SC.task==task),:]

# no directionality so check if roi pair is not reversed
colsr = ['_'.join(s.split('_')[::-1]) for s in cols]
cols = list(set.intersection(set(cols+colsr), set(df_SC2.columns)))

x = df_SC2[cols].reset_index(drop=True)
# y = df_EC_AC[['learning']].reset_index(drop=True)
y = df_SC2.loc[:,'learning'].reset_index(drop=True)

groups = list(df_SC2.AG)

lm = linear_model.LinearRegression(fit_intercept=True)
m = lm.fit(x, y)

feature_names = x.columns
model_coefficients = m.coef_

coefficients_df = pd.DataFrame(data = model_coefficients, 
                              index = feature_names, 
                              columns = ['beta'])

logo=LeaveOneGroupOut()
scaler = StandardScaler()
yhat = np.zeros(len(y))
for train, test in logo.split(x, y, groups):
    # standardize x[train]
    x_std = scaler.fit_transform(x.loc[train,:])
    y_std = scaler.fit_transform(np.array(y[train]).reshape(-1, 1))
    m = lm.fit(X=x_std, y=y_std.flatten())
    yhat[test] = m.predict(x.values[test, :])
    
df_all = pd.concat([df_all,
                    pd.DataFrame({'y':y,'yhat':yhat,'modality':'SC', 'task':task})])

df_all.to_csv(os.path.join(learning_dir,'task-{}_desc-LOGO_df.tsv'.format(task)), index=False, sep='\t')

