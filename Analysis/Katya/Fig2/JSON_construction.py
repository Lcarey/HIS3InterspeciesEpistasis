
# coding: utf-8

# In[1]:

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import os
import distance
import json
pd.set_option('display.max_columns', 100)


# In[2]:

data = {}
clean_data = {}

for f in os.listdir('/home/katya/local/HIS3InterspeciesEpistasis/Data/'):
    if 'csv' in f:
        data[f[:-16]] = pd.DataFrame.from_csv('/home/katya/local/HIS3InterspeciesEpistasis/Data/' + f, sep = '\t')
        clean_data[f[:-16]] = data[f[:-16]][(data[f[:-16]].nonsense == 0) & (data[f[:-16]].middle == 1)]


# In[ ]:

positions = pd.DataFrame.from_csv('/home/katya/local/HIS3InterspeciesEpistasis/Data_Small_Tables/positions.csv', sep = '\t')

for chunk in clean_data:
    print (chunk)
    sqs = list(clean_data[chunk].index)
    
    for ind in range(len(sqs)):
        sqs[ind] = sqs[ind][:int(positions[chunk].ix['len1'])] + sqs[ind][-int(positions[chunk].ix['len2']):]
        
    data = {}
    data['nodes'] = []
    for i in range(len(sqs)):
        data['nodes'].append({'name':sqs[i], 
                              'x':(clean_data[chunk].dist_Scer[i] + np.random.normal(0, .1)), 
                              'y':clean_data[chunk].s[i]*10})
        
    data['connections'] = []
    
    done = []
    for i1 in range(len(sqs)):
        if i1%1000 == 0:
            print (i1)
        done.append(sqs[i1])
        for i2 in range(len(sqs)):
            if sqs[i2] not in done and clean_data[chunk].dist_Scer[i1]!=clean_data[chunk].dist_Scer[i2]             and distance.hamming(sqs[i1],sqs[i2]) == 1:
                if clean_data[chunk].s[i1]>=3 and clean_data[chunk].s[i2]>=3:
                    data['connections'].append({'source':sqs[i1], 'target':sqs[i2], 'color':1})
                else:
                    data['connections'].append({'source':sqs[i1], 'target':sqs[i2], 'color':2})
                
    with open('/home/katya/local/HIS3InterspeciesEpistasis/Analysis/Katya/Fig2/' + chunk + '.json', 'w+') as outfile:  
        json.dump(data, outfile)


# In[ ]:



