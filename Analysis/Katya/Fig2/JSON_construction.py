
# coding: utf-8

# In[78]:

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import os
import distance
from multiprocessing import Pool
import json
from bisect import bisect_left
# %matplotlib inline
pd.set_option('display.max_columns', 100)


# In[201]:

chunks = ['S'+str(i) for i in range(1,13)]
data = {}
clean_data = {}

for f in os.listdir('/home/katya/start/HIS3InterspeciesEpistasis/Data/'):
    if 'csv' in f:
        data[f[:-16]] = pd.DataFrame.from_csv('/home/katya/start/HIS3InterspeciesEpistasis/Data/' + f, sep = '\t')
        clean_data[f[:-16]] = data[f[:-16]][(data[f[:-16]].nonsense == 0) & (data[f[:-16]].middle == 1)]
        clean_data[f[:-16]] = clean_data[f[:-16]].sample(6000)


# In[202]:

color_dict={}
color_dict[11] = '#00AEE9'
color_dict[10] = '#34ADD3'
color_dict[9] = '#70ADC7'
color_dict[8] = '#70ADC7'
color_dict[7] = '#70ADC7'
color_dict[6] = '#9AABB4'
color_dict[5] = '#9AABB4'
color_dict[4] = '#A7A9AC'
color_dict[3] = '#A7A9AC'
color_dict[2] = '#A7A9AC'
color_dict[1] = '#A7A9AC'


# In[203]:

def write_json(chunk):
    print chunk
    
    fitThres = [0.45*(10-i)*0.1 for i in range(11)]
    fitThres.sort()
    
    sqs = list(clean_data[chunk].index)
    
    data = {}
    data['nodes'] = []
    counter = 0
    for i in range(len(sqs)):
        counter+=1
        data['nodes'].append({'name':'%s, distance: %d, fitness: %.2f' % (clean_data[chunk].mut_list[i], clean_data[chunk].dist_Scer[i], clean_data[chunk].s[i]),
                              'x':(clean_data[chunk].dist_Scer[i] + np.random.normal(0, .05)), 
                              'y':-(clean_data[chunk].s[i]*10)})
        
    data['connections'] = []
    
    done = []
    for i1 in range(len(sqs)):
        if i1%1000 == 0:
            print i1
        done.append(sqs[i1])
        for i2 in range(len(sqs)):
            if sqs[i2] not in done and clean_data[chunk].dist_Scer[i1]!=clean_data[chunk].dist_Scer[i2]             and distance.hamming(sqs[i1],sqs[i2]) == 1:
                
                minimum = min(clean_data[chunk].s[i1], clean_data[chunk].s[i2])
                ind = np.searchsorted(fitThres, minimum, side='right')
                data['connections'].append({'source':data['nodes'][i1]['name'], 'target':data['nodes'][i2]['name'], 'color':color_dict[ind], 'min':minimum})
                   
    with open('/home/katya/start/HIS3InterspeciesEpistasis/Analysis/Katya/Fig2/Build/data/' + chunk + '.json', 'w+') as outfile:  
        json.dump(data, outfile)


# In[204]:

pool = Pool()
pool.map(write_json, chunks)


# ***

# In[ ]:



