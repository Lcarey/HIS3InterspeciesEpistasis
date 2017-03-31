
# coding: utf-8

# In[1]:

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import os
import distance
from multiprocessing import Pool
import json


# In[68]:

chunks = ['S'+str(i) for i in range(1,13)]
data = {}
clean_data = {}

for f in os.listdir('/users/fk/eputintseva/projects/HIS/Data/'):
    if 'csv' in f:
        data[f[:-16]] = pd.DataFrame.from_csv('/users/fk/eputintseva/projects/HIS/Data/' + f, sep = '\t')
        clean_data[f[:-16]] = data[f[:-16]][(data[f[:-16]].nonsense == 0) & (data[f[:-16]].middle == 1)]        


# In[69]:

def write_json(chunk):
    print chunk
    
    fitness_threshold_1 = 0.8*clean_data[chunk].s.max()
    fitness_threshold_2 = 0.6*clean_data[chunk].s.max()
    fitness_threshold_3 = 0.4*clean_data[chunk].s.max()
    fitness_threshold_4 = 0.2*clean_data[chunk].s.max()
    
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
                
                if clean_data[chunk].s[i1]>=fitness_threshold_1 and clean_data[chunk].s[i2]>=fitness_threshold_1:
                    data['connections'].append({'source':data['nodes'][i1]['name'], 'target':data['nodes'][i2]['name'], 'color':'#00AEE9'})
                
                elif clean_data[chunk].s[i1]>=fitness_threshold_2 and clean_data[chunk].s[i2]>=fitness_threshold_2:
                    data['connections'].append({'source':data['nodes'][i1]['name'], 'target':data['nodes'][i2]['name'], 'color':'#34ADD3'})
                
                elif clean_data[chunk].s[i1]>=fitness_threshold_3 and clean_data[chunk].s[i2]>=fitness_threshold_3:
                    data['connections'].append({'source':data['nodes'][i1]['name'], 'target':data['nodes'][i2]['name'], 'color':'#70ADC7'})
                
                elif clean_data[chunk].s[i1]>=fitness_threshold_4 and clean_data[chunk].s[i2]>=fitness_threshold_4:
                    data['connections'].append({'source':data['nodes'][i1]['name'], 'target':data['nodes'][i2]['name'], 'color':'#9AABB4'})
                
                else:
                    data['connections'].append({'source':data['nodes'][i1]['name'], 'target':data['nodes'][i2]['name'], 'color':'#A7A9AC'})
                
    with open('/users/fk/eputintseva/projects/HIS/Fig2/Build/data/' + chunk + '.json', 'w+') as outfile:  
        json.dump(data, outfile)


# In[1]:

pool = Pool()
pool.map(write_json, chunks)


# In[ ]:



