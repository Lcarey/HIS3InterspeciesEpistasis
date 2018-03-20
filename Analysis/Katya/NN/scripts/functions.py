from __future__ import print_function
from __future__ import division
from sklearn import cross_validation
from sklearn.preprocessing import MinMaxScaler
import tensorflow as tf
import numpy as np
import itertools
from scipy.stats import gaussian_kde
from sklearn.model_selection import train_test_split
from keras.models import Sequential
from keras.layers import Dense, Activation
from keras.layers.core import Dense, Dropout, Activation
from keras import optimizers,initializers
import pandas as pd
from keras.optimizers import RMSprop
from keras.utils import np_utils
from scipy import stats


def containsMutations(genotype):
    return genotype.split(':')

def makeBinary(unique_mutations, genotype):
    genotypeList = containsMutations(genotype)
    indexList = []
    
    for i in range(len(genotypeList)):
        indexList.append(unique_mutations.index(genotypeList[i]))
    
    line = np.zeros((1,len(unique_mutations)))
    line[:,indexList] = 1.
    
    return line

#############READING DATA#############
######EXTRACTING DATA AND LABELS######

def read_data(chunk):
    input_file = '~/Jupyter/HIS3InterspeciesEpistasis/Analysis/Katya/NN/data/' + chunk + '.txt'
    data = pd.read_table(input_file)
    data.columns = ['mutList', 'fitness', 'aa_seq']
    data.mutList = data.mutList.fillna('')
    unique_mutations = set(':'.join(data.mutList).split(':'))
    unique_mutations = sorted(list(unique_mutations))
    if '' in unique_mutations:
        unique_mutations.remove('')

#     data = data.reindex(np.random.permutation(data.index))
    nn_genotypes_values = np.zeros((len(data), len(unique_mutations)))
    nn_brightness_values = data.fitness.values
    aa_seq = data.aa_seq
    
    for i in range(len(data.mutList)):
        if data.mutList[i] != '':
            nn_genotypes_values[i] = makeBinary(unique_mutations, data.mutList[i])[0]
    
    return nn_genotypes_values, nn_brightness_values, unique_mutations, aa_seq


def makeBinaryDoubles(combinations,genotype):
    mutList = containsMutations(genotype)
    indexList = []
    
    for i in range(len(combinations)):
        if combinations[i][0] in mutList and combinations[i][1] in mutList:
            indexList.append(i)
            
    line = np.zeros((1,len(combinations)))
    line[:,indexList] = 1.
    
    return line   


#############READING DATA#############
######EXTRACTING DATA AND LABELS######
####features are double mutations####

def read_data_doubles(chunk):
    input_file = '~/Jupyter/HIS3InterspeciesEpistasis/Analysis/Katya/NN/data/' + chunk + '.txt'
    data = pd.read_table(input_file)
    data.columns = ['mutList', 'fitness', 'aa_seq']
    data.mutList = data.mutList.fillna('')
    unique_mutations = set(':'.join(data.mutList).split(':'))
    unique_mutations = sorted(list(unique_mutations))
    if '' in unique_mutations:
        unique_mutations.remove('')

    combinations = []

    for sub in itertools.combinations(unique_mutations, 2):
        combinations.append(sub)        
        
    data = data.reindex(np.random.permutation(data.index))
    nn_genotypes_values = np.zeros((len(data), len(combinations)))
    nn_brightness_values = data.fitness.values
        
    for i in range(len(data.mutList)):
        if data.mutList[i] != '':
            nn_genotypes_values[i] = makeBinaryDoubles(combinations, data.mutList[i])[0]

    nn_brightness_values = (nn_brightness_values - min(nn_brightness_values)) / max(
        nn_brightness_values - min(nn_brightness_values))

#     nn_genotypes_values = nn_genotypes_values[:, nn_genotypes_values.sum(axis=0) != 0]
    
    return nn_genotypes_values, nn_brightness_values, combinations


def makeBinaryTriples(combinations, genotype):
    mutList = containsMutations(genotype)
    indexList = []
    
    for i in range(len(combinations)):
        if combinations[i][0] in mutList and combinations[i][1] in mutList and combinations[i][2] in mutList:
            indexList.append(i)
            
    line = np.zeros((1,len(combinations)))
    line[:,indexList] = 1.
    
    return line   


#############READING DATA#############
######EXTRACTING DATA AND LABELS######
####features are triple mutations####

def read_data_triples(chunk):
    input_file = '~/Jupyter/HIS3InterspeciesEpistasis/Analysis/Katya/NN/data/' + chunk + '.txt'
    data = pd.read_table(input_file)
    data.columns = ['mutList', 'fitness', 'aa_seq']
    data.mutList = data.mutList.fillna('')
    unique_mutations = set(':'.join(data.mutList).split(':'))
    unique_mutations = sorted(list(unique_mutations))
    if '' in unique_mutations:
        unique_mutations.remove('')
    
    combinations = []

    for sub in itertools.combinations(unique_mutations, 3):
        combinations.append(sub)
    
    data = data.reindex(np.random.permutation(data.index))
    nn_genotypes_values = np.zeros((len(data), len(combinations)))
    nn_brightness_values = data.fitness.values
        
    for i in range(len(data.mutList)):
        if data.mutList[i] != '':
            nn_genotypes_values[i] = makeBinaryTriples(combinations, data.mutList[i])[0]

    nn_brightness_values = (nn_brightness_values - min(nn_brightness_values)) / max(
        nn_brightness_values - min(nn_brightness_values))
    
    nn_genotypes_values = nn_genotypes_values[:, nn_genotypes_values.sum(axis=0) != 0]
    
    return nn_genotypes_values, nn_brightness_values
    
    
#############PLOTTING FITNESS POTENTIAL VS FITNESS#############

def fitness_potential_plotting(weights_set, num=1):
    plt.figure(figsize=(4*5+3,3*5))
    plt.suptitle('Fitness potential vs Fitness', size=30, x=0.8, y=1.6)
    count=1
    for chunk in chunks:
        plt.subplot(3,4,count)
        plt.subplots_adjust(top = 1.5,right=1.5)
        plt.title(chunk, fontsize=20)

        plt.plot(fitness_potential[chunk][:,weights_set][:10000], true[chunk][:10000], 'ok', alpha = 0.01)
        
        plt.plot(fitness_potential[chunk][:,weights_set][:10000], predicted[chunk][:10000], '.', c='#36D1C4', alpha = 0.03)
        
        if num == 2:
            plt.plot(fitness_potential[chunk][:,weights_set+1][:10000], predicted[chunk][:10000], '.', c='#FF006C', alpha = 0.01)

        plt.grid(True, ls='--', lw=0.5, dash_capstyle = 'round', c='gray')
        plt.xlabel('Fitness potential', fontsize=15)
        plt.ylabel('Observed values', fontsize=15)
        count+=1
        
        
chunks = [('S'+str(x)) for x in range(1,13)]

#############EXTRACTING UNIQUE MUTATIONS FOR EACH OF THE SEGMENTS#############

unique_mutations = {}
for chunk in chunks:
    input_file = '~/Jupyter/HIS3InterspeciesEpistasis/Analysis/Katya/NN/data/' + chunk + '.txt'
    df = pd.read_table(input_file)
    df.columns = ['mutList', 'fitness', 'aa_seq']
    df.mutList = df.mutList.fillna('')
    unique_mutations[chunk] = set(':'.join(df.mutList).split(':'))
    unique_mutations[chunk] = sorted(list(unique_mutations[chunk] ))
    if '' in unique_mutations[chunk]:
        unique_mutations[chunk] .remove('')
        
        
def disectSeq (sq):
    aaList = []
    counter = 0
    for a in sq:
        aaList.append(str(counter)+a)
        counter+=1
    return ':'.join(aaList)


def read_data_all_positions(chunk):
    
    input_file = '~/Jupyter/HIS3InterspeciesEpistasis/Analysis/Katya/NN/data/' + chunk + '.txt'
    data = pd.read_table(input_file)
    data.columns = ['mutList', 'fitness', 'aa_seq']
    
    mut_list = list(data.mutList)
    genotypeList = [disectSeq(sq) for sq in data.aa_seq]
    data.mutList += ':'
    data.mutList += pd.Series(genotypeList)
    
    data.mutList = data.mutList.fillna('')
    
    unique_mutations = set(':'.join(data.mutList).split(':'))
    unique_mutations = sorted(list(unique_mutations))
    if '' in unique_mutations:
        unique_mutations.remove('')

#     data = data.reindex(np.random.permutation(data.index))
    nn_genotypes_values = np.zeros((len(data), len(unique_mutations)))
    nn_brightness_values = data.fitness.values
    aa_seq = data.aa_seq
    
    for i in range(len(data.mutList)):
        if data.mutList[i] != '':
            nn_genotypes_values[i] = makeBinary(unique_mutations, data.mutList[i])[0]
    
    return nn_genotypes_values, nn_brightness_values, unique_mutations, aa_seq, data.mutList