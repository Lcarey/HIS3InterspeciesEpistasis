from functions import *
from sklearn.metrics import r2_score, mean_squared_error
from scipy.stats import spearmanr, pearsonr
from optparse import OptionParser
import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
import seaborn as sns
from keras.callbacks import EarlyStopping

#min_max_scaler = MinMaxScaler()
n_iter = 100
parser = OptionParser()
parser.add_option("-c", "--chunk", type="string",
                  help="Select the chunk to process",
                  dest="chunk")


(options, args) = parser.parse_args()

chunk = options.chunk
print (chunk)

data, labels, unique_mutations[chunk], aa_seq, mut_list = read_data_all_positions(chunk)

print 'Splitting the data'
x_train, x_valid, y_train, y_valid = train_test_split(data, labels, test_size = 0.05)

n_neurons = []

mse_val = []
mse_train = []
r_val=[]
r_train=[]

r2_weights = []

for i in range(1,11):
    
    print '\nNumber of weights combinations = ', i
    temp_mse_train_list=[]
    temp_mse_val_list=[]
    temp_r_train_list=[]
    temp_r_val_list=[]
    temp_weights_r2={}
    it=0
    loop_count=5
    
    while it<loop_count and loop_count<100:
        print it
        model = Sequential()

        model.add(Dense(i,input_dim=data.shape[1],activation='sigmoid',kernel_initializer='glorot_normal'))
        model.add(Dense(20,activation='sigmoid'))
        model.add(Dense(1,activation='relu'))

        opt = optimizers.RMSprop(lr=0.01, rho=0.9, epsilon=1e-08, decay=0.0)
        
        early_stopping_monitor=EarlyStopping(patience=10)
        
        model.compile(optimizer=opt,
                      loss='mean_squared_error')

        hist = model.fit(x_train, y_train, validation_data=[x_valid, y_valid], 
                                epochs=n_iter, batch_size=500, shuffle=True, callbacks=[early_stopping_monitor],verbose=1)

        proba = model.predict_proba(x_valid, batch_size=500,verbose=0)
        predicted_val = proba.flatten()

        proba = model.predict_proba(x_train, batch_size=500,verbose=0)
        predicted_train = proba.flatten()
        
        weights = model.layers[0].get_weights()[0]
        
        temp_mse_val = mean_squared_error(y_valid,predicted_val)
        temp_mse_train = mean_squared_error(y_train,predicted_train)
        temp_r_val = pearsonr(y_valid,predicted_val)[0]
        temp_r_train = pearsonr(y_train,predicted_train)[0]
        
        it+=1
        
        #Sanity checks
        if temp_mse_val<0.1:
            temp_mse_val_list.append(temp_mse_val)
            temp_mse_train_list.append(temp_mse_train)
            temp_r_val_list.append(temp_r_val)
            temp_r_train_list.append(temp_r_train)
            
            for combination in list(itertools.combinations([x for x in range(i)], 2)):
                if combination in temp_weights_r2:
                    temp_weights_r2[combination].extend([spearmanr(weights[:,combination[0]],weights[:,combination[1]])])
                else:
                    temp_weights_r2[combination] = [spearmanr(weights[:,combination[0]],weights[:,combination[1]])]
                    
        else:
            loop_count+=1
                    
                    
    n_neurons.append(i)
    mse_val.append(temp_mse_val_list)
    mse_train.append(temp_mse_train_list)
    r_val.append(temp_r_val_list)
    r_train.append(temp_r_train_list)
    
    if i>1:
        r2_weights.append([np.median(temp_weights_r2[x]) for x in temp_weights_r2])
        
plt.figure(figsize=[6,6])
for i,r_list in enumerate(r_val):
    plt.plot([i+1]*len(r_list),[float(x)**2 for x in r_list],'o',alpha=0.5,color='#283149',ms=5)
    plt.plot([i+0.5,i+1.5],[np.median([float(x)**2 for x in r_list])]*2,'-',color='#DA0463')
plt.grid('--k',lw=0.5)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.savefig('/nfs/scistore08/kondrgrp/eputints/Jupyter/HIS3InterspeciesEpistasis/Analysis/Katya/NN/tmp/r2_'+chunk+'.pdf')