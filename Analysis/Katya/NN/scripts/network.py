from functions import *
from sklearn.metrics import r2_score, mean_squared_error
from scipy.stats import spearmanr, pearsonr
from optparse import OptionParser
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from keras.callbacks import EarlyStopping


parser = OptionParser()
parser.add_option("-c", "--chunk", type="string",
                  help="Select the chunk to process",
                  dest="chunk")
parser.add_option("-n", "--epochs", type="int",
                  help="Number of epochs to perform the iterations for",
                  dest="n_epochs")

(options, args) = parser.parse_args()

#min_max_scaler = MinMaxScaler()
n_iter = options.n_epochs

chunk = options.chunk
print (chunk)

data, labels, unique_mutations[chunk], aa_seq, mut_list = read_data_all_positions(chunk)

print 'Splitting the data'
x_train, x_valid, y_train, y_valid = train_test_split(data, labels, test_size = 0.01)

n_neurons = []

mse_val = []
mse_train = []
r_val=[]
r_train=[]

r2_weights = []

for i in range(1,22,4):
    
    print '\nNumber of weights combinations = ', i
    temp_mse_train_list=[]
    temp_mse_val_list=[]
    temp_r_train_list=[]
    temp_r_val_list=[]
    temp_weights_r2={}
    it=0
    loop_count=10
    
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
        temp_r_val = r2_score(y_valid,predicted_val)
        temp_r_train = r2_score(y_train,predicted_train)
        
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

plt.figure(figsize=[6,7])
for i,r_list in enumerate(r_val):
    plt.plot([i*4+1]*len(r_list),r_list,'o',alpha=0.5,color='#283149')
    plt.plot([i*4+0.5,i*4+1.5],[np.median(r_list)]*2,'-',lw=0.9,color='#DA0463')
plt.grid('--k',lw=0.5)
plt.title(chunk,fontsize=15)
plt.ylabel('R2 between predicted and true values',fontsize=13)
plt.xlabel('Number of linear combinations of substitutions',fontsize=13)
plt.savefig('/nfs/scistore08/kondrgrp/eputints/Jupyter/HIS3InterspeciesEpistasis/Analysis/Katya/NN/complexity/20_iterations/r_'+chunk+'.pdf')

plt.figure(figsize=[6,7])
for i,mse_list in enumerate(mse_val):
    plt.plot([i*4+1]*len(mse_list),mse_list,'o',alpha=0.5,color='#283149')
    plt.plot([i*4+0.5,i*4+1.5],[np.median(mse_list)]*2,'-',lw=0.9,color='#DA0463')
plt.grid('--k',lw=0.5)
plt.title(chunk,fontsize=15)
plt.ylabel('MSE',fontsize=13)
plt.xlabel('Number of linear combinations of substitutions',fontsize=13)
plt.savefig('/nfs/scistore08/kondrgrp/eputints/Jupyter/HIS3InterspeciesEpistasis/Analysis/Katya/NN/complexity/20_iterations/mse_'+chunk+'.pdf')

plt.figure(figsize=[6,7])
for i,r2_list in enumerate(r2_weights):
    plt.plot([i*4+1]*len(r2_list),r2_list,'o',alpha=0.5,color='#283149')
    plt.plot([i*4+0.5,i*4+1.5],[np.median(r2_list)]*2,'-',lw=0.9,color='#DA0463')
plt.ylim(-1,1)
plt.title(chunk,fontsize=15)
plt.grid('--k',lw=0.5)
plt.ylabel('Median Spearman R between sets of weights',fontsize=13)
plt.xlabel('Number of linear combinations of substitutions',fontsize=13)
plt.savefig('/nfs/scistore08/kondrgrp/eputints/Jupyter/HIS3InterspeciesEpistasis/Analysis/Katya/NN/complexity/20_iterations/r2_weights_'+chunk+'.pdf')
