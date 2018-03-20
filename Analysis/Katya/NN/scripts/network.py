from functions import *
from sklearn.metrics import r2_score, mean_squared_error
from scipy.stats import spearmanr
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
#scaledLabels = min_max_scaler.fit_transform(labels.reshape(-1,1))

print 'Splitting the data'
x_train, x_valid, y_train, y_valid = train_test_split(data, labels, test_size = 0.01)

n_neurons = []

r2_val = []
mse_val = []

r2_train = []
mse_train = []

r2_weights = []

for i in range(1,21):
    
    print '\nNumber of weights combinations = ', i
    temp_r2_train_list=[]
    temp_r2_val_list=[]
    temp_mse_train_list=[]
    temp_mse_val_list=[]
    temp_weights_r2={}
    it=0
    loop_count=5
    
    while it<loop_count and loop_count<100:
        print it
        model = Sequential()

        model.add(Dense(i,input_dim=data.shape[1],activation='relu',kernel_initializer='glorot_normal'))
        model.add(Dense(40,activation='relu'))
        model.add(Dense(20,activation='relu'))
        model.add(Dense(1,activation='relu'))

        opt = optimizers.RMSprop(lr=0.01, rho=0.9, epsilon=1e-08, decay=0.0)
        
        early_stopping_monitor=EarlyStopping(patience=3)
        
        model.compile(optimizer=opt,
                      loss='mean_squared_error')

        hist = model.fit(x_train, y_train, validation_data=[x_valid, y_valid], 
                                epochs=n_iter, batch_size=500, shuffle=True, callbacks=[early_stopping_monitor],verbose=0)

        proba = model.predict_proba(x_valid, batch_size=500,verbose=0)
        #predicted_val = min_max_scaler.inverse_transform(proba)
        predicted_val = proba.flatten()

        proba = model.predict_proba(x_train, batch_size=500,verbose=0)
        #predicted_train = min_max_scaler.inverse_transform(proba)
        predicted_train = proba.flatten()
        
        weights = model.layers[0].get_weights()[0]
        
        temp_r2_val = r2_score(y_valid,predicted_val)
        temp_mse_val = mean_squared_error(y_valid,predicted_val)
        temp_r2_train = r2_score(y_train,predicted_train)
        temp_mse_train = mean_squared_error(y_train,predicted_train)
        
        it+=1
        
        #Sanity checks
        if temp_r2_val>-1:
            temp_r2_val_list.append(temp_r2_val)
        else:
            loop_count+=1
            
        if temp_r2_train>-1:
            temp_r2_train_list.append(temp_r2_train)
            
        if temp_mse_val<0.1:
            temp_mse_val_list.append(temp_mse_val)
        
        if temp_mse_train<0.1:
            temp_mse_train_list.append(temp_mse_train)    
            
            for combination in list(itertools.combinations([x for x in range(i)], 2)):
                if combination in temp_weights_r2:
                    temp_weights_r2[combination].extend([spearmanr(weights[:,combination[0]],weights[:,combination[1]])])
                else:
                    temp_weights_r2[combination] = [spearmanr(weights[:,combination[0]],weights[:,combination[1]])]
                    
                    
    n_neurons.append(i)
    r2_val.append(np.median(temp_r2_val_list))
    mse_val.append(np.median(temp_mse_val_list))
    r2_train.append(np.median(temp_r2_train_list))
    mse_train.append(np.median(temp_mse_train_list))
    print temp_weights_r2
    if i>1:
        r2_weights.append([np.median(temp_weights_r2[x]) for x in temp_weights_r2])

plt.figure(figsize=[10,8])
plt.plot(n_neurons,r2_val,'-ok',label='Validation',lw=0.5)
plt.plot(n_neurons,r2_train,'--or',label='Train',lw=0.5)
plt.ylabel('Median R2')
plt.xlabel('# neurons first layer')
plt.legend()
plt.savefig('/nfs/scistore08/kondrgrp/eputints/Jupyter/HIS3InterspeciesEpistasis/Analysis/Katya/NN/complexity/20_iterations/r2_'+chunk+'.pdf')

plt.figure(figsize=[10,8])
plt.plot(n_neurons,mse_val,'-ok',label='Validation',lw=0.5)
plt.plot(n_neurons,mse_train,'--or',label='Train',lw=0.5)
plt.ylabel('Median MSE')
plt.xlabel('# neurons first layer')
plt.legend()
plt.savefig('/nfs/scistore08/kondrgrp/eputints/Jupyter/HIS3InterspeciesEpistasis/Analysis/Katya/NN/complexity/20_iterations/mse_'+chunk+'.pdf')

plt.figure(figsize=[10,8])
for i,r2_list in enumerate(r2_weights):
    plt.plot([i+2]*len(r2_list),r2_list,'o',alpha=0.7)
plt.ylabel('Median R2 between sets of weights')
plt.xlabel('# neurons first layer (= number of sets of weights)')
plt.savefig('/nfs/scistore08/kondrgrp/eputints/Jupyter/HIS3InterspeciesEpistasis/Analysis/Katya/NN/complexity/20_iterations/r2_weights_'+chunk+'.pdf')
