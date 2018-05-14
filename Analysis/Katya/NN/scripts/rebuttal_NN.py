from functions import *
from sklearn.metrics import r2_score, mean_squared_error
from scipy.stats import spearmanr, pearsonr
from optparse import OptionParser
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns
from keras.callbacks import EarlyStopping
from collections import OrderedDict


######FUNCTIONS######


def disectSeq(sq):
    aaList = []
    counter = 0
    for a in sq:
        aaList.append(str(counter)+a)
        counter+=1
    return ':'.join(aaList)

def make_feature_matrix(genotype_list,unique_mutations):
    feature_matrix = np.zeros((len(genotype_list), len(unique_mutations)))
    
    for i in range(len(genotype_list)):
        if genotype_list[i] != '':
            feature_matrix[i] = makeBinary(unique_mutations, genotype_list[i])[0]
            
    return feature_matrix

def train_NN(n_iter, labels, data, f, patience):
    
    min_max_scaler = MinMaxScaler()
    scaledLabels = min_max_scaler.fit_transform(labels.reshape(-1,1))

    x_train, x_valid, y_train, y_valid = train_test_split(data, scaledLabels, test_size = 0.01)

    model = Sequential()

    model.add(Dense(1,input_dim=data.shape[1],kernel_initializer='glorot_normal',activation=f))
    model.add(Dense(50,activation=f))
    model.add(Dense(100,activation=f))
    model.add(Dense(1,activation=f))

    opt = optimizers.RMSprop(lr=0.001)

    early_stopping_monitor=EarlyStopping(patience=patience)

    model.compile(optimizer=opt,
                  loss='mean_squared_error')

    hist = model.fit(x_train, y_train, validation_data=[x_valid, y_valid], 
                            epochs=n_iter, batch_size=500, shuffle=True, callbacks=[early_stopping_monitor],verbose=1)

    proba = model.predict_proba(data, batch_size=500)

    weights = model.layers[0].get_weights()[0]
    biases = model.layers[0].get_weights()[1]

    labels
    predicted = min_max_scaler.inverse_transform(proba)
    predicted = predicted.flatten()
    fitness_potential = data.dot(weights) + biases

    return weights, predicted, fitness_potential

def separate_predicted_weights_by_aa(ref_weights,unique_mutations,weights):
    predicted_weights_per_aa = OrderedDict()
    for a in ref_weights.index:
        predicted_weights_per_aa[a] = []
        for i,aa in enumerate(unique_mutations):
            if a in aa:
                predicted_weights_per_aa[a].extend(weights[i])
                
    return predicted_weights_per_aa

def plot_model(labels, predicted, n, function):
    inds = np.random.choice([x for x in range(len(labels))],n)

    plt.figure(figsize=[10,12])
    r=r2_score(labels,predicted)
    X = fitness_potential[inds]
    Y = labels[inds]
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.plot(X, Y, 'o', alpha = 0.1, color='#F9A828')
    plt.plot(X, predicted[inds], 'o', c='k', alpha = 0.1)
    plt.grid(True, ls='--', lw=0.5, dash_capstyle = 'round', c='gray')
    plt.title('%s, r2=%.2f' % (function, r),fontsize=20)
    plt.xlabel('FP',fontsize=15)
    plt.ylabel('Fitness',fontsize=15)
    plt.savefig('../rebuttal/'+function+'_predicted_model.pdf')
    
def plot_weights_correlations(ref_weights,predicted_weights_per_aa,function):
    plt.figure(figsize=[10,8])
    plt.plot(ref_weights, [np.median(predicted_weights_per_aa[x]) for x in predicted_weights_per_aa], 'ok')
    plt.grid(lw=0.5)
    plt.xlabel('True weights',fontsize=15)
    plt.ylabel('Predicted weights',fontsize=15)
    r2=pearsonr(ref_weights.FP, [np.median(predicted_weights_per_aa[x]) for x in predicted_weights_per_aa])[0]**2
    plt.title('%s, r2=%.2f' % (function, r2),fontsize=20)
    plt.savefig('../rebuttal/'+function+'_weight_corr.pdf')



######ACTIONS######

parser = OptionParser()
parser.add_option("-f", "--function_n", type="int",
                  help="Name the number of function you want to use",
                  dest="function_n")

(options, args) = parser.parse_args()
function_list=['fitness_linear', 'fitness_logsig','fitness_sin','fitness_sin1','fitness_sin2']
function = function_list[options.function_n]
print function

ref_weights=pd.DataFrame.from_csv('../../../Lucas/Revisions/Simulations_NN_Detects_Linear_Nonmonotonic_fitness_functions__AAstates.tab', sep='\t')
lucas=pd.DataFrame.from_csv('../../../Lucas/Revisions/Simulations_NN_Detects_Linear_Nonmonotonic_fitness_functions.tab', sep='\t')

genotype_list = [disectSeq(sq) for sq in lucas.index]
unique_mutations = set(':'.join(genotype_list).split(':'))
unique_mutations = sorted(list(unique_mutations))

data = make_feature_matrix(genotype_list,unique_mutations)

labels = np.array(lucas[function])
f='sigmoid'
patience=100
n_iter=1000
weights, predicted, fitness_potential = train_NN(n_iter, labels, data, f, patience)
predicted_weights_per_aa = separate_predicted_weights_by_aa(ref_weights,unique_mutations,weights)

plot_model(labels, predicted, 5000, function)
plot_weights_correlations(ref_weights,predicted_weights_per_aa,function)
    