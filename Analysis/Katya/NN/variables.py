import matplotlib
matplotlib.use('Agg')
import argparse
# import boto3

# s3 = boto3.resource('s3')

parser = argparse.ArgumentParser(description='''

*** Only 3-layered nets are currently supported ***

Parameters file has to contain 8 parameters -- each on a new line.
Each line has to start with a parameter name, followed by tab and the parameter meaning.


The parameters are:

1) input_file (Path to the input file)

2) output_folder (Path to the output folder)

3) learning_rate (Learning rate)

4) batch_size (Batch size)

5) number_of_iterations (Number of iterations)

6) reshuffling_frequency (Once per how many iterations would you like to reshuffle your data?)

7) optimizer (Optimizer method,
please choose one from https://www.tensorflow.org/versions/r0.9/api_docs/python/train.html#optimizers)

8) net_structure (Net structure)


Net structure has to be formatted the following way:

net_structure	1,tf.tanh	3,tf.tanh	1,tf.tanh

Each layer of the network has to be separated by a tab.
The layer description consists of 2 parameters, separated by comma:

-- the number of neurons in the layer
-- the neuron output function

It is possible to choose from a variety of output functions, listed here:

https://www.tensorflow.org/versions/r0.9/api_docs/python/nn.html#activation-functions

!! It is very important to precisely follow the punctuation !!
Otherwise, the code will give an error.

''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('parameters_file')

arguments = parser.parse_args()
net_structure = {}

# Parsing the parameters document (details are in help above)
for line in open(arguments.parameters_file).readlines():
    if 'INPUT' in line.upper():
        input_file = str(line.rstrip('\n').split('\t')[1])
    if 'OUTPUT' in line.upper():
        output_folder = str(line.rstrip('\n').split('\t')[1])
    if 'LEARNING RATE' in line.upper() or 'LEARNING_RATE' in line.upper():
        learning_rate = float(line.rstrip('\n').split('\t')[1])
    if 'BATCH SIZE' in line.upper() or 'BATCH_SIZE' in line.upper():
        batch_size = int(line.rstrip('\n').split('\t')[1])
    if 'ITERATIONS' in line.upper():
        max_epochs = int(line.rstrip('\n').split('\t')[1])
    if 'RESHUFFLING' in line.upper():
        reshuffling_frequency = int(line.rstrip('\n').split('\t')[1])
    if 'OPTIMIZE' in line.upper():
        optimizer_method = str(line.rstrip('\n').split('\t')[1])
    if 'STRUCTURE' in line.upper():
        counter = 1
        for i in line.split('\t')[1:]:
            net_structure['layer' + str(counter)] = i.split(',')
            counter += 1

if output_folder[-1] != '/':
    output_folder += '/'

# Help dummy files that will be generated locally, before sending to s3. They'll be overwritten each iteration.
cost_stats_file = output_folder + 'cost_stats.txt'
temp_fig_file = output_folder + 'temp_fig.png'
temp_test_fig_file = output_folder + 'temp_test_fig.png'
temp_model_ckpt_file = output_folder + 'model.ckpt'