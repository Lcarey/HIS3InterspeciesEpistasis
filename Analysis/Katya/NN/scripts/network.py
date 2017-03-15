import matplotlib
matplotlib.use('Agg')
from classes import *
from variables import *
import pandas as pd
import tensorflow as tf
import numpy as np
from matplotlib import pyplot as plt
from sklearn.metrics import explained_variance_score

print('Folder name for this run is %s' % output_folder)
print('# learning rate = %s' % learning_rate)
print('# number of iterations = %s' % max_epochs)
print('# batch size = %s' % batch_size)
print('# reshuffling frequency = %s' % reshuffling_frequency)

# SESSION

# 1) Setting the data parameters. Building up the neural network.
data = Data(input_file, batch_size)
net = TFNet(net_structure, data, optimizer_method, learning_rate, batch_size, cost_stats_file)

# 2) Running the session.
with tf.Session(config=tf.ConfigProto(log_device_placement=True)) as sess:
    # Creating the file with cost statistics, writing initial parameters.
    cost_stats = open(net.cost_stats_file, 'w+')
    cost_stats.write('# net')
    for i in net_structure:
        cost_stats.write('_%s' % (net_structure[i][0]))
    cost_stats.write('\n# learning rate = %s\n' % learning_rate + \
                     '# max iteration limit = %s\n' % max_epochs + \
                     '# batch size = %s\n' % batch_size)
    cost_stats.write('iteration,cost\n')

    # Initializing variables
    sess.run(net.init)

    # Initiating the session run for a specified number of iterations.
    for e in range(max_epochs):
        for batch, batch_brightness in data.batches:
            sess.run(net.optimizer, feed_dict={data.nn_genotypes: batch, data.nn_brightness: batch_brightness})

        # Write down the outputs every 10th iteration.
        if e % 10 == 0:

            # Extracting net cost function output.
            to_plot_predicted = np.zeros(data.batch_number * batch_size)
            test_figure_name = output_folder + 'figures/test_iteration_%05d' % e
            figure_name = output_folder + 'figures/net'
            for i in net_structure:
                figure_name += '_%s' % (net_structure[i][0])
            figure_name += '_iteration_%05d' % e

            costs = 0

            # train data prediction
            train_score = []

            for index, (batch, batch_brightness) in enumerate(data.batches):
                cost_value, l3_value = sess.run([net.cost, net.output['layer3']],
                                                feed_dict={data.nn_genotypes: batch,
                                                           data.nn_brightness: batch_brightness})
                costs += cost_value

                to_plot_predicted[(index * batch_size):((index + 1) * batch_size)] = l3_value.reshape(batch_size)
                train_score.append(
                    explained_variance_score(l3_value.reshape(batch_size), batch_brightness.reshape(batch_size)))

            costs /= data.batch_number

            train_score = np.median(train_score)
            print('Train explained variance score is %.2f' % train_score)

            # Plotting observed versus predicted brightness. Saving the plot locally to a temp_fig_file.
            print('Iteration %s: cost=%.7f' % (e, costs))

            cost_stats.write('%s,%s\n' % (e, costs))
            density_plot(data.to_plot_observed, to_plot_predicted, e, costs, train_score)
            plt.savefig(temp_fig_file)
            plt.close('all')

            # Saving ckpt file and sending it and figure file to s3://landscapes-tensorflow.
            net.saver.save(sess, temp_model_ckpt_file)

            with open(temp_model_ckpt_file + '.meta', 'rb') as ckpt:
                s3.Bucket('landscapes-tensorflow').put_object(Key=figure_name + '.ckpt', Body=ckpt)
            with open(temp_fig_file, 'rb') as figure_file:
                s3.Bucket('landscapes-tensorflow').put_object(Key=figure_name + '.png', Body=figure_file)

            # Saving layer1 inputs and layer3 outputs and sending those to s3://landscapes-tensorflow.
            layer1_inputs = np.zeros(data.batch_number * batch_size)
            layer3_outputs = np.zeros(data.batch_number * batch_size)

            for index, (batch, batch_brightness) in enumerate(data.batches):
                l1_values, l3_values = sess.run([net.input['layer1'], net.input['layer3']],
                                                feed_dict={data.nn_genotypes: batch,
                                                           data.nn_brightness: batch_brightness})
                layer1_inputs[(index * batch_size):((index + 1) * batch_size)] = l1_values.reshape(batch_size)
                layer3_outputs[(index * batch_size):((index + 1) * batch_size)] = l3_values.reshape(batch_size)

            neuronal_values = pd.DataFrame()
            neuronal_values['layer1_inputs'] = layer1_inputs
            neuronal_values['layer3_outputs'] = layer3_outputs
            #
            neuronal_values_filename = output_folder + 'neuronal_values_iteration_%s.csv' % e
            neuronal_values.to_csv(neuronal_values_filename, index=False)
            with open(neuronal_values_filename, 'rb') as f:
                s3.Bucket('landscapes-tensorflow').put_object(Key=neuronal_values_filename, Body=f)

        # Reshuffling the data with the specified reshuffling frequency.
        if e % reshuffling_frequency == 0:

            # Saving parameters before reshuffling
            cost_stats.close()
            with open(cost_stats_file, 'rb') as cost_stats:
                s3.Bucket('landscapes-tensorflow').put_object(Key=cost_stats_file, Body=cost_stats)
            cost_stats = open(cost_stats_file, 'a')

            weights = sess.run(
                [net.weights['layer1'], net.biases['layer1'], net.weights['layer2'], net.biases['layer2'],
                 net.weights['layer3'], net.biases['layer3']],
                feed_dict={data.nn_genotypes: batch,
                           data.nn_brightness: batch_brightness})

            mutations_weights = pd.DataFrame()
            mutations_weights["mutation"] = data.unique_mutations
            mutations_weights["weight"] = weights[0].reshape(len(data.unique_mutations))
            mutations_weights_filename = output_folder + 'unique_mutations_scores_iteration_%s.csv' % e
            mutations_weights.to_csv(mutations_weights_filename, index=False)
            with open(mutations_weights_filename, 'rb') as mutations_weights:
                s3.Bucket('landscapes-tensorflow').put_object(Key=mutations_weights_filename,
                                                              Body=mutations_weights)

            if e != 0:

                # reshuffling data
                data.reshuffle()

                # test data prediction
                for index, (batch, batch_brightness) in enumerate(data.test_batches):
                    test_l3_value = sess.run([net.output['layer3']],
                                             feed_dict={data.nn_genotypes: batch,
                                                        data.nn_brightness: batch_brightness})

                    test_values = test_l3_value

                    test_score = explained_variance_score(test_values[0].reshape(batch_size),
                                                          batch_brightness.reshape(batch_size))

                print('Test explained variance score is %.2f' % test_score)
                density_plot(data.to_plot_observed_test, np.array(test_l3_value).reshape(batch_size), e, costs,
                             test_score)
                plt.savefig(temp_test_fig_file)
                plt.close('all')
                with open(temp_test_fig_file, 'rb') as figure_file:
                    s3.Bucket('landscapes-tensorflow').put_object(Key=test_figure_name + '.png', Body=figure_file)

    cost_stats.close();

    # sending cost_stats file to s3://landscapes-tensorflow
    with open(cost_stats_file, 'rb') as cost_stats:
        s3.Bucket('landscapes-tensorflow').put_object(Key=cost_stats_file, Body=cost_stats)
