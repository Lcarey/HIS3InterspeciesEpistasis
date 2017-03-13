import matplotlib

matplotlib.use('Agg')
import tensorflow as tf
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import gaussian_kde


# Format data, split into feature matrix and target value vector
def format_data(data):
    # shuffling rows in the data df
    unique_mutations = data.unique_mutations
    data = data.data.reindex(np.random.permutation(data.data.index))

    # formatting data for the nn input
    print('Normalizing data...')
    nn_genotypes_values = np.zeros((len(data), len(unique_mutations)))
    nn_brightness_values = data.medianBrightness.values
    for i in range(len(unique_mutations)):
        nn_genotypes_values[:, i] = data.aaMutations.str.contains(unique_mutations[i]).astype(np.float32)

    nn_brightness_values = (nn_brightness_values - min(nn_brightness_values)) / max(
        nn_brightness_values - min(nn_brightness_values)) * 2 - 1

    return nn_genotypes_values, nn_brightness_values


# Split data into test and train sets randomly
def split(data):
    length = len(data.data)
    test_id = np.random.randint(0, length, data.batch_size)
    train_id = np.array([x for x in np.arange(length) if x not in test_id])

    genotypes_test, genotypes_train = data.nn_genotypes_values[test_id], data.nn_genotypes_values[train_id]
    brightness_test, brightness_train = data.nn_brightness_values[test_id], data.nn_brightness_values[train_id]

    return genotypes_test, brightness_test, genotypes_train, brightness_train


# Shuffle train set, split train set into batches. Get train batches and test batches
def get_batches(data):
    length = len(data.nn_brightness_train)
    batches = []
    test_batches = []
    ids = np.array([x for x in np.arange(length)])
    to_plot_observed = []

    for i in range(data.batch_number):
        current_ids = np.random.choice(ids, data.batch_size, replace=False)
        ids = np.array([x for x in ids if x not in current_ids])

        current_batch = data.nn_genotypes_train[current_ids].reshape(data.batch_size, 1, len(data.unique_mutations))
        current_batch_brightness = data.nn_brightness_train[current_ids].reshape(data.batch_size, 1, 1)
        to_plot_observed.append(data.nn_brightness_train[current_ids])
        batches.append((current_batch, current_batch_brightness))

    test_batch = data.nn_genotypes_test.reshape(data.batch_size, 1, len(data.unique_mutations))
    test_batch_brightness = data.nn_brightness_test.reshape(data.batch_size, 1, 1)
    test_batches.append((test_batch, test_batch_brightness))
    to_plot_observed_test = data.nn_brightness_test.reshape(data.batch_size)
    to_plot_observed = np.array(to_plot_observed).reshape(data.batch_number * data.batch_size)

    return batches, test_batches, to_plot_observed, to_plot_observed_test


# Broadcast a tensor (used before multiplication with weights)
def broadcast(tensor, batch_size):
    return tf.tile(tensor, (batch_size, 1, 1))


# Plot the observed values versus predicted using density plot
def density_plot(x, y, iteration_number, costs, test_score):
    ''' x = observed, y = predicted '''
    x = x[(~np.isnan(x)) & (~np.isnan(y))]
    y = y[(~np.isnan(x)) & (~np.isnan(y))]

    # Calculate the point density
    xy = np.vstack([x, y])
    z = gaussian_kde(xy)(xy)

    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]

    # formatting
    plt.figure(figsize=[6, 4])
    plt.scatter(x, y, c=z, s=3, edgecolor='', cmap='viridis_r')
    plt.xlim(-1, 1)
    plt.ylim(-1, 1)
    plt.xlabel('Observed brightness')
    plt.ylabel('Predicted brightness')
    plt.title('Iteration %s: cost=%.7f, EVS=%.2f' % (iteration_number, costs, test_score))
    plt.tick_params(axis="both", which="both", bottom="off", top="off",
                    labelbottom="on", left="off", right="off", labelleft="on")
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)
    plt.gca().spines["bottom"].set_color('gray')
    plt.gca().spines["left"].set_color('gray')
    plt.gca().xaxis.grid(True)
    plt.gca().yaxis.grid(True)


# Negative control -- getting random batches, brightness not corresponding to the features
def get_random_batches(data):
    length = len(data.nn_brightness_train)
    batches = []
    test_batches = []
    ids = np.array([x for x in np.arange(length)])
    to_plot_observed = []

    for i in range(data.batch_number):
        current_batch = data.nn_genotypes_train[np.random.choice(ids, data.batch_size, replace=False)].reshape(
            data.batch_size, 1, len(data.unique_mutations))
        current_batch_brightness = data.nn_brightness_train[
            np.random.choice(ids, data.batch_size, replace=False)].reshape(data.batch_size, 1, 1)
        to_plot_observed.append(data.nn_brightness_train[np.random.choice(ids, data.batch_size, replace=False)])
        batches.append((current_batch, current_batch_brightness))

    test_batch = data.nn_genotypes_test.reshape(data.batch_size, 1, len(data.unique_mutations))
    test_batch_brightness = data.nn_brightness_test.reshape(data.batch_size, 1, 1)
    test_batches.append((test_batch, test_batch_brightness))
    to_plot_observed_test = data.nn_brightness_test.reshape(data.batch_size)
    to_plot_observed = np.array(to_plot_observed).reshape(data.batch_number * data.batch_size)

    return batches, test_batches, to_plot_observed, to_plot_observed_test
