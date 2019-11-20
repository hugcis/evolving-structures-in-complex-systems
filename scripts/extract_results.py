"""
Module for extracting the output data created by the automaton binary upon
running.
"""
import os
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import average_precision_score
from DataItem import DataItem

ST = 3
DATA_DIR = os.path.join(os.getcwd(), "data_2d_{}/".format(ST))
LABELS_DIR = os.path.join(os.getcwd(), "data/")

def split_data(data_item):
    """ Split a string supposed to represent an item and ignore if wrong
    format.
    """
    left_right = data_item.strip().split(':')
    if len(left_right) != 2:
        print("Wrong data format: {}".format(left_right))
    return left_right

def get_data_item_or_add(results_dic, name, n_hid, epochs, horizon, timesteps):
    """ Return or create a new DataItem in `results_dic` with the corresponding
    metadata.
    """
    if name not in results_dic:
        results_dic[name] = []

    found = False
    for item in results_dic[name]:
        if item.is_metadata(n_hid, epochs, horizon, timesteps):
            found = True
            return item

    if not found:
        results_dic[name].append(
            DataItem(n_hid, epochs, horizon, timesteps))

    return results_dic[name][-1]

def process_nn_file(i, results_dic):
    """ Process a file produced when computing the neural network score """
    with open(os.path.join(os.path.join(DATA_DIR, 'nn/'), i)) as nn_file:
        body = nn_file.read().strip()
        if not body:
            return

        name = i.split('.')[0].strip('nn')

        for data_item in body.split('\n'):
            wrong_format = False
            try:
                left_right = split_data(data_item)
            except ValueError:
                wrong_format = True

            if not wrong_format:
                metadata = [int(i) for i in left_right[0].split(' ')]
                data = [float(i) for i in left_right[1].split(' ')]

                data_item = get_data_item_or_add(results_dic,
                                                 name,
                                                 metadata[0],
                                                 metadata[1],
                                                 metadata[2],
                                                 metadata[3])

                data_item.train, data_item.test = data[0], data[1]


        if len(left_right) != 2:
            return


def extract_labels(results_dic):
    """ Extract IDs and labels from the labels file. """
    results_dic_dataset = {}
    labels = pd.read_csv(
        os.path.join(LABELS_DIR, "data_labels.csv"), header=None).values

    interest = dict(labels)

    for key in interest:
        results_dic_dataset[key] = results_dic[key]

    return results_dic_dataset, interest


def get_options(x):
    return set(map(str, (itertools.chain(*x.values()))))


def get_res(results_dic, horizon, ts, epochs, n_hid, function):
    available_opts = get_options(results_dic)
    if DataItem.make_repr(n_hid, epochs, horizon, ts) not in available_opts:
        raise ValueError(
            "No data available for this combination of values: {}".format(
                DataItem.make_repr(n_hid, epochs, horizon, ts)))
    res = ([], [])
    for rule in results_dic:
        for item in results_dic[rule]:
            # Take first that matches
            result = function(item)
            if (item.is_metadata(n_hid, epochs, horizon, ts) and
                    result is not None):
                res[0].append(rule)
                res[1].append(result)
                continue
    res = (np.array(res[0]), np.array(res[1]))
    return res


def get_results(results_dic, horizon, ts,
                epochs=30, n_hid=100, function=(lambda x, y: x/y)):
    return get_res(results_dic, horizon, ts, epochs, n_hid,
                   lambda item: function(item.train, item.test))

def update_lookup(ent_file, score_lookup, name):
    """ Parse the lookup table score files and add them. """

    splitted_line = ent_file.readline().strip().split('    ')
    content = list(map(float, splitted_line))
    score_lookup["score_2_300"][name] = content[0]/content[1]
    score_lookup["score_2_50"][name] = content[3]/content[4]
    score_lookup["score_2_5"][name] = content[6]/content[7]

    splitted_line = ent_file.readline().strip().split('    ')
    content = list(map(float, splitted_line))
    score_lookup["score_1_300"][name] = content[0]/content[1]
    score_lookup["score_1_50"][name] = content[3]/content[4]
    score_lookup["score_1_5"][name] = content[6]/content[7]


def get_interest(names, interest):
    return np.array([interest[n] for n in names if n in interest])


def get_precision(results_dic_dataset, interest, neigh, ts, n_hid,
                  fun=lambda x, y: y/x):
    names, score = get_results(results_dic_dataset, neigh, ts, n_hid=n_hid,
                               function=fun)
    inter = get_interest(names, interest)
    return average_precision_score(inter, score)


def make_plots(results_dic_dataset, interest, function, save=None):
    fig, axs = plt.subplots(1, 3, sharey=True, figsize=(10, 3))
    rnge = range(1, 6)

    for time, axe in zip([5, 50, 300], axs):
        res = [get_precision(results_dic_dataset,
                             interest,
                             i, time, 10, function) for i in rnge]

        axe.plot(rnge, res, label="10 hidden units")

        # Format for a Latex table
        print(" & ".join(map('{:.3f}'.format, res)))

        axe.plot(rnge, [
            get_precision(results_dic_dataset, interest, i, time, 20,
                          function) for i in rnge
        ], label="20 hidden units")

        axe.set_title("Average precision {} timesteps".format(time))

    axs[0].legend()
    axs[0].set_xlabel("Radius")
    axs[0].set_ylabel("Average precision")

    fig.tight_layout()

    if save is not None:
        fig.savefig(save)

def shuffle_names(results_dic_dataset):
    rand_idx = np.arange(len(
        get_results(results_dic_dataset, 1, 50,
                    n_hid=20, function=lambda x, y: x/y)[1]
    ))
    np.random.shuffle(rand_idx)
    names, _ = get_results(results_dic_dataset, 1, 50, n_hid=20,
                           function=lambda x, y: x/y)
    rand_names = names[rand_idx]

    return rand_names


def print_accuracies(results_dic_dataset, interest, rand_names):
    # Enumerate timesteps
    for timesteps in [5, 50, 300]:
        # Enumerate sizes of hidden layer
        for hid in [10, 20]:
            accuracies = []
            print(timesteps, hid)

            # Enumerate radii
            for radius in range(1, 6):
                comp = get_results(results_dic_dataset, radius, timesteps,
                                   n_hid=hid, function=lambda x, y: x)[1]
                names, score = get_results(results_dic_dataset, radius,
                                           timesteps, n_hid=hid,
                                           function=lambda x, y: x/y)

                inter = get_interest(names, interest)

                data = sorted(list(zip(comp, score, names, inter)),
                              key=lambda x: -x[1])

                train = [i for i in data
                         if i[2] in rand_names[:int(.7*len(data))]]
                test = [i for i in data
                        if i[2] in rand_names[int(.7*len(data)):]]

                sorted_train = sorted(train, key=lambda x: -x[1])

                thresh = sorted_train[
                    np.argmin([sum([(v[3] - 1)**2 for v in sorted_train[:i]]) +
                               sum([(v[3])**2 for v in sorted_train[i:]])
                               for i in range(len(train))])
                ][1]

                positive_proportion = np.sum([v[3] for v in test])/len(test)
                accuracy = sum([(int(v[3] == 1)
                                 if v[1] >= thresh
                                 else int(v[3] == 0))
                                for v in test])/len(test)
                baseline_accuracy = sum([1
                                         if (v[3] == 0)
                                         else 0
                                         for v in test])/len(test)
                print(
                    ("{} | Accuracy: {:.3f}, Baseline: {:.2f}, " +
                     "Prop: {:.2f}").format(radius, accuracy,
                                            baseline_accuracy,
                                            positive_proportion))

                accuracies.append(accuracy)

def main():
    """ Create and process items from the data directory. """
    results_dic = {}

    # List all files for which we have data
    file_list = [t for t in os.listdir(os.path.join(DATA_DIR, 'nn/'))
                 if t.endswith('dat') and t.startswith("nn")]

    # Process NN files
    for i in file_list:
        process_nn_file(i, results_dic)

    results_dic_dataset, interest = extract_labels(results_dic)

    score_lookup = {
        'score_2_300': {},
        'score_2_50': {},
        'score_2_5': {},
        'score_1_300': {},
        'score_1_50': {},
        'score_1_5': {}
    }
    compressed_length = {}

    for name in results_dic_dataset:
        with open(os.path.join(os.path.join(DATA_DIR, "ent/"),
                               "ent" + name + ".dat")) as ent_file:
            update_lookup(ent_file, score_lookup, name)

        with open(os.path.join(os.path.join(DATA_DIR, "out/"),
                               "out" + name + ".dat")) as out_file:

            content = out_file.read().strip().split('\n')[-1]
            if len(content.split('    ')) >= 5:
                compressed_length[name] = int(content.split('    ')[-1])

if __name__ == "__main__":
    main()
