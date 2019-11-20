"""
Module for extracting the output data created by the automaton binary upon
running.
"""
import os
import pandas as pd
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
