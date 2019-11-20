""" A module that generates the results found in the paper.
"""
import os
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score

COLS = ["ID", "score", "interesting", "compressed_len", "lookup_table"]
DATA_DIR = os.path.join(os.getcwd(), "data/")

def read_data(name):
    """ Read data from the CSV """
    dataset = pd.read_csv(os.path.join(DATA_DIR, name),
                          index_col=0,
                          dtype=str)
    for i in COLS:
        if i not in dataset.columns:
            raise ValueError(
                "Column {} not found in dataset {}".format(i, name))

    dataset = dataset[COLS].astype({
        "ID": str,
        "score": float,
        "interesting": int,
        "compressed_len": int,
        "lookup_table": float
    })
    return dataset


def compute_average_precision(data, y_true):
    """ Compute average precision of each of the datasets """
    avg_prec = {}
    for i in data.columns[1:]:
        avg_prec[i] = average_precision_score(y_true, data[i].values)
    return avg_prec

def compute_accuracy_with_column(col, train_data, test_data):
    """ Compute the accuracy of a step classifier that chooses a thresholds
    minimizing MSE.
    """
    sorted_train = sorted(train_data.values, key=lambda x: -x[col])

    # Find the threshold values that minimizes MSE
    thresh = sorted_train[
        np.argmin([
            sum([(v[2] - 1)**2 for v in sorted_train[:i]]) +
            sum([(v[2])**2 for v in sorted_train[i:]])
            for i in range(len(train_data))
        ])][col]

    test_accuracy = sum([(int(v[2] == 1)
                          if v[col] >= thresh
                          else int(v[2] == 0))
                         for v in test_data.values])/len(test_data)
    return test_accuracy


def compute_accuracy_with_double_threshold(train_data, test_data):
    """ Compute the accuracy of a classifier that chooses a pair a thresholds
    minimizing MSE when all examples in between are positive.
    """
    sorted_train_comp = sorted(train_data.values, key=lambda x: -x[3])
    losses = []

    # Compute losses for all possible pairs of threshold
    for i in range(len(sorted_train_comp) - 1):
        for j in range(i + 1, len(sorted_train_comp) + 1):
            loss = sum([(1 - v[2])**2 for v in sorted_train_comp[i:j]])
            loss += sum([(v[2])**2 for v in sorted_train_comp[:i]])
            loss += sum([(v[2])**2 for v in sorted_train_comp[j:]])
            losses.append((i, j, loss))

    # Find corresponding thresholds
    idx0, idx1, _ = sorted(losses, key=lambda x: x[2])[0]
    tr0 = sorted_train_comp[idx0][3]
    tr1 = sorted_train_comp[idx1][3]

    accuracy = sum([(int(v[2] == 1)
                     if (v[3] <= tr0 and v[3] >= tr1)
                     else int(v[2] == 0))
                    for v in test_data.values])/len(test_data)

    return accuracy


def main():
    """ Compute accuracies """
    train_data = read_data("train_data.csv")
    test_data = read_data("test_data.csv")

    # Proportion of positive elements
    pos_proportion = np.sum(test_data["interesting"].sum())/len(test_data)

    # Baseline accuracy if everything is set to 0
    baseline_accuracy = sum([(v[2] == 0)
                             for v in test_data.values])/len(test_data)

    # Compute accuracy with the score
    score_accuracy = compute_accuracy_with_column(1, train_data, test_data)

    # Compute accuracy with the score
    lookup_accuracy = compute_accuracy_with_column(4, train_data, test_data)

    # Double threshold accuracy
    dbl_thresh_accuracy = compute_accuracy_with_double_threshold(train_data,
                                                                 test_data)


    print("Proportion of positives: {:.1f}%".format(100*pos_proportion))
    print("Accuracy baseline: {:.1f}%".format(100*baseline_accuracy))
    print("Double threshold accuracy: {:.1f}%".format(100*dbl_thresh_accuracy))
    print("Lookup accuracy: {:.1f}%".format(100*lookup_accuracy))
    print("Score accuracy: {:.1f}%".format(100*score_accuracy))

if __name__ == "__main__":
    main()
