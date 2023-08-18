#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 10:18:11 2020

@author: meyerct6

Code to build Booleabayes rules for a given network and dataset
Use cross validation and caluculate AUC
"""
# =============================================================================
# If running from command line, Read in the arguments
# =============================================================================
import sys
import numpy as np
if False:
    if len(sys.argv) < 5 or len(sys.argv) > 8:
        raise OSError(
            'Wrong number of arguments.  Must pass <directory> <network file path from parent directory> <data file path from parent directory> <cellID-table> <node_normalization (optional default=.5) node_threshold (optional default=.1)>')
    dir_prefix = sys.argv[1]
    network_path = sys.argv[2]
    data_path = sys.argv[3]
    cellID_table = sys.argv[4]
    node_normalization = .5 if len(sys.argv) >= 5 else float(sys.argv[5])
    node_threshold = .1 if len(sys.argv) == 7 else float(sys.argv[6])
else:  # Use the example data
    dir_prefix = '/Users/sarahmaddox/Dropbox (Vanderbilt)/Parthenon_Pancreas_project/BooleaBayes/'
    network_path = 'Network/TF-Lit-Network/test_network.csv'
    data_path = 'magic_imputed_genes/magicTF_train.csv'
    cellID_table = 'cellID-lookuptable.csv'
    node_normalization = 0.3
    node_threshold = 0.1
    brcd = str(300)

if dir_prefix[-1] != os.sep:
    dir_prefix = dir_prefix + os.sep
if not network_path.endswith('.csv') or not os.path.isfile(dir_prefix + network_path):
    raise Exception('Network path must be a .csv file.  Check file name and location')
if not data_path.endswith('.csv') or not os.path.isfile(dir_prefix + data_path):
    raise Exception('data path must be a .csv file.  Check file name and location')
if cellID_table is not None:
    if not cellID_table.endswith('.csv') or not os.path.isfile(dir_prefix + cellID_table):
        raise Exception('CellID path must be a .csv file.  Check file name and location')

# =============================================================================
# Import packages.
# See spec-file.txt for packages
# To recreate the environment : conda create --name booleabayes --file spec-file.txt
# =============================================================================
import graph_utils_update, graph_fit
import numpy as np
import os
import os.path as op
import pandas as pd
import resource
import time
import sklearn.model_selection as ms
from matplotlib import pyplot as plt
from graph_fit import reorder_binary_decision_tree
import seaborn as sns
import os
from scipy import stats
from matplotlib.patches import Patch

### test accuracy of predictions vs. actual expression
customPalette = sns.color_palette('tab10')


# =============================================================================
# Functions used for calculating AUC
# =============================================================================
# run in a loop for each gene
def parent_heatmap(data, regulators_dict, gene):

    regulators = [i for i in regulators_dict[gene]]
    n = len(regulators)

    # This is the distribution of how much each sample reflects/constrains each leaf of the Binary Decision Diagram
    heat = np.ones((data.shape[0], 2 ** n))
    for leaf in range(2 ** n):
        binary = graph_utils_update.idx2binary(leaf, len(regulators))
        binary = [{'0': False, '1': True}[i] for i in binary]
        # binary becomes a list of lists of T and Fs to represent each column
        for i, idx in enumerate(data.index):
            # for each row in data column...
            # grab that row (df) and the expression value for the current node (left side of rule plot) (val)
            df = data.loc[idx]
            val = np.float(data.loc[idx, gene])
            for col, on in enumerate(binary):

                # for each regulator in each column in decision tree...
                regulator = regulators[col]
                # if that regulator is on in the decision tree, multiply the weight in the heatmap for that
                # row of data and column of tree with a weight that = probability that that node is on in the data
                # df(regulator) = expression value of regulator in data for that row
                # multiply for each regulator (parent TF) in leaf
                if on:
                    heat[i, leaf] *= np.float(df[regulator])
                else:
                    heat[i, leaf] *= 1 - np.float(df[regulator])

    regulator_order = [i for i in regulators]

    return heat, regulator_order


def r2(x, y):
    return stats.pearsonr(x, y)[0] ** 2

def plot_accuracy(data, g, regulators_dict, rules, phenotypes=None, plot_clusters=False, dir_prefix=None,
                  clusters=None, save_plots=None, plot=False, save_df=False, ):
    try:
        os.mkdir(
            op.join(f'{dir_prefix}', f"{save_plots}"))
    except FileExistsError:
        pass

    h, order = parent_heatmap(data, regulators_dict, g)
    # print("Order",order)
    # print(f"Regulators_dict[{g}]", regulators_dict[g])
    # importance_order = reorder_binary_decision_tree(order, regulators_dict[g])
    rule = rules[g]
    # dot product of weights of test sample and rule will give the predicted value for that sample for that TF
    predicted = (np.dot(h, rule))
    p = pd.DataFrame(predicted, columns=['predicted'], index=data.index)
    p['actual'] = data[g]
    if save_df == True:
        p.to_csv(f'{dir_prefix}/{save_plots}/{g}_validation.csv')
    if plot == True:
        if plot_clusters == True:
            plt.figure()
            predicted = pd.DataFrame(predicted, index=data.index, columns=['predicted'])

            for i in set(clusters['class']):
                clines = data.loc[clusters.loc[clusters['class'] == i].index].index
                sns.scatterplot(x=data.loc[clines][g], y=predicted.loc[clines]['predicted'],
                                label=phenotypes[int(i - 1)])
            plt.xlabel("Actual Normalized Expression")
            plt.ylabel("Predicted Expression from Rule")
            legend_elements = []

            for i, j in enumerate(phenotypes):
                legend_elements.append(Patch(facecolor=customPalette[i], label=j))

            plt.legend(handles=legend_elements, loc='best')
            plt.title(str(g))
            plt.savefig(
                f'{dir_prefix}/{save_plots}/{g}_{save_plots}.pdf')
            plt.close()
        else:
            plt.figure()
            sns.regplot(x=data[g], y=predicted)
            plt.xlabel("Actual Normalized Expression")
            plt.ylabel("Predicted Expression from Rule")
            plt.title(str(g))

            if r2(data[g], predicted) == 0:
                plt.title(str(g))
            else:
                plt.title(str(g) + "\n" + str(round(r2(data[g], predicted), 2)))
            plt.savefig(f'{dir_prefix}/{save_plots}/{g}_{save_plots}.pdf')
            plt.show()
            plt.close()
    return p

def roc(validation, g, n_thresholds, save_plots, plot=False, save=False, dir_prefix=None):
    tprs = []
    fprs = []
    for i in np.linspace(0, 1, n_thresholds, endpoint=False):
        p, r = calc_roc(validation, i)
        tprs.append(p)
        fprs.append(r)
    # area = auc(fprs, tprs) #### AUC function wasn't working... replace with np.trapz
    area = np.abs(np.trapz(x = fprs, y = tprs))
    if plot == True:
        fig = plt.figure()
        ax = plt.subplot()
        plt.plot(fprs, tprs, '-', marker = 'o')
        plt.title(g + " ROC Curve" + "\n AUC: " + str(round(area,3)))
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        # ax.plot([0,1], [1,1], ls="--", c=".3")
        # ax.plot([1,1], [0,1], ls="--", c=".2")
        ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3")

        plt.ylabel("True Positive Rate")
        plt.xlabel("False Positive Rate")
        if save == True:
            plt.savefig(f'{dir_prefix}/{save_plots}/{g}_roc.pdf')
            # plt.close()
        else:
            plt.show()
    return tprs, fprs, area

def calc_roc(validation, threshold):
    # P: True positive over predicted condition positive (of the ones predicted positive, how many are actually
    # positive?)

    # R: True positive over all condition positive (of the actually positive, how many are predicted to be positive?)
    predicted = validation.loc[validation['predicted'] > threshold]
    actual = validation.loc[validation['actual'] > 0.5]
    predicted_neg = validation.loc[validation['predicted'] <= threshold]
    actual_neg = validation.loc[validation['actual'] <= 0.5]
    true_positive = len(set(actual.index).intersection(set(predicted.index)))
    false_positive = len(set(actual_neg.index).intersection(set(predicted.index)))
    true_negative = len(set(actual_neg.index).intersection(set(predicted_neg.index)))
    if len(actual.index.values) == 0 or len(actual_neg.index.values) == 0:
        return -1, -1
    else:
        # print((true_positive+true_negative)/(len(validation)))
        tpr = true_positive / len(actual.index)
        fpr = false_positive / len(actual_neg.index)
        return tpr, fpr

def auc(fpr, tpr): #this function is broken for some reason UGH
    # fpr is x axis, tpr is y axis
    print("Calculating area under discrete ROC curve")
    area = 0
    i_old, j_old = 0, 0
    for c, i in enumerate(fpr):
        j = tpr[c]
        if c == 0:
            i_old = i
            j_old = j
        else:
            area += np.abs(i - i_old) * j_old + .5 * np.abs(i - i_old) * np.abs(j - j_old)
            i_old = i
            j_old = j
    return area


# =============================================================================
# Start the timer and generate a barcode identifier for this job
# =============================================================================
t1 = time.time()
# brcd = str(np.random.randint(0, 10001))
# Test if the barcode is in use.  If it is generate a new one.
while True:
    if os.path.isdir(dir_prefix + brcd):
        brcd = str(np.random.randint(0, 10001))
    else:
        os.mkdir(dir_prefix + brcd)
        break
# Set random seed for reproducibility
np.random.seed(1)
print(brcd)
# =============================================================================
# Load the network
# =============================================================================
graph, vertex_dict = graph_utils_update.load_network(op.join(dir_prefix, f'{network_path}'), remove_sinks=False,
                                                     remove_sources=False)
v_names = graph.vertex_properties['name']  # Map of vertex -> name (basically the inverse of vertex_dict)
nodes = sorted(vertex_dict.keys())
print('Reading in data')
data = graph_utils_update.load_data(op.join(dir_prefix, f"{data_path}"), nodes, norm=node_normalization, delimiter=',',
                                    log1p=False, transpose=True, sample_order=False, fillna = 0)

# =============================================================================
# Load the cell clusters
# =============================================================================
print('Reading in cell cluster labels')
if cellID_table is not None:
    clusters = pd.read_csv(op.join(dir_prefix, f"{cellID_table}"), index_col=0, header=0, delimiter=',')
    # clusters.columns = ["cell.stage", "class", "cell.protocol"]
    clusters.columns = ["class", "pheno_color","Tuft_score","nonNE_score","NE_score","NEv2_score","NEv1_score","treatment"]

    if sum(np.in1d(data.index, clusters.index)) != len(data):
        raise Exception('The cells in the expression data do not match the cell label data')
else:
    clusters = pd.DataFrame([0]*len(data.index), index = data.index, columns=['class'])

# =============================================================================
# Binarize data and save to a file
# =============================================================================
print('Binarizing data')
binarized_data = graph_utils_update.binarize_data(data, phenotype_labels=clusters)

with open(dir_prefix + brcd + os.sep + 'binarized_data_' + brcd + '.csv', 'w+') as outfile:
    for k in binarized_data.keys():
        outfile.write(f"{k}: {binarized_data[k]}\n")

# =============================================================================
# Generate and save rules with 5-fold cross validation
# Do density dependent sampling
# =============================================================================
n = len(nodes)
n_splits = 5
kf = ms.StratifiedKFold(n_splits=n_splits)
test_set = 'validation_set'
os.mkdir(f"{dir_prefix}{brcd}/{test_set}/")

i = 0  # Index of fold
aucs = pd.DataFrame(index=nodes, columns=[str(i) for i in range(n_splits)])
for train_index, test_index in kf.split(data.index, clusters.loc[data.index, 'class']):
    print(f'Generating Rules for K-Fold {i}')

    try:
        os.mkdir(f"{dir_prefix}{brcd}/{test_set}/{i}/")
    except FileExistsError:
        pass

    # save rules for each barcode and cross-validation fold
    rules, regulators_dict = graph_fit.get_rules(data.iloc[train_index], vertex_dict,
                                                 directory=dir_prefix + brcd + os.sep + "rules_" + brcd + '_' + str(i),
                                                 plot=False, threshold=node_threshold)
    graph_fit.save_rules(rules, regulators_dict, fname=f"{dir_prefix}{brcd}/{test_set}/{i}/rules_{brcd}_{i}.txt")


    # make one file of tprs (true positive rate) and fprs for each cross-validation fold


    outfile = open(f"{dir_prefix}{brcd}/{test_set}/{i}/tprs_fprs_{brcd}_{i}.csv", 'w+')
    ind = [x for x in np.linspace(0, 1, 50)]
    tpr_all = pd.DataFrame(index=ind)
    fpr_all = pd.DataFrame(index=ind)
    area_all = []

    outfile.write(f",,")
    for j in ind:
        outfile.write(str(j)+',')
    outfile.write('\n')
    for g in nodes:
        print(g)

        validation = plot_accuracy(data.iloc[test_index], g, regulators_dict, rules, save_plots=i,
                                   plot=False, plot_clusters=False, save_df=True,
                                   dir_prefix=dir_prefix + brcd + os.sep + str(test_set) + os.sep)
        tprs, fprs, area = roc(validation, g, n_thresholds=50, save_plots=i, plot=True, save=True,
                               dir_prefix=dir_prefix + brcd + os.sep + str(test_set) + os.sep)
        tpr_all[g] = tprs
        fpr_all[g] = fprs
        outfile.write(f"{g},tprs,{tprs}\n")
        outfile.write(f"{g},fprs,{fprs}\n")
        area_all.append(area)
    outfile.close()

    # save AUC values by gene for each cross-validation fold and barcode
    outfile = open(f"{dir_prefix}{brcd}/{test_set}/{i}/auc_{brcd}_{i}.csv", 'w+')
    for n, a in enumerate(area_all):
        outfile.write(f"{nodes[n]},{a} \n")
        aucs.loc[nodes[n]][str(i)] = a
    outfile.close()
    i += 1
aucs.to_csv(f"{dir_prefix}{brcd}/{test_set}/aucs.csv")

# =============================================================================
# Write out information about the this job
# =============================================================================
# Append the results to a MasterResults file
T = {}
t2 = time.time()
T['time'] = (t2 - t1) / 60.
# How much memory did I use?   Only can use on linux platform
if os.name == 'posix':
    T['memory_Mb'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000.
else:
    T['memory_Mb'] = np.nan
T['barcode'] = brcd
T['dir_prefix'] = dir_prefix
T['network_path'] = network_path
T['data_path'] = data_path
T['cellID_table'] = cellID_table
T['node_normalization'] = node_normalization
T['node_threshold'] = node_threshold

T = pd.DataFrame([T])
if not os.path.isfile(dir_prefix + 'Job_specs.csv'):
    T.to_csv(dir_prefix + 'Job_specs.csv')
else:
    with open(dir_prefix + 'Job_specs.csv', 'a') as f:
        T.to_csv(f, header=False)