import pandas as pd
import os
# from utils.config import TROPPO_RESULTS_PATH, PROJECT_PATH
import matplotlib.pyplot as plt
# import prince
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.manifold import TSNE
import sys
import seaborn as sns
import numpy as np
import seaborn as sns


def get_res_summary(troppo_res):
    output_file = os.path.join(TROPPO_RESULTS_PATH, troppo_res)
    df = pd.read_csv(output_file, index_col=0)
    samples = list(df.index)
    number_of_reacs = {}
    list_of_reacs = {}
    number_of_unique_reactions = {}

    for ind in df.index:
        row = df.loc[ind, :]
        reactions = list(row[row].index)

        number_of_reacs[ind] = len(reactions)
        list_of_reacs[ind] = reactions

    print('Number of reaction in each model')
    for k in number_of_reacs:
        print(k, ':', number_of_reacs[k], 'reactions')

    full_list = []
    for v in list_of_reacs.values():
        full_list.extend(v)

    unique_reactions = [r for r in full_list if full_list.count(r) == 1]

    for ur in unique_reactions:
        for sample in list_of_reacs:
            if sample not in number_of_unique_reactions:
                number_of_unique_reactions[sample] = 0
            reacs = list_of_reacs[sample]
            if ur in reacs:
                number_of_unique_reactions[sample] += 1

    print('Number of unique reactions in each model')
    for k in number_of_unique_reactions:
        print(k, ':', number_of_unique_reactions[k], 'unique reactions')

    pairwise_comparison = {}

    for i in range(0, len(samples) - 1):
        for j in range(i + 1, len(samples)):
            model1 = list_of_reacs[samples[i]]
            model2 = list_of_reacs[samples[j]]
            in_comon = set(model1).intersection(model2)

            key = (samples[i], samples[j])

            pairwise_comparison[key] = in_comon

    print('Number of reactions in common between two models')
    for k in pairwise_comparison:
        print(k, ':', len(pairwise_comparison[k]), 'reactions in common')


# get_res_summary('vvinif2021_v801_GSE36128_fastcore_default.csv')
# get_res_summary('vvinif2021_v801_GSE36128_fastcore_Local2_0_1_4.csv')

def create_binary_matrix(data: pd.DataFrame, path: str):
    output_file = os.path.join(path)
    # df = pd.read_csv(output_file, index_col=0)
    # samples = [x.replace('_t10', '') for x in df.index]
    # df.index = samples
    df = data.astype(int)
    df_filename = output_file.replace('.csv', '_bin.csv')
    df.to_csv(df_filename)


# create_binary_matrix('vvinif2021_v801_GSE36128_fastcore_Local1_0_1_4.csv')

def get_number_of_reactions(data):
    list_of_reacs = {}
    number_of_reacs = {}

    for ind in data.index:
        row = data.loc[ind]
        reactions = row[row == True].index
        reactions_nodrains = [r for r in reactions if not r.startswith('EX_')]
        number_of_reacs[ind] = len(reactions_nodrains)
        list_of_reacs[ind] = reactions_nodrains

    df_total = pd.DataFrame.from_dict(number_of_reacs, orient='index', columns=['Number of reactions'])
    return df_total, list_of_reacs


def get_number_of_reactions_by_group(data, groups):
    list_of_reacs_tissue = {}
    number_of_reacs = {}

    for ind in data.index:
        row = data.loc[ind, :]
        tissue = groups[ind]
        reactions = list(row[row].index)
        if tissue not in list_of_reacs_tissue:
            list_of_reacs_tissue[tissue] = reactions
        else:
            for r in reactions:
                if r not in list_of_reacs_tissue[tissue]:
                    list_of_reacs_tissue[tissue].append(r)

    for t in list_of_reacs_tissue:
        number_of_reacs[t] = len(list_of_reacs_tissue[t])

    df_total = pd.DataFrame.from_dict(number_of_reacs, orient='index', columns=['Number of reactions'])
    return df_total, list_of_reacs_tissue


def get_number_of_unique_reactions(list_of_reacs):
    number_of_unique_reactions = {}
    full_list = []
    for v in list_of_reacs.values():
        full_list.extend(v)

    unique_reactions = [r for r in full_list if full_list.count(r) == 1]
    unique_reactions_noe = [x for x in unique_reactions if not x.startswith('e-')]

    for ur in unique_reactions_noe:
        for sample in list_of_reacs:
            if sample not in number_of_unique_reactions:
                number_of_unique_reactions[sample] = 0
            reacs = list_of_reacs[sample]
            if ur in reacs:
                number_of_unique_reactions[sample] += 1

    new_df = pd.DataFrame.from_dict(number_of_unique_reactions, orient='index', columns=['Unique reactions'])
    return new_df


# troppo_res = 'vvinif2021_v801_GSE36128_fastcore_Local2_0_1_4.csv'
# output_file = os.path.join(TROPPO_RESULTS_PATH, troppo_res)
# df = pd.read_csv(output_file, index_col=0)
# samples = list(df.index)
# samples = [x.replace('_Local2_0_1_4_fastcore_t0_35', '').replace('vvinif2021_v801_', '') for x in samples]
# df.index = samples
#
# tissues = {'leaf_FS': 'leaf', 'leaf_S': 'leaf', 'leaf_Y': 'leaf',
#            'stem_green': 'stem', 'stem_woody': 'stem',
#            'berry_PFS': 'berry', 'berry_R': 'berry'}
#
# df_number, list_reactions = get_number_of_reactions_by_tissue(df, tissues)
# unique = get_number_of_unique_reactions(list_reactions)
# print(unique)


def get_pair_comparison(list_of_reacs, outputname, samples):
    pairwise_comparison = {}

    df_pair = pd.DataFrame(columns=samples, index=samples)

    for i in range(0, len(samples) - 1):
        line = []
        totali = len(list_of_reacs[samples[i]])
        df_pair.loc[samples[i], samples[i]] = totali

        for j in range(i + 1, len(samples)):
            model1 = list_of_reacs[samples[i]]
            model2 = list_of_reacs[samples[j]]
            in_common = set(model1).intersection(model2)

            key = (samples[i], samples[j])
            pairwise_comparison[key] = in_common

            df_pair.loc[samples[i], samples[j]] = len(in_common)
            df_pair.loc[samples[j], samples[i]] = len(in_common)

            totalj = len(list_of_reacs[samples[j]])
            df_pair.loc[samples[j], samples[j]] = totalj

    df_pair.to_csv(os.path.join(TROPPO_RESULTS_PATH, outputname))
    return df_pair


def run_truncatedsvd(n_components: int, data: pd.DataFrame):
    svd = TruncatedSVD(n_components=n_components)

    sv = svd.fit_transform(data)

    columns = [f'PC {i + 1}' for i in range(n_components)]

    df_svd = pd.DataFrame(data=sv, index=data.index, columns=columns)
    explained_variance = svd.explained_variance_ratio_

    return df_svd, explained_variance


def run_pca(n_components: int, data: pd.DataFrame):
    pca = PCA(n_components=n_components)

    pc = pca.fit_transform(data)

    columns = [f'PC {i + 1}' for i in range(n_components)]

    df_pca = pd.DataFrame(data=pc, index=data.index, columns=columns)
    explained_variance = pca.explained_variance_ratio_

    names = pd.DataFrame(abs(pca.components_), columns=data.columns, index=columns)

    loadings = pd.DataFrame(
        data=pca.components_.T * np.sqrt(pca.explained_variance_),
        columns=columns,
        index=data.columns
    )

    return df_pca, explained_variance, names, loadings


def run_tsne(n_components: int, data: pd.DataFrame):
    tsne = TSNE(n_components, random_state=42, perplexity=30)
    tsne_result = tsne.fit_transform(data)

    columns = [f'tsne {i + 1}' for i in range(n_components)]

    df_tsne = pd.DataFrame(data=tsne_result, index=data.index, columns=columns)

    return df_tsne


def run_mca(n_components: int, data: pd.DataFrame):
    data = data.astype(str)
    mca = prince.MCA(n_components=n_components)
    mc = mca.fit(data)
    df_mca = mc.transform(data)
    df_mca.columns = [f'PC {i + 1}' for i in range(n_components)]
    explained_inertia = mca.explained_inertia_
    return df_mca, explained_inertia


def plot_pca(data: pd.DataFrame, explained_variance: list, c1: str, c2: str, title: str, name_fig: str):
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_title(title, fontsize=20)
    x_label = f'{c1} ({round(explained_variance[0] * 100, 2)} %)'
    y_label = f'{c2} ({round(explained_variance[1] * 100, 2)} %)'
    ax.set_xlabel(x_label, fontsize=15)
    ax.set_ylabel(y_label, fontsize=15)

    labels = set(data.loc[:, 'factor'])

    if len(labels) == 3:
        labels_dic = {2012: 'g', 2013: 'b', 2014: 'r'}
    elif len(labels) == 4:
        labels_dic = {'stem': 'b', 'leaf': 'g', 'berry_green': 'm', 'berry_mature': 'r'}
    else:
        labels_dic = {'mature': 'r', 'green': 'g'}

    for label, color in labels_dic.items():
        mask = data.loc[:, 'factor'] == label

        pc1 = data.loc[mask, c1]
        pc2 = data.loc[mask, c2]
        ax.scatter(pc1, pc2, c=color, s=60)

        # names = pc1.index
        # for pc1_pt, pc2_pt, annotation in zip(pc1, pc2, names):
        #     ax.annotate(annotation, (pc1_pt, pc2_pt))

    ax.legend(labels_dic.keys(), loc='best')
    ax.grid()
    fig_path = os.path.join(name_fig + '.png')
    plt.savefig(fig_path)
    plt.show()


def plot_tsne(data: pd.DataFrame, name_fig: str, title: str):
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    sns.set_style("darkgrid")

    sns.scatterplot(x='tsne 1', y='tsne 2', hue='factor', data=data, ax=ax, s=60, palette=dict(stem="#1f77b4", leaf="#ff7f0e", berry_green="#2ca02c",
                                                                                               berry_mature="#d62728"))

    #
    # for label, color in labels_dic.items():
    #     mask = data.loc[:, 'factor'] == label
    #
    #     pc1 = data.loc[mask, ['tsne 1', 'tsne 2']]
    #     # pc2 = data.loc[mask, 'tsne 2']
    #     sns.scatterplot(pc1, c=color, s=60, ax=ax)

    # lim = (data.loc[:, ['tsne 1', 'tsne 2']].to_numpy().min() - 3,
    #        data.loc[:, ['tsne 1', 'tsne 2']].to_numpy().max() + 3)

    # ax.set_xlim(lim)
    # ax.set_ylim(lim)

    # ax.set_xticks(range(-9, 17, 3))
    # ax.set_yticks(range(-15, 17, 3))

    # ax.set_aspect('equal')
    ax.set_title(title, fontsize=20)
    ax.legend(loc='upper left', borderaxespad=0.0)
    fig_path = os.path.join(name_fig + '.png')
    plt.savefig(fig_path)
    plt.show()


def biplot(data, comps):
    xs = data.iloc[:, 0]
    ys = data.iloc[:, 1]
    n = comps.shape[0]
    scalex = 1.0 / (xs.max() - xs.min())
    scaley = 1.0 / (ys.max() - ys.min())
    plt.scatter(xs * scalex, ys * scaley)

    for i in range(n):
        plt.arrow(0, 0, comps[i, 0], comps[i, 1], color='r', alpha=0.5)
        # if labels is None:
        #     plt.text(comps[i, 0] * 1.15, comps[i, 1] * 1.15, "Var" + str(i + 1), color='g', ha='center', va='center')
        # else:
        #     plt.text(comps[i, 0] * 1.15, comps[i, 1] * 1.15, comps[i], color='g', ha='center', va='center')
    plt.xlim(-1, 1)
    plt.ylim(-1, 1)
    plt.xlabel("PC{}".format(1))
    plt.ylabel("PC{}".format(2))
    plt.grid()
    plt.show()


PROJECT_PATH = 'C:/Users/BiSBII/Documents/plantdb'


def get_dfa_reactions():
    folder = os.path.join(PROJECT_PATH, 'reconstruction_results', 'vvinif2023', 'results_troppo', 'RNAseq',
                          'dfa')

    files = [f for f in os.listdir(folder) if f.endswith('reaction_result_all.csv')]

    all_reactions = []

    for f in files:
        df = pd.read_csv(os.path.join(folder, f))
        reactions = df['Reaction'].tolist()
        for r in reactions:
            if r not in all_reactions:
                all_reactions.append(r)

    return all_reactions
