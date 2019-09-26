from sys import setrecursionlimit

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

from functions.create_html_summary import create_html_summary

from constants import TOTAL_LIST_PATH, SUMMARY


def make_summary():
    df = pd.read_csv(TOTAL_LIST_PATH, index_col='gene').T

    n_genomes = df.shape[0]
    n_clusters = df.shape[1]

    # Plot of number of clusters by size threshold
    n_good_cols = {}
    col_sums = np.array([df[col].sum() for col in tqdm(df.columns)])
    for t in tqdm(np.arange(100)):
        good_cols = (col_sums > t).sum()
        n_good_cols[t] = good_cols

    plt.plot(list(n_good_cols.keys()), list(n_good_cols.values()), 'r-')
    plt.title('Number of clusters by size threshold')
    plt.xlabel('Size threshold')
    plt.ylabel('Number of clusters')
    plt.grid()
    plt.savefig('./plots/n_clusters_by_threshold.png', bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()

    # Plot distribution of genes in clusters
    genes_in_cluster = {}

    for col in tqdm(df.columns):
        number_of_genes = (df[col] > 0).sum()
        genes_in_cluster[col] = number_of_genes
    n_genes = np.array(list(genes_in_cluster.values()))

    plt.figure(figsize=(10, 4))
    plt.subplot(121)
    sns.distplot(n_genes, kde=False)
    plt.title('Distribution of the sizes of clusters')
    plt.xlabel('Number of genes')
    plt.ylabel('Number of clusters')

    plt.subplot(122)
    sns.distplot(n_genes[n_genes > 5], kde=False)
    plt.title('Distribution of the sizes of clusters with size > 5')
    plt.xlabel('Number of genes')
    plt.ylabel('Number of clusters')
    plt.savefig('./plots/clusters_sizes.png', bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()

    good_cols = np.array([col for col in df.columns if (df[col] > 0).sum() > 5])
    n_big_clusters = len(good_cols)

    # Clean df
    df = df[good_cols]

    # Plot matrix of genes presence
    plt.figure(figsize=(25, 5), dpi=300)
    sns.heatmap((df > 0), cbar=False)  # 1 if genome has a gene, 0 otherwise
    plt.savefig('./plots/genes_presence.png', bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()

    # Plot distribution of number of clusters in genomes
    genes_in_genome = np.array([(df.T[col] > 0).sum() for col in df.T.columns])  # TODO: this seems to be stupid
    sns.distplot(genes_in_genome, kde=False)
    plt.title('Distribution of number of genes in genomes')
    plt.xlabel('Number of clusters')
    plt.ylabel('Number of genomes')
    plt.savefig('./plots/clusters_in_genome.png', bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()

    setrecursionlimit(2500)  # This can kill your RAM a bit, but clustering of huge arrays is hard with default limit
    try:
        # Plot clustermap
        plt.figure(figsize=(25, 5), dpi=300)
        sns.clustermap((df > 0), metric='cityblock', cbar=False)
        plt.savefig('./plots/genes_presence_clustering.png', bbox_inches='tight')
    except RecursionError:
        print('Clustering was not performed because of the recursion error')

    create_html_summary(n_genomes, n_clusters, n_big_clusters)
