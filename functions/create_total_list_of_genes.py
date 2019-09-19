from collections import Counter, defaultdict
import json

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm

from functions.perform_cd_hit_clustering import perform_cd_hit_clustering
from constants import THRESHOLD, CLUSTERS_INFO, TEMP_PATH, TOTAL_LIST_PATH, TOTAL_LIST_LOG


def create_total_list_of_genes():  # TODO: document

    log = open(TOTAL_LIST_LOG, 'w', buffering=1)

    log.write('Performing CD-hit clustering\n')
    clusters_df = perform_cd_hit_clustering(THRESHOLD)
    log.write('Done. Shape of a dataframe: {}'.format(clusters_df.shape))

    clusters_sizes = clusters_df['cluster'].value_counts()

    size_threshold = int(0.02*clusters_sizes[-1])  # 2% smallest clusters
    filtered_clusters = {cluster: clusters_sizes[cluster] for cluster in clusters_sizes
                         if clusters_sizes[cluster] > size_threshold}

    # Save distribution of clusters sizes
    plt.figure(figsize=(10, 3))

    plt.subplot(121)
    sns.distplot(clusters_sizes, kde=False)
    plt.title('Distribution of sizes of {} clusters'.format(len(clusters_sizes)))
    plt.xlabel('Size of a cluster')
    plt.ylabel('Number of clusters')

    plt.subplot(122)
    sns.distplot(list(filtered_clusters.values()), kde=False)
    plt.title('Clusters, which size > {}, sizes distribution'.format(size_threshold))
    plt.xlabel('Size of a cluster')
    plt.ylabel('Number of clusters')

    plt.savefig('./plots/clusters_sizes.png')

    clusters_info = {}
    for cluster in filtered_clusters.keys():
        all_info = Counter(clusters_df[clusters_df['cluster'] == cluster]['info'])
        cluster_info = all_info.most_common(1)[0][0]  # The name of most common element
        clusters_info[cluster] = cluster_info

    with open(CLUSTERS_INFO, 'w') as f:
        json.dump(clusters_info, f)

    genomes = list(set(clusters_df['genome']))
    COLS = ['gene'] + genomes

    total_list = pd.DataFrame(data=None, columns=COLS)
    temp_df = pd.DataFrame(data=None, columns=COLS)

    log.write('Creating a total list')
    for i, cluster in enumerate(tqdm(filtered_clusters.keys())):
        cluster_genomes = np.array(clusters_df[clusters_df['cluster'] == cluster]['genome'])
        genes_in_genome = defaultdict(int)

        for genome in genomes:
            genes_in_genome[genome] += (cluster_genomes == genome).sum()

        row = pd.DataFrame(data=[[cluster] + list(genes_in_genome.values())], columns=COLS)
        temp_df = pd.concat([temp_df, row], axis=0, ignore_index=True, sort=False)

        # Concat each 500 rows to total list, so we don't recreate it every iteration
        if i % 500 == 0:
            temp_df.to_csv(TEMP_PATH, index=False)
            total_list = pd.concat([total_list, temp_df], axis=0, ignore_index=True, sort=False)
            temp_df = pd.DataFrame(data=None, columns=COLS)

    total_list = pd.concat([total_list, temp_df], axis=0, ignore_index=True, sort=False)
    total_list.to_csv(TOTAL_LIST_PATH, index=False)
    log.write('Done. Shape of a dataframe: {}'.format(total_list.shape))

    log.close()
