from collections import Counter, defaultdict
import json

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm

from functions.perform_cd_hit_clustering import perform_cd_hit_clustering
from constants import THRESHOLD, CLUSTERS_INFO, TEMP_PATH, TOTAL_LIST_PATH, TOTAL_LIST_LOG, CLUSTERS_DF_PATH


def create_total_list_of_genes():  # TODO: document

    log = open(TOTAL_LIST_LOG, 'w', buffering=1)

    log.write('Performing CD-hit clustering\n')
    clusters_df = perform_cd_hit_clustering(THRESHOLD)
    log.write('Done. Shape of a dataframe: {}\n'.format(clusters_df.shape))
    clusters_df.to_csv(CLUSTERS_DF_PATH, index=False)

    clusters = set(clusters_df['cluster'])

    clusters_info = {}
    for cluster in clusters:
        all_info = Counter(clusters_df[clusters_df['cluster'] == cluster]['info'])
        clusters_info[cluster] = all_info

    with open(CLUSTERS_INFO, 'w') as f:  # Write all the info about clusters in .json file
        json.dump(clusters_info, f)

    genomes = list(set(clusters_df['genome']))
    COLS = ['gene'] + genomes

    total_list = pd.DataFrame(data=None, columns=COLS)
    temp_df = pd.DataFrame(data=None, columns=COLS)

    log.write('Creating a total list\n')
    for i, cluster in enumerate(tqdm(clusters)):
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
    log.write('Done. Shape of the total list: {}\n'.format(total_list.shape))

    log.close()
