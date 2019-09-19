import numpy as np
import pandas as pd


def read_cd_hit_output(path):
    """Create and return table from cd-hit output file.

    Parameters
    ----------
    path : str
        Path to .clstr file with cd-hit output

    Returns
    -------
    clusters : pandas.DataFrame
        DataFrame with following columns:
        1. 'cluster' — cluster of a sequence
        2. 'representative' — boolean variable, which indicates, if the given sequence is representative in the cluster
        3. 'gene' — name of the sequence from fasta header
        4. 'length' — length of the gene's sequence
        5. 'identity' — identity of a sequence to representative sequence of its cluster
    """
    clusters = []

    with open(path) as f:
        file = f.read()

    for line in file.rstrip().split('\n'):
        if line.startswith('>Cluster'):
            current_cluster = line[1:]
            continue

        length = line[line.find('\t') + 1: line.find('aa')]
        gene_id = line[line.find('>') + 1: line.find('...')]
        is_representative = line.endswith('*')
        identity = 100 if is_representative else line[line.rfind(' ') + 1: -1]

        clusters.append((current_cluster, is_representative, gene_id, length, identity))

    clusters = pd.DataFrame(data=clusters, columns=['cluster', 'representative', 'gene', 'length', 'identity'])
    clusters['identity'] = clusters['identity'].astype(np.float32)
    clusters['length'] = clusters['length'].astype(np.int32)

    return clusters
