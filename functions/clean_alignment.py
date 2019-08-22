import pandas as pd
import numpy as np


def clean_alignment(df, identity_threshold=50):
    """Create data frame with best aligned genes only.

    Best aligned gene is defined as gene, which identity is better than a given threshold and
    which start is closer to a database gene, that start of all other genes in alignment.
    This is controversal. Maybe we should develop better criterion. Also, we should remember, that a chromosome
    is circular. A lot of things TODO better here!

    Parameters
    ----------
    df : pandas.DataFrame
        Data frame in a format returned by map_genome().
    identity_threshold : int, optional
        Threshold for `df['identity']`.

    Returns
    -------
    clean_df : pandas.DataFrame
        Data frame in the same format with only the best aligned genes.
    not_aligned : set
        Set of excess genes, which didn't find corresponding gene in database.
    """
    db_genes = set(df['db_gene'])
    clean_df = pd.DataFrame(data=None, columns=df.columns)
    not_aligned = set()

    for db_gene in db_genes:
        rows = df[(df['db_gene'] == db_gene) & (df['identity'].astype(np.float64) > identity_threshold)]

        if rows.shape[0] == 0:
            continue

        min_delta_starts = min(abs(rows['db_gene_start'] - rows['query_start']))

        best_alignment = rows[abs(rows['db_gene_start'] - rows['query_start']) == min_delta_starts]
        clean_df = pd.concat((clean_df, best_alignment))

        worse_alignments = rows[(rows['db_gene_start'] - rows['query_start']) != min_delta_starts]
        # Choose only those genes, which are not aligned yet
        worse_alignments = worse_alignments[worse_alignments['query_gene'].isin(clean_df['query_gene']) == False]
        not_aligned.update(worse_alignments['query_gene'])

        # Remove best_alignment from not_aligned
        for el in best_alignment[best_alignment['query_gene'].isin(not_aligned)]['query_gene']:
            if el in not_aligned:  # For some reason it can still be False
                not_aligned.remove(el)

    return clean_df, not_aligned
