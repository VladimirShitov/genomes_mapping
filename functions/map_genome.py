import pandas as pd
import numpy as np
from Bio.Blast.Applications import NcbiblastpCommandline

import constants

DB_PATH = constants.DB_PATH
THRESHOLD = constants.THRESHOLD


def map_genome(genome, reference, save_alignment=True):
    """Run BLAST search and return dataframe with well aligned genes.

    Parameters
    ----------
    genome : Genome
        An object of Genome() class, which `genes` will be searched against reference genome
    reference : Genome
        An object of Genome() class, which `genes` are needed to get exact positions of database genes.
    save_alignment : bool
        If True, good alignments will be saved to alignment_df and return. If False, return None instead of alignment_df

    Returns
    -------
    alignment_df : pandas.DataFrame
        Data frame with following columns:
        [0]  query_gene - id of a gene from `genome`, which had a good enough BLAST hit
        [1]  db_gene - id of a gene from `reference`, on which query_gene has aligned
        [2]  identity - % of identity between query_gene and db_gene
        [3]  alignment_length
        [4]  mismatches
        [5]  gap_opens
        [6]  query_start - start position of a query_gene
        [7]  query_end - end position of a query_gene
        [8]  db_gene_start
        [9]  db_gene_end
        [10] E_value
        [11] bit_score
        [12] query_genome - assembly number of a genome
        [13] db_genome
        If `save_alignment` is False, alignment_df will not be created and None will be returned.
    not_aligned_genes : set
        Set with genes, which didn't have good enough BLAST hit.
    """
    mapping_list = []
    not_aligned_genes = set()

    cline = NcbiblastpCommandline(query=genome.files['protein.faa'],
                                  db='{}cur_db'.format(DB_PATH),
                                  evalue=0.001,
                                  outfmt=6,  # tab-separated
                                  num_threads=7,
                                  max_target_seqs=3000)
    stdout, stderr = cline()
    found = set()

    if len(stdout) > 1:
        results = stdout.rstrip().split('\n')

        for result in results:
            cols = result.split('\t')
            gene = cols[0]

            if float(cols[2]) > THRESHOLD:
                found.add(gene)

                if save_alignment:
                    query_gene = genome.get_gene_by_id(cols[0])
                    db_gene = reference.get_gene_by_id(cols[1])

                    cols[6] = query_gene.start
                    cols[7] = query_gene.end
                    cols[8] = db_gene.start
                    cols[9] = db_gene.end
                    cols.append(query_gene.genome)
                    cols.append(db_gene.genome)

                    mapping_list.append(cols)
    else:
        raise Exception('Blast failed')

    for gene in genome.genes:
        if gene.id not in found:
            not_aligned_genes.add(gene.id)

    if save_alignment:
        alignment_df = pd.DataFrame(mapping_list, columns=['query_gene', 'db_gene', 'identity',
                                                           'alignment_length', 'mismatches', 'gap_opens',
                                                           'query_start', 'query_end', 'db_gene_start',
                                                           'db_gene_end', 'E_value', 'bit_score',
                                                           'query_genome', 'db_genome'])
        
        alignment_df['identity'] = alignment_df['identity'].astype(np.float32)
        alignment_df['query_start'] = alignment_df['query_start'].astype(np.int32)
        alignment_df['query_end'] = alignment_df['query_end'].astype(np.int32)
        alignment_df['db_gene_start'] = alignment_df['db_gene_start'].astype(np.int32)
        alignment_df['db_gene_end'] = alignment_df['db_gene_end'].astype(np.int32)
        
    else:
        alignment_df = None

    return alignment_df, not_aligned_genes
