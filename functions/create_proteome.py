import os

from Genome import Genome

from constants import PROTEOME_PATH, GENOMES_DIR


def create_proteome():
    """Create one file with protein sequences from all the genomes in data folder"""

    # Get all directories in GENOMES_DIR
    folders = list(filter(lambda x: os.path.isdir(GENOMES_DIR + x), os.listdir(GENOMES_DIR)))

    with open(PROTEOME_PATH, 'a') as f:

        for folder in folders:
            prefix = folder + '_'

            genome = Genome(GENOMES_DIR+folder, prefix)
            genome.read_protein_faa()
            genome.read_feature_table()
            genome.set_gene_positions()

            for gene in genome.genes:
                f.write(gene.get_fasta())
