import os
from tqdm import tqdm  # Progress bar

from Genome import Genome
from functions.map_genome import map_genome
import constants

DB_PATH = constants.DB_PATH
GENOMES_DIR = constants.GENOMES_DIR
LOG_PATH = constants.ALIGNMENT_LOG


def align_to_database():

    reference = Genome(folder_path=DB_PATH, files_prefix='')
    reference.read_protein_faa(filename=DB_PATH+'current_db.faa')
    reference.set_gene_positions(from_fasta=True)

    # Get all directories in GENOMES_DIR
    folders = list(filter(lambda x: os.path.isdir(GENOMES_DIR + x), os.listdir(GENOMES_DIR)))
    log = open(LOG_PATH, 'w', buffering=1)

    for folder in tqdm(folders):

        prefix = folder + '_'

        log.write('Working with {}\n'.format(GENOMES_DIR + folder))
        genome = Genome(folder_path=GENOMES_DIR+folder+'/', files_prefix=prefix)
        genome.read_protein_faa()
        genome.read_feature_table()
        genome.set_gene_positions()

        log.write('Created genome\n')

        if genome.raised_errors:
            log.write('ERROR\n{}\n'.format(genome.log))
            log.close()
            raise Exception(genome.log)

        log.write('Mapping genomes...\n')
        alignment_df, not_aligned = map_genome(genome, reference, save_alignment=True)
        log.write('Mapping is done. df.shape: {},  not aligned genes: {}\n'.format(alignment_df.shape,
                                                                                   not_aligned))

        alignment_df.to_csv(GENOMES_DIR + folder + '/db_alignment.csv', index=False)

    log.close()
