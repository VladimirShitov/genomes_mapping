import os

import matplotlib.pyplot as plt
from Bio import Entrez
from tqdm import tqdm

from Blast_Database import BlastDatabase
from Genome import Genome
from functions.map_genome import map_genome

Entrez.email = 'vladimirs@intern.bii.a-star.edu.sg'

GENOMES_DIR = '/home/vladimir/Documents/E_coli/seqs_mapping/'
TEMP_PATH = '/home/vladimir/Documents/E_coli/_temp_2/temp'
LOG_PATH = '/home/vladimir/Documents/E_coli/all_genomes_mapping_2.log'
DB_PATH = '/home/vladimir/Documents/E_coli/blastdb_2/'
REFERENCE = 'GCF_000005845.2_ASM584v2'

# Get all directories in GENOMES_DIR
folders = list(filter(lambda x: os.path.isdir(GENOMES_DIR + x), os.listdir(GENOMES_DIR)))

# Make REFERENCE the first element in the list
folders.sort(key=lambda x: x == REFERENCE, reverse=True)

log = open(LOG_PATH, 'w', buffering=1)
db_sizes = []

blast_db = BlastDatabase()

for i, folder in enumerate(tqdm(folders)):

    try:
        prefix = folder + '_'
        has_feature_table = (prefix + 'feature_table.txt') in os.listdir(GENOMES_DIR + folder)
        has_protein_faa = (prefix + 'protein.faa') in os.listdir(GENOMES_DIR + folder)

        # Delete folder with genome if it doen't have feature_table or protein_faa
        if not (has_feature_table and has_protein_faa):
            log.write('Deleting {}\n'.format(GENOMES_DIR + folder))
            os.system('rm -rf {}'.format(GENOMES_DIR + folder))

        else:
            log.write('Working with {}\n'.format(GENOMES_DIR + folder))
            genome = Genome(folder_path=GENOMES_DIR + folder + '/', files_prefix=prefix)
            genome.read_protein_faa()
            genome.read_feature_table()
            genome.set_gene_positions()

            log.write('Created genome\n')

            if genome.raised_errors:
                log.write('ERROR\n{}\n'.format(genome.log))
                log.close()
                raise Exception(genome.log)

            if folder == REFERENCE:  # The very first run

                blast_db.feed_genes(genome.genes)
                blast_db.create_db(append=False, path=DB_PATH)

                reference = genome

            log.write('Mapping genomes...\n')
            alignment_list, not_aligned = map_genome(genome, reference)
            log.write('Mapping is done. Shape of dataframe: {}, not aligned genes: {}\n'.format(alignment_list.shape,
                                                                                                len(not_aligned)))

            # alignment_list.to_csv(GENOMES_DIR + folder + '/alignment_list.csv', index=False)

            blast_db.feed_genes([genome[gene] for gene in not_aligned])
            log.write('Updated database. Size now: {}\n'.format(len(blast_db.genes)))

            # We have to add new genes to reference, so we could always know their positions.
            reference.genes.extend([genome[gene] for gene in not_aligned])

            if not_aligned:
                blast_db.create_db(append=True, path=DB_PATH)

            db_sizes.append(len(blast_db.genes))

            log.write('End of the iteration\n\n')

            if i % 10 == 5 or i == len(folders) - 1:
                plt.plot(db_sizes)
                plt.title('Size of database')
                plt.show()

    except Exception as e:
        log.write(str(e))
        log.write('\n')
        continue
log.close()
