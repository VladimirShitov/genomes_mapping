import os

from tqdm import tqdm

from Genome import Genome

from constants import PROTEOME_PATH, GENOMES_DIR, GENES_INFO_PATH, PROTEOME_LOG


def create_proteome():
    """Create one file with protein sequences from all the genomes in data folder"""

    log = open(PROTEOME_LOG, 'w', buffering=1)

    # Get all directories in GENOMES_DIR
    folders = list(filter(lambda x: os.path.isdir(GENOMES_DIR + x), os.listdir(GENOMES_DIR)))

    with open(PROTEOME_PATH, 'a') as proteome:
        with open(GENES_INFO_PATH, 'a') as genes_info:
            genes_info.write('gene\tinfo\n')

            for folder in tqdm(folders):
                log.write('Working with folder {}\n'.format(folder))

                prefix = folder + '_'
                has_feature_table = (prefix + 'feature_table.txt') in os.listdir(GENOMES_DIR + folder)
                has_protein_faa = (prefix + 'protein.faa') in os.listdir(GENOMES_DIR + folder)

                # Delete folder with genome if it doen't have feature_table or protein_faa
                if not (has_feature_table and has_protein_faa):
                    log.write('Deleting {}\n'.format(GENOMES_DIR + folder))
                    os.system('rm -rf {}'.format(GENOMES_DIR + folder))
                    continue

                genome = Genome(GENOMES_DIR+folder+'/', prefix)
                genome.read_protein_faa()
                genome.read_feature_table()
                genome.set_gene_positions()

                if genome.raised_errors:
                    log.write('Error:\n{}\n'.format(genome.log))
                    log.write('Genome: {}\n'.format(str(genome)))

                for gene in genome.genes:
                    proteome.write(gene.get_fasta(join_genome_to_name=True))
                    genes_info.write('{gene}[genome:{genome}]\t{info}\n'.format(gene=gene.id,
                                                                                genome=gene.genome,
                                                                                info=gene.info))
    log.close()
