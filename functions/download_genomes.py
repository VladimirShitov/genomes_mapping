import os
import subprocess as sp

import numpy as np
import pandas as pd
from tqdm import tqdm  # Progress bar

import constants

FTP_LINK = constants.FTP_LINK
ORGANISM = constants.ORGANISM
LOG_PATH = constants.DOWNLOADING_LOG


def download_genomes():

    log = open(LOG_PATH, 'w', buffering=1)

    os.chdir('./data')

    os.system('wget {}'.format(FTP_LINK + ORGANISM + '/assembly_summary.txt'))
    assembly = pd.read_csv('assembly_summary.txt', sep='\t', header=1)
    assembly = assembly[assembly['assembly_level'].isin(['Complete Genome', 'Chromosome'])]

    is_downloaded = np.zeros(assembly.shape[0], dtype=bool)

    for i, link in enumerate(tqdm(assembly['ftp_path'])):

        if link == '-' or is_downloaded[i]:
            log.write('No link found\n')
            log.write('End of the iteration. {}/{} done\n\n'.format(i + 1, assembly.shape[0]))
            continue

        try:
            log.write('Working with {}\n'.format(link))

            name = link.split('/')[-1]
            os.makedirs(name, exist_ok=True)
            os.chdir(name)
            log.write('Working in dir {}\n'.format(os.getcwd()))

            for file in ('protein.faa.gz', 'feature_table.txt.gz', 'genomic.fna.gz', 'protein.gpff.gz'):
                file_link = link + '/' + name + '_' + file
                log.write('Trying to download {}\n'.format(file_link))

                wget = sp.run(['wget', file_link], cwd=os.getcwd(), capture_output=True)
                if wget.returncode:
                    log.write('Error\n')
                else:
                    log.write('Success\n')
                    is_downloaded[i] = 1

                log.write('Trying to extract archive\n')
                os.system('gunzip ' + name + '_' + file)

            os.chdir('..')  # Out of genome directory

        except Exception:
            log.write('Error occured\n')
        finally:
            log.write('End of the iteration. {}/{} done\n\n'.format(i + 1, assembly.shape[0]))

    os.chdir('..')  # Out of data directory
