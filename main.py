from functions.set_enviroment import set_enviroment
from functions.download_genomes import download_genomes
from functions.create_proteome import create_proteome
from functions.create_total_list_of_genes import create_total_list_of_genes
from functions.make_summary import make_summary

from constants import LOG_PATH, DOWNLOADING_LOG, TOTAL_LIST_LOG, PROTEOME_LOG

set_enviroment()
log = open(LOG_PATH, 'w', buffering=1)

try:
    log.write('Downloading genomes. More details in {}\n'.format(DOWNLOADING_LOG))
    download_genomes()

    log.write('Trying to create proteome. More information in {}\n'.format(PROTEOME_LOG))
    create_proteome()

    log.write('Trying to create total list of genes. More information in {}\n'.format(TOTAL_LIST_LOG))
    create_total_list_of_genes()

    log.write('Making summary\n')
    make_summary()
    log.write('All done!')
finally:
    log.close()
