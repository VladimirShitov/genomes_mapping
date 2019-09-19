from functions.set_enviroment import set_enviroment
from functions.download_genomes import download_genomes
from functions.create_proteome import create_proteome
from functions.create_total_list_of_genes import create_total_list_of_genes

from constants import LOG_PATH, DOWNLOADING_LOG, THRESHOLD

set_enviroment()
log = open(LOG_PATH, 'w', buffering=1)

try:
    log.write('Downloading genomes. More details in {}\n'.format(DOWNLOADING_LOG))
    download_genomes()
except Exception as e:
    log.write('Error happened:\n{}\n'.format(e))

try:
    log.write('Trying to create proteome\n')
    create_proteome()
except Exception as e:
    log.write('Error happened:\n{}\n'.format(e))

try:
    log.write('Trying to create total list of genes\n')
    create_total_list_of_genes()
except Exception as e:
    log.write('Error happened:\n{}\n'.format(e))


log.close()
