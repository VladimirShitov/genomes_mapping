from functions.set_enviroment import set_enviroment
from functions.download_genomes import download_genomes
from functions.create_proteome import create_proteome
from functions.perform_cd_hit_clustering import perform_cd_hit_clustering

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
    log.write('Trying to perform CD-hit clustering\n')
    perform_cd_hit_clustering(THRESHOLD)
except Exception as e:
    log.write('Error happened:\n{}\n'.format(e))


log.close()
