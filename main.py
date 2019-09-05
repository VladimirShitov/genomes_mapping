from functions.set_enviroment import set_enviroment
from functions.download_genomes import download_genomes
from functions.create_db import create_database
from functions.align_to_db import align_to_database
import constants

set_enviroment()
log = open(constants.LOG_PATH, 'w', buffering=1)

try:
    log.write('Downloading genomes. More details in {}\n'.format(constants.DOWNLOADING_LOG))
    download_genomes()
except Exception as e:
    log.write('Error happened:\n{}\n'.format(e))

try:
    log.write('Creating database. More details in {}\n'.format(constants.DB_LOG))
    create_database()
except Exception as e:
    log.write('Error happened:\n{}\n'.format(e))

try:
    log.write('Aligning genes to database. More details in {}\n'.format(constants.ALIGNMENT_LOG))
    align_to_database()
except Exception as e:
    log.write('Error happened:\n{}\n'.format(e))

# build_graph()

log.close()
