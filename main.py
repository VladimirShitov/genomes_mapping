from functions.set_enviroment import set_enviroment
from functions.download_genomes import download_genomes
from functions.create_db import create_database
from functions.align_to_db import align_to_database

set_enviroment()
download_genomes()
create_database()
align_to_database()
# build_graph()