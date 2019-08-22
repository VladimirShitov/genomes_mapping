import os

import matplotlib.pyplot as plt
from Bio import Entrez
from tqdm import tqdm  # progress bar

import constants
from functions.create_db import create_database
from functions.align_to_db import align_to_database

Entrez.email = constants.EMAIL

GENOMES_DIR = constants.GENOMES_DIR
TEMP_PATH = constants.TEMP_PATH
LOG_PATH = constants.LOG_PATH
DB_PATH = constants.DB_PATH
REFERENCE = constants.REFERENCE

align_to_database()
