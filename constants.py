FTP_LINK = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/'
ORGANISM = 'Escherichia_coli'  # TODO: make a parameter to generalize the pipeline
REFERENCE = 'GCF_000005845.2_ASM584v2'

DB_PATH = './blast_db/'
GENOMES_DIR = './data/'
TEMP_PATH = './temp'
CD_HIT_PATH = 'cd-hit'
PROTEOME_PATH = './data/proteome.faa'

LOG_PATH = './logs/pipeline.log'
DOWNLOADING_LOG = './logs/downloading.log'
DB_LOG = './logs/database.log'
ALIGNMENT_LOG = './logs/alignment.log'

THRESHOLD = 50


EMAIL = 'vladimirs@intern.bii.a-star.edu.sg'
