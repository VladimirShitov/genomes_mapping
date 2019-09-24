EMAIL = 'vladimirs@intern.bii.a-star.edu.sg'

FTP_LINK = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/'

CD_HIT_PATH = 'cd-hit'  # path to a tool iteself. If it is not in $PATH variable, please, set it

ORGANISM = 'Escherichia_coli'  # TODO: make a parameter to generalize the pipeline
REFERENCE = 'GCF_000005845.2_ASM584v2'

DB_PATH = './blast_db/'
GENOMES_DIR = './data/'
TEMP_PATH = './temp'
PROTEOME_PATH = './results/proteome.faa'
GENES_INFO_PATH = './results/genes_info.tsv'
CLUSTERS_INFO = './results/cluster_info.json'
TOTAL_LIST_PATH = './results/total_list.csv'  # TODO: maybe it is better to create a separate folder for results
CLUSTERS_DF_PATH = './results/clusters.csv'
SUMMARY = './results/summary.txt'

LOG_PATH = './logs/pipeline.log'
DOWNLOADING_LOG = './logs/downloading.log'
DB_LOG = './logs/database.log'
ALIGNMENT_LOG = './logs/alignment.log'
TOTAL_LIST_LOG = './logs/total_list.log'
PROTEOME_LOG = './logs/proteome.log'

THRESHOLD = 90
