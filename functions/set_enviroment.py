import os


def set_enviroment():  # TODO: document

    os.makedirs('blast_db', exist_ok=True)
    os.makedirs('data', exist_ok=True)
    os.makedirs('plots', exist_ok=True)
    os.makedirs('logs', exist_ok=True)
