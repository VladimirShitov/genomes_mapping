import os

from functions.read_cd_hit_output import read_cd_hit_output

from constants import CD_HIT_PATH, PROTEOME_PATH


def perform_cd_hit_clustering(threshold):
    os.system(CD_HIT_PATH
              + ' -i ' + PROTEOME_PATH
              + ' -o ./cd_hit/cd_hit_output'
              + ' -n 5 -d 0 -T 7 -M 4096'
              + ' -c {}'.format(threshold / 100))

    return read_cd_hit_output('./cd_hit/cd_hit_output.clstr')
