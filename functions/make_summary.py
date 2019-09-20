import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from constants import TOTAL_LIST_PATH


def make_summary():
    df = pd.read_csv(TOTAL_LIST_PATH)
