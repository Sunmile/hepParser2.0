import numpy as np
import pandas as pd
from sklearn.metrics import r2_score
from sklearn.linear_model import Lasso
from src.theory_sp import *


def generate_feature(the_SP_path, the_HP, max_ion, is_HNa):
    the_spectra = get_the_sp(the_SP_path, the_HP[3], the_HP[2], 1, max_ion, is_HNa=is_HNa)







