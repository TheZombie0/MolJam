import os
import warnings
from functools import partial
from math import pi
from multiprocessing import Pool, cpu_count

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, Lipinski
from rdkit.rdBase import DisableLog

warnings.filterwarnings('ignore')
DisableLog('rdApp.warning')
DisableLog('rdApp.error')

