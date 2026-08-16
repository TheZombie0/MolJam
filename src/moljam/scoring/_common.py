import os
import time
import warnings
from collections import defaultdict
from functools import partial
from multiprocessing import Pool, cpu_count

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.rdBase import DisableLog

from .chem import (
    calculate_mol_fingerprint,
    calculate_mol_scaffold,
    calculate_qed_value,
    check_mol_chirality,
    process_mol_for_consistency,
    process_single_smiles,
)

warnings.filterwarnings('ignore')
DisableLog('rdApp.warning')
