"""
_summary_
"""

# system dependecies
import os
from pathlib import Path
import sys
from typing import Union

# job stuff
import logging
from joblib import delayed, Parallel


# library dependencies
import pandas as pd
from tqdm import tqdm

# biopython
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# pyhmmer
import pyhmmer