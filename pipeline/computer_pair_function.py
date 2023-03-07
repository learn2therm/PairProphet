"""
This script parses hmmer output and tests if the paired protein are from a similar family

The packages you need to run this script are the following:
- biopython
- HMMER (http://hmmer.org/documentation.html)
- pandas

You also need to have:
- pfam db locally
- protein db
"""

# system dependecies
from pathlib import Path
import subprocess



# library dependencies
from collections import defaultdict
import pandas as pd



## biopython
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SearchIO


# local dependencies/utils

## Paths
PFAM_PATH = Path("/Users/humoodalanzi/pfam/Pfam-A.hmm")
ID_DB_PATH = Path("/Users/humoodalanzi/pfam/proteins_id.zip")