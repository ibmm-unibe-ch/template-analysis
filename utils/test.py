from pathlib import Path
import argparse
import re
from prody import parsePDB
import pandas as pd
import pickle
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from read_files import read_ost_scores

if __name__ == "__main__":
    read_ost_scores() 