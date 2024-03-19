import argparse

from pdb_manipulation import flatten_pdb

parser = argparse.ArgumentParser()
parser.add_argument("--startpdb", required=True)
parser.add_argument("--targetpdb", required=True)
parser.add_argument("--pca", required=True)
args = vars(parser.parse_args())
flatten_pdb(args["startpdb"], args["targetpdb"], int(args["pca"]))
