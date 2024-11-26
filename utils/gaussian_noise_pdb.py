import argparse
from pdb_manipulation import noise_gaussian_pdb

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--startpdb", required=True)
    parser.add_argument("--targetpdb", required=True)
    parser.add_argument("--std", required=True)
    args = vars(parser.parse_args())
    noise_gaussian_pdb(args["startpdb"], args["targetpdb"], float(args["std"]))
