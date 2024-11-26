import argparse
from pdb_manipulation import get_plddt

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--startpdb", required=True)
    args = vars(parser.parse_args())
    get_plddt(args["startpdb"])
