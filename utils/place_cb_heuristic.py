import argparse
from pdb_manipulation import select_bb_c_beta

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--startpdb", required=True)
    parser.add_argument("--targetpdb", required=True)
    args = vars(parser.parse_args())
    select_bb_c_beta(args["startpdb"], args["targetpdb"], only_backbone=True)
