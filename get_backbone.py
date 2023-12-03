import argparse

from pdb_manipulation import select_backbone

parser = argparse.ArgumentParser()
parser.add_argument("--startpdb", required=True)
parser.add_argument("--targetpdb", required=True)
args = vars(parser.parse_args())
select_backbone(args["startpdb"], args["targetpdb"])
