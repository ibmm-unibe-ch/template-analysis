import argparse
import pickle

from protein_learning.assessment.sidechain import assess_sidechains, summarize

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--predictedpdb", required=True)
    parser.add_argument("--targetpdb", required=True)
    parser.add_argument("--outputpath", required=True)
    args = vars(parser.parse_args())
    res_level_stats = assess_sidechains(
        args["predictedpdb"], args["targetpdb"], steric_tol_fracs=[1, 0.9, 0.8]
    )
    with open(args["outputpath"], "wb") as outfile:
        pickle.dump(res_level_stats, outfile)
