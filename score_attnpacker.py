import argparse
import json

from protein_learning.assessment.sidechain import assess_sidechains, summarize

parser = argparse.ArgumentParser()
parser.add_argument("--predictedpdb", required=True)
parser.add_argument("--targetpdb", required=True)
parser.add_argument("--outputpath", required=True)
args = vars(parser.parse_args())
res_level_stats = assess_sidechains(
    args["predictedpdb"], args["targetpdb"], steric_tol_fracs=[1, 0.9, 0.8]
)
target_level_stats = summarize(res_level_stats)
dikt = {
    "rmsd": target_level_stats["rmsd"].item(),
    "mae_1": target_level_stats["mean_mae"][0].item(),
    "mae_2": target_level_stats["mean_mae"][1].item(),
    "mae_3": target_level_stats["mean_mae"][2].item(),
    "mae_4": target_level_stats["mean_mae"][3].item(),
    "acc": target_level_stats["mae_sr"].item(),
}
with open(args["outputpath"], "w") as outfile:
    json.dump(dikt, outfile)
