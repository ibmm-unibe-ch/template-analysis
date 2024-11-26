from pathlib import Path
import argparse
import re
from prody import parsePDB
import pandas as pd
import pickle
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from read_files import read_ost_scores

one_to_three = {
    "A": "Ala",
    "B": "Asx",
    "C": "Cys",
    "D": "Asp",
    "E": "Glu",
    "F": "Phe",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "K": "Lys",
    "L": "Leu",
    "M": "Met",
    "N": "Asn",
    "P": "Pro",
    "Q": "Gln",
    "R": "Arg",
    "S": "Ser",
    "T": "Thr",
    "U": "Sec",
    "V": "Val",
    "W": "Trp",
    "X": "Xaa",
    "Y": "Tyr",
    "Z": "Glx",
}


def get_plddt(pdb_path: str):
    cas = parsePDB(pdb_path).select("name CA")
    aui = []
    for ca in cas:
        curr = {
            "af_pLDDT": float(ca.getBeta()),
            "af_resnum": ca.getResnum(),
            "af_resname": ca.getResname().upper(),
        }
        aui.append(curr)
    return pd.DataFrame(aui)


def get_sasa(pdb_path: str):
    structure = PDBParser().get_structure("pdb", pdb_path)
    model = structure[0]
    dssp = list(DSSP(model, pdb_path, dssp="mkdssp"))
    res_name = [one_to_three[lis[1]].upper() for lis in dssp]
    sasa_score = [lis[3] for lis in dssp]
    return pd.DataFrame(
        zip(res_name, sasa_score), columns=["sasa_resname", "sasa_sasa"]
    )


def get_sec_struc(pdb_path: str):
    structure = PDBParser().get_structure("pdb", pdb_path)
    model = structure[0]
    dssp = list(DSSP(model, pdb_path, dssp="mkdssp"))
    res_name = [one_to_three[lis[1]].upper() for lis in dssp]
    sec_struc = [lis[2] for lis in dssp]
    return pd.DataFrame(
        zip(res_name, sec_struc), columns=["struc_resname", "struc_struc"]
    )


def get_scores(score_path: str):
    with open(score_path, "rb") as outfile:
        score_dict = pickle.load(outfile)
    rmsd_score = score_dict["rmsd"]["per-prot"]
    rmsd_score_central = score_dict["rmsd"]["per-prot-central"]
    rmsd_score_surface = score_dict["rmsd"]["per-prot-surface"]
    rmsd_per_res = score_dict["rmsd"]["per-res"]
    mae_score = score_dict["dihedral"]["mae"]
    mae_mask = score_dict["dihedral"]["mask"]
    seq = score_dict["sequence"]
    aui = []
    for it, char in enumerate(seq):
        mae_1, mae_2, mae_3, mae_4 = [
            abs(mae) if mae_mask[it][it2] else None
            for it2, mae in enumerate(mae_score[it])
        ]
        curr = {
            "att_char": one_to_three[char].upper(),
            "att_rmsd": rmsd_per_res[it],
            "att_mae1": mae_1,
            "att_mae2": mae_2,
            "att_mae3": mae_3,
            "att_mae4": mae_4,
        }
        aui.append(curr)
    return pd.DataFrame(aui), rmsd_score, rmsd_score_central, rmsd_score_surface


def make_sasa_df():
    for method in [
        "TEMPLATE_AF",
        "AF_FULL",
        "ATTNPACK",
        "ATTNPACK_AF",
        "BB_CB",
        "CB_HEURISTIC",
        "FASPR_AF",
        "PSF_AF",
    ]:
        print(f"Starting with {method}")
        out = []
        for dataset in ["CASP13", "CASP14"]:
            print(f"Started with {dataset}")
            method_path = Path(f"{dataset}_{method}")
            for target_path in method_path.glob("*"):
                target = target_path.name
                print(f"Starts with {target}")
                best_pdb_list = list(target_path.glob("*rank_001*6217.pdb"))
                best_pdb = (
                    best_pdb_list[0]
                    if len(best_pdb_list) > 0
                    else list(target_path.glob("*.pdb"))[0]
                )
                curr_sasa_df = get_sasa(f"{dataset}/{target}.pdb")
                curr_plddt = get_plddt(str(best_pdb))["af_pLDDT"]
                curr_lddt = read_ost_scores(target_path / "ost_scores.json", True)[
                    "local_lddt"
                ]
                curr_out = pd.concat([curr_sasa_df, curr_plddt, curr_lddt], axis=1)
                out.append(curr_out)
            pd.concat(out, axis=0, ignore_index=True).to_csv(
                f"correlations/corr_{dataset}_{method}.csv"
            )


def compile_scores(parent_path: Path, correct_pdb: Path = None):
    best_pdb_list = list(parent_path.glob("*rank_001*6217.pdb"))
    best_pdb = (
        best_pdb_list[0]
        if len(list(parent_path.glob("*rank_001*6217.pdb"))) > 0
        else list(parent_path.glob("*.pdb"))[0]
    )
    correct_pdb = best_pdb if correct_pdb is None else correct_pdb
    pLDDT_df = get_plddt(str(best_pdb))
    scores_df, rmsd, rmsd_central, rmsd_surface = get_scores(
        parent_path / "attn_scores.json"
    )
    sasa_df = get_sasa(correct_pdb)
    ost = read_ost_scores(parent_path / "ost_scores.json", True)
    output = {
        "rmsd": rmsd,
        "rmsd_central": rmsd_central,
        "rmsd_surface": rmsd_surface,
        "lddt": ost["lddt"],
        "ost_rmsd": ost["ost_rmsd"],
        "tm_score": ost["tm_score"],
    }
    concat_df = pd.concat([pLDDT_df, scores_df, sasa_df, ost["local_lddt"]], axis=1)
    other_columns = [
        "af_pLDDT",
        "att_rmsd",
        "att_mae1",
        "att_mae2",
        "att_mae3",
        "att_mae4",
        "local_lddt",
    ]
    for other in other_columns:
        output[other + "_mean"] = concat_df[other].mean()
        output[other + "_x_sasa"] = concat_df[other].corr(
            concat_df["sasa_sasa"], "pearson"
        )
    for other in other_columns[1:]:
        output[other + "_x_pLDDT"] = concat_df[other].corr(
            concat_df["af_pLDDT"], "pearson"
        )
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--parentpath", required=True)
    parser.add_argument("--correctpath", required=True)
    args = vars(parser.parse_args())
    parent_path = Path(args["parentpath"])
    template = parent_path.stem
    dataset_method = parent_path.parent.stem
    dataset, method = re.split(r"_", dataset_method, 1)
    score_df = {"method": method, "dataset": dataset, "template": template}
    score_df = score_df | compile_scores(parent_path, args["correctpath"])
    with open(parent_path / "scores.pkl", "wb") as output_path:
        pickle.dump(score_df, output_path)
