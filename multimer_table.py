import json
import os
from pathlib import Path
from statistics import mean

import pandas as pd


# fix 7cg5
def read_scores(scores_path: Path) -> dict:
    with open(scores_path, "r") as file:
        data = json.load(file)
    output = {
        "rmsd": data.get("rmsd", None),
        "tm_score": data.get("tm_score", None),
        "lddt": data.get("lddt", None),
    }
    if "dockq_scores" in data:
        output["dockq_scores"] = data["dockq_scores"][0] if data["dockq_scores"] else None
    return output


def get_omega_result(parent_path: Path, identifier: str) -> dict:
    output = read_scores(parent_path / f"{identifier}_omega" / "scores.json")
    basic_info = {"identifier": identifier, "tool": "OmegaFold"}
    return output | basic_info


def get_esm_result(parent_path: Path, identifier: str) -> dict:
    output = read_scores(parent_path / f"{identifier}_esm" / "scores.json")
    basic_info = {"identifier": identifier, "tool": "ESM"}
    return output | basic_info


def read_colabfold_scores(colabfold_path: Path) -> dict:
    with open(colabfold_path, "r") as file:
        data = json.load(file)
    output = {
        "plddt": mean(data["plddt"]),
        "max_pae": data["max_pae"],
        "ptm": data["ptm"],
    }
    return output


def get_colabfold_result(
    parent_path: Path, identifier: str, msa_mode: str, template_id: str
) -> dict:
    folder_path = parent_path / f"{identifier}_{msa_mode}_{template_id}"
    file_names = [
        filename
        for filename in os.listdir(folder_path)
        if filename.startswith(f"{identifier}_scores_rank_001_")
    ]
    if len(file_names) == 0:
        return None
    file_name = file_names[0]
    colabfold_results = read_colabfold_scores(folder_path / file_name)
    basic_info = {
        "identifier": identifier,
        "msa_mode": msa_mode,
        "template_id": template_id,
        "tool": "ColabFold",
    }
    score = read_scores(folder_path / "scores.json")
    return basic_info | colabfold_results | score


def get_results_identifier(parent_path: Path, identifier: str) -> dict:
    output = [get_omega_result(parent_path, identifier)]
    output.append(get_esm_result(parent_path, identifier))
    for (msa_mode, template_id) in [
        (msa_mode, template_id)
        for msa_mode in ["single", "full"]
        for template_id in ["none", "omega", "truth", "comb", "esm"]
    ]:
        result = get_colabfold_result(parent_path, identifier, msa_mode, template_id)
        if result:
            output.append(result)
    result = get_colabfold_result(parent_path, identifier, "single", "af")
    if result:
        output.append(result)
    return output


def get_multimer_esm_result(parent_path: Path, identifier: str):
    output = read_scores(parent_path / f"{identifier}_fullesm_full" / "scores.json")
    basic_info = {"identifier": identifier, "tool": "ESM", "parts": None}
    return output | basic_info


def get_multimer_colabfold_result(parent_path: Path, identifier: str, parts_id: str):
    folder_path = parent_path / f"{identifier}_{parts_id}_full"
    file_names = [
        filename
        for filename in os.listdir(folder_path)
        if filename.startswith(f"{identifier}_scores_rank_001_")
    ]
    if len(file_names) == 0:
        return None
    file_name = file_names[0]
    colabfold_results = read_colabfold_scores(folder_path / file_name)
    basic_info = {"identifier": identifier, "parts": parts_id, "tool": "ColabFold"}
    score = read_scores(folder_path / "scores.json")
    return basic_info | colabfold_results | score


def get_results_identifier_multimer(parent_path, identifier):
    output = [get_multimer_esm_result(parent_path, identifier)]
    for parts_id in ["af", "comb", "fullaf", "omega", "truth", "esm", "msaaf"]:
        result = get_multimer_colabfold_result(parent_path, identifier, parts_id)
        if result:
            output.append(result)
    return output


def get_results_folder(folder_path: Path):
    output = []
    for file in os.listdir(folder_path):
        if not str(file).endswith(".fasta") or str(file).endswith("_single.fasta"):
            continue
        identifier = str(file)[:4]
        output = output + get_results_identifier_multimer(folder_path, identifier)
    dataframe = pd.DataFrame(output)
    return dataframe
    #dataframe.to_csv(folder_path / "table_multimer.csv")


if __name__ == "__main__":
    get_results_folder(Path("casp_oligomers"))
