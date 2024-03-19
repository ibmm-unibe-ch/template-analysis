import json
import os
from pathlib import Path
from statistics import mean
from typing import Optional

import pandas as pd


def read_ost_scores(path: Path, local_lddt: bool = False) -> dict:
    with open(path, "r") as json_file:
        read_json = json.load(json_file)
    if read_json["status"] == "FAILURE":
        return None
    output = {
        "lddt": read_json.get("lddt", None),
        "ost_rmsd": read_json.get("rmsd", None),
        "tm_score": read_json.get("tm_score", None),
    }
    if "dockq_scores" in read_json and len(read_json["dockq_scores"]):
        output["dockq_scores"] = read_json["dockq_scores"][0]
        output["irmsd"] = read_json["irmsd"][0]
    if local_lddt:
        local_lddt = read_json.get("local_lddt", None).values()
        output["local_lddt"] = pd.Series(local_lddt, name="local_lddt")
    return output


def read_ppi_scores(path: Path) -> dict:
    if not path.exists():
        return None
    with open(path, "r") as json_file:
        read_json = json.load(json_file)
    if read_json["status"] == "FAILURE":
        return None
    output = {"ppi-rmsd": read_json.get("rmsd", None)}
    if "dockq_scores" in read_json and len(read_json["dockq_scores"]):
        output["dockq_scores"] = read_json["dockq_scores"][0]
        output["irmsd"] = read_json["irmsd"][0]
    return output


def get_omega_result(parent_path: Path, identifier: str) -> dict:
    output = read_ost_scores(parent_path / f"{identifier}_omega" / "scores.json")
    basic_info = {"identifier": identifier, "tool": "OmegaFold"}
    if output:
        return output | basic_info
    return None


def get_esm_result(parent_path: Path, identifier: str) -> dict:
    output = read_ost_scores(parent_path / f"{identifier}_esm" / "scores.json")
    basic_info = {"identifier": identifier, "tool": "ESM"}
    if output:
        return output | basic_info
    return None


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
    score = read_ost_scores(folder_path / "scores.json")
    if score and colabfold_results:
        return basic_info | colabfold_results | score
    return None


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
    output = read_ost_scores(parent_path / f"{identifier}_esm" / "scores.json")
    basic_info = {"identifier": identifier, "tool": "ESM", "parts": None}
    if output:
        return output | basic_info
    return None


def get_multimer_colabfold_result(parent_path: Path, identifier: str, parts_id: str):
    folder_path = parent_path / f"{identifier}_full_{parts_id}"
    if not folder_path.exists():
        return None
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
    score = read_ost_scores(folder_path / "scores.json")
    if colabfold_results and score:
        return basic_info | colabfold_results | score
    return None


def get_ppi_colabfold_result(parent_path: Path, identifier: str, parts_id: str):
    folder_path = Path(str(parent_path) + "_" + parts_id + "_full") / identifier
    print(folder_path)
    if not folder_path.exists():
        return None
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
    score = read_ost_scores(folder_path / "scores.json")
    ppi_scores = read_ppi_scores(folder_path / "scores_backbone.json")
    if ppi_scores:
        score = score | ppi_scores
    if colabfold_results and score:
        return basic_info | colabfold_results | score
    return None


def get_results_identifier_multimer(parent_path, identifier):
    output = [get_multimer_esm_result(parent_path, identifier)]
    for parts_id in ["af", "comb", "fullaf", "omega", "truth", "esm", "msaaf"]:
        result = get_multimer_colabfold_result(parent_path, identifier, parts_id)
        if result:
            output.append(result)
    return output


def get_results_identifier_ppi(parent_path, identifier):
    # output = [get_multimer_esm_result(parent_path, identifier)]
    output = []
    for parts_id in ["af", "comb", "fullaf", "omega", "truth", "esm", "msaaf"]:
        result = get_ppi_colabfold_result(parent_path, identifier, parts_id)
        if result:
            output.append(result)
    return output


def get_results_identifer_monomer(parent_path, identifier):
    output = [get_multimer_esm_result(parent_path, identifier)]
    for parts_id in ["af", "comb", "fullaf", "omega", "truth", "esm", "msaaf"]:
        result = get_multimer_colabfold_result(parent_path, identifier, parts_id)
        if result:
            output.append(result)
    return output


def get_results_folder(folder_path: Path, save_path: Optional[Path] = None):
    output = []
    for file in os.listdir(folder_path):
        if not str(file).endswith(".fasta") or str(file).endswith("_single.fasta"):
            continue
        identifier = str(file)[:4]
        print(identifier)
        output = output + get_results_identifier_multimer(folder_path, identifier)
    print(output)
    dataframe = pd.DataFrame([output_item for output_item in output if output_item])
    if save_path:
        dataframe.to_csv(save_path)
    else:
        return dataframe


def get_results_folder_ppi(folder_path: Path, save_path: Optional[Path] = None):
    output = []
    for file in os.listdir(folder_path):
        if not str(file).endswith(".fasta") or str(file).endswith("_single.fasta"):
            continue
        identifier = str(file)[:4]
        print(identifier)
        output = output + get_results_identifier_ppi(folder_path, identifier)
    print(output)
    dataframe = pd.DataFrame([output_item for output_item in output if output_item])
    if save_path:
        dataframe.to_csv(save_path)
    else:
        return dataframe


def get_refinement_results_folder_monomer(
    folder_path: Path, save_path: Optional[Path] = None
):
    output = []
    for file in os.listdir(folder_path):
        if not str(file).endswith(".fasta") or str(file).endswith("_single.fasta"):
            continue
        identifier = str(file)[:4]
        print(identifier)
        output = output + get_results_identifier(folder_path, identifier)
    dataframe = pd.DataFrame([output_item for output_item in output if output_item])
    if save_path:
        dataframe.to_csv(save_path)
    else:
        return dataframe
