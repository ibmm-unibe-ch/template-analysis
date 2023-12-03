import argparse
import os
from pathlib import Path

import dill
import pandas as pd
from biopandas.pdb import PandasPdb

three_to_one = {
    "CYS": "C",
    "ASP": "D",
    "SER": "S",
    "GLN": "Q",
    "LYS": "K",
    "ILE": "I",
    "PRO": "P",
    "THR": "T",
    "PHE": "F",
    "ASN": "N",
    "GLY": "G",
    "HIS": "H",
    "LEU": "L",
    "ARG": "R",
    "TRP": "W",
    "ALA": "A",
    "VAL": "V",
    "GLU": "E",
    "TYR": "Y",
    "MET": "M",
}


def get_identifier(identifier):
    # magic numbers are magically working for the 100 test set
    additive = int(identifier[11]) + 4 if identifier[11].isnumeric() else 0
    upd = str(
        (
            int(identifier[0])
            + int(identifier[8]) * 9
            + int(identifier[10]) * 7
            + additive
        )
        % 10
    )
    return upd + str(identifier[1:4])


def write_file(file_path, input):
    os.makedirs(Path(file_path).parent, exist_ok=True)
    with open(file_path, "w") as f:
        f.write("\n".join(input))


def parse_dill(fp, parent_path):
    with open(fp, "rb") as f:
        data = dill.load(f)
    pdb_A, code_A = parse_df(data[1], "A")
    pdb_B, code_B = parse_df(data[2], "B")
    path = Path(fp)
    parent_path = Path(parent_path)
    identifier = get_identifier(path.name)
    other_identifer = str((int(identifier[0]) + 1) % 10) + identifier[1:]
    pdb_combined = pd.concat([pdb_A, pdb_B])
    write_pdb(pdb_A, parent_path / f"{identifier}_truth" / f"{identifier}.pdb")
    write_pdb(pdb_B, parent_path / f"{identifier}_truth" / f"{other_identifer}.pdb")
    write_pdb(pdb_combined, parent_path / f"{identifier}" / f"{identifier}.pdb")
    fasta_single = f"{code_A}:{code_B}"
    write_file(
        parent_path / f"{identifier}_single.fasta",
        [f">{identifier}", code_A, f">{other_identifer}", code_B],
    )
    write_file(parent_path / f"{identifier}.fasta", [f">{identifier}", fasta_single])


def parse_df(df, chain_id):
    """
    Parse PDB DataFrame
    """
    # extract dataframe values
    df["residue"] = df["residue"].str.extract("(\d+)").astype(int)
    df = df.rename(
        columns={
            "aid": "atom_number",
            "resname": "residue_name",
            "residue": "residue_number",
            "x": "x_coord",
            "y": "y_coord",
            "z": "z_coord",
            "element": "element_symbol",
        },
        errors="raise",
    )
    df["record_name"] = "ATOM"
    df["occupancy"] = 1
    df["b_factor"] = 0
    df["charge"] = None
    df["line_idx"] = df["atom_number"]
    df["blank_1"] = ""
    df["alt_loc"] = ""
    df["blank_2"] = ""
    df["insertion"] = ""
    df["blank_3"] = ""
    df["blank_4"] = ""
    df["segment_id"] = ""
    df["chain_id"] = chain_id
    df = df[
        [
            "record_name",
            "atom_number",
            "blank_1",
            "atom_name",
            "alt_loc",
            "residue_name",
            "blank_2",
            "chain_id",
            "residue_number",
            "insertion",
            "blank_3",
            "x_coord",
            "y_coord",
            "z_coord",
            "occupancy",
            "b_factor",
            "blank_4",
            "segment_id",
            "element_symbol",
            "charge",
            "line_idx",
        ]
    ]
    residues = df.groupby(["residue_number"]).agg("first").sort_index()["residue_name"]
    code = "".join([three_to_one[g] for g in residues])
    return df, code


def write_pdb(df, output_path):
    os.makedirs(Path(output_path).parent, exist_ok=True)
    pdb1 = PandasPdb()
    pdb1._df["ATOM"] = df
    pdb1.to_pdb(str(output_path))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input")
    parser.add_argument("--output")
    args = vars(parser.parse_args())
    input_path = args["input"]
    output_path = args["output"]
    if (input_path is None) or (output_path is None):
        dips100 = pd.read_csv("data_file_100_test.csv")["path"]
        for dip in dips100:
            print(dip)
            parse_dill(f"DIPS/{dip}", "dips100")
    else:
        parse_dill(input_path, output_path)
