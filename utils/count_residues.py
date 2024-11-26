import argparse
from prody import parsePDB

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputpdb", required=True)
    parser.add_argument("--rfformat")
    args = vars(parser.parse_args())
    atoms = parsePDB(args["inputpdb"])
    num_residues = len(set(atoms.getResnums()))
    output = (
        "" + "contigmap.contigs=[" + f"{num_residues}-{num_residues}]" + ""
        if args["rfformat"]
        else num_residues
    )
    print(output)
