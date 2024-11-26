import os
from pathlib import Path
import numpy as np
from prody import calcTransformation, parsePDB, writePDB
from scipy.spatial.transform import Rotation as R
from sklearn.decomposition import PCA
from sklearn.preprocessing import RobustScaler

default_angles = {"1": [60, 14, 42], "2": [83, -100, 17], "3": [66, -14, 132]}


def get_plddt(input_file):
    atoms = parsePDB(input_file)
    plddt = atoms.getBetas()
    print(plddt)
    return atoms


def select_residues(atoms, start=None, end=None):
    start = 0 if start is None else start
    end = max(atoms.getResnums()) if end is None else end
    before = (
        atoms.select(f"resnum < {start}").toAtomGroup()
        if not atoms.select(f"resnum < {start}") is None
        else None
    )
    middle = (
        atoms.select(f"resnum {start}:{end}").toAtomGroup()
        if not atoms.select(f"resnum {start}:{end}") is None
        else None
    )
    after = (
        atoms.select(f"resnum >= {end}").toAtomGroup()
        if not atoms.select(f"resnum >= {end}") is None
        else None
    )
    return before, middle, after


def put_together(parts):
    atoms = []
    for part in parts:
        if not part is None:
            coords = part if isinstance(part, np.ndarray) else part.getCoords()
            atoms.append(coords)
    return np.vstack(atoms)


def flatten_atoms(atoms, n_components=2):
    scaler = RobustScaler()
    coords = atoms.getCoords()
    standardised_coords = scaler.fit_transform(coords)
    pca = PCA(n_components=n_components)
    flattened_atoms = pca.fit_transform(standardised_coords)
    updated_coords = scaler.inverse_transform(pca.inverse_transform(flattened_atoms))
    atoms.setCoords(updated_coords)
    return atoms


def flatten_pdb(
    input_file: str, output_file: str, n_components: int, start=None, end=None
):
    atoms = parsePDB(input_file)
    before, middle, after = select_residues(atoms, start, end)
    flattened_middle = flatten_atoms(middle, n_components=n_components)
    flattened_atoms = put_together((before, flattened_middle, after))
    atoms.setCoords(flattened_atoms)
    os.makedirs(Path(output_file).parent, exist_ok=True)
    writePDB(output_file, atoms, renumber=False)


def noise_gaussian_atoms(atoms, std=1):
    noise = np.random.normal(scale=std, size=atoms.getCoords().shape)
    result = atoms.getCoords() + noise
    return result


def noise_gaussian_pdb(
    input_file: str, output_file: str, std: float = 1, start=None, end=None
):
    atoms = parsePDB(input_file)
    before, middle, after = select_residues(atoms, start, end)
    noised_middle = noise_gaussian_atoms(middle, std=std)
    noised_atoms = put_together((before, noised_middle, after))
    atoms.setCoords(noised_atoms)
    os.makedirs(Path(output_file).parent, exist_ok=True)
    writePDB(output_file, atoms, renumber=False)


def randomize_atoms(atoms, radius=20):
    coords = atoms.getCoords()
    height = coords.shape[0]
    randomized_atoms = 2 * radius * np.random.random((height, 3)) - radius
    atoms.setCoords(randomized_atoms)
    return atoms


def randomize_pdb(input_file: str, output_file: str, radius=20, start=None, end=None):
    atoms = parsePDB(input_file)
    before, middle, after = select_residues(atoms, start, end)
    randomized_middle = randomize_atoms(middle, radius=radius)
    randomized_atoms = put_together((before, randomized_middle, after))
    atoms.setCoords(randomized_atoms)
    os.makedirs(Path(output_file).parent, exist_ok=True)
    writePDB(output_file, atoms, renumber=False)


def rotate_atoms(atoms, angles=[0, 0, 0]):
    coords = atoms.getCoords()
    angles = default_angles[angles] if isinstance(angles, str) else angles
    r = R.from_euler("xyz", angles, degrees=True)
    atoms.setCoords(r.apply(coords))
    return atoms


def rotate_pdb(input_file: str, output_file: str, angles=[90, 60, -180]):
    atoms = parsePDB(input_file)
    rotated_atoms = rotate_atoms(atoms, angles)
    os.makedirs(Path(output_file).parent, exist_ok=True)
    writePDB(output_file, rotated_atoms, renumber=False)


def select_backbone(input_file: str, output_file: str):
    atoms = parsePDB(input_file)
    atoms = atoms.select("backbone")
    os.makedirs(Path(output_file).parent, exist_ok=True)
    writePDB(output_file, atoms, renumber=False)


def place_c_beta(c, n, ca, length=1.522, angle=1.927, dihedral=-2.143):
    """
    Function used to add C-Beta to glycine residues
    input: 3 coords (c,n,ca), (L)ength, (A)ngle, and (D)ihedral
    output: C-Beta coords
    https://github.com/jproney/AF2Rank/blob/master/test_templates.py#L164C5-L164C5
    """
    N = lambda x: x / np.sqrt(np.square(x).sum(-1, keepdims=True) + 1e-8)
    nca = N(n - ca)
    n = N(np.cross(n - c, nca))
    m = [nca, np.cross(n, nca), n]
    d = [
        length * np.cos(angle),
        length * np.sin(angle) * np.cos(dihedral),
        -length * np.sin(angle) * np.sin(dihedral),
    ]
    return ca + sum([m * d for m, d in zip(m, d)])


def format_c_beta(ca_atom, placed_c_beta_coordinates: np.array):
    atom_identifier = int(ca_atom.getResindices()[0]) + 1
    resnum = ca_atom.getResnums()[0]
    return f"ATOM    {str(atom_identifier):>3}  CB  {ca_atom.getResnames()[0]} {ca_atom.getChids()[0]} {str(resnum):>3}    {placed_c_beta_coordinates[0]:>8.3f}{placed_c_beta_coordinates[1]:>8.3f}{placed_c_beta_coordinates[2]:>8.3f}  1.00 54.89           C  "


def generate_c_beta(atoms, resnum, length=1.522, angle=1.927, dihedral=-2.143):
    c = atoms.select(f"resnum {resnum} name C")
    n = atoms.select(f"resnum {resnum} name N")
    ca = atoms.select(f"resnum {resnum} name CA")
    placed_c_beta_coordinates = place_c_beta(
        c.getCoords()[0],
        n.getCoords()[0],
        ca.getCoords()[0],
        length=length,
        angle=angle,
        dihedral=dihedral,
    )
    return format_c_beta(c, placed_c_beta_coordinates)


def select_bb_c_beta(
    input_file: str,
    output_file: str,
    replacement_seq=None,
    only_backbone=False,
    beta_mode=None,
):
    atoms = parsePDB(input_file)
    output = []
    if replacement_seq and len(replacement_seq) == len(set(atoms.getResnums())):
        for iterator, res_num in enumerate(set(atoms.getResnums())):
            curr_residue = atoms.select(f"resnum {res_num}")
            curr_residue.setResnames(replacement_seq[iterator])
    elif replacement_seq:
        print("===== Replacement seq and current protein have different length.")
    if only_backbone:
        atoms = atoms.select("backbone")
        for resnum in set(atoms.getResnums()):
            if not (atoms.select(f"resnum {resnum}").getResnames()[0] == "GLY"):
                output.append(generate_c_beta(atoms, resnum))
    else:
        atoms = atoms.select("backbone or name CB")
    if beta_mode and beta_mode == "generate":
        glycine_resnums = set(atoms.select(f"resname GLY").getResnums())
        for glycine_resnum in glycine_resnums:
            if atoms.select(f"resnum {glycine_resnum} and name CB") is None:
                output.append(generate_c_beta(atoms, glycine_resnum))
    elif beta_mode and beta_mode == "dummy":
        atoms.select("name CB").setCoords([0.001, 0.000, 0.000])
    writePDB(output_file, atoms, renumber=False)
    with open(output_file, "a") as file:
        for line in output:
            file.write(line + "\n")


def average_files(input_file1: str, input_file2: str, output_file: str):
    atoms1 = parsePDB(input_file1)
    atoms2 = parsePDB(input_file2)
    transformation = calcTransformation(atoms1, atoms2)
    atoms1 = transformation.apply(atoms1)
    coords1 = atoms1.getCoords()
    coords2 = atoms2.getCoords()
    avg = (coords1 + coords2) / 2
    atoms1.setCoords(avg)
    os.makedirs(Path(output_file).parent, exist_ok=True)
    writePDB(output_file, atoms1, renumber=False)
