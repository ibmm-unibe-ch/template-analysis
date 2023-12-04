import os
from pathlib import Path

import numpy as np
from prody import calcTransformation, parsePDB, writePDB
from scipy.spatial.transform import Rotation as R
from sklearn.decomposition import PCA
from sklearn.preprocessing import RobustScaler

default_angles = {"1": [60, 14, 42], "2": [83, -100, 17], "3": [66, -14, 132]}


def select_residues(atoms, start=None, end=None):
    start = 0 if start is None else start
    end = max(atoms.getResnums()) if end is None else end
    print(start)
    print(end)
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
        print("no")
        if not part is None:
            print("hey")
            print(part.getCoords().shape)
            atoms.append(part.getCoords())
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
    writePDB(output_file, atoms)


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
    print(randomized_atoms.shape)
    atoms.setCoords(randomized_atoms)
    os.makedirs(Path(output_file).parent, exist_ok=True)
    writePDB(output_file, atoms)


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
    writePDB(output_file, rotated_atoms)


def select_backbone(input_file: str, output_file: str):
    atoms = parsePDB(input_file)
    atoms = atoms.select("backbone")
    os.makedirs(Path(output_file).parent, exist_ok=True)
    writePDB(output_file, atoms)


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
    writePDB(output_file, atoms1)
