import numpy as np
import argparse
from function import function_shift, function_shift_dipoles
from utils import GetPrevAndNext, Normalized
from molecule import Molecule

approximated_negative_charge_shift = function_shift("Neggrid.txt")
approximated_positive_charge_shift = function_shift("Posgrid.txt")
approximated_dipole_shift = function_shift_dipoles()

parser = argparse.ArgumentParser(description=" This program calculate max lambda for chrompohore ")
parser.add_argument("-pdbfilename", default="1M0L.pdb",
                    type=str,
                    help="Name of file in format: X.pdb")
parser.add_argument("-resid", default=210,
                    type=int,
                    help="Aminoacide id in pdb file")
parser.add_argument("-chrid", default=212,
                    type=int,
                    help="Chromophore id in pdb")
args = parser.parse_args()

PDB_FILENAME = args.pdbfilename  # имя PDB-структуры
RESID = args.resid  # номер аминокислоты в PDB, указываемый пользователем
CHROMOPHORE_RESID = args.chrid  # номер хромофора в PDB

# CHROMOPHORE
CHROMOPHORE_CONSID_ATOMS = ["C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12", "C13", "C14", "C15", "NZ"]
CHROMOPHORE_CONSID_ATOMS_A = ["C3"] + CHROMOPHORE_CONSID_ATOMS + ["HZ1"]
CYLINDER_1 = ["C4", "C5", "C6"]
CYLINDER_2 = ["C7", "C8", "C9", "C10", "C11"]
CYLINDER_3 = ["C11", "C12", "C13", "C14", "C15", "NZ"]
RING_ATOM = "C6"
# CHROMOPHORE

# RUN
molecule = Molecule()
molecule.ReadFromPdb(PDB_FILENAME)

# AXIS CYLINDERS
atom1, atom2 = molecule.GetAtomByRes(CHROMOPHORE_RESID, "C4"), molecule.GetAtomByRes(CHROMOPHORE_RESID, "C6")
AXIS_CYLINDER_1 = atom2.Coordinate - atom1.Coordinate

atom1, atom2 = molecule.GetAtomByRes(CHROMOPHORE_RESID, "C7"), molecule.GetAtomByRes(CHROMOPHORE_RESID, "C11")
AXIS_CYLINDER_2 = atom2.Coordinate - atom1.Coordinate

atom1, atom2 = molecule.GetAtomByRes(CHROMOPHORE_RESID, "C11"), molecule.GetAtomByRes(CHROMOPHORE_RESID, "NZ")
AXIS_CYLINDER_3 = atom2.Coordinate - atom1.Coordinate
# AXIS CYLINDERS

CALC_TYPE = molecule.GetAminoacidType(RESID)

resid_position = None
dipole_vec = None
negative_charge_pos = None
charge_sign = ''

if CALC_TYPE == "C":
    charge_sign, resid_position = molecule.GetChargePosition(RESID)
elif CALC_TYPE == "P":
    resid_position, dipole_vec, negative_charge_pos, _ = molecule.GetDipolePosition(RESID)
elif CALC_TYPE == "O":
    print("lambda max 0 nm, neutral non-polar residue")
    exit()

chrNearestAtom = None
atom_best, v_best, distance_best = None, None, None

for atom_name in CHROMOPHORE_CONSID_ATOMS:
    atom = molecule.GetAtomByRes(CHROMOPHORE_RESID, atom_name)
    v = resid_position - atom.Coordinate
    distance = np.linalg.norm(v)
    if distance_best is None or distance < distance_best:
        v_best = np.array(v)
        atom_best = atom
        chrNearestAtom = atom_name
        distance_best = distance

if chrNearestAtom is None:
    print("Cannot define atom in chromophore")
    exit()
else:
    if chrNearestAtom == "NZ":
        chrNearestAtom = 16
    else:
        chrNearestAtom = int(chrNearestAtom.replace("C", ""))

v = v_best
atom1 = atom_best
atom_name2, atom_name3 = GetPrevAndNext(atom1.Name, CHROMOPHORE_CONSID_ATOMS_A)
atom2 = molecule.GetAtomByRes(CHROMOPHORE_RESID, atom_name2)
atom3 = molecule.GetAtomByRes(CHROMOPHORE_RESID, atom_name3)

x1 = atom1.Coordinate
x2 = atom2.Coordinate
x3 = atom3.Coordinate

angle = None

if CALC_TYPE == "C":
    if atom1.Name in CYLINDER_1:
        AXIS_CYLINDER = AXIS_CYLINDER_1
    elif atom1.Name in CYLINDER_3:
        AXIS_CYLINDER = AXIS_CYLINDER_3
    elif atom1.Name in CYLINDER_2:
        AXIS_CYLINDER = AXIS_CYLINDER_2
    else:
        assert False, "Logical error"

    n1 = Normalized(np.cross(x2 - x1, x3 - x1))
    n2 = Normalized(AXIS_CYLINDER)
    n3 = Normalized(np.cross(n1, n2))

    n = Normalized(np.cross(n1, n3))
    u = Normalized(v - v.dot(n) * n)

    angle = np.arccos(u.dot(n1))

    idx = molecule.GetAtomIdx(atom1)
    print(f"{CALC_TYPE} {charge_sign} {idx} {distance_best:.5f} {angle:.5f}")
    if charge_sign == "+":
        lambda_max_val = approximated_dipole_shift
    else:
        lambda_max_val = approximated_negative_charge_shift(chrNearestAtom, distance_best, angle)
    print(f"lambda max {lambda_max_val:.6f} nm")

if CALC_TYPE == "P":
    n = Normalized(np.cross(x2 - x1, x3 - x1))
    u = Normalized(v)

    sign = 1 if n.dot(dipole_vec) > 0 else -1

    angle = np.arccos(sign * u.dot(n))

    vec_to_ring = molecule.GetAtomByRes(CHROMOPHORE_RESID, RING_ATOM).Coordinate - negative_charge_pos
    orientation = "L" if np.dot(dipole_vec, vec_to_ring) > 0 else "R"

    idx = molecule.GetAtomIdx(atom1)
    lambda_max_val = approximated_dipole_shift(chrNearestAtom, distance_best, angle, orientation)

    print(f"{CALC_TYPE} {idx} {orientation} {distance_best:.5f} {angle:.5f}")
    print(f"lambda max {lambda_max_val:.6f} nm")
