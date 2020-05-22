import numpy as np
from Atom import Atom


class Molecule:
    CHARGED_AMINOACIDS = ["ARG", "HIS", "LYS", "ASP", "GLU"]
    POLAR_AMINOACIDS = ["SER", "THR", "ASN", "GLN", "CYS", "PRO", "TYR", "TRP", "GLH", "ASH"]

    def __init__(self):
        self._Atoms = dict()
        self._Residues = dict()

    @property
    def Atoms(self):
        return self._Atoms

    @property
    def Residues(self):
        return self._Residues

    def GetAtomByRes(self, resid, name) -> 'Atom':
        return self._Residues[resid][name]

    def GetAtomIdx(self, atom):
        for k, v in self._Atoms.items():
            if atom == v:
                return k

    def ReadFromPdb(self, fname):
        self._Residues = dict()
        self._Atoms = dict()

        with open(fname, "r") as file:
            for line in file.readlines():
                data_type = line[0:6].strip()
                if data_type not in ['ATOM', 'HETATM']:
                    continue
                atom = Atom()
                atom.DataType = line[0:6].strip()
                atom.Name = line[12:16].strip()
                atom.AltLoc = line[16].strip()
                atom.ResName = line[17:20].strip()
                atom.ChainId = line[21].strip()
                atom.ResSeq = int(line[22:26])
                atom.ResCode = line[26].strip()
                atom.Coordinate = np.array(list(map(float, [line[30:38], line[38:46], line[46:54]])))
                atom.Occup = 0.0  # float(line[54:60])
                atom.Tempfac = 0.0  # float(line[60:66])
                atom.Element = atom.Name[0]  # line[76:78].strip()
                num = int(line[6:11])

                if atom.ResSeq not in self._Residues:
                    self._Residues[atom.ResSeq] = dict()

                self._Residues[atom.ResSeq][atom.Name] = atom
                self._Atoms[num] = atom

    def SaveToPdb(self, fname):
        with open(fname, "w") as f:
            for (idx, atom) in self._Atoms.items():
                line = '{a0:<6}{a1:>5}{s}{a2:>4}{a3:>1}{a4:>3}{s}{a5:>1}' \
                       '{a6:>4}{a7:<1}{s:>3}{a8[0]:>8.3f}{a8[1]:>8.3f}{a8[2]:>8.3f}' \
                       '{a9:>6.2f}{a10:>6.2f}{s:>11}{a11:<2}\n'.format(
                    a0=atom.DataType, a1=idx, a2=atom.Name, a3=atom.AltLoc,
                    a4=atom.ResName, a5=atom.ChainId, a6=atom.ResSeq,
                    a7=atom.ResCode, a8=atom.Coordinate, a9=atom.Occup,
                    a10=atom.TempFac, a11=atom.Element, s=' '
                )
                f.write(line)

    def GetChargePosition(self, resid):
        position: np.ndarray

        residue = self._Residues[resid]
        resname = residue['C'].ResName

        if resname == "ARG":
            position = residue["NH1"].Coordinate
        elif resname == "HIS":
            position = residue["ND1"].Coordinate
        elif resname == "LYS":
            position = residue["NZ"].Coordinate
        elif resname == "ASP" or resname == "ASH":
            position = 0.5 * (residue["OD1"].Coordinate + residue["OD2"].Coordinate)
        elif resname == "GLU" or resname == "GLH":
            position = 0.5 * (residue["OE1"].Coordinate + residue["OE2"].Coordinate)
        else:
            assert False, f"{resname} not charged"

        return position

    def GetDipolePosition(self, resid):
        positive_charge: str
        negative_charge: str

        residue = self._Residues[resid]
        resname = residue['C'].ResName

        if resname == "SER":
            positive_charge = "HG"
            negative_charge = "OG"
        elif resname == "THR":
            positive_charge = "HG1"
            negative_charge = "OG1"
        elif resname == "ASN":
            positive_charge = "CG"
            negative_charge = "OD1"
        elif resname == "GLN":
            positive_charge = "CD"
            negative_charge = "OE1"
        elif resname == "CYS":
            positive_charge = "HG"
            negative_charge = "SG"
        elif resname == "TYR":
            positive_charge = "HH"
            negative_charge = "OH"
        elif resname == "TRP":
            positive_charge = "HE1"
            negative_charge = "NE1"
        elif resname == "PRO":
            positive_charge = "NH"
            negative_charge = "N"
        elif resname == "GLH":
            positive_charge = "HE2"
            negative_charge = "OE2"
        elif resname == "ASH":
            positive_charge = "HD2"
            negative_charge = "OD2"
        else:
            assert False, f"Logical error"

        positive_charge_pos = residue[positive_charge].Coordinate
        negative_charge_pos = residue[negative_charge].Coordinate

        dipole_vec = positive_charge_pos - negative_charge_pos
        dipole_pos = 0.5 * (positive_charge_pos + negative_charge_pos)

        return dipole_pos, dipole_vec, negative_charge_pos, positive_charge_pos

    def GetAminoacidType(self, resid):
        resname = self._Residues[resid]['C'].ResName
        if resname in self.POLAR_AMINOACIDS:
            return "P"  # Polar
        elif resname in self.CHARGED_AMINOACIDS:
            return "C"  # Charged
        else:
            return "O"  # Other
