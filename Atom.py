import numpy as np


class Atom:
    def __init__(self):
        self._DataType = 'ATOM'  # "ATOM"/"HETATM"
        self._Name = ''  # Atom name
        self._AltLoc = ''  # Alternate location indicator.
        self._ResName = ''  # Residue name
        self._ChainId = 'A'  # Chain identifier
        self._ResSeq = 0  # Residue sequence number
        self._ResCode = ''  # Code for insertions of residues
        self._Coordinate = np.array([0, 0, 0])  # (X,Y,Z)
        self._Occup = 0.0  # Occupancy
        self._TempFac = 0.0  # Temperature factor
        self._Element = 'Xx'  # Element symbol
        self._Charge = 0.0  # Atom charge

    def __eq__(self, other: 'Atom'):
        if self._Name != other.Name:
            return False
        if self._ResName != other.ResName:
            return False
        if self._ChainId != other.ChainId:
            return False
        if self._ResSeq != other.ResSeq:
            return False
        return True

    @property
    def DataType(self):
        return self._DataType

    @DataType.setter
    def DataType(self, value):
        self._DataType = value

    @property
    def Name(self):
        return self._Name

    @Name.setter
    def Name(self, value):
        self._Name = value

    @property
    def AltLoc(self):
        return self._AltLoc

    @AltLoc.setter
    def AltLoc(self, value):
        self._AltLoc = value

    @property
    def ResName(self):
        return self._ResName

    @ResName.setter
    def ResName(self, value):
        self._ResName = value

    @property
    def ChainId(self):
        return self._ChainId

    @ChainId.setter
    def ChainId(self, value):
        self._ChainId = value

    @property
    def ResSeq(self):
        return self._ResSeq

    @ResSeq.setter
    def ResSeq(self, value):
        self._ResSeq = value

    @property
    def Coordinate(self):
        return np.array(self._Coordinate)

    @Coordinate.setter
    def Coordinate(self, value):
        self._Coordinate = np.array(value)

    @property
    def Occup(self):
        return self._Occup

    @Occup.setter
    def Occup(self, value):
        self._Occup = value

    @property
    def TempFac(self):
        return self._TempFac

    @TempFac.setter
    def TempFac(self, value):
        self._TempFac = value

    @property
    def Element(self):
        return self._Element

    @Element.setter
    def Element(self, value):
        self._Element = value

    @property
    def Charge(self):
        return self._Charge

    @Charge.setter
    def Charge(self, value):
        self._Charge = value
