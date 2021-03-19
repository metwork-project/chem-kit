import numpy as np
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole

AROMATICITY_MODEL = Chem.rdmolops.AromaticityModel.AROMATICITY_MDL


class Molecule:
    """
    Usage:
         
        
        ```python
        from chem_kit import Molecule
        mol = Molecule("CC")
        ```
    """

    _AROMATICITY_MODEL = AROMATICITY_MODEL

    def __init__(
        self,
        smiles: str,
        aromaticity: str = None,
        preserve_H: bool = False,
    ):
        """
        Args:
            smiles: SMILES of the molecule to create.
            aromaticity:
                Aromacity model to apply.

                Possible values are :

                - custom
                - default
                - mdl
                - rdkit
                - simple

                If None, default is "mdl".
            preserve_H: Preserve H that are explicit in smiles.

        """
        self._set_aromaticity_model(aromaticity)
        self._set_rdkit(smiles, preserve_H=preserve_H)
        self._smiles = Chem.MolToSmiles(self.rdkit)

    def __eq__(self, other):
        if not isinstance(other, Molecule):
            return NotImplemented

        return self.smiles == other.smiles

    @classmethod
    def from_rdkit(cls, rdkit):
        smiles = Chem.MolToSmiles(rdkit)
        return cls(smiles)

    def _set_aromaticity_model(self, aromaticity):
        if aromaticity:
            self._AROMATICITY_MODEL = getattr(
                Chem.rdmolops.AromaticityModel,
                "AROMATICITY_{}".format(aromaticity.upper()),
            )

    def _set_rdkit(self, smiles, preserve_H):
        if preserve_H:
            mol = self._mol_from_smiles_preserve_H(smiles)
        else:
            mol = Chem.MolFromSmiles(smiles)
        Chem.Kekulize(mol, True)
        Chem.SetAromaticity(mol, self._AROMATICITY_MODEL)
        self._rdkit = mol

    @staticmethod
    def _mol_from_smiles_preserve_H(smiles):
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        sanitize_params = [
            # "ADJUSTHS",
            # "ALL" ,
            "CLEANUP",
            "CLEANUPCHIRALITY",
            "FINDRADICALS",
            "KEKULIZE",
            # "NONE",
            "PROPERTIES",
            "SETAROMATICITY",
            "SETCONJUGATION",
            "SETHYBRIDIZATION",
            "SYMMRINGS",
        ]
        arr = np.array(
            [
                getattr(Chem.rdmolops.SanitizeFlags, "SANITIZE_" + param)
                for param in sanitize_params
            ]
        )
        sanitize_param = int(np.bitwise_or.reduce(arr))
        Chem.SanitizeMol(mol, sanitize_param)
        return mol

    @property
    def smiles(self):
        return self._smiles

    @property
    def smiles_kekulized(self):
        return Chem.MolToSmiles(self.rdkit_kekulized)

    @property
    def rdkit(self):
        return self._rdkit.__copy__()

    @property
    def rdkit_kekulized(self):
        mol = self.rdkit
        Chem.Kekulize(mol, clearAromaticFlags=True)
        return mol

    def _repr_svg_(self):
        # Modified to highlight bonds, from :
        # https://github.com/rdkit/rdkit/blob/bc4d5478478542045a0f89e28a5da8b6d1f0e5d0/rdkit/Chem/Draw/IPythonConsole.py#L112
        mol = self.rdkit
        return Chem.Draw._moltoSVG(
            mol,
            IPythonConsole.molSize,
            getattr(mol, "__sssAtoms", []),
            "",
            IPythonConsole.kekulizeStructures,
            drawOptions=IPythonConsole.drawOptions,
            highlightBonds=getattr(mol, "__sssBonds", []),
        )

    def get_substructure_match(self, smarts):
        substructure = Chem.MolFromSmarts(smarts)
        mol = self.rdkit
        matches = mol.GetSubstructMatches(substructure)
        highlight_atoms = getattr(mol, "__sssAtoms", [])
        highlight_bonds = [
            bond.GetIdx()
            for bond in mol.GetBonds()
            if self._is_highlight_bond(bond, highlight_atoms)
        ]
        setattr(self._rdkit, "__sssAtoms", highlight_atoms)
        setattr(self._rdkit, "__sssBonds", highlight_bonds)
        return matches

    def _is_highlight_bond(self, bond, highlight_atoms):
        return set(self._get_bond_atoms(bond)).issubset(set(highlight_atoms))

    @staticmethod
    def _get_bond_atoms(bond):
        return [bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()]

    def highlight_aromatics(self):
        mol = self.rdkit
        ri = mol.GetRingInfo()
        aromatic_bonds = []
        aromatic_atoms = []
        for ring_bonds in ri.BondRings():
            for bond_idx in ring_bonds:
                bond = mol.GetBondWithIdx(bond_idx)
                if mol.GetBondWithIdx(bond_idx).GetIsAromatic():
                    aromatic_bonds.append(bond_idx)
                    aromatic_atoms += self._get_bond_atoms(bond)
        setattr(self._rdkit, "__sssAtoms", aromatic_atoms)
        setattr(self._rdkit, "__sssBonds", aromatic_bonds)
        return self

    def iter_atoms(self):
        return self.rdkit.GetAtoms()

    def get_ring_systems(self, includeSpiro=False):
        return get_ring_systems(self.rdkit, includeSpiro=includeSpiro)


def get_ring_systems(mol, includeSpiro=False):
    ring_info = mol.GetRingInfo()
    systems = []
    for ring in ring_info.AtomRings():
        ringAts = set(ring)
        nSystems = []
        for system in systems:
            nInCommon = len(ringAts.intersection(system))
            if nInCommon and (includeSpiro or nInCommon > 1):
                ringAts = ringAts.union(system)
            else:
                nSystems.append(system)
        nSystems.append(ringAts)
        systems = nSystems
    return systems