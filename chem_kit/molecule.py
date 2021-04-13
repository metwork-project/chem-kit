from __future__ import annotations
from typing import List, Tuple
import re
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Draw import IPythonConsole

AROMATICITY_MODEL = Chem.rdmolops.AromaticityModel.AROMATICITY_MDL


class Molecule:
    """
    Instance of molecule using RDKit `Mol` class.

    It's instantiate by providing SMILES but it's also possible to use `from_rdkit()` class method
    to instantiate from RDKit `Mol` instance.

    !!! example
        ```python
        from chem_kit import Molecule
        mol = Molecule("CCO")
        ```
    """

    def __init__(
        self,
        smiles: str,
        aromaticity: str = "mdl",
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

    @property
    def smiles(self) -> str:
        """Cannonical SMILES of the molecule"""
        return self._smiles

    @property
    def rdkit(self) -> Chem.Mol:
        """`rdkit.Chem.Mol` instance of the molecule. Use it to directly apply RDKit methods"""
        return self._rdkit.__copy__()

    @property
    def smiles_kekulized(self):
        """Kekulized SMILES of the molecule"""
        return Chem.MolToSmiles(self._kekulize_rdkit())

    @property
    def mass(self):
        return Descriptors.MolWt(self.rdkit)

    def _kekulize_rdkit(self):
        mol = self.rdkit
        Chem.Kekulize(mol, clearAromaticFlags=True)
        return mol

    def __eq__(self, other):
        """
        Comparaison betwwen two Molecule instances is done
        by comparing cannonical SMILES (`.smiles` property)
        """
        if not isinstance(other, Molecule):
            return NotImplemented

        return self.smiles == other.smiles

    @classmethod
    def from_rdkit(cls, rdkit: Chem.Mol) -> Molecule:
        """
        Instanciate by providing a RDKit instance rather than a smiles.

        Args:
            rdkit: `rdkit.Chem.Mol` instance

        Returns:
            Molecule instance
        """
        rdkit = cls.resolve_charge(rdkit)
        smiles = Chem.MolToSmiles(rdkit)
        return cls(smiles)

    @staticmethod
    def resolve_charge(rdkit):
        map_charge = re.findall(r"\[#\d+\&(\++):(\d+)]", Chem.MolToSmarts(rdkit))
        map_charge = {int(map_num): len(charge) for charge, map_num in map_charge}
        for atom in rdkit.GetAtoms():
            charge = map_charge.get(atom.GetAtomMapNum())
            if charge:
                atom.SetFormalCharge(charge)
        return rdkit

    def _set_aromaticity_model(self, aromaticity):
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

    def _repr_svg_(self):
        """
        Modified to highlight bonds, from :
        https://github.com/rdkit/rdkit/blob/bc4d5478478542045a0f89e28a5da8b6d1f0e5d0/rdkit/Chem/Draw/IPythonConsole.py#L112
        """
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

    def get_substructure_match(self, smarts: str) -> Tuple[tuple]:
        """
        Apply RDKit `GetSubstructMatches` method.

        !!! note
            It highlight bonds as it seems to not working with the python API of RDKit.

        Args:
            smarts: the substructure to find.

        Returns:
            Tuple of atom ids for each match
        """
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
        """Highlight aromatics rings"""
        mol = self.rdkit
        ring_info = mol.GetRingInfo()
        aromatic_bonds = []
        aromatic_atoms = []
        for ring_bonds in ring_info.BondRings():
            for bond_idx in ring_bonds:
                bond = mol.GetBondWithIdx(bond_idx)
                if mol.GetBondWithIdx(bond_idx).GetIsAromatic():
                    aromatic_bonds.append(bond_idx)
                    aromatic_atoms += self._get_bond_atoms(bond)
        setattr(self._rdkit, "__sssAtoms", aromatic_atoms)
        setattr(self._rdkit, "__sssBonds", aromatic_bonds)
        return self

    def get_atoms(self):
        """Yield every atom in the molecule"""
        return self.rdkit.GetAtoms()

    def get_ring_systems(self, includeSpiro: bool = False) -> List[set]:
        """
        Get all ring systems.

        !!! note
            From [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#count-ring-systems)

        Returns:
            List of set if atom id for each ring systems
        """
        return self._get_ring_systems(self.rdkit, includeSpiro=includeSpiro)

    @staticmethod
    def _get_ring_systems(mol, includeSpiro: bool = False):
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
