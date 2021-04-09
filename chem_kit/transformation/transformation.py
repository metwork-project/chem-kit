from __future__ import annotations
import re
from typing import List
from rdkit import Chem
from rdkit.Chem import AllChem
from chem_kit.molecule import Molecule
from ..chem_doodle import ChemDoodle
from .generator import SmartsFromSmiles, SmartsFromMetWorkV0
from .simplifier import TransformationSimplifier, SimplifierParams


class Transformation:
    """
    Instance of (bio-)transformation, using RDKit `Reaction` class.

    It's instantiate by providing SMARTS
    but it's also possible to `from_smiles()`
    and `from_metwork_v0()` class method.

    !!! example
        ```python
        from chem_kit import Transformation
        tsf = Transformation("[#6:1]-[#8:2]-[#1:3]>>[#6:1]-[#8:2]-[#6:3](-[#1])(-[#1])-[#1]")
        ```
    """

    def __init__(self, smarts: str):
        """
        Args:
            smarts: SMARTS of the transformation to create
        """
        self._set_from_smarts(smarts)

    def _set_from_smarts(self, smarts):
        self._smarts = smarts
        self._rdkit = AllChem.ReactionFromSmarts(smarts)

    @classmethod
    def from_smiles(cls, reactant: str, product: str) -> Transformation:
        """
        Args:
            reactant: SMILES of reactant
            product: SMILES of product

        Returns:
            Transformation instance
        """
        return cls(SmartsFromSmiles(reactant, product))

    @classmethod
    def from_metwork_v0(cls, smarts: str) -> Transformation:
        """
        Args:
            smarts: smarts from MetWork v0 reaction database

        Returns:
            Transformation instance
        """
        return cls(SmartsFromMetWorkV0(smarts))

    @property
    def smarts(self) -> str:
        """The SMARTS of the transformation."""
        return self._smarts

    @property
    def rdkit(self) -> Chem.rdChemReactions.ChemicalReaction:
        """
        The `Chem.rdChemReactions.ChemicalReaction` instance
        associate with the transformation.
        """
        return self._rdkit

    @property
    def reactant(self):
        return Molecule.from_rdkit(self.rdkit.GetReactants()[0])

    @property
    def product(self):
        return Molecule.from_rdkit(self.rdkit.GetProducts()[0])

    @property
    def mass_delta(self):
        return self.product.mass - self.reactant.mass

    @property
    def chemdoodle_json(self) -> dict:
        """
        [ChemDoodle JSON](https://web.chemdoodle.com/docs/chemdoodle-json-format)
        format of the transformation.
        """
        return ChemDoodle.react_to_json(self.rdkit)

    def _repr_png_(self):
        return self.rdkit._repr_png_()

    def _repr_svg_(self):
        return self.rdkit._repr_svg_()

    def run(self, smiles: str):
        """
        Run the transformation for a reactant.

        Args:
            smiles: SMILES of the reactant
        """
        react = Molecule(smiles).rdkit
        if re.search(r"(?:#1|H)", self.smarts):
            react = Chem.AddHs(react)
        products = self.rdkit.RunReactant(react, 0)
        res = []
        for r in products:
            try:
                res.append(Molecule.from_rdkit(r[0]))
            except Exception as ex:
                pass
                # print(type(ex))
        return res

    def simplify(self, **params) -> List[Transformation]:
        """
        Simplify the transformation by removing non changing atoms that are not
        'connected' to the reaction site (i.e. atoms that changes).

        Args:
            params: See [SimplifierParams](transformation.md#simplifierparams)

        Returns:
            List of simplified transformations.
        """
        results = TransformationSimplifier(self.smarts, **params)
        return [
            self.__class__(smarts)._reindex_mapping()
            for smarts in results.simplified_smarts
        ]

    def reverse(self) -> Transformation:
        """
        Reverse the way of the transformation,
        i.e. reactant become product and product become reactant.
        """
        smarts = self.smarts.split(">>")
        smarts.reverse()
        self._set_from_smarts(">>".join(smarts))
        return self

    def _reindex_mapping(self):
        mapping_origin = {
            atom.GetAtomMapNum()
            for react in self.rdkit.GetReactants()
            for atom in react.GetAtoms()
        } - {0}
        mapping_origin = list(mapping_origin)
        mapping_origin.sort()
        mapping_change = {origin: new + 1 for new, origin in enumerate(mapping_origin)}
        smarts = ">>".join(
            [
                ".".join(
                    [
                        Chem.MolToSmarts(self._reindex_mol(mol, mapping_change))
                        for mol in mols
                    ]
                )
                for mols in [self.rdkit.GetReactants(), self.rdkit.GetProducts()]
            ]
        )
        self._set_from_smarts(smarts)
        return self

    @staticmethod
    def _reindex_mol(mol, mapping_change):
        for atom in mol.GetAtoms():
            map_num = atom.GetAtomMapNum()
            if map_num:
                atom.SetAtomMapNum(mapping_change.get(map_num, 0))
        return mol
