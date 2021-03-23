import re
from rdkit import Chem
from rdkit.Chem import AllChem
from chem_kit.molecule import Molecule
from ..chem_doodle import ChemDoodle
from .generator import SmartsFromSmiles, SmartsFromMetWorkV0
from .simplifier import TransformationSimplifier


class Transformation:
    """
    Instance of (bio-)transformation, a.k.a Reaction in RDKit.

    It's instantiate by providing SMARTS but it's also possible to `from_smiles()` and `from_metwork_v0()` class method.

    !!! example
        ```python
        from chem_kit import Transformation
        tsf = Transformation("[#6:1]-[#8:2]-[#1:3]>>[#6:1]-[#8:2]-[#6:3](-[#1])(-[#1])-[#1]")
        ```
    """

    def __init__(self, smarts):
        self._set_from_smarts(smarts)

    def _set_from_smarts(self, smarts):
        self._smarts = smarts
        self._rdkit = AllChem.ReactionFromSmarts(smarts)

    @classmethod
    def from_smiles(cls, reactant: str, product: str) -> "Transformation":
        """
        Args:
            reactant: SMILES of reactant
            product: SMILES of product
        """
        return cls(SmartsFromSmiles(reactant, product))

    @classmethod
    def from_metwork_v0(cls, smarts: str) -> "Transformation":
        """
        Args:
            smarts: smarts from MetWork v0 reaction database
        """
        return cls(SmartsFromMetWorkV0(smarts))

    @property
    def smarts(self) -> str:
        return self._smarts

    @property
    def rdkit(self) -> Chem.rdChemReactions.ChemicalReaction:
        return self._rdkit

    @property
    def chemdoodle_json(self) -> dict:
        return ChemDoodle.react_to_json(self.rdkit)

    def _repr_png_(self):
        return self.rdkit._repr_png_()

    def _repr_svg_(self):
        return self.rdkit._repr_svg_()

    def run(self, smiles: str):
        """
        Args:
            smiles: SMILES of reactant.
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

    def simplify(self, **params):
        results = TransformationSimplifier(self.smarts, **params)
        return [
            self.__class__(smarts)._reindex_mapping()
            for smarts in results.simplified_smarts
        ]

    def reverse(self):
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
