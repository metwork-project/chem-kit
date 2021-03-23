import re
from rdkit import Chem
from rdkit.Chem import AllChem
from chem_kit.molecule import Molecule
from ..chem_doodle import ChemDoodle
from .generator import SmartsFromSmiles, SmartsFromMetWorkV0
from .simplifier import TransformationSimplifier

class Transformation:
    def __init__(self, smarts):
        self.set_from_smarts(smarts)

    def set_from_smarts(self, smarts):
        self._smarts = smarts
        self._rdkit = AllChem.ReactionFromSmarts(smarts)

    @classmethod
    def from_smiles(cls, reactant, product):
        return cls(SmartsFromSmiles(reactant, product))

    @classmethod
    def from_metwork_v0(cls, smarts):
        return cls(SmartsFromMetWorkV0(smarts))

    @property
    def smarts(self):
        return self._smarts

    @property
    def rdkit(self):
        return self._rdkit

    @property
    def mdl_rxn(self):
        return Chem.rdChemReactions.ReactionToRxnBlock(self.rdkit)

    @property
    def chemdoodle_json(self):
        return ChemDoodle.react_to_json(self.rdkit)

    def _repr_png_(self):
        return self.rdkit._repr_png_()

    def _repr_svg_(self):
        return self.rdkit._repr_svg_()

    def run(self, smiles):
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
            self.__class__(smarts).reindex_mapping()
            for smarts in results.simplified_smarts
        ]

    def reverse(self):
        smarts = self.smarts.split(">>")
        smarts.reverse()
        self.set_from_smarts(">>".join(smarts))
        return self

    def reindex_mapping(self):
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
                        Chem.MolToSmarts(self.reindex_mol(mol, mapping_change))
                        for mol in mols
                    ]
                )
                for mols in [self.rdkit.GetReactants(), self.rdkit.GetProducts()]
            ]
        )
        self.set_from_smarts(smarts)
        return self

    @staticmethod
    def reindex_mol(mol, mapping_change):
        for atom in mol.GetAtoms():
            map_num = atom.GetAtomMapNum()
            if map_num:
                atom.SetAtomMapNum(mapping_change.get(map_num, 0))
        return mol


