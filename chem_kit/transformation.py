import re
from rdkit import Chem
from rdkit.Chem import rdFMCS, AllChem
from chem_kit.molecule import Molecule, AROMATICITY_MODEL
from chem_kit.transformation_reductor import TransformationReductor
from .chem_doodle import ChemDoodle

class Transformation:
    def __init__(self, smarts):
        self.set_from_smarts(smarts)

    def set_from_smarts(self, smarts):
        self._smarts = smarts
        self._rdkit = AllChem.ReactionFromSmarts(smarts)

    @classmethod
    def from_smiles(cls, reactant, product):
        return TransformationFromSmiles(reactant, product)

    @classmethod
    def from_metwork_v0(cls, smarts):
        return TransformationFromMetWorkV0(smarts)

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
        reductor = TransformationReductor(self.smarts, **params)
        return [
            self.__class__(smarts).reindex_mapping()
            for smarts in reductor.reduced_smarts
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


class TranformationGenerator:

    EXPLICIT_H_MAP_FROM = 1000

    def __new__(cls, *arg, **kwargs):
        smarts = cls.get_smarts(*arg, **kwargs)
        return Transformation(smarts)

    @classmethod
    def get_smarts(cls, *source):
        mols = cls.list_mols(*source)
        mols = [Chem.AddHs(mol) for mol in mols]
        matches = cls.get_common_atoms(*mols)
        cls.map_matches(matches, mols)
        h_mapped = cls.get_h_mapped(mols)
        mols = [cls.remove_h_non_mapped(mol) for mol in mols]
        mols = cls.explicit_h_when_mapped_with_h(mols, h_mapped)
        mols_smarts = [Chem.MolToSmarts(mol) for mol in mols]
        return ">>".join(mols_smarts)

    @staticmethod
    def get_h_mapped(mols):
        h_mapped = []
        for mol in mols:
            for atom in mol.GetAtoms():
                map_num = atom.GetAtomMapNum()
                if atom.GetAtomicNum() == 1 and map_num != 0:
                    h_mapped.append(map_num)
        return h_mapped

    @staticmethod
    def explicit_h_when_mapped_with_h(mols, h_mapped):
        mols_res = []
        for mol in mols:
            mol = Chem.MolFromSmarts(Chem.MolToSmarts(mol))
            mol.UpdatePropertyCache(strict=False)
            idxs = [
                atom.GetIdx()
                for atom in mol.GetAtoms()
                if (atom.GetAtomMapNum() in h_mapped)
            ]
            if idxs:
                mol = Chem.AddHs(mol, onlyOnAtoms=idxs)
            mols_res.append(mol)
        return mols_res

    @staticmethod
    def list_mols(*smarts):
        mols = [Chem.MolFromSmarts(sm) for sm in smarts]
        for mol in mols:
            Chem.SanitizeMol(mol)
        return mols

    @classmethod
    def get_common_atoms(cls, *mols):
        mcs = rdFMCS.FindMCS(
            mols,
            atomCompare=rdFMCS.AtomCompare.CompareAny,
            bondCompare=rdFMCS.BondCompare.CompareAny,
            # ringCompare=Chem.rdFMCS.RingCompare.StrictRingFusion,
            # ringCompare=Chem.rdFMCS.RingCompare.PermissiveRingFusion,
            # ringMatchesRingOnly=True,
            # completeRingsOnly=True,
        ).queryMol
        return [mol.GetSubstructMatch(mcs) for mol in mols]

    @staticmethod
    def map_matches(matches, mols):
        map_count = 1
        for atoms_idx in zip(*matches):
            atoms = [mol.GetAtomWithIdx(idx) for mol, idx in zip(mols, atoms_idx)]
            if set(atom.GetSymbol() for atom in atoms) != {"H"}:
                [atom.SetAtomMapNum(map_count) for atom in atoms]
                map_count += 1

    @classmethod
    def remove_h_non_mapped(cls, mol):
        rwmol = Chem.RWMol(mol)
        atoms_to_remove = []
        for atom in rwmol.GetAtoms():
            if atom.GetAtomMapNum() == 0 and atom.GetSymbol() == "H":
                atoms_to_remove.append(atom.GetIdx())
            if atom.GetAtomMapNum() > cls.EXPLICIT_H_MAP_FROM:
                atom.SetAtomMapNum(0)
        atoms_to_remove.sort(reverse=True)
        for idx in atoms_to_remove:
            rwmol.RemoveAtom(idx)
        return rwmol


class TransformationFromSmiles(TranformationGenerator):
    @classmethod
    def list_mols(cls, *smiles):
        mols = [Molecule(sm, preserve_H=True).rdkit for sm in smiles]
        for mol in mols:
            h_map_idx = cls.EXPLICIT_H_MAP_FROM + 1
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() == 1:
                    atom.SetAtomMapNum(h_map_idx)
                    h_map_idx += 1
        return mols


class TransformationFromMetWorkV0(TranformationGenerator):
    @classmethod
    def get_smarts(cls, smarts):
        smarts = smarts.replace("=,:", "=")
        smarts = smarts.split(">>")
        mols = [Chem.MolFromSmarts(sm) for sm in smarts]
        mols = [cls.unset_mapping(mol) for mol in mols]
        for mol in mols:
            Chem.Kekulize(mol)
            Chem.SetAromaticity(mol, AROMATICITY_MODEL)
        return super().get_smarts(*mols)

    @staticmethod
    def list_mols(*source):
        return source

    @staticmethod
    def unset_mapping(mol):
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        return mol
