from rdkit import Chem
from rdkit.Chem import rdFMCS
from chem_kit.molecule import Molecule, AROMATICITY_MODEL


class SmartsGenerator:

    EXPLICIT_H_MAP_FROM = 1000

    def __new__(cls, *arg, **kwargs):
        smarts = cls.get_smarts(*arg, **kwargs)
        return smarts

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
            # atomCompare=rdFMCS.AtomCompare.CompareElements,
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


class SmartsFromSmiles(SmartsGenerator):
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


class SmartsFromMetWorkV0(SmartsGenerator):
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
