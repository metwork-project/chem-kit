import re
from pydantic import BaseModel, Field
from rdkit import Chem
from rdkit.Chem import rdFMCS
from chem_kit.molecule import Molecule


class SimplifierParams(BaseModel):
    hetero_atoms: bool = Field(
        True,
        description="Include hetero atoms close to transformation site",
    )
    aromatic: bool = Field(
        True,
        description="Include aromatic cycles close to transformation site",
    )
    conjugated: bool = Field(
        False,
        description="Include conjugated bonds close to transformation site",
    )
    cycles: bool = Field(
        False,
        description="Include cycles close to transformation site",
    )


class TransformationSimplifier:

    simplified_smarts = None

    def __init__(self, smarts, **params):
        self.params = SimplifierParams(**params)
        self._full_smarts = smarts
        self.simplify_smarts()

    @property
    def full_smarts(self):
        return self._full_smarts

    def simplify_smarts(self):
        mols = [Chem.MolFromSmarts(sm) for sm in self.full_smarts.split(">>")]
        mcs = rdFMCS.FindMCS(mols).queryMol

        _ = [mol.GetSubstructMatch(mcs) for mol in mols]

        mol_simplifiers = [MoleculeSimplifier(mol, self.params) for mol in mols]

        for ms in mol_simplifiers:
            ms.init_keep()

        self.keep_for_bond_diff(mol_simplifiers)

        for ms in mol_simplifiers:
            ms.propagate_map_to_keep()

        self.keep_for_bond_diff(mol_simplifiers)

        self.join_mapped_atoms(mol_simplifiers)

        smarts = [Chem.MolToSmarts(ms.mol) for ms in mol_simplifiers]
        splits = [sm.split(".") for sm in smarts]

        splits = [
            [
                {
                    "smarts": sm,
                    "mapped": set(re.findall(r"\[#\d+(?:\&(?:\++|-))?:(\d+)\]", sm)),
                }
                for sm in split
            ]
            for split in splits
        ]
        smarts = []
        for dic1 in splits[0]:
            for dic2 in splits[1]:
                if dic1["mapped"] == dic2["mapped"]:
                    smarts.append(">>".join([dic1["smarts"], dic2["smarts"]]))

        smarts = list(set(smarts) - {">>"})

        self.simplified_smarts = smarts

    def join_mapped_atoms(self, mol_simplifiers):
        map_to_keep = set([])
        for ms in mol_simplifiers:
            map_to_keep = map_to_keep | ms.map_to_keep
        for ms in mol_simplifiers:
            ms.map_to_keep = map_to_keep
            ms.remove_atoms()

    @staticmethod
    def keep_for_bond_diff(mol_simplifiers):
        map_to_keep = set()
        bond_types = [ms.get_bond_types() for ms in mol_simplifiers]
        bond_types_keys = [set(bond_type.keys()) for bond_type in bond_types]
        bond_types_keys = set().union(*bond_types_keys)
        for key in bond_types_keys:
            values = set([bond_type.get(key, tuple()) for bond_type in bond_types])
            if len(values) > 1:
                map_to_keep = map_to_keep | set(key)
        for ms in mol_simplifiers:
            ms.map_to_keep = ms.map_to_keep | map_to_keep


class MoleculeSimplifier:
    def __init__(self, mol, params: SimplifierParams):
        self.params = params
        mol = Molecule.resolve_charge(mol)
        self.mol = mol
        self.mol_smiles = Molecule(Chem.MolToSmiles(mol), preserve_H=True)
        self.map_to_keep = set()

    def append_map_to_keep(self, value):
        self.map_to_keep = self.map_to_keep | {value}

    def init_keep(self):
        matched_atoms = getattr(self.mol, "__sssAtoms")
        all_idx = {atom.GetIdx() for atom in self.mol.GetAtoms()}
        idx_to_keep = all_idx - set(matched_atoms)
        self.map_to_keep = {
            self.mol.GetAtomWithIdx(idx).GetAtomMapNum() for idx in idx_to_keep
        } | {0}

    def propagate_map_to_keep(self):
        self.propagate_hetero(to_connector=True)
        if self.params.cycles:
            self.propagate_ring()
            self.propagate_hetero(to_connector=False)

    def propagate_hetero(self, to_connector):
        mol_smiles = self.mol_smiles.rdkit
        for idx in self.idx_to_keep_from_map(mol_smiles):
            atom = mol_smiles.GetAtomWithIdx(idx)
            visitor = AtomVisitor(self, atom, params=self.params)
            visitor.propagate(to_connector=to_connector)

    def propagate_ring(self):
        mol_smiles = self.mol_smiles
        systems = mol_smiles.get_ring_systems()
        idx_to_keep = set([])

        for atom in mol_smiles.get_atoms():
            if atom.GetAtomMapNum() in self.map_to_keep:
                for system in systems:
                    if atom.GetIdx() in system:
                        idx_to_keep = idx_to_keep | system

        map_to_keep = self.map_to_keep
        systems_to_keep = {
            mol_smiles.rdkit.GetAtomWithIdx(idx).GetAtomMapNum() for idx in idx_to_keep
        }
        map_to_keep = map_to_keep | systems_to_keep
        self.map_to_keep = map_to_keep

    def get_idx_to_remove(self):
        mol = self.mol
        idx_to_keep = self.idx_to_keep_from_map(mol)
        all_idx = {atom.GetIdx() for atom in mol.GetAtoms()}
        idx_to_remove = list(all_idx - idx_to_keep)
        idx_to_remove.sort(reverse=True)
        return idx_to_remove

    def remove_atoms(self):
        idx_to_remove = self.get_idx_to_remove()
        rwmol = Chem.RWMol(self.mol)
        for idx in idx_to_remove:
            rwmol.RemoveAtom(idx)
        self.mol = rwmol

    def idx_to_keep_from_map(self, mol):
        return {
            atom.GetIdx()
            for atom in mol.GetAtoms()
            if atom.GetAtomMapNum() in self.map_to_keep
        }

    def get_bond_types(self):
        bond_types = dict()
        for bond in self.mol.GetBonds():
            mapped_atoms = [
                bond.GetBeginAtom().GetAtomMapNum(),
                bond.GetEndAtom().GetAtomMapNum(),
            ]
            mapped_atoms.sort()

            if all(mapped_atoms):
                bond_types[tuple(mapped_atoms)] = bond.GetBondType()
        return bond_types


class AtomVisitor:
    def __init__(self, mol_simplifier, atom, params: SimplifierParams):
        self.params = params
        self.mol_simplifier = mol_simplifier
        self.atom = atom

    @property
    def is_hetero(self):
        return self.atom.GetAtomicNum() != 6

    @property
    def is_aromatic(self):
        return self.atom.GetIsAromatic()

    @property
    def map_num(self):
        return self.atom.GetAtomMapNum()

    @property
    def is_keeped(self):
        return self.map_num in self.mol_simplifier.map_to_keep

    def to_propagate_from(self, from_visitor, bond, to_connector):
        if from_visitor.is_keeped and self.is_keeped:
            return False, False
        if from_visitor.is_aromatic and self.is_aromatic and self.params.aromatic:
            return True, False
        if bond.GetIsConjugated() and self.params.conjugated:
            return True, to_connector
        if self.is_hetero and self.params.hetero_atoms:
            return True, to_connector
        if from_visitor.is_hetero and self.params.hetero_atoms and to_connector:
            return True, self.is_hetero
        return to_connector, False

    def set_keep(self):
        self.mol_simplifier.append_map_to_keep(self.map_num)

    def propagate(self, to_connector):
        for atom, bond in self.connected_atoms():
            visitor = self.__class__(self.mol_simplifier, atom, params=self.params)
            propagate, to_connector = visitor.to_propagate_from(
                self, bond, to_connector=to_connector
            )
            if propagate:
                visitor.set_keep()
                visitor.propagate(to_connector=to_connector)

    def connected_atoms(self):
        for bond in self.atom.GetBonds():
            for atom in (bond.GetBeginAtom(), bond.GetEndAtom()):
                if atom.GetIdx() != self.atom.GetIdx():
                    yield atom, bond
