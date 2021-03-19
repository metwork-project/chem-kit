from chem_kit import Molecule
from rdkit import Chem

SMILES = "OC1=CC2=C(NC=C2)C=C1"
CANONICAL_SMILES = "Oc1ccc2c(c1)C=CN2"
OTHER_SMILES = "CCCCO"


def test_mol_from_smiles():
    mol = Molecule(CANONICAL_SMILES)
    assert mol.smiles == CANONICAL_SMILES
    mol = Molecule(SMILES)
    assert mol.smiles == CANONICAL_SMILES


def test_equal():
    assert Molecule(CANONICAL_SMILES) == Molecule(SMILES)
    assert Molecule(CANONICAL_SMILES) != Molecule(OTHER_SMILES)
    assert Molecule(SMILES) != SMILES


def test_mol_from_rdkit():
    mol_rdkit = Chem.MolFromSmiles(SMILES)
    mol = Molecule.from_rdkit(mol_rdkit)
    assert mol.smiles == CANONICAL_SMILES


def test_aromaticity():
    mol_default = Molecule(CANONICAL_SMILES)
    mol_aromat_simple = Molecule(CANONICAL_SMILES, aromaticity="simple")
    assert mol_default != mol_aromat_simple


def test_smiles_kekulized():
    mol = Molecule(CANONICAL_SMILES)
    assert mol.smiles_kekulized == "OC1=CC=C2NC=CC2=C1"


def test__repr_svg_(shared_datadir):
    svg = Molecule(CANONICAL_SMILES)._repr_svg_()
    assert svg == (shared_datadir / "mol.svg").read_text()


def test_get_substructure_match():
    mol = Molecule("OC1=CC2=C(NC=C2)C=C1")
    mol.get_substructure_match("c1ccccc1")

    assert mol._rdkit.__sssAtoms == [1, 2, 3, 4, 8, 9]
    assert mol._rdkit.__sssBonds == [1, 2, 3, 7, 8, 9]


def test_highlight_aromatics():
    mol = Molecule("OC1=CC2=C(NC=C2)C=C1")
    mol.highlight_aromatics()
    assert mol._rdkit.__sssAtoms == [9, 1, 8, 9, 4, 8, 3, 4, 2, 3, 1, 2, 3, 4]
    assert mol._rdkit.__sssBonds == [9, 8, 7, 3, 2, 1, 3]


def test_preserve_H():
    smiles = "[H]C([H])([H])Oc1ccccc1"
    mol = Molecule(smiles)
    assert mol.smiles == "COc1ccccc1"
    mol = Molecule(smiles, preserve_H=True)
    assert mol.smiles == smiles