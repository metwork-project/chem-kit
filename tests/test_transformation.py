import pytest
from rdkit.Chem import AllChem
from chem_kit.transformation import Transformation, TranformationGenerator

SMARTS = "[#8:1]-[#6:2]1:[#6:3]:[#6:5]:[#6:8]:[#6:6]:[#6:4]:1-[#1:7]>>[#8:1]-[#6:2]1:[#6:3]:[#6:5]:[#6:8]:[#6:6]:[#6:4]:1-[#35:7]"


def test_from_smarts():
    tsf = Transformation(SMARTS)
    assert tsf.smarts == SMARTS
    tsf = TranformationGenerator(
        "[#8]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1",
        "[#8]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#35]",
    )
    assert tsf.smarts == SMARTS


def test_from_smiles():
    tsf = Transformation.from_smiles("OC1=CC=CC=C1", "OC1=CC=CC=C1Br")
    assert tsf.smarts == SMARTS


# TODO: More test data
def test_from_metwork_v0():
    tsf = Transformation.from_metwork_v0(
        "[#6:1](=,:[#8:2])(-[#6:3])-[#6:4]>>[#6:4]-[#6:1](-[#8:2]-[H])(-[H])-[#6:3]"
    )
    assert (
        tsf.smarts == "[#6:1](=[#8:2])(-[#6:3])-[#6:4]>>[#6:3]-[#6:1](-[#8:2])-[#6:4]"
    )
    tsf = Transformation.from_metwork_v0(
        "[#6:1](-[#7:2]-[H:3])-[#6:4]-[#6:5]>>[#7:2](-[#6:3]1-[#6](-[#6](-[#6](-[#6](-[#8]-1)-[#6](=,:[#8])-[#8])-[#8])-[#8])-[#8])-[#6:1]-[#6:4]-[#6:5]"
    )
    assert (
        tsf.smarts
        == "[#6:3](-[#7:1]-[#1:2])-[#6:4]-[#6:5]>>[#7:1](-[#6:2]1(-[#6](-[#6](-[#6](-[#6](-[#8]-1)-[#6](=[#8])-[#8])-[#8])-[#8])-[#8])-[H])-[#6:3]-[#6:4]-[#6:5]"
    )


def test_rdkit():
    tsf = Transformation(SMARTS)
    assert AllChem.ReactionToSmarts(tsf.rdkit) == AllChem.ReactionToSmarts(
        AllChem.ReactionFromSmarts(SMARTS)
    )


def test_repr_png(shared_datadir):
    tsf = Transformation(SMARTS)
    assert (shared_datadir / "tsf.png").read_bytes() == tsf._repr_png_()


#  TODO: Get repr_svg
def test_repr_svg():
    tsf = Transformation(SMARTS)
    assert tsf._repr_svg_() is None


def test_run():
    tsf = Transformation.from_smiles("N1C=CC2=C1C=CC=C2", "OC1=CC2=C(NC=C2)C=C1")
    products = tsf.run("N1C=CC2=C1C=CC=C2")
    assert [p.smiles for p in products] == ["Oc1ccc2c(c1)C=CN2"]

    tsf = Transformation.from_smiles("N1C=CC2=C1C=CC=C2", "OC1=CC2=C(NC=C2)C=C1")
    products = tsf.run("N1C=CC2=C1C=CC=C2")
    assert [p.smiles for p in products] == ["Oc1ccc2c(c1)C=CN2"]

    tsf = Transformation.from_smiles("OC1=CC=CC=C1", "OC1=CC=CC=C1Br")
    products = tsf.run("CCC1=CC=CC=C1O")
    assert [p.smiles for p in products] == ["CCc1cccc(Br)c1O"]
    products = tsf.run("CCC1=CC=CC(CC)=C1O")
    assert not products


def test_preserve_explicit_H():
    tsf = Transformation.from_smiles("[H]C([H])([H])OC1=CC=CC=C1", "OC1=CC=CC=C1")
    assert (
        tsf.smarts
        == "[#1]-[#6:3](-[#1])(-[#1])-[#8:1]-[#6:2]1:[#6:4]:[#6:6]:[#6:8]:[#6:7]:[#6:5]:1>>[#8:1](-[#6:2]1:[#6:4]:[#6:6]:[#6:8]:[#6:7]:[#6:5]:1)-[#1:3]"
    )


def test_reverse():
    tsf = Transformation(
        "[#6:1]1:[#6:2]:[#6:4]:[#6:6](-[#8:7]-[#1:8]):[#6:5]:[#6:3]:1>>[#6:8]-[#8:7]-[#6:6]1:[#6:4]:[#6:2]:[#6:1]:[#6:3]:[#6:5]:1"
    )
    tsf.reverse()
    assert (
        tsf.smarts
        == "[#6:8]-[#8:7]-[#6:6]1:[#6:4]:[#6:2]:[#6:1]:[#6:3]:[#6:5]:1>>[#6:1]1:[#6:2]:[#6:4]:[#6:6](-[#8:7]-[#1:8]):[#6:5]:[#6:3]:1"
    )


def test_reindex_mapping():
    tsf = Transformation(
        "[#6:11]-[#8:10]-[#6:9]1:[#6:7]:[#6:5]:[#6:4]:[#6:6]:[#6:8]:1>>[#6:4]1:[#6:5]:[#6:7]:[#6:9](-[#8:10]-[#1:11]):[#6:8]:[#6:6]:1"
    )
    tsf.reindex_mapping()
    assert (
        tsf.smarts
        == "[#6:8]-[#8:7]-[#6:6]1:[#6:4]:[#6:2]:[#6:1]:[#6:3]:[#6:5]:1>>[#6:1]1:[#6:2]:[#6:4]:[#6:6](-[#8:7]-[#1:8]):[#6:5]:[#6:3]:1"
    )


def test_simplify():
    tsf = Transformation.from_smiles(
        "OCC1OC(Oc2cc(O)c3c(occc3=O)c2O)C(O)C(O)C1O",
        "CC(=O)OCC1OC(Oc2cc(O)c3c(occc3=O)c2O)C(O)C(O)C1O",
    )
    assert (
        tsf.simplify()[0].smarts
        == "[#8:1](-[#6:2])-[#1:3]>>[#6]-[#6:3](=[#8])-[#8:1]-[#6:2]"
    )


@pytest.mark.parametrize(
    "smiles_set,smarts_simplified",
    (
        (
            ("CCC1=CC=C(O)C=C1", "CCC1=CC=C(O)C(O)=C1"),
            "[#6:1]1:[#6:2]:[#6:4]:[#6:6](-[#8:8]):[#6:5](:[#6:3]:1)-[#1:7]>>[#6:1]1:[#6:2]:[#6:4]:[#6:6](-[#8:8]):[#6:5](-[#8:7]-[#1]):[#6:3]:1",
        ),
        (
            (
                "CC(=O)OCC1OC(OC2C(Oc3cc(O)c4c(=O)cc(-c5ccc(O)cc5)oc4c3O)OC(CO)C(O)C2O)C(O)C(O)C1O",
                "CC(=O)OCC1OC(OC2C(Oc3cc(O)c4c(=O)cc(-c5ccc(O)c(O)c5)oc4c3O)OC(CO)C(O)C2O)C(O)C(O)C1O",
            ),
            "[#6:1]1:[#6:2]:[#6:4]:[#6:6](-[#8:8]):[#6:5](:[#6:3]:1)-[#1:7]>>[#6:1]1:[#6:2]:[#6:4]:[#6:6](-[#8:8]):[#6:5](-[#8:7]-[#1]):[#6:3]:1",
        ),
        (
            (
                "Oc1ccc(CCC(=O)c2ccc(O)cc2O)cc1",
                "O=C(CCc1ccc(O)cc1)c1ccc(O)cc1OC1OC(CO)C(O)C(O)C1O",
            ),
            "[#6:1]1:[#6:2]:[#6:4]:[#6:7](-[#8:9]):[#6:5]:[#6:3]:1-[#8:6]-[#1:8]>>[#6:1]1:[#6:2]:[#6:4]:[#6:7](-[#8:9]):[#6:5]:[#6:3]:1-[#8:6]-[#6:8]1(-[#8]-[#6](-[#6]-[#8])-[#6](-[#8])-[#6](-[#8])-[#6]-1-[#8])-[#1]",
        ),
        (
            (
                r"COc1ccc(\C=C\N2C(=O)N\C(=C/CCNC([NH3+])=N)C2=O)cc1",
                r"COc1ccc(/C=C/N2C(=O)N/C(=C\CCNC(=N)[NH3+])C2=O)cc1Br",
            ),
            "[#8:1]-[#6:2]1:[#6:3]:[#6:5]:[#6:8]:[#6:6]:[#6:4]:1-[#1:7]>>[#8:1]-[#6:2]1:[#6:3]:[#6:5]:[#6:8]:[#6:6]:[#6:4]:1-[#35:7]",
        ),
    ),
)
def test_simplify_keep_aromatization(smiles_set, smarts_simplified):
    tsf = Transformation.from_smiles(*smiles_set)
    assert tsf.simplify()[0].smarts == smarts_simplified


@pytest.mark.parametrize(
    "smarts_metwork,smarts_simplified",
    (
        (
            r"[#6:4]1-[#6:3]2(-[#6:8](-[#6:7]-[#6:6]-[#6:5]-1)/[#6:12]=,:[#6:13]\[#6:14]#[#6:15])-[#6:2]-[#6:1]-[#6:11]-[#6:10]-[#7:9]-2>>[#6:4]1-[#6:3]2(-[#6:8](-[#6:7]-[#6:6]-[#6:5]-1)/[#6:12]=,:[#6:13]\[#6:14]=,:[#6:15])-[#6:2]-[#6:1]-[#6:11]-[#6:10]-[#7:9]-2",
            "[#6:1]-[#6:2]#[#6:3]>>[#6:1]-[#6:2]=[#6:3]",
        ),
        (
            "[#6:4]1:[#6:3](:[#6:8](:[#6:7]:[#6:6]:[#6:5]:1)-[#8:9])-[#8:2]-[#6:1]>>[#6:4]1:[#6:3]2:[#6:8](:[#6:7]:[#6:6]:[#6:5]:1)-[#8:9]-[#6:1]-[#8:2]-2",
            "[#6:1]1:[#6:2](:[#6:4](:[#6:7]:[#6:6]:[#6:3]:1)-[#8:8])-[#8:5]-[#6:9]>>[#6:1]1:[#6:2]2:[#6:4](:[#6:7]:[#6:6]:[#6:3]:1)-[#8:8]-[#6:9]-[#8:5]-2",
        ),
    ),
)
def test_simplify_keep_bond_diff(smarts_metwork, smarts_simplified):
    tsf = Transformation.from_metwork_v0(smarts_metwork)
    assert tsf.simplify()[0].smarts == smarts_simplified


@pytest.mark.parametrize(
    "smarts_metwork,smarts_simplified",
    (
        (
            "[#6:1]1-[#6:2]2-[#6:3](-[#7:4]-[#6](-[#6:5]-1)-[#6:6]-[#6:7]-[#6:8])-[#6:9]-[#6:10]-[#6:11]-[#6:12]-2-[#6:13]>>[#6:1]1-[#6:2]2-[#6:3](-[#7:4]-[#6](-[#6:5]-1)-[#6:6]-[#6:7]-[#6:8])-[#6:9]-[#6:10]-[#6:11](-[#6:12]-2-[#6:13])-[#8]",
            "[#6:1]1-[#6:2]2-[#6:4](-[#7:7]-[#6:6]-[#6:3]-1)-[#6:8]-[#6:10]-[#6:9](-[#6:5]-2)-[#1:11]>>[#6:1]1-[#6:2]2-[#6:4](-[#7:7]-[#6:6]-[#6:3]-1)-[#6:8]-[#6:10]-[#6:9](-[#6:5]-2)-[#8:11]-[#1]",
        ),
    ),
)
def test_simplify_keep_cycles(smarts_metwork, smarts_simplified):
    tsf = Transformation.from_metwork_v0(smarts_metwork)
    assert tsf.simplify(cycles=True)[0].smarts == smarts_simplified


@pytest.mark.parametrize(
    "smiles_set,smarts_simplified",
    (
        (
            (
                "CC(=O)OCC1OC(OC2C(Oc3cc(O)c4c(=O)cc(-c5ccc(O)cc5)oc4c3O)OC(CO)C(O)C2O)C(O)C(O)C1O",
                "CC(=O)OCC1OC(OC2C(Oc3cc(O)c4c(=O)cc(-c5ccc(O)c(O)c5)oc4c3O)OC(COC(C)=O)C(O)C2O)C(O)C(O)C1O",
            ),
            [
                "[#6:1]1:[#6:2]:[#6:4]:[#6:6](-[#8:8]):[#6:5](:[#6:3]:1)-[#1:7]>>[#6:1]1:[#6:2]:[#6:4]:[#6:6](-[#8:8]):[#6:5](-[#8:7]-[#1]):[#6:3]:1",
                "[#6:1]-[#8:2]-[#1:3]>>[#6:1]-[#8:2]-[#6:3](-[#6])=[#8]",
            ],
        ),
        (
            (
                "COc1ccc(-c2cc(=O)c3c(O)cc(OC4OC(COC(C)=O)C(O)C(O)C4OC4OC(COC(C)=O)C(O)C(O)C4O)c(O)c3o2)cc1",
                "CC(=O)OCC1OC(OC2C(Oc3cc(O)c4c(=O)cc(-c5ccc(O)cc5)oc4c3O)OC(CO)C(O)C2O)C(O)C(O)C1O",
            ),
            [
                "[#6:8](-[#8:7]-[#6:6]1:[#6:4]:[#6:2]:[#6:1]:[#6:3]:[#6:5]:1)(-[#1])(-[#1])-[#1]>>[#6:1]1:[#6:2]:[#6:4]:[#6:6](-[#8:7]-[#1:8]):[#6:5]:[#6:3]:1",
                "[#6:1]-[#8:2]-[#6:3](-[#6])=[#8]>>[#6:1]-[#8:2]-[#1:3]",
            ],
        ),
    ),
)
def test_simplify_to_many_transformations(smiles_set, smarts_simplified):
    tsf = Transformation.from_smiles(*smiles_set)
    assert set([simp.smarts for simp in tsf.simplify()]) == set(smarts_simplified)


@pytest.mark.parametrize(
    "hetero_atoms,smarts",
    (
        (False, "[#8:1]-[#1:2]>>[#6:2](-[#8:1])(-[#1])(-[#1])-[#1]"),
        (
            True,
            "[#8:1](-[#6:2]-[#8:4])-[#1:3]>>[#6:3](-[#8:1]-[#6:2]-[#8:4])(-[#1])(-[#1])-[#1]",
        ),
    ),
)
def test_simplify_keep_hetero_atoms(hetero_atoms, smarts):
    smiles = [r"OCO", r"COCO"]
    tsf = Transformation.from_smiles(*smiles)
    simp = tsf.simplify(hetero_atoms=hetero_atoms)[0]
    assert simp.smarts == smarts


@pytest.mark.parametrize(
    "conjugated,smarts",
    (
        (
            False,
            "[#6:1]1:[#6:2]:[#6:4](-[#8:6]-[#1:8]):[#6:7]:[#6:5]:[#6:3]:1>>[#6:1]1:[#6:2]:[#6:4](-[#8:6]-[#6:8](-[#1])(-[#1])-[#1]):[#6:7]:[#6:5]:[#6:3]:1",
        ),
        (
            True,
            "[#6:1]=[#6:2]-[#6:3]1:[#6:4]:[#6:6](-[#8:8]-[#1:10]):[#6:9]:[#6:7]:[#6:5]:1>>[#6:1]=[#6:2]-[#6:3]1:[#6:4]:[#6:6](-[#8:8]-[#6:10](-[#1])(-[#1])-[#1]):[#6:9]:[#6:7]:[#6:5]:1",
        ),
    ),
)
def test_simplify_keep_conjugated(conjugated, smarts):
    smiles = [r"CC\C=C\C1=CC(O)=CC=C1", r"CC\C=C\C1=CC(OC)=CC=C1"]
    tsf = Transformation.from_smiles(*smiles)
    simp = tsf.simplify(conjugated=conjugated)[0]
    assert simp.smarts == smarts
