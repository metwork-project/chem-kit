import re
from rdkit import Chem


class ChemDoodle:
    """
    Utils to use [ChemDoolde JSON format](https://web.chemdoodle.com/docs/chemdoodle-json-format)
    """

    @classmethod
    def mol_rdkit_to_json(cls, mr, begin_id={"a": 0, "b": 0}):

        json_res = {"a": [], "b": []}
        Chem.rdDepictor.Compute2DCoords(mr)
        ZOOM = 20
        positions = mr.GetConformer().GetPositions()
        BOND_TYPE_DIC = {
            Chem.rdchem.BondType.SINGLE: 1,
            Chem.rdchem.BondType.DOUBLE: 2,
            Chem.rdchem.BondType.TRIPLE: 3,
        }
        chiral_bonds = {}
        CHIRAL_CONFIG = {
            Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW: {
                0: "protruding",
                2: "recessed",
            },
            Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW: {0: "recessed", 2: "protruding"},
        }

        for i, a in enumerate(mr.GetAtoms()):
            mol_json = {
                "i": "a{0}".format(i + begin_id["a"]),
                "x": ZOOM * positions[i][0],
                "y": ZOOM * positions[i][1],
            }
            symbols = re.findall(r"\[([^:]*?)(?:\:\d)?\]", a.GetSmarts())
            if len(symbols) > 0:
                symbols = symbols[0].split(",")
            if len(symbols) > 1:
                v = []
                if "*" in symbols:
                    v = ["a"]
                else:
                    for sy in symbols:
                        find_at = re.findall("#(\d)", sy)
                        if len(find_at) == 1:
                            try:
                                at_value = int(find_at[0])
                            except:
                                at_value = find_at[0]
                            v.append(Chem.Atom(at_value).GetSymbol())
                        else:
                            v.append(sy)
                mol_json["q"] = {"as": {"v": v, "n": False}}
            else:
                symbol = a.GetSymbol()
                if symbol != "C":
                    mol_json["l"] = symbol
                charge = a.GetFormalCharge()
                if charge != 0:
                    mol_json["c"] = charge
                elif a.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
                    ct = a.GetChiralTag()
                    for i, b in enumerate(a.GetBonds()):
                        if i in CHIRAL_CONFIG[ct]:
                            chiral_bonds[b.GetIdx()] = CHIRAL_CONFIG[ct][i]
            json_res["a"].append(mol_json)

        try:
            # Force atom to be aromatic if is in ring
            for a in mr.GetAtoms():
                a.SetIsAromatic(a.IsInRing())
            # Kekulize mol
            Chem.rdmolops.SanitizeMol(mr)
            Chem.Kekulize(mr, True)

        except:
            # Transform aromatic bonds in kekulize form
            # This should not be used ...
            arom_bond_type = Chem.rdchem.BondType.SINGLE

            def propagate_kekulize(atom, arom_bond_type):
                for b in atom.GetBonds():
                    if b.GetBondType() == Chem.rdchem.BondType.AROMATIC:
                        if arom_bond_type == Chem.rdchem.BondType.SINGLE:
                            arom_bond_type = Chem.rdchem.BondType.DOUBLE
                        else:
                            arom_bond_type = Chem.rdchem.BondType.SINGLE
                        b.SetBondType(arom_bond_type)
                        if atom.GetIdx() != b.GetBeginAtomIdx():
                            target = b.GetBeginAtom()
                        else:
                            target = b.GetEndAtom()
                        propagate_kekulize(target, arom_bond_type)

            for b in mr.GetBonds():
                if b.GetBondType() == Chem.rdchem.BondType.AROMATIC:
                    b.SetBondType(arom_bond_type)
                    propagate_kekulize(b.GetBeginAtom(), arom_bond_type)

        for i, b in enumerate(mr.GetBonds()):
            b_id = b.GetIdx()
            bond_json = {
                "i": "b{0}".format(i + begin_id["b"]),
                "b": b.GetBeginAtomIdx(),
                "e": b.GetEndAtomIdx(),
            }
            bond_type = b.GetBondType()
            if bond_type != Chem.rdchem.BondType.SINGLE:
                bond_json["o"] = BOND_TYPE_DIC[bond_type]
            if b_id in chiral_bonds:
                bond_json["s"] = chiral_bonds[b_id]
            json_res["b"].append(bond_json)
        return json_res

    @classmethod
    def react_to_json(cls, reaction):
        PADDING = 80
        ARROW_LENGTH = 60
        json_res = {
            "m": [],
            "s": [],
        }
        x_bound = 0
        atom_maping = {}
        begin_id = {"a": 0, "b": 0}

        def append_mol(m, mol_type):
            m_json = cls.mol_rdkit_to_json(m, begin_id)
            x_max = 0
            x_min = 0
            for a in m_json["a"]:
                x_min = min(x_min, a["x"])
            for a in m_json["a"]:
                a["x"] += x_bound - x_min
                x_max = max(x_max, a["x"])
            for i, a in enumerate(m.GetAtoms()):
                map_num = a.GetAtomMapNum()
                if map_num > 0:
                    atom_id = "a{0}".format(i + begin_id["a"])
                    if mol_type == "reactant":
                        atom_maping[map_num] = [atom_id, None]
                    elif mol_type == "product":
                        atom_maping[map_num][1] = atom_id
            json_res["m"].append(m_json)
            begin_id["a"] += len(m_json["a"])
            begin_id["b"] += len(m_json["b"])

            return x_max

        # Reactants
        for m in reaction.GetReactants():
            x_bound += append_mol(m, "reactant") + PADDING

        # Arrow
        x_bound -= PADDING / 2
        json_res["s"].append(
            {
                "i": "s0",
                "t": "Line",
                "x1": x_bound,
                "y1": 0,
                "x2": x_bound + ARROW_LENGTH,
                "y2": 0,
                "a": "synthetic",
            }
        )
        x_bound += ARROW_LENGTH + PADDING

        # Products
        for m in reaction.GetProducts():
            x_bound += append_mol(m, "product")

        # AtomMapping
        for i, map_ in atom_maping.items():
            json_res["s"].append(
                {
                    "i": "s{0}".format(i),
                    "t": "AtomMapping",
                    "a1": map_[0],
                    "a2": map_[1],
                }
            )

        return json_res
