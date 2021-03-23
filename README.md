# ChemKit

<p align="center">
<a href="https://pypi.org/project/chem-kit" target="_blank">
    <img src="https://img.shields.io/pypi/v/chem-kit?color=%2334D058&label=pypi%20package" alt="Package version">
</a>
</p>

---

**Documentation**: <a href="http://chem-kit.metwork.science/" target="_blank">http://chem-kit.metwork.science/</a>

**Source Code**: <a href="https://github.com/YannBeauxis/chem-kit" target="_blank">https://github.com/YannBeauxis/chem-kit</a>

---

ChemKit is a chemical toolbox based on [RDKit](https://www.rdkit.org/) with currently 2 main purposes :

- Facilitate the usage of the [RDKIt Python API](https://www.rdkit.org/docs/api-docs.html)
 with some more easy to use classes that can occasionally fix some bug (especially with Jupyter rendering).

- Provide custom method for the [MetWork](http://www.metwork.science) project

## Usage

### Manipulate Molecules

```python
    from chem_kit import Molecule
    mol = Molecule("CCO")
```

### Manipulate Transformation

```python
    from chem_kit import Transformation
    tsf = Transformation("[#6:1]-[#8:2]-[#1:3]>>[#6:1]-[#8:2]-[#6:3](-[#1])(-[#1])-[#1]")
```

More examples with [Jupyter notebook](usage_example)
