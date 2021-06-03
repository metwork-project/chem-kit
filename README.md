# ChemKit

<p align="center">
<a href="https://pypi.org/project/chem-kit" target="_blank">
    <img src="https://img.shields.io/pypi/v/chem-kit?color=%2334D058&label=pypi%20package" alt="Package version">
</a>
</p>

---

**Documentation**: <a href="http://chem-kit.metwork.science/" target="_blank">http://chem-kit.metwork.science/</a>

**Source Code**: <a href="https://github.com/metwork-project/chem-kit" target="_blank">https://github.com/metwork-project/chem-kit</a>

---

ChemKit is a chemical toolbox based on [RDKit](https://www.rdkit.org/) with currently 2 main purposes :

- Facilitate the usage of the [RDKIt Python API](https://www.rdkit.org/docs/api-docs.html)
 with some more easy to use classes that can occasionally fix some bug (especially with Jupyter rendering).

- Provide tailored methods for the [MetWork](http://www.metwork.science) project

## Usage

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

More examples with [Jupyter notebook](http://chem-kit.metwork.science/jupyter_example/)

## Install

Like RDKit, ChemKit needs [Conda](https://docs.conda.io) :

```bash
conda env create -f conda-env.yml
conda activate chem_kit
```

To manage other required Python packages, 
the better way is to use [Poetry](https://python-poetry.org) on top of Conda :

```bash
cd /path/to/chem_kit
poetry install
```

> Poetry manipulate Python packages directly on Conda env.
