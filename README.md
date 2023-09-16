<!-- These are examples of badges you might want to add to your README:
     please update the URLs accordingly

[![Built Status](https://api.cirrus-ci.com/github/<USER>/nuclear.svg?branch=main)](https://cirrus-ci.com/github/<USER>/nuclear)
[![ReadTheDocs](https://readthedocs.org/projects/nuclear/badge/?version=latest)](https://nuclear.readthedocs.io/en/stable/)
[![Coveralls](https://img.shields.io/coveralls/github/<USER>/nuclear/main.svg)](https://coveralls.io/r/<USER>/nuclear)

[![Conda-Forge](https://img.shields.io/conda/vn/conda-forge/nuclear.svg)](https://anaconda.org/conda-forge/nuclear)
[![Monthly Downloads](https://pepy.tech/badge/nuclear/month)](https://pepy.tech/project/nuclear)
[![Twitter](https://img.shields.io/twitter/url/http/shields.io.svg?style=social&label=Twitter)](https://twitter.com/nuclear)
-->

# NUCLEAR
NUCLEAR (*NUCLeotide AssembleR*) is a Python package designed primarily to generate oligonucleotide sequences from one or several mono-nucleotide explorations returned by the MCSS docking method.

In a single job, users can request one of two exclusive explorations: either
a molecular hotspots search (to gain insights into the most accessible regions of the receptor) or an oligonucleotide search. In the latter case, the oligo-nucleotide sequence and the receptor region to consider can also be specified.


# Installation
Create and activate a new virtual environment with your preferred tool (here we use conda), and then install NUCLEAR from PyPI using pip.

```
conda create -n nuclear
source activate nuclear
pip install -r requirements.txt nuklear
```

After installation, you should be able to run NUCLEAR as follows:

```$ nuclear path-to-config-file.cfg```


# Documentation
The extensive documentation of the project will be soon available.


# Changelog
All notable changes to this project will be documented [in this file](https://github.com/rglez/nuclear/blob/main/CHANGELOG.md).


# Contributing
Please read the [contributor's guide](https://github.com/rglez/nuclear/blob/main/CONTRIBUTING.md) if you are interested in advancing the project.

### Note
This project has been set up using PyScaffold 4.5. For details and usage
information on PyScaffold see https://pyscaffold.org/.
