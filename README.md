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
