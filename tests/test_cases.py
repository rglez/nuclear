# Created by roy.gonzalez-aleman at 17/09/2023
import filecmp
import os
import shutil
import tarfile
from os.path import join, basename
from subprocess import run

import pytest

from . import helpers as hp

# >>>> directories declaration <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# tests_dir = '/home/roy.gonzalez-aleman/RoyHub/nuclear/tests'
tests_dir = os.path.abspath(os.path.dirname(__file__))
examples_dir = join(tests_dir, 'examples')
results_dir = join(examples_dir, 'results')
gs_basename = 'results-GS-001'
gs_dir = os.path.join(examples_dir, gs_basename)
gs_tar = join(examples_dir, f'{gs_basename}.tar')

# >>>> clean tests data <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
hp.remove_dir(gs_dir)
hp.remove_dir(results_dir)

# >>>> run nuclear on each cfg case <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
cfgs = list(hp.finder('*.cfg', examples_dir))
os.chdir(examples_dir)
nuclear_py = shutil.which('nuclear')


def test_nuclear_is_callable():
    assert nuclear_py


runs = [run(f'{nuclear_py} {cfg}', shell=True) for cfg in cfgs]


@pytest.fixture
def gold_standard_data():
    with tarfile.open(gs_tar) as tar:
        tar.extractall(examples_dir)

    gs_data = dict()
    for case in os.listdir(gs_dir):
        case_files = {basename(x): x for x in
                      hp.finder('*', join(gs_dir, case))}
        gs_data.update({case: case_files})
    return gs_data


@pytest.fixture
def target_data():
    results_data = dict()
    for case in os.listdir(results_dir):
        case_files = {basename(x): x for x in
                      hp.finder('*', join(results_dir, case))}
        results_data.update({case: case_files})
    return results_data


def test_same_number_of_cases(target_data, gold_standard_data):
    assert set(target_data.keys()) == set(gold_standard_data.keys())


def test_same_number_of_files_produced(target_data, gold_standard_data):
    for key in gold_standard_data:

        # test produced files are the same
        assert key in target_data

        # test files identity
        gs_files = gold_standard_data[key]
        target_files = target_data[key]
        for gs_file in gs_files:
            gs_file_path = gs_files[gs_file]
            target_file_path = target_files[gs_file]
            if 'nuclear.cfg' not in gs_file_path and 'html' not in gs_file_path:
                assert filecmp.cmp(gs_file_path, target_file_path,
                                   shallow=False), \
                    f'\n{target_file_path} offending identity comparison.'

    # cleaning after testing
    # hp.remove_dir(gs_dir)
    # hp.remove_dir(results_dir)
