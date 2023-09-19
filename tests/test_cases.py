# Created by roy.gonzalez-aleman at 17/09/2023
import filecmp
import pickle
from os.path import join, dirname, abspath

import pytest


def unpickle_from_file(file_name):
    """ Unserialize a **pickle** file.

    Args:
        file_name (str): file to unserialize.
    Returns:
        (object): an unserialized object.
    """
    with open(file_name, 'rb') as file:
        data = pickle.load(file)
    return data


tests_dir = abspath(dirname(__file__))


@pytest.fixture()
def reference_and_target_fixture(request):
    cwd = dirname(request.module.__file__)
    examples_dir = join(cwd, 'examples')

    reference_data_path = join(examples_dir, 'reference_data.pick')
    reference_data = unpickle_from_file(reference_data_path)

    target_data_path = join(examples_dir, 'target_data.pick')
    target_data = unpickle_from_file(target_data_path)

    return reference_data, target_data


def test_cases_correspondence(reference_and_target_fixture):
    reference_data, target_data = reference_and_target_fixture
    assert set(target_data.keys()) == set(reference_data.keys())


def test_files_correspondence(reference_and_target_fixture):
    reference_data, target_data = reference_and_target_fixture

    for key in reference_data:
        # assert produced files are the same
        assert key in target_data


def test_files_identity(reference_and_target_fixture):
    reference_data, target_data = reference_and_target_fixture
    for key in reference_data:
        # assert all files are identical
        gs_files = reference_data[key]
        target_files = target_data[key]
        for gs_file in gs_files:
            gs_file_path = gs_files[gs_file]
            target_file_path = target_files[gs_file]
            if 'nuclear.cfg' not in gs_file_path and 'html' not in gs_file_path:
                assert filecmp.cmp(gs_file_path, target_file_path,
                                   shallow=False), \
                    f'\n{target_file_path} offending identity comparison.'
