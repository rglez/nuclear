# Created by roy.gonzalez-aleman at 17/09/2023
import pickle
from os.path import join, dirname, getsize

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
    # assert same data will be compared
    reference_data, target_data = reference_and_target_fixture
    assert set(target_data.keys()) == set(reference_data.keys())


def test_files_correspondence(reference_and_target_fixture):
    # assert produced files are the same
    reference_data, target_data = reference_and_target_fixture
    for key in reference_data:
        assert key in target_data


def test_files_identity(reference_and_target_fixture):
    # assert all files are identical
    reference_data, target_data = reference_and_target_fixture

    for key in reference_data:
        gs_files = reference_data[key]
        target_files = target_data[key]

        for gs_file in gs_files:
            gs_file_path = gs_files[gs_file]
            target_file_path = target_files[gs_file]
            sms1 = f'File {target_file_path} has different size than reference.'
            sms2 = f'File {target_file_path} is not identical to reference.'

            if ('nuclear.cfg' not in gs_file_path) and (
                    'html' not in gs_file_path):
                assert getsize(gs_file_path) == getsize(target_file_path), sms1
                gs_opened = open(gs_file_path, 'rb').read()
                target_opened = open(target_file_path, 'rb').read()
                assert gs_opened == target_opened, sms2
