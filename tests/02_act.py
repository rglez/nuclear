# Created by roy.gonzalez-aleman at 19/09/2023
import fnmatch
import os
import pickle
import sys
import tarfile
from os.path import basename, join, abspath, dirname


def finder(pattern, root=os.curdir):
    """
    Find all files matching **pattern** inside **root** directory

    Args:
        pattern (str): specified pattern
        root (str)   : path where to begin the recursive search
    """
    for path, dirs, files in os.walk(os.path.abspath(root), followlinks=True):
        for filename in fnmatch.filter(files, pattern):
            yield os.path.join(path, filename)


def pickle_to_file(data, file_name):
    """ Serialize data using **pickle**.

    Args:
        data (object)  : any serializable object.
        file_name (str): name of the **pickle** file to be created.
    Returns:
        (str): file_name
    """
    with open(file_name, 'wb') as file:
        pickle.dump(data, file)
    return file_name


root_dir = dirname(abspath(__file__))
examples_dir = join(root_dir, 'examples')
results_dir = join(root_dir, 'examples', 'results')
gs_dir = join(root_dir, 'examples', 'results-GS-001')
gs_tar = join(root_dir, 'examples', 'results-GS-001.tar')

print(root_dir)
assert os.path.exists(gs_tar)
assert os.path.exists(examples_dir)
assert os.path.exists(results_dir)

# ____ extract this_version data ______________________________________________
results_data = dict()
for case in os.listdir(results_dir):
    case_files = {basename(x): x for x in finder('*', join(results_dir, case))}
    results_data.update({case: case_files})
pickle_to_file(results_data, 'target_data.pick')

# ____ extract gold standard data _____________________________________________
with tarfile.open(gs_tar) as tar:
    tar.extractall(examples_dir)
gs_data = dict()
for case in os.listdir(gs_dir):
    case_files = {basename(x): x for x in finder('*', join(gs_dir, case))}
    gs_data.update({case: case_files})
pickle_to_file(gs_data, 'reference_data.pick')
