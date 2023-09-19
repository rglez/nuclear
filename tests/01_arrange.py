# Created by roy.gonzalez-aleman at 19/09/2023
import fnmatch
import os
import pickle
import shutil
import tarfile
from os.path import join, basename, abspath, dirname
from subprocess import run


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


def remove_dir(path_to_dir):
    try:
        shutil.rmtree(path_to_dir)
    except FileNotFoundError:
        pass


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


# ____ directories declaration ________________________________________________
tests_dir = abspath(dirname(__file__))
# tests_dir = abspath('/home/roy.gonzalez-aleman/RoyHub/nuclear/tests')
examples_dir = join(tests_dir, 'examples')
results_dir = join(examples_dir, 'results')
gs_basename = 'results-GS-001'
gs_dir = join(examples_dir, gs_basename)
gs_tar = join(examples_dir, f'{gs_basename}.tar')

# ____ clean tests data _______________________________________________________
remove_dir(gs_dir)
remove_dir(results_dir)

# ____ run nuclear on each cfg case ___________________________________________
cfgs = list(finder('*.cfg', examples_dir))
os.chdir(examples_dir)
nuclear_py = shutil.which('nuclear')
runs = [run(f'{nuclear_py} {cfg}', shell=True) for cfg in cfgs]

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
