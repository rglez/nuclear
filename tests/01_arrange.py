# Created by roy.gonzalez-aleman at 19/09/2023
import fnmatch
import os
import shutil
from os.path import join, abspath, dirname


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

for cfg in cfgs:
    # process = subprocess.Popen([nuclear_py, cfg], stdout=subprocess.PIPE,
    #                            stderr=subprocess.PIPE)
    # output, error = process.communicate()
    # print(output)
    # if process.returncode != 0:
    #     raise Exception(
    #         f"File handling failed {process.returncode} {output} {error}")
    os.system(f'{nuclear_py} {cfg}')
