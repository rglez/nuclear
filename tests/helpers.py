# Created by roy.gonzalez-aleman at 18/09/2023
import fnmatch
import os
import shutil


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
