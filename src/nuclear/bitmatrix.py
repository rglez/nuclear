#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  9 22:00:32 2021

@author: roy.gonzalez-aleman
"""
import platform
import itertools as it
import multiprocessing as mp
from multiprocessing import Pool

import numpy as np
import pandas as pd
import scipy.spatial.ckdtree as ckd
from bitarray import bitarray as ba


if platform.system() != 'Linux':
    mp.set_start_method('fork')

def from_array_to_bitarray(array, N):
    zero_arr = np.zeros(N, dtype=bool)
    zero_arr[array] = 1
    bitarr = ba()
    bitarr.pack(zero_arr.tobytes())
    return bitarr


def get_pose_clashes(chunk, whole_kd, dist, poses_idx, N):
    chunk_arrays = []
    for fragment_kd in chunk:
        contact_list = fragment_kd.query_ball_tree(whole_kd, r=dist)
        contact_conc = pd.unique(np.fromiter(it.chain.from_iterable(contact_list), dtype=int))
        contact_real = pd.unique(poses_idx[contact_conc])
        bitarray = from_array_to_bitarray(contact_real, N)
        chunk_arrays.append(bitarray)
    return chunk_arrays


# def get_pose_clashes2(fragment_kd, whole_kd, dist, poses_idx, N):
#     contact_list = fragment_kd.query_ball_tree(whole_kd, r=dist)
#     contact_conc = pd.unique(np.fromiter(it.chain.from_iterable(contact_list), dtype=int))
#     contact_real = pd.unique(poses_idx[contact_conc])
#     bitarray = from_array_to_bitarray(contact_real, N)
#     return bitarray


def split_list(list_, N):
    return np.array_split(list_, N)


class builder:
    def __init__(self, args, sampling):
        self.args = args
        self.sampling = sampling

    # @profile
    def get_clash_matrix(self):
        ncores = self.args.ncores
        pool = Pool(processes=ncores)

        # create a kdtree with all fragment coordinates merged
        dist = self.args.clash_dist
        poses_idx = self.sampling.poses_idx
        poses_all = self.sampling.poses_all
        N = len(self.sampling.poses_all)
        tmp1 = [self.sampling.poses_all[x].data for x in range(N)]
        tmp2 = np.concatenate(tmp1)
        tmp1 = None
        whole_kd = ckd.cKDTree(tmp2)
        tmp2 = None
        # compute clashes of every fragment in the whole_kd
        chunks = split_list(poses_all, ncores)
        chunks_idx = split_list(range(N), ncores)
        inputs = [[chunk, whole_kd, dist, poses_idx, N] for chunk in chunks]
        results = pool.starmap(get_pose_clashes, inputs)
        clash_matrix = {chunks_idx[i][idx]: y for i, x in enumerate(results)
                        for idx, y in enumerate(x)}
        self.clash_matrix = clash_matrix

    # def get_clash_matrix2(self):
    #     # create a kdtree with all fragment coordinates merged
    #     dist = self.args.clash_dist
    #     poses_idx = self.sampling.poses_idx
    #     poses_all = self.sampling.poses_all
    #     N = len(self.sampling.poses_all)
    #     tmp1 = [self.sampling.poses_all[x].data for x in range(N)]
    #     tmp2 = np.concatenate(tmp1)
    #     tmp1 = None
    #     whole_kd = ckd.cKDTree(tmp2)
    #     tmp2 = None
    #     # compute clashes of every fragment in the whole_kd
    #     results = [get_pose_clashes2(fragment_kd, whole_kd, dist, poses_idx, N)
    #                         for i, fragment_kd in enumerate(poses_all)]
    #     clash_matrix = {i: x for i, x in enumerate(results)}
    #     self.clash_matrix = clash_matrix

    def get_link_matrix(self):
        # reduce linker candidates by O3'->C5' distance
        N = len(self.sampling.poses_all)
        C5_kdtree = self.sampling.C5_tree
        O3_kdtree = self.sampling.O3_tree
        # compute clashes of every fragment in the whole_kd
        link_matrix = dict()
        for i, fragment_kd in enumerate(self.sampling.poses_all):
            # under_max_O3C5 = O3_kdtree.query_ball_point(C5_kdtree.data[i], self.args.max_dist_O3C5)
            under_max_O3C5 = C5_kdtree.query_ball_point(O3_kdtree.data[i], self.args.max_dist_O3C5)
            bit_candidates = from_array_to_bitarray(np.fromiter(under_max_O3C5, dtype=int), N)
            link_matrix[i] = bit_candidates & (~self.clash_matrix[i])
        self.link_matrix = link_matrix


# =============================================================================
# Demonstration
# =============================================================================
# import time
# from confread import parser as config_parser
# from sampling import parser as sampling_parser
# config = '/home/roy.gonzalez-aleman/rprojects/nuclear/examples/example.cfg'
# args = config_parser(config)
# sampling = sampling_parser(args)
# self = builder(args, sampling)

# start = time.time()
# self.get_clash_matrix()
# matrix = self.clash_matrix
# print(time.time() - start)

# start = time.time()
# self.get_clash_matrix2()
# matrix2 = self.clash_matrix
# print(time.time() - start)
