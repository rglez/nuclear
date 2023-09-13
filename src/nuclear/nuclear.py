#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  9 14:01:12 2021

@author: roy.gonzalez-aleman
"""
import multiprocessing as mp
import sys
import time

import hotspots as hs
import oligogen
from bitmatrix import builder
from confread import parser as config_parser
from minimizer import minimizer
from oligogen import generator
from sampling import parser as sampling_parser

try:
    mp.set_start_method('fork')
except RuntimeError:
    pass

start = time.time()

# >>>> 1. read arguments from the config file
if len(sys.argv) != 2:
    raise ValueError('\nNUCLEAR syntax is: nuclear.py path-to-config-file.cfg')
config = sys.argv[1]
# config = '/home/roy.gonzalez-aleman/RoyHub/Chapter4/hotspots_search/' \
#          '1SGZ_whole_stdRXN300.cfg'
args = config_parser(config)
print('\nElapsed time until reading arguments: {:8.3f} secs'.format(
    time.time() - start))

# >>>> 2. parse the MCSS exploration
sampling = sampling_parser(args)
print('\nElapsed time until parsing MCSS explorations: {:8.3f} secs'.format(
    time.time() - start))

# >>>> 3. Decide main workflow to follow: NUCLEAR or SPOTIFY?
if 'spots_params' in args.sections:
    # .... starting a spotify job
    spotify = hs.spotify(args, sampling)

    # .... retrieve clusters
    spotify.get_supersets()
    spotify.find_spots()
    spotify.get_crdspots()

    # .... analyses
    spotify.get_dataframes()
    spotify.get_overall_best_rank()
    spotify.get_overall_counts()
    spotify.get_countingby_tops()
    spotify.get_axis_views()
    spotify.get_fingerlog()
    spotify.get_fingerlog2()
    spotify.create_reps()

else:
    # .... starting a nuclear job
    matrices = builder(args, sampling)
    matrices.get_clash_matrix()
    matrices.get_link_matrix()
    print(
        '\nElapsed time until building of binary matrices: {:8.3f} secs'.format(
            time.time() - start))

    # search for oligonucleotides
    oligos = generator(args, sampling, matrices)
    top_N_pool, counts, leaders, clust_array, passed = oligos.get_sequences()
    oligos.write_sequences(top_N_pool)
    oligos.write_log(top_N_pool)
    print('\n{} sequences explored in {:8.3f} secs.'
          .format(counts, time.time() - start))
    print('\nNumber of heap entries was {}\n'.format(passed))

    # path counting analysis
    kmers_dict = oligogen.ksplit_top_N_pool(top_N_pool, 1)
    oligos.write_log_popularity(kmers_dict, top_N_pool)

    # produce the minimization inputs ?
    if args.minimize:
        minim = minimizer(args)
        minim.read_and_split_miniprot()
        minim.write_minim_inp()

print('\nNormal termination. Total runtime: {:8.3f} secs'.format(
    time.time() - start))
