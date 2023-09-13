#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 08:35:30 2021

@author: roy.gonzalez-aleman
"""
import heapq
import itertools as it
import os
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
from bitarray import bitarray as ba


def get_tanimoto(A, B):
    intersection = len(A & B)
    if not intersection:
        return 1
    return 1 - (intersection / len(A | B))


def split_list(list_, N):
    return np.array_split(list_, N)


def get_indices(L):
    return it.combinations(range(L + 1), 2)


def ksplit_sequence(sequence, k, indices):
    strseq = [str(x) for x in sequence]
    subchains = [' '.join(strseq[s:e]) for s, e in indices if (e - s >= k)]
    return subchains


def ksplit_top_N_pool(top_N_pool, k):
    countlist = defaultdict(list)
    for size in top_N_pool:
        if size >= k:
            sequences = [top_N_pool[size][x][2]
                         for x in range(len(top_N_pool[size]))]
            for sequence in sequences:
                indices = get_indices(size)
                kmers = ksplit_sequence(sequence, k, indices)
                for kmer in kmers:
                    countlist[len(kmer.split())].append(kmer)
    return countlist


class generator:
    def __init__(self, args, sampling, matrices):
        self.args = args
        self.sampling = sampling
        self.matrices = matrices
        self.refseq = self.args.seq_to_search

    def get_sequences(self):
        # try to explore only nodes with at least 1 degree
        degrees = [self.matrices.link_matrix[x].count()
                   for x in self.matrices.link_matrix]
        degrees = np.fromiter(degrees, dtype=int)
        to_explore = degrees.nonzero()[0]
        # 0th-index of to_explore must match the self.refseq constraints
        if self.refseq != 'all':
            if self.refseq[0] != '_':
                reduced = np.fromiter(
                    [to_explore[i] for i, x in enumerate(
                        self.sampling.poses_resnames[to_explore])
                     if self.refseq[0] == x], dtype=int)
            else:
                reduced = to_explore
        else:
            reduced = to_explore

        # restrict to best scored top_X solutions
        counts = 0
        passed = 0
        top_X_pool = defaultdict(list)
        self.interval = range(self.args.seq_min_size,
                              self.args.seq_max_size + 1)

        for replica in reduced:
            # start the DFS binary search of sequences
            top_X_pool, count, npass = self.grow_replica(replica, top_X_pool)
            counts += count
            passed += npass

        zones = self.sampling.zones
        leaders = []
        clust_array = np.full(len(zones), -1, dtype=int)
        return top_X_pool, counts, leaders, clust_array, passed

    # @profile
    def grow_replica(self, index, top_X_pool):
        one = ba('1')
        top_X = self.args.top_X
        scores = self.sampling.scores_all
        link_matrix = self.matrices.link_matrix
        clash_matrix = self.matrices.clash_matrix
        # initialize search variables
        paths = [self.matrices.link_matrix[index].copy()]
        current_seq = [index]
        clashes = [self.matrices.clash_matrix[index].copy()]
        score_seq = [scores[index]]
        counter = 0
        passed = 0

        if self.refseq != 'all':
            refseq = self.args.seq_to_search.copy()
            resnames = self.sampling.poses_resnames

        while paths:
            try:
                if self.refseq == 'all':
                    # 1. binary DFS exploration of the dinucleotide graph
                    current_idx = next(paths[-1].itersearch(one))
                    paths[-1][current_idx] = False

                else:
                    # 1. binary DFS exploration of the dinucleotide graph
                    tmp_idx = paths[-1].find(one)
                    if tmp_idx == -1:
                        del paths[-1]
                        del current_seq[-1]
                        del clashes[-1]
                        del score_seq[-1]
                        continue
                    else:
                        paths[-1][tmp_idx] = False
                        position = len(current_seq)
                        if refseq[position] == '_'\
                                or (refseq[position] == resnames[tmp_idx]):
                            current_idx = tmp_idx
                        else:
                            continue
            # 3. when the branch search is exhausted, stop DFS and:
            except Exception:
                # 3.3 close paths started from current node
                del paths[-1]
                del current_seq[-1]
                del clashes[-1]
                del score_seq[-1]
                continue

            # 2. grow if no clashes
            # update the nuclear_sequence
            current_seq.append(current_idx)
            score_seq.append(score_seq[-1] + scores[current_idx])
            seq_len = len(current_seq)
            # update feasible links of current_idx
            current_bit = link_matrix[current_idx]
            paths.append((current_bit ^ clashes[-1]) & current_bit)
            # erase from current expansion, the branches that lead to clashes
            clashes.append(clashes[-1] | clash_matrix[current_idx])
            if seq_len in self.interval:
                counter += 1
                seq_score = score_seq[-1] / seq_len
                # scorings = self.sampling.scores_all[current_seq]
                # if any(scorings >= 0):
                #     geometric = seq_score
                # else:
                #     geometric = -gm(-scorings)
                # check if pool is not full [free entry]
                if len(top_X_pool[seq_len]) < top_X:
                    heapq.heappush(
                        top_X_pool[seq_len],
                        (-seq_score, seq_len, current_seq.copy()))
                        # (-geometric, seq_len, current_seq.copy()))
                    passed += 1
                else:
                    if seq_score < -top_X_pool[seq_len][0][0]:
                        # if geometric < -top_X_pool[seq_len][0][0]:
                        heapq.heappushpop(
                            top_X_pool[seq_len],
                            (-seq_score, seq_len, current_seq.copy()))
                            # (-geometric, seq_len, current_seq.copy()))
                        passed += 1
                    else:
                        # pass
                        del paths[-1]
                        del current_seq[-1]
                        del clashes[-1]
                        del score_seq[-1]
        return top_X_pool, counter, passed

    def write_sequences(self, top_X_pool):
        # create SEQUENCES subdirs
        for subdir in top_X_pool:
            os.makedirs(os.path.join(self.args.out_seqdir, str(subdir)))

        # init constants
        for k, heap in enumerate(top_X_pool):
            prefix = len(str(self.args.top_X))
            pnames = self.sampling.O1PO2P_names

            outnames = []
            top_X = heapq.nlargest(self.args.top_X, top_X_pool[heap])

            for r, s in enumerate(top_X):
                name = os.path.join(self.args.out_seqdir,
                                    str(heap),
                                    'seq-{}nt-{}.pdb'
                                    .format(heap, str(r + 1).zfill(prefix)))
                sequence_info = top_X[r]
                remark = 'REMARK  {:>6} {:>8.5f}  {:>4}  {:<80}\n'\
                    .format(r, -sequence_info[0], sequence_info[1],
                            str(sequence_info[2]))
                outnames.append((name, sequence_info, remark))

            string = '{:<6}{:>5}{:>1}{:>4}{:>1}{:<4}{:>1}{:>4}{:>1}{:>3}' \
                     '{:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}{:>6}{:<4}' \
                     '{:>2}{:>2}\n'
            # write one nuclear_sequence per pdb file
            for name, sequence_info, remark in outnames:
                with open(name, 'wt') as outseq:
                    num = it.count(1)
                    outseq.write(remark)
                    for i, nucl in enumerate(sequence_info[2], 1):
                        coords = self.sampling.poses_all[nucl].data
                        names = self.sampling.names_dict[self.sampling
                                                         .poses_str_all[nucl]]
                        score = self.sampling.scores_all[nucl]
                        resname = self.sampling.poses_resnames[nucl]
                        for idx, coor in enumerate(coords):
                            outseq.write(string.format(
                                'ATOM',
                                next(num),
                                '',
                                names[idx],
                                '',
                                resname,
                                '',
                                i,
                                '',
                                '',
                                coor[0],
                                coor[1],
                                coor[2],
                                score,
                                score,
                                '',
                                'SEQ',
                                names[idx][0],
                                ''
                                ))

                        for c, coorp in enumerate(
                                self.sampling.O1PO2P_coords[nucl]):
                            outseq.write(string.format(
                                'ATOM',
                                next(num),
                                '',
                                pnames[c],
                                '',
                                resname,
                                '',
                                i,
                                '',
                                '',
                                coorp[0],
                                coorp[1],
                                coorp[2],
                                score,
                                score,
                                '',
                                'SEQ',
                                names[idx][0],
                                ''
                                ))

    def write_log(self, top_X_pool):
        for i, heap in enumerate(top_X_pool):
            outfile = os.path.join(self.args.out_topXdir,
                                   'topX_sequences_{}nt.log'.format(heap))
            top_X = heapq.nlargest(len(top_X_pool[heap]), top_X_pool[heap])
            maxim = len(str(len(top_X)))
            outdf = pd.DataFrame(columns=['TopSeq', 'SeqScore', 'SeqSize',
                                          'SeqNucs', 'NucScore', 'TopNuc'])
            outdf['TopSeq'] = range(1, len(top_X) + 1)
            outdf['SeqScore'] = [-x[0] for x in top_X]
            outdf['SeqSize'] = [x[1] for x in top_X]
            resnames = []
            scores = []
            nucranks = []
            ncontacts = []
            distances = []

            for tup in top_X:
                seqnames = []
                seqscores = []
                seqtops = []
                ncont = []
                dists = []
                for idx, num in enumerate(tup[2]):
                    seqnames.append(self.sampling.poses_resnames[num])
                    seqscores.append(self.sampling.scores_all[num])
                    seqtops.append(self.sampling.babel[num])
                    ncont.append(len(self.sampling.zones[num]))
                    if idx != 0:
                        x = tup[2][idx]
                        y = tup[2][idx - 1]
                        str_x = self.sampling.poses_str_all[x]
                        str_y = self.sampling.poses_str_all[y]
                        char1 = self.sampling.names_dict[str_x]
                        char2 = self.sampling.names_dict[str_y]
                        c5_index = int(np.where(char1.find("C5'") == 0)[0])
                        o3_index = int(np.where(char2.find("O3'") == 0)[0])
                        c5_coords = self.sampling.poses_all[x].data[c5_index]
                        o3_coords = self.sampling.poses_all[y].data[o3_index]
                        dist = np.linalg.norm(c5_coords - o3_coords)
                        assert dist <= self.args.max_dist_O3C5
                        dists.append('{:3.2f}'.format(round(dist, 2)))
                resnames.append(seqnames)
                scores.append(seqscores)
                nucranks.append(seqtops)
                ncontacts.append(sum(ncont))
                distances.append(dists)

            outdf['SeqNucs'] = [' '.join(x) for x in resnames]
            outdf['NucScore'] = [' '.join(['{:>6.2f}'.format(round(y, 2))
                                           for y in x]) for x in scores]
            outdf['TopNuc'] = [' '.join(['{:>{}}'.format(str(y), maxim)
                                         for y in x]) for x in nucranks]
            outdf['nContacts'] = ncontacts
            outdf['InterDists'] = [' '.join(x) for x in distances]
            with open(outfile, 'wt') as out:
                outdf.to_string(out, index=False)

    def write_log_popularity(self, kmers_dict, top_X_pool):
        ksize = []
        totalCounts = []
        distinctCounts = []
        uniqueCounts = []

        for size in kmers_dict:
            outfile = os.path.join(self.args.out_popdir,
                                   'seq_popularity_{}nt.log'.format(size))
            top_X = kmers_dict[size]

            ksize.append(size)
            totalCounts.append(len(top_X))
            distinctCounts.append(len(set(top_X)))
            uniqueCounts.append(
                (np.asarray(list(Counter(top_X).values())) == 1).sum())

            maxim = len(str(len(top_X)))
            outdf = pd.DataFrame()
            counter = Counter(top_X).most_common()

            outdf['TopSeq'] = range(1, len(counter) + 1)

            seqcounts = [tup[1] for tup in counter]
            outdf['SeqCounts'] = seqcounts

            seqscores = [np.average(
                self.sampling.scores_all[list(map(int, string[0].split()))])
                for string in counter]
            outdf['SeqScore'] = seqscores

            TopN_seqs = set([' '.join(map(str, x[2]))
                             for x in top_X_pool[size]])
            isin_TopN = ['yes' if (x[0] in TopN_seqs) else 'no'
                         for x in counter]
            outdf['inTopN'] = isin_TopN

            maxval = counter[0][1]
            popularities = ['{:>2.4f}'.format(tup[1] / maxval)
                            for tup in counter]
            outdf['SeqPopularity'] = popularities

            outdf['SeqSize'] = [len(x[0].split()) for x in counter]

            resnames = []
            scores = []
            nucranks = []
            for string_ in counter:
                seqnames = []
                seqscores = []
                seqtops = []
                for num in map(int, string_[0].split()):
                    seqnames.append(self.sampling.poses_resnames[num])
                    seqscores.append(self.sampling.scores_all[num])
                    seqtops.append(self.sampling.babel[num])
                resnames.append(seqnames)
                scores.append(seqscores)
                nucranks.append(seqtops)
            outdf['SeqNucs'] = [' '.join(x) for x in resnames]
            outdf['NucScore'] = [' '.join(['{:>6.2f}'.format(round(y, 2))
                                           for y in x]) for x in scores]
            outdf['TopNuc'] = [' '.join(['{:>{}}'.format(str(y), maxim)
                                         for y in x]) for x in nucranks]
            with open(outfile, 'wt') as out:
                outdf.to_string(out, index=False)

        outfile2 = os.path.join(self.args.out_popdir, 'seq_popularity_all.log')
        outdf2 = pd.DataFrame()
        outdf2['ksize'] = ksize
        outdf2['tot_counts'] = totalCounts
        outdf2['distinct_counts'] = distinctCounts
        outdf2['unique_counts'] = uniqueCounts
        with open(outfile2, 'wt') as out:
            outdf2.to_string(out, index=False)


# =============================================================================
# Demonstration
# =============================================================================
# from confread import parser as config_parser
# from sampling import parser as sampling_parser
# from bitmatrix import builder
#
# args = config_parser('/home/roy.gonzalez-aleman/rprojects/nuclear/examples/'
#                     '1.0.0/2xnr_0_6_preliminar.cfg')
# sampling = sampling_parser(args)
# matrices = builder(args, sampling)
# matrices.get_clash_matrix()
# matrices.get_link_matrix()
# self = generator(args, sampling, matrices)
