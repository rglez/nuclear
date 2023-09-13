#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  9 14:01:12 2021

@author: roy.gonzalez-aleman
"""
import heapq
import itertools as it
import os
import pickle
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import plotly.express as px


class spotify:
    def __init__(self, args, sampling):
        self.args = args
        self.sampling = sampling
        self.fingerprints = np.asarray(self.sampling.zones_int)
        self.scores = self.sampling.scores_all

    def get_supersets(self):
        # order fingerprints by ascending score
        # ordered_idx = (-self.scores).argsort().tolist()
        sizes = np.fromiter([len(x) for x in self.fingerprints], dtype=int)
        ordered_idx2 = (sizes).argsort().tolist()
        # unclustered = set(ordered_idx)
        unclustered = set(ordered_idx2)

        # #### resid-level subset clustering
        # first pass detecting leaders
        clusters = []
        # while ordered_idx:
        while ordered_idx2:
            # lowest_score = ordered_idx.pop()
            biggest_size = ordered_idx2.pop()
            # if lowest_score in unclustered:
            if biggest_size in unclustered:
                # subset = [lowest_score]
                subset = [biggest_size]
                # unclustered.remove(lowest_score)
                unclustered.remove(biggest_size)
                # reference = self.fingerprints[lowest_score]
                reference = self.fingerprints[biggest_size]
                for target in unclustered:
                    if self.fingerprints[target].issubset(reference):
                        subset.append(target)
                clusters.append(subset)
                print('Identified Cluster-{} with {} elements.'.format(len(clusters), len(subset)))
                unclustered = unclustered.difference(set(subset))
            else:
                continue
        # reorganize by score
        # get densities
        densities = []
        minima = []
        for cluster in clusters:
            seed = cluster[0]
            fingerlenght = len(self.fingerprints[seed])
            clustersize = len(cluster)
            densities.append(clustersize / fingerlenght)
            minima.append(self.scores[cluster].min())
        max_den = max(densities)
        min_den = min(densities)
        maxmin = max_den - min_den
        normalized_densities = [(x - min_den) / maxmin for x in densities]
        self.densities = normalized_densities
        self.clusters = clusters
        self.minima = minima

    def find_spots(self):
        iterable = []
        iterable_densities = []
        iterable_minima = []
        for i, cluster in enumerate(self.clusters):
            cluster_density = self.densities[i]
            if cluster_density > self.args.density_cut:
                iterable.append(cluster)
                iterable_densities.append(cluster_density)
                iterable_minima.append(self.minima[i])
        assert len(iterable) > 0, 'The density_cut specified is too high to produce hotspots.'

        # reorganize by score
        order = np.argsort(iterable_minima)
        self.iterable = []
        self.iterable_densities = []
        self.iterable_minima = []
        for idx in order:
            self.iterable.append(iterable[idx])
            self.iterable_densities.append(iterable_densities[idx])
            self.iterable_minima.append(iterable_minima[idx])

        # merge to best_ranked using Tanimoto
        leaders = [self.fingerprints[x][0] for x in self.iterable]
        dejavu = set()
        for acceptor, leader1 in enumerate(leaders):
            for donor, leader2 in enumerate(leaders):
                if acceptor < donor:
                    tanimoto = get_tanimoto(leader1, leader2)
                    if (tanimoto <= self.args.merge_cut) and (donor not in dejavu):
                        self.iterable[acceptor].extend(self.iterable[donor])
                        dejavu.add(donor)
                        print('{} <-- {}'.format(acceptor, donor))
        [self.iterable.pop(x) for x in sorted(list(dejavu))[::-1]]
        [self.iterable_densities.pop(x) for x in sorted(list(dejavu))[::-1]]
        [self.iterable_minima.pop(x) for x in sorted(list(dejavu))[::-1]]

    def get_crdspots(self):
        rank_dict = defaultdict(list)
        spots = self.iterable
        for spot in spots:
            crds = self.sampling.poses_str_all[spot]
            ranks = self.sampling.babel[spot]
            for i, crd in enumerate(crds):
                rank_dict[crd].append(ranks[i])
        rank_dict2 = {}
        for crd in rank_dict:
            crd_path = os.path.join(self.args.input_dir, crd)
            rank_dict2.update({crd_path: np.asarray(rank_dict[crd])})

        file_name = os.path.join(self.args.output_dir, "crdspots.pick")
        with open(file_name, 'wb') as file:
            pickle.dump(rank_dict2, file)

    def get_dataframes(self):
        codes = sorted([get_code(x) for x in set(self.sampling.poses_str_all)])
        sites = list(range(len(self.iterable)))
        self.df = pd.DataFrame(columns=sites, index=codes, dtype=float)
        self.counts = self.df.copy()

        for i, cluster in enumerate(self.iterable):
            strings = self.sampling.poses_str_all[cluster]
            ranking = self.sampling.babel[cluster]

            for x, y in Counter(strings).items():
                self.counts.loc[get_code(x), i] = y

            queue = list(zip(strings, ranking))
            heapq.heapify(queue)
            dejavu = set()
            while queue:
                name, rank = heapq.heappop(queue)
                if name not in dejavu:
                    self.df.loc[get_code(name), i] = rank
                    dejavu.add(name)
        self.best_rank = self.df.copy()

    def get_overall_best_rank(self):
        df = self.best_rank[self.best_rank <= 10000]
        fig1 = px.imshow(df,
                         color_continuous_scale='Rainbow_r',
                         title='Best Top-X Ranked Group per Spot',
                         labels={'x': 'Spot-ID', 'y': 'Group', 'color': 'BestRank'})
        outname = os.path.join(self.args.output_dir, "overall_best_ranks.html")
        fig1.write_html(outname)
        df.to_csv(os.path.join(self.args.TSV, "overall_best_ranks.tsv"), sep='\t')

    def get_overall_counts(self):
        fig2 = px.imshow(self.counts, color_continuous_scale='Rainbow_r',
                         title='Count of Groups per Spot',
                         labels={'x': 'Spot-ID', 'y': 'Group', 'color': 'Count'})
        outname = os.path.join(self.args.output_dir, "overall_counts.html")
        fig2.write_html(outname)
        self.counts.to_csv(os.path.join(self.args.TSV, "overall_counts.tsv"), sep='\t')

    def get_countingby_tops(self):
        self.tops = [1000, 500, 100, 50, 10, 5]
        self.spotview = []
        self.groupview = []
        for top in self.tops:
            top_df = self.df.copy()
            top_df[:] = np.nan
            for i, cluster in enumerate(self.iterable):
                strings = self.sampling.poses_str_all[cluster]
                ranking = self.sampling.babel[cluster]
                uniques = np.unique(strings)
                dictRank = {}
                for unique in uniques:
                    ranks = ranking[strings == unique]
                    dictRank.update({get_code(unique): ranks})
                for code in dictRank:
                    top_df.loc[code, i] = (dictRank[code] < top).sum()
            # spotview
            self.spotview.append((top_df > 0).sum(axis=0))
            # groupview
            self.groupview.append((top_df > 0).sum(axis=1))
            # wholeview
            fig3 = px.imshow(top_df, color_continuous_scale='Rainbow_r',
                             title='Count of Groups under Top-{}'.format(top),
                             labels={'x': 'Spot-ID', 'y': 'Group', 'color': 'Count'})
            outname = os.path.join(self.args.output_dir, "counting_by_top_{}.html".format(top))
            fig3.write_html(outname)
            top_df.to_csv(os.path.join(self.args.TSV, "counting_by_top_{}.tsv".format(top)), sep='\t')

    def get_axis_views(self):
        gv = pd.DataFrame([tuple(x) for x in self.groupview],
                          index=['top-{}'.format(x) for x in self.tops])
        fig4 = px.bar(gv.T, barmode='overlay', opacity=1,
                      title='Distribution of spots per group and Top-X',
                      labels={'x': 'Group-ID', 'y': 'No. of Spots'})
        outname = os.path.join(self.args.output_dir, "overall_groupview.html")
        fig4.write_html(outname)
        gv.T.to_csv(os.path.join(self.args.TSV, "overall_groupview.tsv"), sep='\t')

        sv = pd.DataFrame(self.spotview,
                          index=['top-{}'.format(x) for x in self.tops])
        fig5 = px.bar(sv.T, barmode='overlay', opacity=1,
                      title='Distribution of groups per spot and Top-X',
                      labels={'x': 'Spot-ID', 'y': 'No. of Groups'})
        outname = os.path.join(self.args.output_dir, "overall_spotview.html")
        fig5.write_html(outname)
        sv.T.to_csv(os.path.join(self.args.TSV, "overall_spotview.tsv"), sep='\t')

    def get_fingerlog(self, tops=[5, 10, 50, 100, 500, 1000]):
        fingers = [self.fingerprints[itb] for itb in self.iterable]
        self.residues = ['residue ' + ' '.join([str(x) for x in sorted(list(frozenset.union(*finger)))]) for finger in fingers]

        spots_3D = pd.DataFrame()
        spots_3D['ID'] = range(len(self.iterable))
        spots_3D['nGroupsTotal'] = (self.counts > 0).sum(axis=0)
        spots_3D['PopTotal'] = self.counts.sum(axis=0)
        spots_3D['BestRank'] = self.best_rank.min(axis=0)
        spots_3D['BestScore'] = self.iterable_minima
        spots_3D['Density'] = self.iterable_densities
        spots_3D['Residues'] = self.residues
        fig3d1 = px.scatter_3d(spots_3D, x='nGroupsTotal',
                               y='PopTotal',
                               z='BestRank', color='Density',
                               hover_name=spots_3D.ID)
        fig3d1.write_html(os.path.join(self.args.output_dir, 'spots_3D_Total.html'))
        for top in tops:
            topN_idx = [(self.sampling.babel[x] < top) for i, x in enumerate(self.iterable)]
            spots_3D['PopTop-{}'.format(top)] = [x.sum() for x in topN_idx]
            spots_3D['nGroupsTop-{}'.format(top)] = [
                np.unique(self.sampling.poses_str_all[x][topN_idx[i]]).size
                for i, x in enumerate(self.iterable)]
            # 3D graphics
            fig3d2 = px.scatter_3d(spots_3D, x='nGroupsTop-{}'.format(top),
                                   y='PopTop-{}'.format(top),
                                   z='BestRank', color='Density',
                                   hover_name=spots_3D.ID)
            fig3d2.write_html(os.path.join(self.args.output_dir, 'spots_3D_Top-{}.html'.format(top)))

        # fingerprints file
        outname = 'fingerprints'
        spots_3D.to_csv(os.path.join(self.args.TSV, "{}.tsv".format(outname)), sep='\t')
        spots_3D.to_string(os.path.join(self.args.output_dir, "{}.log".format(outname)),
                           index=False, justify='center')

    def get_fingerlog2(self, tops=[5, 10, 50, 100, 500, 1000]):
        fingers = [self.fingerprints[itb] for itb in self.iterable]
        indices = [frozenset.union(*finger) for finger in fingers]
        pp = self.args.prot_parsed
        babelDict = dict(set(zip(pp.getResindices(), pp.getResnums())))
        pdbnums = ['resid ' + ' '.join([str(f) for f in sorted([babelDict[y] for y in list(x)])]) for x in indices]

        spots_3D = pd.DataFrame()
        spots_3D['ID'] = range(len(self.iterable))
        spots_3D['nGroupsTotal'] = (self.counts > 0).sum(axis=0)
        spots_3D['PopTotal'] = self.counts.sum(axis=0)
        spots_3D['BestRank'] = self.best_rank.min(axis=0)
        spots_3D['BestScore'] = self.iterable_minima
        spots_3D['Density'] = self.iterable_densities
        spots_3D['Residues'] = pdbnums
        for top in tops:
            topN_idx = [(self.sampling.babel[x] < top) for i, x in enumerate(self.iterable)]
            spots_3D['PopTop-{}'.format(top)] = [x.sum() for x in topN_idx]
            spots_3D['nGroupsTop-{}'.format(top)] = [
                np.unique(self.sampling.poses_str_all[x][topN_idx[i]]).size
                for i, x in enumerate(self.iterable)]
        # fingerprints file
        outname = 'fingerprints2'
        spots_3D.to_csv(os.path.join(self.args.TSV, "{}.tsv".format(outname)), sep='\t')
        spots_3D.to_string(os.path.join(self.args.output_dir, "{}.log".format(outname)),
                           index=False, justify='center')

    def create_reps(self):
        colors = it.cycle(range(32))
        outname = os.path.join(self.args.output_dir, "representations2.vmd")
        with open(outname, 'wt') as ln:
            ln.write('mol new {}\n'.format(self.args.prot_coords))
            ln.write('animate style Loop\n')
            ln.write('menu graphics off\n')
            ln.write('menu graphics on\n')
            ln.write('mol modstyle 0 0 NewCartoon 0.300000 10.000000 4.100000 0\n')
            ln.write('mol color Name\n')
            for i, spot in enumerate(self.residues, 1):
                ln.write('mol addrep 0\n')
                ln.write('mol modselect {} 0 {}\n'.format(i, spot))
                ln.write('mol modcolor {} 0 ColorID {}\n'.format(i, next(colors)))
                ln.write('mol representation Lines 2.0\n')
            for i, rep in enumerate(self.residues, 1):
                ln.write('mol showrep 0 {} 0\n'.format(i))


def get_tanimoto(A, B):
    intersection = len(A & B)
    if not intersection:
        return 1
    return 1 - (intersection / len(A | B))


def get_code(string):
    return '-'.join(string.split('_')[2:])


def get_outliers_1X(x_arr, y_arr):

    # Average magnitudes ------------------------------------------------------
    y_ave = y_arr.mean()

    # flores-garza candidates ordered by gamma --------------------------------
    # candidates_idx = np.where(y_arr > y_ave)[0]
    candidates_idx = range(y_arr.size)
    candidates = y_arr[candidates_idx]
    order = (-candidates).argsort()
    ordered = candidates[order]

    # flores-garza distances and retrieval by last gap procedure --------------
    d_i = abs(np.diff(ordered))
    if d_i.size == 0:
        raise ValueError('\n\n>>> Algorithm Inconsistency\nThe autodetection'
                         ' of cluster centers can not proceed for the'
                         ' specified argument values. Please set -auto_centers'
                         ' to False')
    d_ave = d_i.mean()
    nnodes = 2
    nodes_by_level = []
    nodes_by_index = []
    while nnodes > 1:
        center_gaps = order[np.where(d_i > d_ave)[0]]
        if any(center_gaps):
            nnodes = center_gaps.size
            last_center_gap = center_gaps[-1]
            last_center_idx = order.tolist().index(last_center_gap)
            nodes = order[:last_center_idx + 1]
            nodes_by_level.append(candidates[nodes])
            nodes_by_index.append(nodes)
            d_ave = d_i[np.where(d_i > d_ave)[0]].mean()
        else:
            break
    return nodes_by_index, nodes_by_level


# from confread import parser as config_parser
# from sampling import parser as sampling_parser
# config = '/home/roy.gonzalez-aleman/Dr.Alzheimer/03-period/thesis/nuclear_jobs/hotspots_search/5MCQ_X_hotspots_std.cfg'
# args = config_parser(config)
# sampling = sampling_parser(args)
# self = spotify(args, sampling)

# # retrieve clusters
# self.get_supersets()
# self.find_spots()

# self.get_dataframes()
# self.get_overall_best_rank()
# self.get_overall_counts()
# self.get_countingby_tops()
# self.get_axis_views()
# self.get_fingerlog()
# self.create_reps()
