#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  9 17:27:59 2021

@author: roy.gonzalez-aleman
"""

import itertools as it
import os
import sys

import numpy as np
import scipy.spatial.ckdtree as ckd
from parmed.charmm import CharmmCrdFile

root_modules = os.sep.join(__file__.split(os.sep)[:-3])
sys.path.insert(0, root_modules)
# from central.formats.crd import crd


def get_tanimoto(A, B):
    intersection = len(A & B)
    if not intersection:
        return 0
    return intersection / len(A | B)


def clustering_orthogonal(zones, cutoff):
    members = range(len(zones))
    available = set(members)
    ordered = list(members[::-1])
    clust_array = np.full(len(zones), -1, dtype=int)
    id_clust = it.count(0)
    leaders = []
    while available:
        # take the lowest score as seed
        minim = ordered.pop()
        if minim not in available:
            continue
        available.remove(minim)
        # add all similar frames to seed's cluster
        minim_cluster = [minim]
        leaders.append(minim)
        zone_seed = zones[minim]
        for index in available:
            if get_tanimoto(zone_seed, zones[index]) > cutoff:
                minim_cluster.append(index)
        # update
        clust_array[minim_cluster] = next(id_clust)
        available -= set(minim_cluster)
    return leaders, clust_array


def calc_rmsd(A, B, k):
    return np.sqrt(((A - B) ** 2).sum() * k)


def reduce_poses_by_rmsd(poses, cutoff, k):
    members = range(len(poses))
    available = set(members)
    ordered = list(members[::-1])
    clust_array = np.full(len(poses), -1, dtype=int)
    id_clust = it.count(0)
    leaders = []
    while available:
        # take the lowest score as seed
        minim = ordered.pop()
        if minim not in available:
            continue
        available.remove(minim)
        # add all similar frames to seed's cluster
        minim_cluster = [minim]
        leaders.append(minim)
        pose_seed = poses[minim]
        for index in available:
            # if get_tanimoto(pose_seed, poses[index]) > cutoff:
            if calc_rmsd(pose_seed, poses[index], k) <= cutoff:
                minim_cluster.append(index)
        # update
        clust_array[minim_cluster] = next(id_clust)
        available -= set(minim_cluster)
    return leaders, clust_array


def get_mcss_minimal_box(kdtrees, outdir):
    coordinates = np.concatenate([x.data for x in kdtrees])
    xmin, ymin, zmin = coordinates.min(axis=0)
    xmax, ymax, zmax = coordinates.max(axis=0)
    dx = xmax - xmin
    dy = ymax - ymin
    dz = zmax - zmin
    outname = os.path.join(outdir, 'mcss_box.log')
    with open(outname, 'wt') as mcss:
        mcss.write('Box size in x, y, z is: {:3.2f} {:3.2f} {:3.2f}\n'.format(dx, dy, dz))
        mcss.write('{:3.2f}\n{:3.2f}\n{:3.2f}\n{:3.2f}\n{:3.2f}\n{:3.2f}\n'.format(xmin, ymin, zmin, xmax, ymax, zmax))

def get_path_explored(prot_coords, zones_int, outdir):
    outname = os.path.join(outdir, 'path_explored.vmd')
    with open(outname, 'wt') as pts:
        pts.write('mol new {}\n'.format(os.path.abspath(prot_coords)))
        pts.write('animate style Loop\n')
        pts.write('menu graphics off\n')
        pts.write('menu graphics on\n')
        pts.write('mol modstyle 0 0 NewCartoon 0.300000 10.000000 4.100000 0\n')
        # pts.write('mol color Name\n')
        pts.write('mol modcolor 0 0 ColorID 0\n')

        pts.write('mol addrep 0\n')
        union = set()
        [union.update(x) for x in zones_int]
        pts.write('mol modselect 1 0 residue {}\n'.format(' '.join([str(x) for x in union])))
        pts.write('mol modcolor 1 0 ColorID 1\n')
        pts.write('mol representation Lines 2.0\n')


class parser2:
    def __init__(self, args):
        self.args = args
        self.prot_atom_names = np.char.asarray(self.args.prot_parsed.getNames())
        self.prot_pose = self.args.prot_parsed.getCoords()
        self.prot_kdtree = ckd.cKDTree(self.prot_pose)

        # containers of the all-groups sampling
        names_dict = dict()
        poses_all = []
        scores_all = []
        poses_str_all = []
        poses_resnames_all = []
        babel = []
        C5_coords = []
        O3_coords = []
        poxygens_coords = []
        zones = []
        zones_int_red = []
        zones_str_red = []

        N = len(self.args.crd_list)
        for i, (file, maxscore, nposes) in enumerate(sorted(self.args.crd_list), 1):
            # parse the whole crd_file
            print('Parsing {} from {}: {}'.format(i, N, file))

            crd_object = crd(os.path.join(self.args.input_dir, file),
                             self.args.frag_sel, top=self.args.top)
            crd_object.parse_crd_lines()
            sel_atoms = crd_object.sel_atoms
            frag_resname = crd_object.resname
            scores = crd_object.weights
            poses = crd_object.poses
            oxygens_coords = crd_object.oxygens_coords
            oxygens_names = crd_object.oxygens_names

            frag_name = os.path.basename(file).split('.')[0]
            names_dict.update({frag_name: sel_atoms})

            # get all contact zones and their size
            zones_str, zones_int, zones_both = self.get_zones(
                poses, sel_atoms)
            zones_sizes = [len(x) for x in zones_int]

            if 'search_params' in args.cfg.sections():
                # path_to_search reduction
                sele_path = set(
                    [i for i, x in enumerate(zones_int)
                     if len(frozenset.intersection(self.args.searchpath, x)) > 0])
            else:
                sele_path = set([i for i, x in enumerate(zones_int)])

            # percentile reduction
            num_contacts = 0
            sele_size = set([i for i, x in enumerate(zones_sizes) if x > num_contacts])

            # combining reductions
            sele_idxs = sorted(list(set.intersection(sele_path, sele_size)))

            # # Np reduction
            # if (nposes == 'N') and (maxscore == 'N'):
            #     Np = len(poses)
            # elif maxscore != 'N':
            #     Np = np.nonzero(scores >= float(maxscore))[0][0]
            # elif nposes != 'N':
            #     Np = int(nposes)
            # reduced_idx = [x for x in sele_idxs if x < Np]
            reduced_idx = [x for x in sele_idxs]

            # reducing the dimension of all objects
            sele_scores = scores[reduced_idx]
            sele_poses = [poses[x] for x in reduced_idx]
            if oxygens_coords.size > 0:
                sele_oxygens_coords = oxygens_coords[reduced_idx]

            sele_zones_int = [zones_int[x] for x in reduced_idx]
            sele_zones_str = [zones_str[x] for x in reduced_idx]
            sele_zones_both = [zones_both[x] for x in reduced_idx]

            # updating containers
            M = len(reduced_idx)
            scores_all.append(sele_scores)
            poses_all.extend(sele_poses)
            poses_str_all.extend(it.repeat(frag_name, M))
            poses_resnames_all.extend(it.repeat(frag_resname, M))
            babel.extend(reduced_idx)
            zones.extend(sele_zones_both)
            zones_int_red.extend(sele_zones_int)
            zones_str_red.extend(sele_zones_str)
            if len(oxygens_coords) > 0:
                poxygens_coords.extend(sele_oxygens_coords)
            for pose in sele_poses:
                C5_coords.extend(pose[np.where(sel_atoms == "C5'")[0]])
                O3_coords.extend(pose[np.where(sel_atoms == "O3'")[0]])

        # declaring class attributes
        if (C5_coords) and (O3_coords):
            self.C5_tree = ckd.cKDTree(C5_coords)
            self.O3_tree = ckd.cKDTree(O3_coords)
        else:
            print('No C5 or O3 coordinates for this job')

        if 'search_params' in args.cfg.sections():
            self.names_dict = names_dict
            self.poses_all = [ckd.cKDTree(pose) for pose in poses_all]
            self.poses_resnames = np.asarray(poses_resnames_all)
            self.poses_idx = np.fromiter([i for i, pose in enumerate(poses_all) for y in pose], dtype=int)
            self.O1PO2P_coords = poxygens_coords
            self.O1PO2P_names = oxygens_names

        self.babel = np.fromiter(babel, dtype=int)
        self.poses_str_all = np.asarray(poses_str_all)
        self.scores_all = np.concatenate(scores_all)
        self.zones = zones
        self.zones_int = zones_int_red
        self.zones_str = zones_str_red

        # get the mcss_box
        get_mcss_minimal_box([ckd.cKDTree(pose) for pose in poses_all], self.args.output_dir)
        # get the check of path_to_search
        get_path_explored(self.args.prot_coords, self.zones_int, self.args.output_dir)

        C5_coords = None
        O3_coords = None
        poses_all = None

    def get_zones(self, frag_poses, sel_atoms):
        cut = self.args.inter_dist
        zones_int = []
        zones_str = []
        zones_both = []
        for pose in frag_poses:
            dist, neighbors = self.prot_kdtree.query(pose.data, k=1,
                                                     distance_upper_bound=cut)
            close = neighbors[dist <= cut]
            frag_contact_names = sel_atoms[dist <= cut]
            prot_contact_names = self.prot_atom_names[close]
            prot_contact_residx = self.args.prot_parsed.getResindices()[close]
            # simplest zones nomenclature only with protein residues
            zones_int.append(frozenset(prot_contact_residx))
            # zones nomenclature 'centered' in the protein
            zone_str = []
            for i, dd in enumerate(prot_contact_residx):
                zone_str.append('{}-{}'.format(dd, prot_contact_names[i]))
            zones_str.append(frozenset(zone_str))
            # zones nomenclature 'englobing' protein and fragment
            zone_both = []
            for i, dd in enumerate(prot_contact_residx):
                zone_both.append('{}-{}_{}'.format(dd, prot_contact_names[i],
                                                   frag_contact_names[i]))
            zones_both.append(frozenset(zone_both))
        return zones_str, zones_int, zones_both


class parser:
    def __init__(self, args):
        self.args = args
        self.prot_atom_names = np.char.asarray(self.args.prot_parsed.getNames())
        self.prot_pose = self.args.prot_parsed.getCoords()
        self.prot_kdtree = ckd.cKDTree(self.prot_pose)

        # containers of the all-groups sampling
        names_dict = dict()
        poses_all = []
        scores_all = []
        poses_str_all = []
        poses_resnames_all = []
        babel = []
        C5_coords = []
        O3_coords = []
        poxygens_coords = []
        zones = []
        zones_int_red = []
        zones_str_red = []

        N = len(self.args.crd_list)
        for i, (file, maxscore, nposes) in enumerate(sorted(self.args.crd_list), 1):
            # parse the whole crd_file
            print('\nParsing {} from {}: {}'.format(i, N, file))
            sel_atoms, frag_resname, scores, poses, oxygens_coords,\
                oxygens_names = self.parse_crd_file(file)
            frag_name = os.path.basename(file).split('.')[0]
            names_dict.update({frag_name: sel_atoms})

            # reduce poses by rmsd orthogonal pre-clustering
            if self.args.pre_clustering:
                print('Reducing geometrical redundancy ...')
                k = (1.0 / poses[0].shape[0])
                seeds, clust_array = reduce_poses_by_rmsd(
                    poses, self.args.rmsd_cut, k)
                poses_reduced = [poses[x] for i, x in enumerate(seeds)]
                print('from {} poses to {} using a {}A cutoff\n'
                      .format(len(poses), len(poses_reduced),
                              self.args.rmsd_cut))
            else:
                seeds = range(len(poses))
                poses_reduced = [poses[x] for i, x in enumerate(seeds)]

            # get all contact zones and their size
            zones_str, zones_int, zones_both = self.get_zones(poses_reduced,
                                                              sel_atoms)
            zones_sizes = [len(x) for x in zones_int]

            if 'search_params' in args.cfg.sections():
                # path_to_search reduction
                sele_path = set(
                    [i for i, x in enumerate(seeds)
                      if len(frozenset.intersection(
                            self.args.searchpath, zones_int[i])) > 0])
            else:
                sele_path = set(range(len(seeds)))

            # percentile reduction
            num_contacts = 0
            sele_size = set([i for i, x in enumerate(seeds)
                              if zones_sizes[i] > num_contacts])

            # combining reductions
            sele_idxs = sorted(list(
                set.intersection(sele_path, sele_size)))

            # Np reduction
            if (nposes == 'N') and (maxscore == 'N'):
                Np = len(poses)
            elif maxscore != 'N':
                Np = np.nonzero(scores >= float(maxscore))[0][0]
            elif nposes != 'N':
                Np = int(nposes)
            reduced_idx = [seeds[x] for x in sele_idxs if seeds[x] <= Np]

            # reducing the dimension of all objects
            sele_scores = scores[reduced_idx]
            sele_poses = [poses[x] for x in reduced_idx]
            if oxygens_coords.size > 0:
                sele_oxygens_coords = oxygens_coords[reduced_idx]

            sele_zones_int = [zones_int[i] for i, x in enumerate(reduced_idx)]
            sele_zones_str = [zones_str[i] for i, x in enumerate(reduced_idx)]
            sele_zones_both = [zones_both[i] for i, x in enumerate(reduced_idx)]

            # updating containers
            M = len(reduced_idx)
            scores_all.append(sele_scores)
            poses_all.extend(sele_poses)
            poses_str_all.extend(it.repeat(frag_name, M))
            poses_resnames_all.extend(it.repeat(frag_resname, M))
            babel.extend(reduced_idx)
            zones.extend(sele_zones_both)
            zones_int_red.extend(sele_zones_int)
            zones_str_red.extend(sele_zones_str)
            if len(oxygens_coords) > 0:
                poxygens_coords.extend(sele_oxygens_coords)
            for pose in sele_poses:
                C5_coords.extend(pose[np.where(sel_atoms == "C5'")[0]])
                O3_coords.extend(pose[np.where(sel_atoms == "O3'")[0]])

        # declaring class attributes
        if (C5_coords) and (O3_coords):
            self.C5_tree = ckd.cKDTree(C5_coords)
            self.O3_tree = ckd.cKDTree(O3_coords)
        else:
            print('No C5 or O3 coordinates for this job')

        if 'search_params' in args.cfg.sections():
            self.names_dict = names_dict
            self.poses_all = [ckd.cKDTree(pose) for pose in poses_all]
            self.poses_resnames = np.asarray(poses_resnames_all)
            self.poses_idx = np.fromiter([i for i, pose in enumerate(poses_all) for y in pose], dtype=int)
            self.O1PO2P_coords = poxygens_coords
            self.O1PO2P_names = oxygens_names

        self.babel = np.fromiter(babel, dtype=int)
        self.poses_str_all = np.asarray(poses_str_all)
        self.scores_all = np.concatenate(scores_all)
        self.zones = zones
        self.zones_int = zones_int_red
        self.zones_str = zones_str_red

        # get the mcss_box
        get_mcss_minimal_box([ckd.cKDTree(pose) for pose in poses_all], self.args.output_dir)
        # get the check of path_to_search
        get_path_explored(self.args.prot_coords, self.zones_int, self.args.output_dir)

        C5_coords = None
        O3_coords = None
        poses_all = None

    def parse_crd_file(self, crd_file):
        frags_raw = CharmmCrdFile(os.path.join(self.args.input_dir, crd_file))
        frags_nframes = len(set(frags_raw.resno))
        frags_natoms = int(frags_raw.natom / frags_nframes)
        frag_resname = frags_raw.resname[0]

        atom_names = np.char.asarray(frags_raw.atname)[:frags_natoms]
        selected_atoms_idx = np.asarray(self.args.frag_selections[self.args.frag_sel](atom_names))
        selected_atoms_idx_all = np.concatenate([selected_atoms_idx + frags_natoms * i
                                                 for i in range(frags_nframes)])

        sel_atoms = atom_names[selected_atoms_idx]
        scores = np.asarray(frags_raw.weighting[::frags_natoms])
        poses = np.split(frags_raw.coords[0][selected_atoms_idx_all], frags_nframes)

        selected_oxygens = np.asarray(self.args.frag_selections['poxygens'](atom_names))

        if any(selected_oxygens):
            oxygens_names = atom_names[selected_oxygens]
            oxygens_index = np.concatenate([selected_oxygens + frags_natoms * i for i in range(frags_nframes)])
            oxygens_index = oxygens_index.reshape(int(oxygens_index.size / 2), 2)
            oxygens_coords = frags_raw.coords[0][oxygens_index]
        else:
            oxygens_coords = np.asarray([])
            oxygens_names = np.asarray([])
        return sel_atoms, frag_resname, scores, poses, oxygens_coords, oxygens_names

    def get_zones2(self, frag_poses, sel_atoms):
        cut = self.args.inter_dist

        zones_int = []
        zones_str = []
        zones_both = []
        for pose in frag_poses:
            contacts = self.prot_kdtree.query_ball_point(pose.data, cut)

            close2 = []
            idxs = []
            boths = []
            for i, contact in enumerate(contacts):
                close2.extend(contact)
                both_list = []
                if contact:
                    idxs.append(i)
                    both = list(it.zip_longest(
                        [sel_atoms[i]], self.prot_atom_names[contact],
                        fillvalue=sel_atoms[i]))
                    both_list.append(both)
                if both_list:
                    boths.extend(both_list)

            prot_contact_names2 = self.prot_atom_names[close2]
            prot_contact_residx2 = self.args.prot_parsed.getResindices()[close2]

            # simplest zones nomenclature only with protein residues
            zones_int.append(frozenset(prot_contact_residx2))

            # zones nomenclature 'centered' in the protein
            zone_str = []
            for i, dd in enumerate(prot_contact_residx2):
                zone_str.append('{}-{}'.format(dd, prot_contact_names2[i]))
            zones_str.append(frozenset(zone_str))

            # zones nomenclature 'englobing' protein and fragment
            zone_both = ['-'.join(y) for x in boths for y in x]
            zones_both.append(frozenset(zone_both))

        return zones_str, zones_int, zones_both

    def get_zones(self, frag_poses, sel_atoms):
        cut = self.args.inter_dist
        zones_int = []
        zones_str = []
        zones_both = []
        for pose in frag_poses:
            dist, neighbors = self.prot_kdtree.query(pose.data, k=1,
                                                     distance_upper_bound=cut)
            close = neighbors[dist <= cut]
            frag_contact_names = sel_atoms[dist <= cut]
            prot_contact_names = self.prot_atom_names[close]
            prot_contact_residx = self.args.prot_parsed.getResindices()[close]
            # simplest zones nomenclature only with protein residues
            zones_int.append(frozenset(prot_contact_residx))
            # zones nomenclature 'centered' in the protein
            zone_str = []
            for i, dd in enumerate(prot_contact_residx):
                zone_str.append('{}-{}'.format(dd, prot_contact_names[i]))
            zones_str.append(frozenset(zone_str))
            # zones nomenclature 'englobing' protein and fragment
            zone_both = []
            for i, dd in enumerate(prot_contact_residx):
                zone_both.append('{}-{}_{}'.format(dd, prot_contact_names[i],
                                                   frag_contact_names[i]))
            zones_both.append(frozenset(zone_both))
        return zones_str, zones_int, zones_both


# from confread import parser as config_parser
# args = config_parser('/home/roy.gonzalez-aleman/RoyHub/Chapter4/nuclear_from_hotspots/1sgz_whole.cfg')
# # # args.top = 5
# self = parser(args)
