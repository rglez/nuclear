#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  8 16:04:51 2021

@author: roy.gonzalez-aleman
"""
import configparser
import os
from os.path import join

import numpy as np
import prody as prd


def assert_existence(path, ext=None):
    existence = os.path.exists(path)
    assert existence, '\nNo such file or directory: %s' % path
    if ext:
        sms = '\n {} must have .{} extension.'.format(path, ext)
        assert os.path.basename(path).split('.')[-1] == ext, sms
    return path


def assert_unexistence(path):
    existence = os.path.exists(path)
    sms = '\nFile or dir_ already exists and will not be overwritten: %s' % path
    assert not existence, sms
    return path


def assert_section(section, sections):
    sms = '\nSection [%s] not declared in the config file. Aborting.' % section
    assert section in sections, sms


def check_param_dtype(parsed, param, dtype):
    sms = '\n{} must be of type {}'.format(param, dtype)
    try:
        parsed = dtype(parsed)
        assert isinstance(parsed, dtype), sms
        return parsed
    except ValueError:
        raise ValueError(sms)


def check_param_interval(param, parsed, brackets, mini, maxi):
    if brackets is not None:
        opening, closing = brackets

        if mini is not None:
            if opening == 'o':
                assert parsed > mini, '\n{} must be > {}.'.format(param, mini)
            elif opening == 'c':
                assert parsed >= mini, '\n{} must be >= {}.'.format(param,
                                                                    mini)

        if maxi is not None:
            if closing == 'c':
                assert parsed <= maxi, '\n{} must be <= {}.'.format(param,
                                                                    maxi)
            elif closing == 'o':
                assert parsed < maxi, '\n{} must be < {}.'.format(param, maxi)
        return parsed
    raise ValueError('\nPlease specify a non-null brackets.')


def split_csv(csv_string):
    return [x.strip() for x in csv_string.split(',')]


def read_string(string):
    return string


def read_cfg(cfg_path):
    assert_existence(cfg_path)
    config = configparser.ConfigParser(allow_no_value=True,
                                       inline_comment_prefixes='#')
    config.optionxform = str
    config.read(cfg_path)
    return config


def check_cfg(param_dict, cfg_object):
    checked_param_names = set()
    parsed_param_names = set()
    checked_dict = {}

    for function in param_dict:
        for cfg_argument in param_dict[function]:
            parsed = function(*param_dict[function][cfg_argument])
            section, param = cfg_argument.split('.')
            checked_param_names.add(param)
            checked_dict.update({param: parsed})

    for x in cfg_object.values():
        for y in x.keys():
            parsed_param_names.add(y)
    extraneous_param = set.difference(parsed_param_names, checked_param_names)

    assert not extraneous_param, '\nFollowing extraneous parameters found in \
the config file: {}'.format(extraneous_param)

    return checked_dict


def check_job_sections(sections):
    sms1 = '\nIf [spots_params] section is specified in the NUCLEAR' \
           ' configuration file, only these sections must appear:' \
           ' {}.'
    sms2 = '\nThe following sections are the only ones you must specify' \
           ' in the NUCLEAR configuration file: {}.'
    sms3 = '\n {} section has not been specified.'

    if 'spots_params' in sections:
        allowed_sections = ['inputs', 'mcss_files', 'spots_params']
        assert len(sections) == len(allowed_sections),\
            sms1.format(allowed_sections)
    else:
        allowed_sections = ['inputs', 'mcss_files', 'zones_params',
                            'search_params', 'minimization']
        assert len(sections) == len(allowed_sections),\
            sms2.format(allowed_sections)
    for section in allowed_sections:
        assert section in sections, sms3.format(section)
    return allowed_sections


class parser:
    def __init__(self, config_path):
        # read the configuration file
        self.cfg_path = config_path
        self.cfg = read_cfg(self.cfg_path)
        self.sections = self.cfg.sections()
        self.allowed_sections = check_job_sections(self.sections)

        # prot/frag selection definitions
        self.patch_atoms = ['C5T', 'O5T', 'O1P', 'O2P', 'H5T1', 'H5T2', 'H5T3']
        self.prot_selections = ['all', 'heavy', 'protein']
        self.frag_selections = {
            'all': lambda x:
            list(range(len(x))),
            'heavy': lambda x:
            [i for i, y in enumerate(x) if not y.startswith('H')],
            'all_nopatch': lambda x:
            [i for i, y in enumerate(x)
             if i not in np.in1d(x, self.patch_atoms).nonzero()[0]],
            'heavy_nopatch': lambda x:
            [i for i, y in enumerate(x)
             if (i not in np.in1d(x, self.patch_atoms).nonzero()[0])
             and (not y.startswith('H'))],
            'poxygens': lambda x:
            [i for i, y in enumerate(x)
             if (i in np.in1d(x, ["O1P", "O2P", "S1'", "S2'"]).nonzero()[0])]}
        self.topo_dict = {''}

        # ==== [inputs] section ===============================================
        # input_dir parameter
        input_dir = self.cfg['inputs']['input_dir']
        self.input_dir = assert_existence(input_dir)

        # output_dir parameter
        output_dir = self.cfg['inputs']['output_dir']
        self.output_dir = assert_unexistence(output_dir)

        # prot_coords parameter
        prot_coords = self.cfg['inputs']['prot_coords']
        self.prot_coords = assert_existence(prot_coords, ext='pdb')

        # ncores parameter
        ncores = self.cfg['inputs']['ncores']
        ncores = check_param_dtype(ncores, 'ncores', int)
        self.ncores = check_param_interval('ncores', ncores, 'oc', 0, None)

        # ==== [mcss_files] section ===========================================
        self.crd_files = list(self.cfg['mcss_files'])
        if len(self.crd_files) < 1:
            raise ValueError('\nYou must specify at least one *.crd'
                             ' MCSS exploration')
        crd_list = [tuple(x.split()) for x in self.crd_files]

        for file, maxscore, nposes in crd_list:
            assert_existence(join(self.input_dir, file))
            if maxscore != 'N':
                check_param_dtype(maxscore, 'maxscore', float)
            if nposes != 'N':
                check_param_dtype(nposes, 'nposes', int)
            if 'N' not in [maxscore, nposes]:
                raise ValueError(
                    '\n{} parsing can lead to ambigous '.format(file) +
                    'selection of poses. Set one of the two last'
                    ' columns to N in the configuration file.')
        self.crd_list = crd_list

        # ==== [spots_params] section =========================================
        if 'spots_params' in self.sections:

            # density_cut parameter
            density_cut = self.cfg['spots_params']['density_cut']
            density_cut = check_param_dtype(density_cut, 'density_cut', float)
            self.density_cut = check_param_interval('density_cut', density_cut,
                                                    'cc', 0, 1)
            # merge_cut parameter
            merge_cut = self.cfg['spots_params']['merge_cut']
            merge_cut = check_param_dtype(merge_cut, 'merge_cut', float)
            self.merge_cut = check_param_interval('merge_cut', merge_cut, 'cc',
                                                  0, 1)

        # ==== [(zones/spots)_params] sections ================================
        if ('zones_params' in self.sections) or\
                ('spots_params' in self.sections):
            # prot_sel parameter
            try:
                prot_sel = self.cfg['spots_params']['prot_sel']
            except KeyError:
                prot_sel = self.cfg['zones_params']['prot_sel']

            if prot_sel in self.prot_selections:
                self.prot_sel = prot_sel
            else:
                raise ValueError('\nSpecified prot_sel must be in {}'
                                 .format(self.prot_selections))

            # parsing protein
            prot_raw = prd.parsePDB(self.prot_coords)
            self.prot_parsed = prot_raw.select(self.prot_sel)

            # frag_sel parameter
            try:
                frag_sel = self.cfg['spots_params']['frag_sel']
            except KeyError:
                frag_sel = self.cfg['zones_params']['frag_sel']

            if frag_sel in self.frag_selections:
                self.frag_sel = frag_sel
            else:
                raise ValueError('\nSpecified frag_sel must be in {}'
                                 .format(self.frag_selections))

            # inter_dist parameter
            try:
                inter_dist = self.cfg['spots_params']['inter_dist']
            except KeyError:
                inter_dist = self.cfg['zones_params']['inter_dist']
            inter_dist = check_param_dtype(inter_dist, 'inter_dist', float)
            self.inter_dist = check_param_interval('inter_dist', inter_dist,
                                                   'co', 0, None)

            # pre_clustering parameter
            try:
                pre_clustering = self.cfg['spots_params']['pre_clustering']
            except KeyError:
                pre_clustering = self.cfg['zones_params']['pre_clustering']
            if pre_clustering == 'True':
                self.pre_clustering = True
            elif pre_clustering == 'False':
                self.pre_clustering = False
            else:
                raise ValueError(
                    '\nSpecified pre_clustering must be True or False')

            # rmsd_cut parameter
            try:
                rmsd_cut = self.cfg['spots_params']['rmsd_cut']
            except KeyError:
                rmsd_cut = self.cfg['zones_params']['rmsd_cut']

            if pre_clustering == 'True':
                rmsd_cut = check_param_dtype(rmsd_cut, 'rmsd_cut', float)
                self.rmsd_cut = check_param_interval('rmsd_cut', rmsd_cut,
                                                     'cc', 0, None)
            else:
                assert rmsd_cut == 'None', '\n Set rmsd_cut to None if' \
                                           ' pre_clustering parameter is False'
                self.rmsd_cut = rmsd_cut

        # ==== [search_params] section ========================================
        if 'search_params' in self.allowed_sections:
            # top_X parameter
            top_X = self.cfg['search_params']['top_X']
            top_X = check_param_dtype(top_X, 'top_X', int)
            self.top_X = check_param_interval('top_X', top_X, 'co', 1, None)

            # seq_min_size parameter
            seq_min_size = self.cfg['search_params']['seq_min_size']
            seq_min_size = check_param_dtype(seq_min_size, 'seq_min_size', int)
            self.seq_min_size = check_param_interval(
                'seq_min_size', seq_min_size, 'oo', 1, None)

            # seq_max_size parameter
            seq_max_size = self.cfg['search_params']['seq_max_size']
            if seq_max_size == "N":
                self.seq_max_size = 1000000000
            else:
                seq_max_size = check_param_dtype(seq_max_size, 'seq_max_size', int)
                self.seq_max_size = check_param_interval(
                    'seq_max_size', seq_max_size, 'cc', self.seq_min_size, None)

            # seq_to_search parameter
            seq_to_search = self.cfg['search_params']['seq_to_search']
            if seq_to_search == 'all':
                self.seq_to_search = 'all'
            else:
                seq_split = seq_to_search.split(':')
                if len(seq_split) < 2:
                    raise ValueError('\nSpecified seq_to_search must contain'
                                     ' at least one separator ":"')
                else:
                    self.seq_to_search = seq_split

            # seq_to_search modifies seq_min/max
            if self.seq_to_search != 'all':
                if (self.seq_min_size != len(self.seq_to_search))\
                        or (self.seq_max_size != len(self.seq_to_search)):
                    raise ValueError(
                        "\nIf seq_to_search was specified, seq_min_size and"
                        " seq_max_size must equal seq_to_search's lenght. ({})"
                        .format(len(self.seq_to_search)))
                else:
                    self.seq_min_size = self.seq_max_size \
                        = len(self.seq_to_search)

            # path_to_search parameter
            path_to_search_str = self.cfg['search_params']['path_to_search']
            self.searchpath = frozenset(
                self.prot_parsed.select(path_to_search_str).getResindices())
            assert len(self.searchpath) > 0,\
                'The specified path_to_search in your prot_coords' \
                ' corresponds to no atoms.'

            # max_dist_O3C5 parameter
            max_dist_O3C5 = self.cfg['search_params']['max_dist_O3C5']
            max_dist_O3C5 = check_param_dtype(max_dist_O3C5, 'max_dist_O3C5',
                                              float)
            self.max_dist_O3C5 = check_param_interval(
                'max_dist_O3C5', max_dist_O3C5, 'cc', 0, None)

            # clash_dist parameter
            clash_dist = self.cfg['search_params']['clash_dist']
            clash_dist = check_param_dtype(clash_dist, 'clash_dist', float)
            self.clash_dist = check_param_interval('clash_dist', clash_dist,
                                                   'oo', 0, self.max_dist_O3C5)

        # ==== [minimization] section =========================================
        if 'minimization' in self.allowed_sections:
            # minimize parameter
            minimize = self.cfg['minimization']['minimize']
            if minimize == 'True':
                self.minimize = True
            elif minimize == 'False':
                self.minimize = False
            else:
                raise ValueError('\nSpecified minimize must be True or False')

            # protna_topol/param parameters
            prot_topol = self.cfg['minimization']['prot_topol']
            prot_param = self.cfg['minimization']['prot_param']

            if self.minimize:
                self.prot_topol = assert_existence(prot_topol)
                self.prot_param = assert_existence(prot_param)
            else:
                sms = "\nIf minimize = False, prot_topol and prot_param" \
                      " must be set to None"
                assert prot_topol == 'None', sms
                self.prot_topol = None
                assert prot_param == 'None', sms
                self.prot_param = None

        # ==== Detecting extraneous parameters ================================
        # NUCLEAR valid parameters
        valid_params = {
            'DEFAULT':
                [],
            'inputs':
                ['input_dir', 'output_dir', 'prot_coords', 'ncores'],
            'mcss_files':
                self.crd_files,
            'zones_params':
                ['prot_sel', 'frag_sel', 'inter_dist', 'pre_clustering',
                 'rmsd_cut'],
            'search_params':
                ['top_X', 'seq_min_size', 'seq_max_size', 'seq_to_search',
                 'path_to_search', 'max_dist_O3C5', 'clash_dist'],
            'minimization':
                ['minimize', 'prot_topol', 'prot_param'],
            'spots_params':
                ['density_cut', 'merge_cut', 'prot_sel', 'frag_sel',
                 'inter_dist', 'pre_clustering', 'rmsd_cut']
        }

        # cfg parsed parameters
        parsed_params = {x.name: list(x.keys())
                         for x in list(self.cfg.values())}

        # check for extraneous sections or parameters
        sms = '\n{} parameter is not recognized by this version of NUCLEAR.'
        for section in parsed_params:
            for param in parsed_params[section]:
                assert param in valid_params[section], sms.format(param)

        # ==== NUCLEAR dir_ hierarchy ==========================================
        os.makedirs(self.output_dir)
        if 'spots_params' in self.sections:
            self.TSV = join(self.output_dir, 'TSV')
            os.makedirs(self.TSV)
        else:
            self.out_seqdir = join(self.output_dir, 'SEQUENCES')
            os.makedirs(self.out_seqdir)

            self.out_topXdir = join(self.output_dir, 'topX_logs')
            os.makedirs(self.out_topXdir)

            self.out_popdir = join(self.output_dir, 'popularity_logs')
            os.makedirs(self.out_popdir)

        # write the configuration file for reproducibility
        with open(join(self.output_dir, 'nuclear.cfg'), 'wt') as ini:
            self.cfg.write(ini)


# self = parser('/home/roy.gonzalez-aleman/rprojects/nuclear/examples/1.0.0/'
#               '2xnr_0_6_preliminar.cfg')
