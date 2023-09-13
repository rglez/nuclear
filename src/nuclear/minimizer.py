#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 14:43:40 2021

@author: roy.gonzalez-aleman
"""
import glob
import itertools as it
import os
from collections import defaultdict
from os.path import join

protocol = '''
\n
!++++++++++++++++++++++
! MINIMIZATION PROTOCOL
!++++++++++++++++++++++

! Atom Definition
! ------------------
DEFINE PROT selec segi PRO* show end
DEFINE RNA selec segi RNZ show end

! nonbonded definition
nbonds elec atom shif rdie vatom vdw vshi vdis cutnb 8.5 ctofnb 7.5 -
  ctonnb 6.5 wmin 0.9 eps 3.0 e14f 0.4 nbxm 5

! Minimization 1 / ATOMS LINKING NUCLEOTIDES
! ---------------
cons fix selec .not. (type P .or. type O3' .or. type C3' .or. type O5' .or. type C5' .or. type O1P .or. type O2P .or. segid WAT*) show end
ENERGY
MINI SD NSTE 1000 NPRI 5 TOLG 10.0 STEP 0.02 TOLS 0.0 TOLENR 0.0
MINI ABNR NSTE 10000 NPRI 5 TOLG 0.1 STEP 0.02 TOLS 0.0 TOLENR 0.0

! Minimization 2 / Phosphodiester skeleton
! ---------------
cons fix selec .not. (type P .or. type C+' .or. type O+' .or. type H* .or. segid WAT*) show end
ENERGY
MINI SD NSTE 1000 NPRI 5 TOLG 10.0 STEP 0.02 TOLS 0.0 TOLENR 0.0
MINI ABNR NSTE 10000 NPRI 5 TOLG 0.1 STEP 0.02 TOLS 0.0 TOLENR 0.0

! Minimization 3 / RNA released
! ---------------
cons fix selec .not. (segid RNZ .or. segid WAT*) show end
ENERGY
MINI SD NSTE 10000 NPRI 5 TOLG 10.0 STEP 0.02 TOLS 0.0 TOLENR 0.0
MINI ABNR NSTE 10000 NPRI 5 TOLG 0.1 STEP 0.02 TOLS 0.0 TOLENR 0.0

! Minimization 4 / RNA & Protein released
! --------------
cons fix selec none end
ENERGY
MINI SD NSTE 10000 NPRI 5 TOLG 10.0 STEP 0.02 TOLS 0.0 TOLENR 0.0
MINI ABNR NSTE 10000 NPRI 5 TOLG 0.1 STEP 0.02 TOLS 0.0 TOLENR 0.0

'''

measures = '''
open write card unit 50 name @Savedata1
write title unit 50
* EINT ERNA ERNA0 ETOT
*

interaction select PROT end select RNA end
SET EINT  ?ener

interaction select RNA end
SET ERNA  ?ener

DELEte atom sele PROT end

MINI SD NSTE 10000 NPRI 5 TOLG 10.0 STEP 0.02 TOLS 0.0 TOLENR 0.0
MINI ABNR NSTE 10000 NPRI 5 TOLG 0.1 STEP 0.02 TOLS 0.0 TOLENR 0.0

interaction select all end
SET ERNA0  ?ener

calc ETOT = @EINT + ( @ERNA - @ERNA0)

write title unit 50
* @EINT @ERNA @ERNA0 @ETOT
*

Close unit 50

stop
'''


class minimizer():
    def __init__(self, args):
        self.args = args
        self.out_seqdirs = {x: join(os.path.abspath(self.args.out_seqdir), x)
                            for x in os.listdir(self.args.out_seqdir)}
        self.sequences = {x: sorted(glob.glob(join(self.args.out_seqdir,
                                                   x + '/*seq*pdb')))
                          for x in self.out_seqdirs}

    def read_and_split_miniprot(self):
        prot_path = self.args.prot_coords
        basename = os.path.basename(prot_path).split('.')[0].lower()
        # read & split all sequences found in ATOM | HETATM records
        with open(prot_path, 'rt') as miniprot:
            segnames = defaultdict(list)
            for line in miniprot:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    segname = line[72:76].lower()
                    segnames[segname].append(line)
        # write separated sequences found in the protein
        self.prot_files = {}
        for out_dir in self.out_seqdirs:
            os.chdir(join(self.args.output_dir,
                          self.out_seqdirs[out_dir]))
            prot_files = {}
            for segname in segnames:
                outbasename = '{}_{}.pdb'.format(basename, segname)
                outfile = join(self.out_seqdirs[out_dir],
                               outbasename)
                prot_files.update({segname: outfile})
                with open(outfile, 'wt') as segment:
                    for line in segnames[segname]:
                        segment.write(line)
                    segment.write('END\n')
            self.prot_files.update({out_dir: prot_files})
        self.minim_dir = self.out_seqdirs

    def write_minim_inp(self):
        for ntdir in self.minim_dir:
            os.chdir(join(self.args.output_dir, self.minim_dir[ntdir]))
            for sequence in self.sequences[ntdir]:
                numid = it.count(10, 10)
                seq_name = os.path.basename(sequence).split('.')[0]
                with open('{}.inp'.format(seq_name), 'wt') as inp:
                    # proteins_path topo and param
                    inp.write('DIMENS CHSIZE 5000000\n')
                    inp.write('BOMBlev -1\n')
                    inp.write('set 0 ./\n\n')
                    inp.write('! Read topology and parameter files\n')
                    topid = next(numid)
                    inp.write('open read card unit {} name "{}"\n'.format(
                        topid, self.args.prot_topol))
                    inp.write('read  rtf card unit {}\n\n'.format(topid))
                    parid = next(numid)
                    inp.write('open read card unit {} name "{}"\n'.format(
                        parid, self.args.prot_param))
                    inp.write('read para card unit {}\n\n'.format(parid))

                    # RNZ segment section
                    cmpx_id = next(numid)
                    inp.write('! Read {}\n'.format(seq_name))
                    inp.write('open read card unit {} name {}\n'.format(cmpx_id, os.path.basename(sequence)))
                    inp.write('read nuclear_sequence pdb unit {}\n'.format(cmpx_id))
                    inp.write('generate RNZ setup warn first none last none\n\n')
                    inp.write('open read card unit {} name {}\n'.format(cmpx_id, os.path.basename(sequence)))
                    inp.write('read coor pdb unit {}\n\n'.format(cmpx_id))

                    # PRO* and WAT* segment section
                    for segname in self.prot_files[ntdir]:
                        inp.write('! Read {}\n'.format(segname.upper()))
                        inp.write('open read card unit {} name {}\n'.format(
                            cmpx_id, os.path.basename(self.prot_files[ntdir][segname])))
                        inp.write('read nuclear_sequence pdb unit {}\n'.format(cmpx_id))
                        # the generate sentence
                        if segname.startswith('pro'):
                            inp.write('generate {}\n\n'.format(segname.upper()))
                        elif segname.startswith('wat'):
                            inp.write('generate {} setup warn noangle nodihedral\n\n'.format(segname.upper()))
                        inp.write('open read card unit {} name {}\n'.format(
                            cmpx_id, os.path.basename(self.prot_files[ntdir][segname])))
                        inp.write('read coor pdb unit {} resid\n\n'.format(cmpx_id))
                    # generate and hbuild section
                    inp.write('ic generate\nic param\nic build\nhbuild sele all end\n!coor print\n\n')
                    # minimization protocol
                    inp.write(protocol)
                    # writing output files
                    inp.write('open write card unit {} name complex_{}_mini.pdb\n'.format(cmpx_id, seq_name))
                    inp.write('write coor pdb unit {}\n'.format(cmpx_id))
                    # charmm internal measures
                    inp.write('! Measures\n')
                    inp.write('! --------\n')
                    inp.write('set Savedata1 = @0/complex_{}_mini_measures.dat\n'.format(seq_name))
                    inp.write(measures)
                    inp.write('stop')


# =============================================================================
#
# =============================================================================
# from confread import parser as config_parser
# args = config_parser('/home/roy.gonzalez-aleman/rprojects/nuclear/examples/'
#                     '1.0.0/2xnr_0_6_preliminar.cfg')
# self = minimizer(args)
