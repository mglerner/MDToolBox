#!/usr/bin/env python

from __future__ import division
import sys, os, argparse
import numpy as np, scipy as sp, pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import mdtraj as md
import MDToolBox as mdt
import MDToolBox.util

from MDToolBox import pygro
parser = argparse.ArgumentParser(description='Box sizes',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('trajectory',
                    help='A single, unfolded trajectory file')
parser.add_argument('-t', '--topology',
                    help='A topology file. This is required if you are not feeding an h5-formatted trajectory. If you have at least version 1.1 of mdtraj, you can pass a psf. Otherwise, you will probably want to pass a pdb')
parser.add_argument('-T','--tableout', default=False,action='store_true',
                    help='Format output for a table',)

args = parser.parse_args()

if not os.path.isfile(args.trajectory):
    sys.exit('Could not find trajectory "{t}"'.format(t=args.trajectory))

if args.topology and not os.path.isfile(args.topology):
    sys.exit('Could not find topology "{t}"'.format(t=args.topology))
if args.topology:
    t = md.load(args.trajectory, top=args.topology)
else:
    t = md.load(args.trajectory)

if 1:
    t.time *= 100
    #for a in ('unitcell_angles', 'unitcell_lengths', 'unitcell_vectors', 'openmm_boxes'):
    #    print('{a} {v}'.format(a=a,v=getattr(t,a)))
    t.xyz *= 10
    t.unitcell_vectors *= 10

shape = t.unitcell_lengths.shape
avg = np.average(t.unitcell_lengths,axis=0)
std = np.std(t.unitcell_lengths,axis=0)

if args.tableout:
    print('| {f:10s} | {a:10.6f} | {smalla:10.6f} | {largea:10.6f} | {stda:10.6f} |'.format(a=avg[0], stda=std[0], smalla=avg[0]-std[0], largea=avg[0]+std[0], f=args.trajectory))
else:
    if avg[0] == avg[1] and std[0] == std[1]:
        print('A=B {a} +/- {sa} small {smalla} large {largea}'.format(a=avg[0], sa=std[0], smalla=avg[0]-std[0], largea=avg[0]+std[0]))
        print('C {c} +/- {sc} small {smallc} large {largec}'.format(c=avg[2], sc=std[2], smallc=avg[2]-std[2], largec=avg[2]+std[2]))
    else:
        print('A {a} +/- {sa} small {smalla} large {largea}'.format(a=avg[0], sa=std[0], smalla=avg[0]-std[0], largea=avg[0]+std[0]))
        print('B {a} +/- {sa} small {smalla} large {largea}'.format(a=avg[1], sa=std[1], smalla=avg[1]-std[1], largea=avg[1]+std[1]))
        print('C {c} +/- {sc} small {smallc} large {largec}'.format(c=avg[2], sc=std[2], smallc=avg[2]-std[2], largec=avg[2]+std[2]))

