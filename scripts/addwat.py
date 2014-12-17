#!/usr/bin/env python

from __future__ import division

import os,sys,time
from subprocess import Popen,PIPE,STDOUT
import pygro
def run(cmd):
    print cmd
    f = file('/tmp/tmp.txt','w')
    p = Popen(cmd, shell=True, stdout=f, stderr=STDOUT)
    print "waiting"
    sts = os.waitpid(p.pid,0)
    return p.returncode

def addwat(f,outgro,outtop,includetop,box,watbox,mol,lipid,numlipid):
    ''' add waters '''
    # Set up the box
    b1,b2,b3 = box
    o = file('vdwradii.dat','w')
    o.write('''; Very approximate VanderWaals radii
; only used for drawing atoms as balls or for calculating atomic overlap.
; longest matches are used
; '???' or '*' matches any residue name
; 'AAA' matches any protein residue name
???  C     0.375
???  F     0.12
???  H     0.04
???  N     0.110
???  O     0.105
???  S     0.16
; MARTINI particles
;W    W     1.35 Not the water! We want this to use the parameters from the box!
DPP  NC3   0.375
DPP  PO4   0.375
DPP  GL1   0.375
DPP  GL2   0.375
DPP  C1A   0.375
DPP  C2A   0.375
DPP  C3A   0.375
DPP  C4A   0.375
DPP  C1B   0.375
DPP  C2B   0.375
DPP  C3B   0.375
DPP  C4B   0.375
; Water charge sites
SOL  MW    0
SOL  LP    0
;???  WH   15.00
;???  HP4  15.00
; Masses for vsite construction
GLY  MN1   0
GLY  MN2   0
ALA  MCB1  0
ALA  MCB2  0
VAL  MCG1  0
VAL  MCG2  0
ILE  MCG1  0
ILE  MCG2  0
ILE  MCD1  0
ILE  MCD2  0
LEU  MCD1  0
LEU  MCD2  0
MET  MCE1  0
MET  MCE2  0
TRP  MTRP1 0
TRP  MTRP2 0
THR  MCG1  0
THR  MCG2  0
LYSH MNZ1  0
LYSH MNZ2  0
''')
    o.close()
    if (b1,b2) != (0.,0.):
        run('editconf -f %s -c -box %f %f %f -o %s'%(f,b1,b2,b3,outgro))
    else:
        run('editconf -f %s -c -d 0 -o %s'%(f,outgro))
        b1,b2,_b3 = [float(i) for i in file(outgro).readlines()[-1].split()]
        b3 = _b3 + b3
        print b1,b2,b3,'yo'
        run('editconf -f %s -c -box %f %f %f -o %s'%(outgro,b1,b2,b3,outgro))
    # Add the waters
    run('genbox -cp %s -cs %s -vdwd 0.19 -o %s'%(outgro,watbox,outgro))
    run('rm vdwradii.dat')
    # find out how man waters we've added
    o = file(outgro).readlines()
    header = o[:2]
    box = [o[-1],]
    w = []
    wf = []
    p = []
    for line in o[2:-1]:
        if pygro.groline2parts(line)[1].strip() == 'WF':
        #if 'WF' in line:
            wf.append(line)
        elif pygro.groline2parts(line)[1].strip() == 'W':
        #elif 'W' in line:
            w.append(line)
        else: 
            p.append(line)
    o = file(outgro,'w')
    o.writelines(header+p+w+wf+box)
    o.close()
    run('editconf -f %s -o %s'%(outgro,outgro)) #renumber
    o = file(outtop,'w')
    topdir,atop = os.path.split(includetop)
    tops = ('martini_v2.1.itp','martini_v2.1_aminoacids.itp','martini_v2.0_ions.itp','martini_v2.0_lipids.itp',atop)
    for top in tops:
        o.write('#include "%s"\n'%os.path.join(topdir,top))
    o.write('''[ system ]
SOLV%s
[ molecules ]
%s 1
%s %s
W %d
WF %d
'''%(mol,
     mol,
     lipid,numlipid,
     len(w),
     len(wf)))
    o.close()
        

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-b','--box',nargs=3,type='float',help="If you specify all 3, they'll be treated as the box. If you specify 0 0 Z, we'll keep X and Y the same, but add your specified value in Z. NOTE: editconf is used to figure out the box. X and Y may differ, as is standard in GROMACS, but is not standard in CHARMM.")
    parser.add_option('-w','--watbox',default='water.gro')
    parser.add_option('-f',default='pdb-CG.pdb',help='Input file, can be .pdb or .gro')
    parser.add_option('-o','--outgro',default='pdb-CG-genbox.gro',help='output file, must be .gro')
    parser.add_option('-t','--outtop',default='protsol.top')
    parser.add_option('-m','--mol',default='DIM')
    parser.add_option('-l','--lipid',default='DPP')
    parser.add_option('-L','--numlipid',default=0,type='int',help='If greater than 0, we will add this many lipids to the topology')
    parser.add_option('-I','--includetop',default='./protein_withele.itp',help='Full path a topology file to include. We will also include the standard MARTINI topology files from the same directory')
    (options,args) = parser.parse_args()

    for f in options.watbox,options.f:
        if not os.path.isfile(f):
            raise IOError('Could not find %s'%f)
    for f in (options.outgro,options.outtop):
        if os.path.isfile(f):
            print "Will overwrite",f
    if options.box is None:
        sys.exit('Please specify box, e.g. 10 10 10')
        
    addwat(f=options.f,
           outgro=options.outgro,
           outtop=options.outtop,
           includetop=options.includetop,
           box=options.box,
           watbox=options.watbox,
           mol=options.mol,
           lipid=options.lipid,
           numlipid=options.numlipid,
           )
