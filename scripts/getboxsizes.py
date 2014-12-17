#!/usr/bin/env python

from __future__ import division
import argparse,sys,os,gzip,glob
import numpy as np
import mdtraj as md

parser = argparse.ArgumentParser()
parser.add_argument('filename',
                    help="e.g. dyn2.out. If this is a directory, we will assume it is full of numbered dynXXX.out files and use them all.")
parser.add_argument('-c', '--nsavc',
                    help="NSAVC from your CHARMM input file. If specified, we will only print out box lines from those steps",
                    default=0,type=int)
parser.add_argument('-p', '--printbox',default=False,action='store_true',
                    help="Print the box information to STDOUT")
parser.add_argument('--traj',
                    default=None,
                    help="If you specify a trajectory and topology, we'll write out .crd files corresponding to the indicated frames")
parser.add_argument('--top',
                    default=None,
                    help="If you specify a trajectory and topology, we'll write out .crd files corresponding to the indicated frames")
parser.add_argument('--outprefix',
                    default='boxes',
                    help="Output files will be called <prefix>min.crd <prefix>max.crd etc.")

args = parser.parse_args()

def savecrd(traj,frameno,fname='test.crd',titleinfo=None):
    """
    Note that we have to explicitly convert from nm to A when saving crds.
    """
    with open(fname,'w') as f:
        atoms = list(traj.top.atoms)
        if titleinfo:
            if not titleinfo[-1] == '\n':
                titleinfo = titleinfo + '\n'
                f.write('* {titleinfo}'.format(**locals()))
        f.write('* Frame {frameno}\n*\n'.format(**locals()))
        #     73568  EXT
        #1234567890123456
        f.write('{n:10}  EXT\n'.format(n=len(atoms)))
        for (atom,(x,y,z)) in zip(atoms,traj.xyz[frameno]):
            #         1         1  DPPC      NC3          -112.2567853263     -115.5683915952       23.1908541326  L         1               0.0000000000
            #     73568     51040  W         W            -102.1950077444       57.1395581102      -24.1123896281  WAT       48992           0.0000000000
            #12345678901234567890  12345678  12345678123456789012345678901234567890123456789012345678901234567890  12345678901234567890
            # PDB format is not good for us, but it's what we have to use as a topology for now.
            # NOT GENERAL!!! The calculation of resi and segi is hardcoded to this simulation.
            if atom.residue.name == 'DPPC':
                resn,resi = atom.residue.name,atom.residue.resSeq
                segn,segi = 'L',resi
            elif atom.residue.name == 'W':
                resn,resi = atom.residue.name,atom.serial - 22528
                segn,segi = 'WAT',resi-2048
            else:
                sys.exit('Unknown residue {r}'.format(r=atom.residue.name))
            f.write('{atomid:10}{resi:10}  {resn:<8}  {name:<8}{x:20.10f}{y:20.10f}{z:20.10f}  {segn:<10}{segi:<10}      0.0000000000\n'.format(atomid=atom.serial,resi=resi,resn=resn,
                                                                                                                                                name=atom.name,
                                                                                                                                                x=10*x,y=10*y,z=10*z,
                                                                                                                                                segn=segn,segi=segi))


if os.path.isdir(args.filename):
    filenames = glob.glob(os.path.join(args.filename,'dyn*.out*'))
    def getnum(fname):
        # call splitext twice in case of .gz
        result = os.path.basename(os.path.splitext(os.path.splitext(fname)[0])[0])
        result = result.replace('dyn','')
        return int(result)
    result = [(getnum(f),f) for f in filenames]
    result.sort()
    filenames = [f for (i,f) in result]
else:
    filenames = [args.filename]

stepsandtimes, sizes = [],[]

for filename in filenames:
    print("Processing {f}".format(f=filename))
    if filename.endswith('.gz'):
        opener = gzip.open
    else:
        opener = open
    with opener(filename,'r') as f:
        for line in f:
            if line.startswith('DYNA>'):
                # Deal with lack of spaces e.g. 10000.00000-325371.15008
                parts = line.replace('-',' -').split()
                stepsandtimes.append((int(parts[1]),float(parts[2])))
            elif line.strip().startswith('DYNA A '):
                parts = line.split()
                sizes.append((float(parts[3]),float(parts[6]),float(parts[9])))
            
if len(stepsandtimes) != len(sizes):
    sys.exit('Found {a} steplines but {b} sizelines'.format(a=len(steplines),b=len(sizelines)))

# Remove duplicates. These come because every NTRFRQ steps, the
# information is printed out twice

_stepsandtimes,_sizes = [],[]
for ((step,time),(a,b,c)) in zip(stepsandtimes,sizes):
    if (not _stepsandtimes) or (step != _stepsandtimes[-1][0]):
        _stepsandtimes.append((step,time))
        _sizes.append((a,b,c))

stepsandtimes, sizes = _stepsandtimes, _sizes


def lcm(n1,n2):
    # make sure n2 is bigger
    if n1 > n2:
        n1,n2 = n2,n1

    for i in range(n2,n1*n2+1):
        if (i % n2 == 0) and (i % n1 == 0):
            return i
    return i


if args.nsavc > 0:
    nprint = stepsandtimes[1][0] - stepsandtimes[0][0]
    stepwant = lcm(nprint,args.nsavc)
    stride = int(stepwant/nprint)
    print "nprint",type(nprint),"stepwant",type(stepwant),"nsavc",type(args.nsavc),"stride",type(stride)
    #print "NPRINT is {nprint} and NSAVC is {args.nsavc} so I will take every {stepwant}th frame".format(**locals())
    stepsandtimes, sizes = stepsandtimes[::stride], sizes[::stride]

print "Processing {nf} frames".format(nf=len(stepsandtimes))

mina,imina = np.inf,np.inf
maxa,imaxa = -np.inf,-np.inf
mida,imida,midadiff = np.inf,np.inf,np.inf
avga,iavga,avgadiff = np.inf,np.inf,np.inf
alla = [a for (a,b,c) in sizes]
allc = [c for (a,b,c) in sizes]
_avga = np.average(alla)
_mida = min(alla) + (max(alla) - min(alla))/2
for (i,(a,b,c)) in enumerate(sizes):
    if a < mina:
        imina,mina = i,a
    if a > maxa:
        imaxa,maxa = i,a
    if abs(a-_avga) < avgadiff:
        iavga,avga,avgadiff = i,a,abs(a-_avga)
    if abs(a-_mida) < midadiff:
        imida,mida,midadiff = i,a,abs(a-_mida)
    

maxtit = "Max box size A=B={m} C={c} with area {apl:.2f} at step {s}, time {t} ps".format(m=maxa,
                                                                                          c=allc[imaxa],
                                                                                          apl=(maxa*maxa)/1024,
                                                                                          s=stepsandtimes[imaxa][0],
                                                                                          t=stepsandtimes[imaxa][1],
                                                                                          )
mintit = "Min box size A=B={m} C={c} with area {apl:.2f} at step {s}, time {t} ps".format(m=mina,
                                                                                          c=allc[imina],
                                                                                          apl=(mina*mina)/1024,
                                                                                          s=stepsandtimes[imina][0],
                                                                                          t=stepsandtimes[imina][1],
                                                                                          )
avgtit = "Avg box size A=B={m} C={c} with area {apl:.2f} at step {s}, time {t} ps".format(m=avga,
                                                                                          c=allc[iavga],
                                                                                          apl=(avga*avga)/1024,
                                                                                          s=stepsandtimes[iavga][0],
                                                                                          t=stepsandtimes[iavga][1],
                                                                                          )
midtit = "Mid box size A=B={m} C={c} with area {apl:.2f} at step {s}, time {t} ps".format(m=mida,
                                                                                          c=allc[imida],
                                                                                          apl=(mida*mida)/1024,
                                                                                          s=stepsandtimes[imida][0],
                                                                                          t=stepsandtimes[imida][1],
                                                                                          )
print maxtit,alla[imaxa]
print mintit,alla[imina]
print avgtit,alla[iavga]
print midtit,alla[imida]
if args.traj and args.top:
    t = md.load(args.traj,top=args.top)
    #savecrd(t,imaxa-1,fname='{prefix}max_pre.crd'.format(prefix=args.outprefix),titleinfo=maxtit)
    savecrd(t,imaxa-1,fname='{prefix}max.crd'.format(prefix=args.outprefix),titleinfo=maxtit)
    #savecrd(t,imaxa+1,fname='{prefix}max_post.crd'.format(prefix=args.outprefix),titleinfo=maxtit)
    savecrd(t,imina-1,fname='{prefix}min.crd'.format(prefix=args.outprefix),titleinfo=mintit)
    savecrd(t,iavga-1,fname='{prefix}avg.crd'.format(prefix=args.outprefix),titleinfo=avgtit)
    savecrd(t,imida-1,fname='{prefix}mid.crd'.format(prefix=args.outprefix),titleinfo=midtit)
    
if args.printbox:                                                                             
    for ((step,time),(a,b,c)) in zip(stepsandtimes,sizes):
        print "{step}\t{time}\t{a}\t{b}\t{c}".format(**locals())

print "I do not know why there is an off-by-one error in the calculation of imaxa imina etc."
