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

parser = argparse.ArgumentParser(description='Calculate diffusion constants, via block averaging.',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-o', '--output-dir',
                    help='Name of output directory; use -O to overwrite',
                    default='diffusion')
parser.add_argument('-O', '--overwrite',
                    help='Overwrite existing files. This means we will recalculate the COM and MSD trajectories.',
                    action='store_true',default=False)
parser.add_argument('trajectory',
                    help='A single, unfolded trajectory file')
parser.add_argument('-t', '--topology',
                    help='A topology file. This is required if you are not feeding an h5-formatted trajectory. If you have at least version 1.1 of mdtraj, you can pass a psf. Otherwise, you will probably want to pass a pdb')

parser.add_argument('-c', '--clean-intermediate-files',
                    help='Clean up intermediate files, e.g. the com trajectories and msds. These will be in either h5 or npy (numpy) format. When present, intermediate files are used in favor of recalculations, so make sure to delete them by hand if you need to.',
                    action='store_true',default=False)

parser.add_argument('-a','--one-atom-name',
                    help='The name of one atom in your residue of interest. Used as an intermediate in making the COM trajectory, this could be fixed.',
                    default='PO4')

parser.add_argument('--do-not-fix-length-for-charmm',
                    dest='fix_length_for_charmm',default=True,action='store_false',
                    help='By default, if we are given a dcd and a psf, or a trj and a psf, we will multiply the times by 100 and multiply the coordinates by 10. At the time of writing, this is necessary because CHARMM trajectories are in ns and A internally, whereas mdtraj expects ps and nm. DO NOT DO GIVE THIS OPTION UNLESS YOU KNOW WHAT YOU ARE DOING')

parser.add_argument('--do-not-fix-time-for-charmm',
                    dest='fix_time_for_charmm',default=True,action='store_false',
                    help='By default, if we are given a dcd and a psf, or a trj and a psf, we will multiply the times by 100 and multiply the coordinates by 10. At the time of writing, this is necessary because CHARMM trajectories are in ns and A internally, whereas mdtraj expects ps and nm. DO NOT DO GIVE THIS OPTION UNLESS YOU KNOW WHAT YOU ARE DOING')

parser.add_argument('-start',
                    help='Skip this many nanoseconds at the start of the trajectory as equilibration.',
                    type=int,default=0)

parser.add_argument('-n', '--num-blocks',
                    help='Number of blocks to use for block averaging',
                    type=int,default=10)

parser.add_argument('--stop',
                    help='Stop after this many nanoseconds. -1 means use the whole trajectory.',
                    default=None,
                    )

parser.add_argument('--frame-spacing',
                    help='Only calculate between frames spaced this far apart. E.g. if your trajectory stores frames every 50 ps, but you just want to calculate every 100 ps, pass in 100ps. You can pass in arguments in ns or ps, but you must specify.',default=None)



args = parser.parse_args()

if not (args.fix_time_for_charmm and args.fix_length_for_charmm):
     print('\n\nWARNING: Do not specify --do-not-fix-XXX-for-charmm UNLESS YOU HAVE TESTED THE RESULTS! Check the box sizes and times below!')

if args.clean_intermediate_files:
    sys.exit('Clean intermediate files not yet implemented')

if not os.path.isfile(args.trajectory):
    sys.exit('Could not find trajectory "{t}"'.format(t=args.trajectory))

if args.topology and not os.path.isfile(args.topology):
    sys.exit('Could not find topology "{t}"'.format(t=args.topology))

#if not args.overwrite:
#    if os.path.isdir(args.output_dir) and os.listdir(args.output_dir):
#        sys.exit('Output directory "{d}" is not empty. Please delete the files within or specify -O/--overwrite'.format(d=args.output_dir))

if os.path.exists(args.output_dir) and not os.path.isdir(args.output_dir):
    sys.exit('"{d}" is a file, not a directory'.format(d=args.output_dir))

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

if args.topology:
    t = md.load(args.trajectory, top=args.topology)
else:
    t = md.load(args.trajectory)

# Ready to go. Set up our filenames

extra_name_part = '_s{start}_e{stop}_fs{fs}_b{nb}'.format(start=args.start,
                                                          stop=args.stop,
                                                          fs=args.frame_spacing,
                                                          nb=args.num_blocks,)

comname = os.path.join(args.output_dir,'com{enp}.h5'.format(enp=extra_name_part))
msdname = os.path.join(args.output_dir,'msd{enp}.npy'.format(enp=extra_name_part))


#print(dir(t))
if args.fix_time_for_charmm:
    t.time *= 100
if args.fix_length_for_charmm:
    #for a in ('unitcell_angles', 'unitcell_lengths', 'unitcell_vectors', 'openmm_boxes'):
    #    print('{a} {v}'.format(a=a,v=getattr(t,a)))
    t.xyz *= 10
    t.unitcell_vectors *= 10

# Throughout this, we enforce the idea that frame 0 is time 0.
t.time -= t.time[0]

print("Successfully loaded trajectory with shape {s}.\nFirst two times {a} ps {b} ps.".format(s=t.xyz.shape,
                                                                                              a=t.time[0],
                                                                                              b=t.time[1]))
print("The first box has A={a} A B={b} A C={c} A".format(a=t.unitcell_lengths[0][0],
                                                         b=t.unitcell_lengths[0][1],
                                                         c=t.unitcell_lengths[0][2],))


def frametime2ps(fs):
    if fs.endswith('ns'):
        fs = 1000 * float(fs[:-2].strip())
    elif fs.endswith('ps'):
        fs = float(fs[:-2].strip())
    else:
        sys.exit('Frame spacing should end with "ns" for nanoseconds or "ps" for picoseconds; you said {f}.'.format(f=fs))
    return fs
def ps2frames(ps,times):
    curr_spacing = times[1] - times[0]
    stride,oops = divmod(ps,curr_spacing)
    if oops:
        sys.exit('Trajectory has frame spacing {c} ps, which could not be converted evenly to {fs} ps'.format(c=curr_spacing,
                                                                                                              ps = ps))
    return stride

    
if args.frame_spacing:
    fs = frametime2ps(args.frame_spacing)
    stride = ps2frames(fs,t.time)
    t.xyz = t.xyz[::stride]
    t.time = t.time[::stride]
    if fs != t.time[1] - t.time[0]:
        sys.exit('I tried to convert to frame spacing {fs} ps with stride {s} but got {b}'.format(fs=fs,
                                                                                                  s=stride,
                                                                                                  b = t.time[1] - t.time[0]))
    t.time = t.time - t.time[0]
                

def dropframes(t,drop,side,droptype='frames'):
    """
    `side`: 'start' means drop this many from the stop. 'end' means
            drop this many from the end. 'stop' means stop after this
            many frames.
    """
    print("Dropframes {d} {s} {df}".format(d=drop,s=side,df=droptype))
    if droptype == 'psns':
        print('    first {f} second {s}'.format(f=frametime2ps(drop),s=ps2frames(frametime2ps(drop),t.time)))
        drop = ps2frames(frametime2ps(drop),t.time)
        print("  converted that to {d} frames".format(d=drop))
    elif droptype == 'frames':
        pass
    else:
        raise Exception('Unknown drop type {t}'.format(t=droptype))
    print("Before skipping {s} frames, we have shapes {s1} {s2} {s3}".format(s=drop,
                                                                            s1 = t.xyz.shape,
                                                                            s2 = t.unitcell_vectors.shape,
                                                                            s3 = t.time.shape))
    if side == 'start':
        t.xyz = t.xyz[drop:]
        t.unitcell_vectors = t.unitcell_vectors[drop:]
    elif side == 'end':
        t.xyz = t.xyz[:-drop]
        t.unitcell_vectors = t.unitcell_vectors[:-drop]
    elif side == 'stop':
        t.xyz = t.xyz[:drop]
        t.unitcell_vectors = t.unitcell_vectors[:drop]
    else:
        sys.exit('Bad call to dropframes, unknown side "{s}"'.format(s=side))
    # No need to care which side of the times we lop off, because we
    # normalize them all to t0 = 0.
    t.time = t.time[:drop]
    t.time -= t.time[0]
    print("After skipping {s} frames, we have shapes {s1} {s2} {s3}".format(s=drop,
                                                                            s1 = t.xyz.shape,
                                                                            s2 = t.unitcell_vectors.shape,
                                                                            s3 = t.time.shape))

if args.start:
    dropframes(t,args.start,side='start',droptype='psns')
if args.stop:
    dropframes(t,args.stop,side='stop',droptype='psns')

chunklen,remaining = divmod(t.xyz.shape[0],args.num_blocks)
if remaining:
    print("Can't evenly divide {f} frames into {b} blocks; Using blocksize {bs}, dropping last {df} frames".format(f=t.xyz.shape[0],
                                                                                                                   b=args.num_blocks,
                                                                                                                   bs=chunklen,
                                                                                                                   df=remaining))
    
    dropframes(t,drop=remaining,side='end',droptype='frames')


startstop = zip([i + 1 for i in range(0,t.xyz.shape[0]+1,chunklen)][:-1],[i + 1 for i in range(0,t.xyz.shape[0]+1,chunklen)][1:])
startstop = zip([i for i in range(0,t.xyz.shape[0]+1,chunklen)][:-1],[i for i in range(0,t.xyz.shape[0]+1,chunklen)][1:])

print('Using the following frame numbers for the blocks:\n{b}'.format(b=startstop))
print('First chunk {c1} second {c2}'.format(c1=t.xyz[startstop[0][0]:startstop[0][1]].shape,
                                            c2=t.xyz[startstop[-1][0]:startstop[-1][1]].shape,
                                            ))


atom_indices = []
for r in t.topology.residues:
    atom_indices.append([a.index for a in r.atoms])

if os.path.isfile(comname) and not args.overwrite:
    com_traj = md.load(comname)
else:
    # Make a phony trajectory, fill it iwth the COM traj
    com_traj = t.atom_slice([a.index for a in t.topology.atoms if a.name == args.one_atom_name])
    assert com_traj.xyz.shape[1] == len(atom_indices)
    for (i,atoms) in enumerate(atom_indices):
        # check each atom slice to make sure it's just a simple range
        assert list(range(atoms[0],atoms[-1]+1)) == atoms
        com_traj.xyz[:,i,:] = np.average(t.xyz[:,atoms[0]:atoms[-1]+1,:],axis=1)
    # Dump the COM traj
    com_traj.save(comname)
print('COM traj has shape {s}'.format(s=com_traj.xyz.shape))
print('COM traj has timestep {s} ps'.format(s=com_traj.timestep))


if os.path.isfile(msdname) and not args.overwrite:
    msds = np.load(msdname)
else:
    msds = []
    print('Writing MSDs')
    for (start,stop) in startstop:
        msd = np.zeros((chunklen,com_traj.xyz.shape[1],com_traj.xyz.shape[2]))
        print "\nNew Chunk {start} --> {stop}".format(start=start,stop=stop)
        p = mdt.util.ProgressBar(msd.shape[1])
        for i in range(msd.shape[1]):
            p.animate(i)
            msd[1:,i,:] = [np.average((com_traj.xyz[start+j:stop,i,:] - com_traj.xyz[start:stop-j,i,:])**2,axis=0) for j in range(1,chunklen)]
        msds.append(msd)
    msds = np.array(msds)
    np.save(msdname,msds)
    print ('\nMSDs have been written')
print 'MSDs has shape {s}'.format(s=msds.shape)

msds_xy = [msd[:,:,0] + msd[:,:,1] for msd in msds]
times = com_traj.time[:chunklen]

print("Before I calculate the Ds, the first two times {a} ps {b} ps.".format(
                                                                             a=t.time[0],
                                                                             b=t.time[1]))
print("  the last two are {a} ps and {b} ps.".format(a=t.time[-2],b=t.time[-1]))
print("And the first box has A={a} A B={b} A C={c} A".format(a=t.unitcell_lengths[0][0],
                                                             b=t.unitcell_lengths[0][1],
                                                             c=t.unitcell_lengths[0][2],))
print("     the last box has A={a} A B={b} A C={c} A".format(a=t.unitcell_lengths[-1][0],
                                                             b=t.unitcell_lengths[-1][1],
                                                             c=t.unitcell_lengths[-1][2],))

###numpts = len(times)
###droppts = int(numpts/10)
#### nm^2/ps --> cm^2/s
###conversion = 1e-7 * 1e-7 / 1e-12
#### A^2/ps --> cm^2/s
###conversion = 1e-8 * 1e-8 / 1e-12
###ds = []
###for (i,msdxy) in enumerate(msds_xy):
###    av_msd = np.average(msdxy,axis=1)
###    m,b = np.polyfit(times[droppts:-droppts],av_msd[droppts:-droppts],1)
###    d = m*conversion/4
###    ds.append(d)
###    #plt.plot(times,np.average(msdxy,axis=1),label='c{i} {d:.2f}'.format(i=i,d=d/1e-7))
###    #plt.plot(times,m*times+b,'k--')
###ds = np.array(ds)
###print('All ds: {ds}'.format(ds=ds))
###print('D = {d:.2f} +/- {s:.2f} x 10^-7 cm^2/s'.format(d=np.average(ds)/1e-7,s=np.std(ds)/1e-7))

def calcds(times,msds_xy,droppct=0):
    numpts = len(times)
    if droppct > 0:
        droppts = int(numpts/droppct)
    else:
        droppts = 0
    # nm^2/ps --> cm^2/s
    conversion = 1e-7 * 1e-7 / 1e-12
    # A^2/ps --> cm^2/s
    conversion = 1e-8 * 1e-8 / 1e-12
    ds = []
    for (i,msdxy) in enumerate(msds_xy):
        av_msd = np.average(msdxy,axis=1)
        if droppts > 0:
            m,b = np.polyfit(times[droppts:-droppts],av_msd[droppts:-droppts],1)
        else:
            m,b = np.polyfit(times,av_msd,1)
        d = m*conversion/4
        ds.append(d)
        #plt.plot(times,np.average(msdxy,axis=1),label='c{i} {d:.2f}'.format(i=i,d=d/1e-7))
        #plt.plot(times,m*times+b,'k--')
    ds = np.array(ds)
    print('Skipping the first and last {p}% of the MSD trajectories'.format(p=droppct))
    print('  All ds: {ds}'.format(ds=ds))
    print('  D = {d:.2f} +/- {s:.2f} x 10^-7 cm^2/s'.format(d=np.average(ds)/1e-7,s=np.std(ds)/1e-7))
    return ds


calcds(times,msds_xy,0)
calcds(times,msds_xy,10)
