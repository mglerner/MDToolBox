#!/usr/bin/env python
"""
Pure Python version of GROMACS utils

Since this is pyGRO, we'll store things internally in GROMACS units, i.e.

length: nm
energy:

"""
from __future__ import division
from __future__ import print_function


__version__ = 0.2
import os,sys,re
import numpy as np
from numpy import random,std,average,array,sqrt
import matplotlib.pyplot as plt
from . import util as U
from .fileio import parts2groline, parts2crdline, groline2parts

class GromacsError(Exception):
    pass


class Conversion(object):
    """Convert from, e.g. GROMACS to CHARMM. Most methods will be static.
    """
    
    def __init__(self, ):
        """Most methods will be static.
        """
        
        pass
    @staticmethod
    def gro2crd(fn,outfn=None,cleanup=True,segmap={},resnmap={},renumber=None,forceext=False,fixlip=0):
        """Converts a .gro file to a CHARMM .crd file, first running
        it through editconf to make sure that it's properly centered.
        
        Arguments:
        - `fn`: input file name.
        - `outfn`: output filename. If None, will be input with extension
                   replaced with .crd
        - `cleanup`: editconf makes a centered file. If cleanup is True, delete that file.
        - `segmap`: a mapping from resns to desired segments
        - `resnmap`: a mapping from existing to desired resnames
        - `renumber`: Force renumbering of atom IDs. If you have more than 100000 lines, 
                      you need to specify this as True or False.

        TODO: allow user to specify cfn. Or, better yet, just center in pure Python.
        """

        cfn = os.path.splitext(fn)[0]+'_extra_centered.gro'
        if outfn is None:
            outfn = os.path.splitext(fn)[0]+'.crd'
        U.run('editconf','-f %s -o %s -center 0 0 0'%(fn,cfn))
        f = file(cfn)
        out = file(outfn,'w')
        title = '* ' + next(f) + '* \n'
        out.write(title)
        nat = next(f)
        #out.write(nat)
        nat = int(nat)
        if forceext:
            out.write('%10i  EXT\n'%nat)
        elif nat < 100000:
            out.write('%5i\n'%nat)
        else:
            out.write('%10i\n'%nat)
        if nat >= 1000000 :
            if renumber is None:
                raise Exception('{fn} had too many lines; you must tell gro2crd whether to explicitly renumber or not'.format(fn=fn))
            elif renumber is False:
                print("NOTE: you are explicitly choosing NOT to renumber long file {fn}".format(fn=fn))
            elif renumber is True:
                print("NOTE: you are explicitly choosing to renumber long file {fn}".format(fn=fn))
            else:
                raise Exception("unknown renumber parameter {r}".format(r=renumber))
        lastrawresi,renumberedresi = 0,0
        for i in range(nat):
            line = next(f)
            if renumber:
                parts = groline2parts(line)
                [resi,resn,atomname,atomnumber],otherparts = parts[:4],parts[4:]
                atomnumber = i + 1
                if resi == lastrawresi:
                    pass
                else:
                    lastrawresi = resi
                    renumberedresi = renumberedresi + 1
                resi = renumberedresi
                parts = [resi,resn,atomname,atomnumber] + otherparts
                _line = parts2crdline(parts,segmap=segmap,resnmap=resnmap,nat={True:10000000,False:nat}[forceext],fixlip=fixlip)
            else:
                _line = parts2crdline(groline2parts(line),segmap=segmap,resnmap=resnmap,fixlip=fixlip)
            out.write(_line)
        line = next(f)
        a,b,c = [float(i) for i in line.strip().split()]
        print("Gromacs listed this box as {a} {b} {c}".format(a=a,b=b,c=c))
        print("A CHARMM line might look like")
        print("crystal define ORTHogonal {a:.4f} {b:.4f} {c:.4f} 90. 90. 90.".format(a=a*10.,b=b*10.,c=c*10.))
        print("{a:.4f} {b:.4f} {c:.4f}".format(a=a*10.,b=b*10.,c=c*10.))
        
        f.close()
        if cleanup:
            U.run('rm',cfn)
        
    @staticmethod
    def crd2gro(fn):
        """Converts a CHARMM .crd file to a GROMACS .gro file, running
        it through editconf afterwards to make sure that it's properly
        centered
        
        Arguments:
         - `fn`: input filename. output will have extension replaced with .crd
        """
        raise(NotImplemented)


class NDX:
    """
    IDs will be internally stored as ints to facilitate comparison.
    They'll be converted to strings upon output.
    """
    @staticmethod
    def idxstogrpstr(idxs,grp):
        """Turns a list of atom id numbers into a group suitable for stamping into an index file.
        For convenience, we sort the idxs
        """
        result = []
        result.append('[ %s ]\n'%grp)
        #   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
        #7267 7268 7269 7270 7271 7272 7273 7274 7275 7276 7277 7278 7279 7280 7281 
        for g in U.splitseq(sorted(idxs),15):
            line = ' '.join(['%4s'%i for i in g])
            line = line + '\n'
            result.append(line)
        return result
    @staticmethod
    def addgrpstofile(grps,inf,outf):
        """Adds a list of groups to an index file"""
        f = file(inf)
        inlines = f.readlines()
        f.close()
        outf = file(outf,'w')
        outf.writelines(inlines)
        for grp in grps:
            outf.writelines(grp)
        outf.close()
    def addgrp(self,idxs,grp):
        print("adding group",grp)
        self.groups[grp] = idxs
    def write(self,fname):
        f = file(fname,'w')
        groups = list(self.groups.keys())
        groups.sort()
        groups.reverse()
        for g in groups:
            f.writelines(self.idxstogrpstr(self.groups[g],g))
        f.close()
    def __init__(self,fname=None):
        """Parse an index file"""
        groups = {}
        grp,ids = None,[]
        if fname:
            for line in open(fname):
                line = line.strip()
                if not line: continue
                if line.startswith('['):
                    if ids:
                        groups[grp] = ids
                        ids=[]
                    grp = line.replace('[','').replace(']','').strip()
                else:
                    ids.extend([int(i) for i in line.split()])
            if ids:
                groups[grp] = ids
        self.groups = groups
    def getcompressedidxs(self,allgroups,group):
        """
        Common use case: we have an xtc file that only contains groups
        allgroups. We want the indexes of e.g. the top leaflet, which
        is defined in a group for the whole system.
        """
        allids = []
        for g in allgroups:
            allids.extend(self.groups[g])
        print(len(allids))
        assert len(allids) == len(set(allids))
        allids.sort()
        return [allids.index(i) + 1 for i in self.groups[group]]
         
                

class MDP:
    def __init__(self):
        self.params = {}
    def write(self,fname):
        f = file(fname,'w')
        for k in self.params:
            try:
                f.write(' '.join([k,'=',]+self.params[k]+['\n',]))
            except TypeError:
                print("Trying to concatenate")
                print([k,'=',])
                print(self.params[k])
                print(['\n',])
                raise
        f.close()
    def genseed(self):
        # extra str/int is to get rid of leading zeroes
        self.params['gen_seed'] = [str(int(str(random.uniform())[2:2+6])),]
    def setsimtime(self,
                   dt,
                   simtime,
                   trrcoordfreq=2,
                   trrvelfreq=2,
                   xtcfreq=0.05,
                   logfreq=0.05,
                   energyfreq=0.05,
                   ):
        ''' everything but dt is specified in ns, and internally converted to ps '''
        self.params['dt'] =        [str(dt),]
        self.params['nsteps'] =    [str(int(1000*simtime/dt)),]
        self.params['nstxout'] =   [str(int(1000*trrcoordfreq/dt)),]
        self.params['nstvout'] =   [str(int(1000*trrvelfreq/dt)),]
        self.params['nstxtcout'] = [str(int(1000*xtcfreq/dt)),]
        self.params['nstlog'] =  [str(int(1000*logfreq/dt)),]
        self.params['nstenergy'] = [str(int(1000*energyfreq/dt)),]
        
        

class MartiniInitialMDMDP(MDP):
    ''' short dt 0.010, 10 ns, generates a random seed. '''
    def __init__(self,tc_grps,pcoupltype='semiisotropic',dt=0.010,simtime=510,**kwargs):
        MDP.__init__(self)
        self.setsimtime(dt=dt,simtime=simtime,**kwargs)
        self.params['integrator']               = 'md'.split()
        self.params['tinit']                    = '0.0'.split()
        self.params['comm-mode']                = 'linear'.split()
        self.params['nstcomm']                  = '1'.split()
        self.params['comm-grps']                = 'Prot SOL'.split()
        self.params['nstfout']                  = '0'.split()
        self.params['xtc_precision']            = '100'.split()
        self.params['xtc-grps']                 = 'System'.split()
        self.params['energygrps']               = 'System'.split()
        self.params['nstlist']                  = '5'.split()
        self.params['ns_type']                  = 'grid'.split()
        self.params['pbc']                      = 'xyz'.split()
        self.params['rlist']                    = '1.2'.split()
        self.params['coulombtype']              = 'Shift '.split()
        self.params['rcoulomb_switch']          = '0.0'.split()
        self.params['rcoulomb']                 = '1.2'.split()
        self.params['epsilon_r']                = '15 '.split()
        self.params['vdw_type']                 = 'Shift '.split()
        self.params['rvdw_switch']              = '0.9'.split()
        self.params['rvdw']                     = '1.2'.split()
        self.params['DispCorr']                 = 'No'.split()
        self.params['fourierspacing']           = '0.12'.split()
        self.params['fourier_nx']               = '10'.split()
        self.params['fourier_ny']               = '10'.split()
        self.params['fourier_nz']               = '10'.split()
        self.params['pme_order']                = '4'.split()
        self.params['ewald_rtol']               = '1e-05'.split()
        self.params['epsilon_surface']          = '0'.split()
        self.params['optimize_fft']             = 'no'.split()
        self.params['tcoupl']                   = 'Nose-Hoover'.split()
        self.params['tc-grps']                  = tc_grps
        self.params['tau_t']                    = ['2' for i in tc_grps]
        self.params['ref_t']                    = ['323' for i in tc_grps]
        self.params['Pcoupl']                   = 'Parrinello-Rahman'.split()
        self.params['Pcoupltype']               = [pcoupltype]
        if pcoupltype == 'semiisotropic':
            self.params['tau_p']                    = '5 5'.split()
            self.params['compressibility']          = '5e-5 5e-5'.split()
            self.params['ref_p']                    = '1.0 1.0'.split()
        elif pcoupltype == 'isotropic':
            self.params['tau_p']                    = '5'.split()
            self.params['compressibility']          = '5e-5'.split()
            self.params['ref_p']                    = '1.0'.split()
        else:
            raise Exception('unknown pressure coupling')
        self.params['annealing']                = ['no' for i in tc_grps]
        self.params['continuation']             = 'no'.split()
        self.params['gen_vel']                  = 'yes'.split()
        self.params['gen_temp']                 = '323'.split()
        self.params['constraints']              = 'none '.split()
        self.params['constraint_algorithm']     = 'Lincs'.split()
        self.params['shake_tol']                = '0.0001'.split()
        self.params['lincs_order']              = '4'.split()
        self.params['lincs_warnangle']          = '30'.split()
        self.params['morse']                    = 'no'.split()
        self.params['disre']                    = 'No'.split()
        self.params['disre_weighting']          = 'Equal'.split()
        self.params['disre_mixed']              = 'no'.split()
        self.params['disre_fc']                 = '1000'.split()
        self.params['disre_tau']                = '1.25'.split()
        self.params['nstdisreout']              = '100'.split()

        self.genseed()


class MartiniContinuationMDMDP(MDP):
    ''' normal dt 0.020, 510 ns, does not generate a random seed. '''
    def __init__(self,tc_grps,pcoupltype='semiisotropic',dt=.010,simtime=510):
        MDP.__init__(self)
        self.setsimtime(dt=dt,simtime=simtime,
                        trrcoordfreq=2,
                        trrvelfreq=2,
                        xtcfreq=.05,
                        logfreq=.05,
                        energyfreq=.05,)
        self.params['integrator']               = 'md'.split()
        self.params['tinit']                    = '0.0'.split()
        self.params['comm-mode']                = 'linear'.split()
        self.params['nstcomm']                  = '1'.split()
        self.params['comm-grps']                = 'System'.split()
        self.params['nstfout']                  = '0'.split()
        self.params['xtc_precision']            = '100'.split()
        self.params['xtc-grps']                 = 'Prot'.split()
        self.params['energygrps']               = 'Prot SOL'.split()
        self.params['nstlist']                  = '5'.split()
        self.params['ns_type']                  = 'grid'.split()
        self.params['pbc']                      = 'xyz'.split()
        self.params['rlist']                    = '1.2'.split()
        self.params['coulombtype']              = 'Shift '.split()
        self.params['rcoulomb_switch']          = '0.0'.split()
        self.params['rcoulomb']                 = '1.2'.split()
        self.params['epsilon_r']                = '15 '.split()
        self.params['vdw_type']                 = 'Shift '.split()
        self.params['rvdw_switch']              = '0.9'.split()
        self.params['rvdw']                     = '1.2'.split()
        self.params['DispCorr']                 = 'No'.split()
        self.params['fourierspacing']           = '0.12'.split()
        self.params['fourier_nx']               = '10'.split()
        self.params['fourier_ny']               = '10'.split()
        self.params['fourier_nz']               = '10'.split()
        self.params['pme_order']                = '4'.split()
        self.params['ewald_rtol']               = '1e-05'.split()
        self.params['epsilon_surface']          = '0'.split()
        self.params['optimize_fft']             = 'no'.split()
        self.params['tcoupl']                   = 'Nose-Hoover'.split()
        self.params['tc-grps']                  = tc_grps
        self.params['tau_t']                    = ['2' for i in tc_grps]
        self.params['ref_t']                    = ['323' for i in tc_grps]
        self.params['Pcoupl']                   = 'Parrinello-Rahman'.split()
        self.params['Pcoupltype']               = [pcoupltype]
        if pcoupltype == 'semiisotropic':
            self.params['tau_p']                    = '5 5'.split()
            self.params['compressibility']          = '5e-5 5e-5'.split()
            self.params['ref_p']                    = '1.0 1.0'.split()
        elif pcoupltype == 'isotropic':
            self.params['tau_p']                    = '5'.split()
            self.params['compressibility']          = '5e-5'.split()
            self.params['ref_p']                    = '1.0'.split()
        else:
            raise Exception('unknown pressure coupling')
        self.params['annealing']                = ['no' for i in tc_grps]
        self.params['continuation']             = 'yes'.split()
        self.params['gen_vel']                  = 'no'.split()
        self.params['constraints']              = 'none '.split()
        self.params['constraint_algorithm']     = 'Lincs'.split()
        self.params['shake_tol']                = '0.0001'.split()
        self.params['lincs_order']              = '4'.split()
        self.params['lincs_warnangle']          = '30'.split()
        self.params['morse']                    = 'no'.split()
        self.params['disre']                    = 'No'.split()
        self.params['disre_weighting']          = 'Equal'.split()
        self.params['disre_mixed']              = 'no'.split()
        self.params['disre_fc']                 = '1000'.split()
        self.params['disre_tau']                = '1.25'.split()
        self.params['nstdisreout']              = '100'.split()

class MartiniBilayerContinuationMDMDP(MartiniContinuationMDMDP):
    ''' uses different COM groups for Prot, SOL '''
    def __init__(self,tc_grps,pcoupltype='semiisotropic',dt=0.010,simtime=510):
        #ASDFASDF MGL HERE
        MartiniContinuationMDMDP.__init__(self,tc_grps=tc_grps,pcoupltype=pcoupltype,dt=dt,simtime=simtime)
        #self.params['comm-grps']                = 'top bottom SOL'.split()
        self.params['comm-grps']                = 'Prot SOL'.split()
        #self.setsimtime(dt=dt,simtime=simtime)

class MartiniBilayerInitialMDMDP(MartiniInitialMDMDP):
    ''' uses different COM groups for Prot, SOL '''
    def __init__(self,tc_grps,pcoupltype='semiisotropic',dt=0.010,simtime=10,):
        MartiniInitialMDMDP.__init__(self,tc_grps=tc_grps,pcoupltype=pcoupltype,dt=dt,simtime=simtime)
        self.params['comm-grps']                = 'Prot SOL'.split()
        #self.setsimtime(dt=dt,simtime=simtime)

class MartiniUnRestrainedMiniMDP(MDP):
    def __init__(self,tc_grps,nsteps='500'):
        MDP.__init__(self)
        self.params['nsteps']                   = nsteps.split()
        self.params['title']                    = 'Yo'.split()
        self.params['cpp']                      = '/lib/cpp'.split()
        self.params['integrator']               = 'steep'.split()
        self.params['tinit']                    = '500.0'.split()
        self.params['dt']                       = '0.010'.split()
        self.params['nstcomm']                  = '0'.split()
        self.params['nstxout']                  = '5000'.split()
        self.params['nstvout']                  = '5000'.split()
        self.params['nstfout']                  = '0'.split()
        self.params['nstlog']                   = '100'.split()
        self.params['nstenergy']                = '100'.split()
        self.params['nstxtcout']                = '1000'.split()
        self.params['xtc_precision']            = '100'.split()
        self.params['xtc-grps']                 = 'System'.split()
        self.params['energygrps']               = 'System'.split()
        self.params['nstlist']                  = '5'.split()
        self.params['ns_type']                  = 'grid'.split()
        self.params['pbc']                      = 'xyz'.split()
        self.params['rlist']                    = '1.2'.split()
        self.params['domain-decomposition']     = 'no'.split()
        self.params['coulombtype']              = 'Shift '.split()
        self.params['rcoulomb_switch']          = '0.0'.split()
        self.params['rcoulomb']                 = '1.2'.split()
        self.params['epsilon_r']                = '15 '.split()
        self.params['vdw_type']                 = 'Shift '.split()
        self.params['rvdw_switch']              = '0.9'.split()
        self.params['rvdw']                     = '1.2'.split()
        self.params['DispCorr']                 = 'No'.split()
        self.params['fourierspacing']           = '0.12'.split()
        self.params['fourier_nx']               = '10'.split()
        self.params['fourier_ny']               = '10'.split()
        self.params['fourier_nz']               = '10'.split()
        self.params['pme_order']                = '4'.split()
        self.params['ewald_rtol']               = '1e-05'.split()
        self.params['epsilon_surface']          = '0'.split()
        self.params['optimize_fft']             = 'no'.split()
        self.params['tcoupl']                   = 'Berendsen'.split()
        self.params['tc-grps']                  = tc_grps
        self.params['tau_t']                    = ['2' for i in tc_grps]
        self.params['ref_t']                    = ['323' for i in tc_grps]
        self.params['Pcoupl']                   = 'berendsen'.split()
        self.params['Pcoupltype']               = 'isotropic'.split()
        self.params['tau_p']                    = '1.20'.split()
        self.params['compressibility']          = '5e-5'.split()
        self.params['ref_p']                    = '1.0 '.split()
        self.params['annealing']                = ['no' for i in tc_grps]
        self.params['gen_vel']                  = 'yes'.split()
        self.params['gen_temp']                 = '300'.split()
        self.params['gen_seed']                 = '473529'.split()
        self.params['constraints']              = 'none '.split()
        self.params['constraint_algorithm']     = 'Lincs'.split()
        self.params['unconstrained_start']      = 'yes'.split()
        self.params['shake_tol']                = '0.0001'.split()
        self.params['lincs_order']              = '4'.split()
        self.params['lincs_warnangle']          = '30'.split()
        self.params['morse']                    = 'no'.split()
        self.params['disre']                    = 'No'.split()
        self.params['disre_weighting']          = 'Equal'.split()
        self.params['disre_mixed']              = 'no'.split()
        self.params['disre_fc']                 = '1000'.split()
        self.params['disre_tau']                = '1.25'.split()
        self.params['nstdisreout']              = '100'.split()
        self.params['free_energy']              = 'no'.split()
        self.params['init_lambda']              = '0'.split()
        self.params['delta_lambda']             = '0'.split()
        self.params['sc-alpha']                 = '0'.split()
        self.params['sc-sigma']                 = '0.3'.split()
        self.params['cos-acceleration']         = '0'.split()


class MartiniRestrainedMiniMDP(MartiniUnRestrainedMiniMDP):
    def __init__(self,tc_grps,nsteps='500'):
        MartiniUnRestrainedMiniMDP.__init__(self,tc_grps=tc_grps,nsteps=nsteps)
        self.params['define']                   = '-DPOSRES'.split()

def isitpsectionline(line):
    line = line.strip().replace('[','[ ').replace(']',' ]')
    parts = line.split()
    if (len(parts) == 3) and (parts[0]+parts[2] == '[]'):
        return parts[1]
    return False
def makeitpsectionline(sec):
    result = '[ %s ]\n'%sec
    if sec == 'header': result = ';' + result
    return result
class ITP:
    @staticmethod
    def bymolecules(f,needtoskipheader=False):
        # We assume the file is only moleculetypes at this point.  If
        # it's not, you'll get all of the leading stuff in the first
        # molecule. This is done so that we can call this after we've
        # already consumed the first 'moleculetype' line.
        if needtoskipheader:
            for line in f:
                s = isitpsectionline(line)
                if s == 'moleculetype': break
        mol = {'moleculetype':[]}
        molpart = 'moleculetype'
        for line in f:
            s = isitpsectionline(line)
            if s == 'moleculetype':
                if mol: yield mol
                mol = {'moleculetype':[]}
                molpart = s
            elif s:
                molpart = s
                mol[molpart] = []
            else:
                mol[molpart].append(line)
        yield mol
class SingleMoleculeITP(ITP):
    """
    At the moment, this only allows one [moleculetype] per ITP
    file. This should obviously be fixed, but you can work around it
    by just making multiple ITP files for now.
    """
    def __init__(self,fname):
        # We use self.orderedparts() to get things in the correct
        # order. In later versions of Python, we could use an ordered
        # dict.
        #

        # NOTE: exclusions not handled properly yet.
        f = file(fname)
        self.parts = {}
        self.setpart('header')
        for line in f:
            s = isitpsectionline(line) 
            if s:
                self.setpart(s)
                print("Now parsing",s,line.strip())
            else:
                self.addline(line)
        lines = [i for i in self.parts['atoms'] if (i.strip() and not i.strip().startswith(';'))]
        self.natoms = len(lines)
        assert self.natoms == int(lines[-1].split()[0])
        self.nres = int(lines[-1].split()[2])
        self.incr(0,0) # normalize the lines

    def orderedparts(self):
        keys = list(self.parts.keys())
        required = 'header moleculetype atoms bonds angles'.split()
        for r in required: 
            assert r in keys
            keys.remove(r)
        atomdefs = 'atomtypes defaults nonbond_params'.split()
        _adefs = []
        for a in atomdefs:
            # If we have these, we're defining new atom types. Those
            # need to come at the beginning, i.e. before we actually
            # try to use the atoms.
            if a in keys:
                keys.remove(a)
                _adefs.append(a)
        specialkeys = required[:1] + _adefs + required[1:] # header still comes first
        keys.sort()
        return specialkeys + keys
            
    def setpart(self,p):
        if p in self.parts:
            print("duplicate for",p)
            sys.exit()
        self.parts[p] = []
        self.currpart = p
    def addline(self,line):
        self.parts[self.currpart].append(line)
    def incr(self,nat,nres):
        """increment all of the atom numbers by nat and all of the residue numbers by nres. None means do not increment"""
        def incrline(line,idxs,incat,incres=None):
            if line.strip().startswith(';'): return line
            if idxs is None: return line
            if not line.strip(): return line
            parts = line.split()
            for idx in idxs:
                try:
                    parts[idx] = str(int(parts[idx]) + incat)
                except IndexError:
                    print("Could not find idx",idx,"in parts",parts)
                    raise
            if incres:
                parts[2] = str(int(parts[2]) + incres)
            return '\t'.join(parts) + '\n'
        # Atoms increment the residue number, atom number, and cgnr.
        # This lists which columns to increment and how much to
        # increment them by. If None is listed, do not in
        incby = {
            'header':([],nat),
            'moleculetype':([],nat),
            'atoms':([0,5],nat,nres),
            'bonds':([0,1],nat),
            'constraints':([0,1],nat),
            'angles':([0,1,2],nat),
            'angle_restraints_z':([0,1],nat),
            'dihedrals':([0,1,2,3],nat),

            # These are atom type definitions, not molecule
            # definitions.
            'atomtypes':(None,None),
            'defaults':(None,None),
            'nonbond_params':(None,None),
            }
        #_parts = []
        for p in self.parts:
            p1 = self.parts[p]
            self.parts[p] = [incrline(line,*incby[p]) for line in self.parts[p]]
            #_parts.append([p0] + [incrline(line,*incby[p0]) for line in p1])
        #self.parts = _parts
        self.natoms = self.natoms + nat
        self.nres = self.nres + nres
    def __add__(self,other):
        nat = self.natoms
        nres = self.nres
        other.incr(nat,nres)
        for p in self.parts:
            if p == 'moleculetype': continue
            self.parts[p] = self.parts[p] + other.parts[p]
        self.natoms = other.natoms
        self.nres = other.nres
        other.incr(-nat,-nres)
        return self
    def setmoltype(self,moltype):
        #print self.parts.keys()
        self.parts['moleculetype'] = ['; molname 	nrexcl\n','\t'.join([moltype,'1']) + '\n','\n']
    def moltype(self):
        lines = [i for i in self.parts['moleculetype'] if (i.strip() and not i.strip().startswith(';'))]
        assert len(lines) == 1
        return lines[0].split()[0]

    def _gettxt(self,includeblanklines):
        result = ''
        for part in self.orderedparts():
            result += makeitpsectionline(part)
            for i in self.parts[part]:
                if (includeblanklines or (i.strip() and not i.strip().startswith(';'))):
                    result += i
            if not includeblanklines:
                result += '\n' # still must have a blank line between
        return result
    def printbasic(self):
        print(self._gettxt(includeblanklines=False))
    def printall(self):
        print(self._gettxt(includeblanklines=True))
    def write(self,fn,includeblanklines=True):
        f = file(fn,'w')
        f.write(self._gettxt(includeblanklines=includeblanklines))
        f.close()
    def __str__(self):
        result = '< ITP[%s]: ' % self.moltype()
        for p in self.parts:
            result += p + '(%s) '%len([i for i in self.parts[p] if (i.strip() and not i.strip().startswith(';'))])
        result += '>'
        return result
    __repr__ = __str__
        
class GeneralITP(ITP):
    """ Can't add or subtract things, but does parse multiple molecules """
    def __init__(self,fname):
        # We use self.orderedparts() to get things in the correct
        # order. In later versions of Python, we could use an ordered
        # dict.
        #

        # NOTE: exclusions not handled properly yet.
        f = file(fname)
        self.fname = fname
        self.parts = {}
        self.molecules = []
        self.setpart('header')
        # process everything up to the molecules
        for line in f:
            s = isitpsectionline(line)
            if s == 'moleculetype':
                break
            if s:
                self.setpart(s)
                print("Now parsing",s,line.strip())
            else:
                self.addline(line)


        if s == 'moleculetype':
            # We can now expect the rest of the file to be individual molecule definitions.
            for mol in ITP.bymolecules(f):
                self.molecules.append(mol)

    def orderedparts(self):
        keys = list(self.parts.keys())
        required = 'header moleculetype atoms bonds angles'.split()
        for r in required: 
            assert r in keys
            keys.remove(r)
        atomdefs = 'atomtypes defaults nonbond_params'.split()
        _adefs = []
        for a in atomdefs:
            # If we have these, we're defining new atom types. Those
            # need to come at the beginning, i.e. before we actually
            # try to use the atoms.
            if a in keys:
                keys.remove(a)
                _adefs.append(a)
        specialkeys = required[:1] + _adefs + required[1:] # header still comes first
        keys.sort()
        return specialkeys + keys
            
    def setpart(self,p):
        if p in self.parts:
            print("duplicate for",p)
            sys.exit()
        self.parts[p] = []
        self.currpart = p
    def addline(self,line):
        self.parts[self.currpart].append(line)

    def _gettxt(self,includeblanklines):
        result = ''
        for part in self.orderedparts():
            result += makeitpsectionline(part)
            for i in self.parts[part]:
                if (includeblanklines or (i.strip() and not i.strip().startswith(';'))):
                    result += i
            if not includeblanklines:
                result += '\n' # still must have a blank line between
        return result
    def printbasic(self):
        print(self._gettxt(includeblanklines=False))
    def printall(self):
        print(self._gettxt(includeblanklines=True))
    def write(self,fn,includeblanklines=True):
        f = file(fn,'w')
        f.write(self._gettxt(includeblanklines=includeblanklines))
        f.close()
    def __str__(self):
        result = '< ITP[%s]: ' % self.fname
        for p in self.parts:
            result += p + '(%s) '%len([i for i in self.parts[p] if (i.strip() and not i.strip().startswith(';'))])
            result += 'molecules(%s) '%len(self.molecules)
        result += '>'
        return result
    __repr__ = __str__
                


#######
# XVG related
#######

def plotdata(x,y=None,yerr=None,numblocks=10,numerrbars=50,colors=('blue','green'),cumavg=True,clear=False,
             plottitle=None,xaxislabel=None,yaxislabel=None,axes=None,grid=False):
    """ Our standard plotting routine.

     - The first argument is the data. We try to extract x and y data from it if possible.
     - Plot 50 error bars if we're given errors.
     - Do block averaging.
     - Plot the cumulative average.
     """
    if axes is None:
        raise Error('Need axes instance')
    txt = ''
    # Figure out x,y,err
    if y is None:
        if len(x.shape) > 1:
            x,y = x[:,0],x[:,1]
        else:
            x,y = array(list(range(len(x)))),x
    if clear: clf()

    annotation_location = (min(x) + (max(x) - min(x))*0.1,min(y) + (max(y) - min(y))*0.9)
    #print annotation_location,max(y)
    axes.plot(x,y,color=colors[0],zorder=10,alpha=0.5)

    if cumavg:
        ya = np.cumsum(y)/np.arange(1,len(y)+1)
        axes.plot(x,ya,'k-',zorder=40)

    if yerr is not None:
        error_step = int(len(x)/numerrbars)
        if error_step == 0: error_step = len(x)
        axes.errorbar(x[::error_step],y[::error_step],yerr[::error_step],color=colors[0],zorder=20)
    blocksize = int(len(x)/numblocks)
    blocksizes = [len(i) for i in U.splitseq(x,blocksize)]
    ravg = average(y)
    rstdev = std(y)
    txt += 'Raw Avg: {ra:.4f} Raw std dev: {rs:.4f}\n'.format(ra=ravg,rs=rstdev)

    bx = [average(i) for i in U.splitseq(x,blocksize)][:numblocks]
    byerr = array([std(i) for i in U.splitseq(y,blocksize)])[:numblocks]
    by = array([average(i) for i in U.splitseq(y,blocksize)])[:numblocks]
    txt += 'Block Avg: %.4f'%average(by)
    txt += ', std err: %.4f, avg err: %.4f'%(std(by)/sqrt(numblocks),average(byerr))
    txt += '\nblocks of length %s'%blocksize
    if len(blocksizes) > numblocks:
        txt += ' discarding %s points at the end'%sum(blocksizes[numblocks:])
    axes.errorbar(bx,by,byerr,elinewidth=20,color=colors[1],barsabove=True,zorder=30)
    if plottitle: axes.set_title(plottitle)
    if xaxislabel: axes.set_xlabel(xaxislabel)
    if yaxislabel: axes.set_ylabel(yaxislabel)
    axes.annotate(txt,annotation_location)
    if grid:
        plt.grid(True)

class XVG:
    def __init__(self,fname,startdir='.'):
        self.fname = os.path.join(startdir,fname)
        self.attributes = {'comments':''}
        self.columns = ['time']
        self.units = {'time':'ns'}
        self.process_file()
        self.startframe = 0
        self.stopframe = None
    def process_file(self):
        print("processing",self.fname)
        f = open(self.fname)
        lines = []
        for line in f:
            if line.startswith('#'):
                self.attributes['comments'] += line
            elif line.startswith('@'):
                self.add_info(line)
            else:
                lines.append([float(i) for i in line.split()])
        f.close()
        if lines and (len(lines[-1]) != len(lines[-2])):
            print("Deleting bad last line :",lines[-1])
            print("Expected something like:",lines[-2])
            lines = lines[:-1]
        self.data = array(lines)
        if len(self.data) > 0:
            self.data[::,0] = self.data[::,0]/1000. # convert to nanoseconds
    def add_info(self,line):
        line = line[1:].strip()
        info_type = line.split()[0]
        rest = line[line.index(info_type)+len(info_type):].strip()
        if info_type in 'title TYPE'.split():
            self.attributes[info_type] = rest.replace('"','').strip()
        elif info_type in 'xaxis yaxis'.split():
            self.attributes[info_type] = rest.replace('label','').replace('"','').strip()
            if info_type == 'xaxis' and self.attributes[info_type] == 'Time (ps)':
                self.attributes[info_type] = 'Time (ns)'
                
            if info_type == 'yaxis':
                self.defaultunit = rest.split('(')[-1].replace(')','').replace('"','').strip()
        elif re.match('s\d+$',info_type):
            #@ s1 legend "Coulomb (SR)"
            rest = rest.split('"')[-2] #Coulomb (SR)
            if '(' in rest:
                label = rest.split('(')[0].strip()
                unit = rest.split('(')[-1].replace(')','').strip()
                self.columns.append(label)
                self.units[label] = unit
            else:
                self.columns.append(rest)
                self.units[rest] = self.defaultunit#'kJ mol\S-1\N'
    def get(self,item):
        try:
            result = self.data[:,self.columns.index(item)]
        except ValueError:
            # Special cases
            if item == 'Area':
                # Don't use self.get here, because it can call through
                # twice, causing you to lop off startframe frames
                # twice
                result = self.data[:,self.columns.index('Box-X')] * self.data[:,self.columns.index('Box-Y')]
                if item not in self.units: self.units[item] = 'nm^2'
            else:
                raise
        if self.stopframe is not None:
            return result[self.startframe:self.stopframe]
        else:
            return result[self.startframe:]
            
    def getxy(self,item):
        if self.stopframe is not None:
            return self.data[self.startframe:self.stopframe,0],self.get(item)
        else:
            return self.data[self.startframe:,0],self.get(item)
    def getKA(self):
        # Ka = A*k*T/sigma_A^2
        #  A is the area per lipid
        #  k is Boltzmann's constant
        #  T is the temperature
        #  sigma^2 is the variance of the area
        
        # We'd like dyn/cm units in the 100s
        # Area comes in nanometers, so multiply by 1e-7 for cgs
        k = 1.3806488e-16 # dyn/cm, cgs
        area = self.get('Area')
        a = average(area) 
        s2 = std(area)**2
        T = average(self.get('Temperature'))
        katot = 1e14 * (a*k*T)/(s2) # in dyn/nm --> dyn/cm
        numblocks = 10
        blocksize = int(len(area)/numblocks)
        blocksizes = [len(i) for i in U.splitseq(area,blocksize)]
        blocks = U.splitseq(area,blocksize)[:numblocks]
        kas = [1e14*(average(block)*k*T)/std(block)**2 for block in blocks]
        return katot,kas
    def printKA(self):
        katot,kas = self.getKA()
        print('K_A = {ka:f}'.format(ka=katot))
        print('K_A blocks: ' + ' '.join(['{ka:f}'.format(ka=ka) for ka in kas]))
        print('K_A block avg: {ka:f} K_A block std dev {ks:f}'.format(ka=average(kas),ks=std(kas)))
    def plot(self,axes,item,color='',truncat=None,
             numblocks=10,deavg=False,grid=False):
        x,y=self.getxy(item)
        if truncat:
            x = x[:truncat]
            y = y[:truncat]
        if deavg:
            print("Subtracting averages y = %s"%(average(y),))
            y = y - average(y)
        plotdata(x,y,numblocks=numblocks,xaxislabel=self.attributes['xaxis'],yaxislabel='%s (%s)'%(item,self.units[item]),plottitle=item+' ('+os.path.split(self.fname)[-1]+')',axes=axes, grid=grid)

    def plotitem(self,idx,axes):
        self.plot(axes,self.columns[idx])
    def plotall(self,figure,ncol=4,tstart=None,skipTime=False):
        cols = self.columns[:]
        if skipTime: cols.remove('time')
        self.plotseveral(figure,cols,ncol)
    def plotseveral(self,figure,groups,ncol,grid=False):
        nrow = np.ceil(len(groups)/ncol)
        for (i,item) in enumerate(groups):
            #print "Plotting",i,"out of",len(groups),"in subplot",(nrow,ncol,i+1)
            axes = figure.add_subplot(nrow,ncol,i+1)
            print("Plotting item",item)
            self.plot(axes,item,grid=grid)
