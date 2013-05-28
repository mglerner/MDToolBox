#!/usr/bin/env python

from __future__ import division
import os,sys
import numpy as np
from numpy import random
import util as U

class GromacsError(StandardError):
    pass

######
## Pure Python version of GROMACS utils
######

class Conversion(object):
    """Convert from, e.g. GROMACS to CHARMM. Most methods will be static.
    """
    
    def __init__(self, ):
        """Most methods will be static.
        """
        
        pass
    @staticmethod
    def gro2crd(fn):
        """Converts a .gro file to a CHARMM .crd file, first running
        it through editconf to make sure that it's properly centered.
        
        Arguments:
        - `fn`: input file name. output will have extension replaced with .crd
        """
        cfn = os.path.splitext(fn)[0]+'_centered.gro'
        ofn = os.path.splitext(fn)[0]+'.crd'
        U.run('editconf','-f %s -o %s -center 0 0 0'%(fn,cfn))
        f = file(cfn)
        out = file(ofn,'w')
        title = '* ' + f.next() + '* \n'
        out.write(title)
        nat = f.next()
        out.write(nat)
        nat = int(nat)
        for i in range(nat):
            line = f.next()
            _line = parts2crdline(groline2parts(line))
            out.write(_line)
        line = f.next()
        print "Gromacs listed this box",[10*float(i) for i in line.strip().split()]
        
        



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
        print "adding group",grp
        self.groups[grp] = idxs
    def write(self,fname):
        f = file(fname,'w')
        groups = self.groups.keys()
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
        print len(allids)
        assert len(allids) == len(set(allids))
        allids.sort()
        return [allids.index(i) + 1 for i in self.groups[group]]
         
                

# Since this is PyGRO, everything will be stored internally in GROMACS units.
# That means multiply by 10 for CHARMM coordinates, etc.
def parts2groline(parts):
    # resi,resn,atomname,atomnumber,x,y,z,vx,vy,vz
    if len(parts) == 10:
        return "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n"%tuple(parts)
    elif len(parts) == 7:
        return "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n"%tuple(parts)
    else:
        raise Exception('Unknown number of parts')
def parts2crdline(parts,segid='L'):
    """
    This will multiply coordinates by 10.
    """
    # resi,resn,atomname,atomnumber,x,y,z,vx,vy,vz
    if len(parts) == 10:
        parts = parts[:7]
    elif len(parts) == 7:
        pass
    else:
        raise Exception('Unknown number of parts')
    resi,resn,atomname,atomnumber,x,y,z = parts
    x,y,z = 10*x,10*y,10*z
    resn,atomname = resn.strip(),atomname.strip()
    '    1    1 DPPC N    -19.21035 -10.18234  21.22432 L    1      0.00000'
    'ATOMNO RESNO   RES  TYPE  X     Y     Z   SEGID RESID Weighting'
    'I5    I5  1X A4 1X A4 F10.5 F10.5 F10.5 1X A4 1X A4 F10.5'
    return '%5i%5i %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4s%10.5f\n'%(atomnumber,resi,resn,atomname,x,y,z,segid,resi,0.0)

#    2    1 DPPC C13  -20.45601 -10.97580  21.35871 L    1      0.00000 (good)
#1111122222 3333 4444555555555566666666667777777777 8888 99990000000000
#    1    1 DPP  NC3   -7.65000  -3.35000  11.48000 L    1      0.00000 (good - us)

def groline2parts(line):
    # resi,resn,atomname,atomnumber,x,y,z,vx,vy,vz
    # %5d%5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f
    #  5  5  5  5  8    8    8    [8   8    8]
    resi,resn,atomname,atomnumber,x,y,z = ( int(line[0:5]),
                                                line[5:10],
                                                line[10:15],
                                            int(line[15:20]),
                                           float(line[20:28]),
                                           float(line[28:36]),
                                           float(line[36:44]),)
    if len(line) > 45:
        vx,vy,vz = (float(line[44:52]),
                    float(line[52:60]),
                    float(line[62:68]),)
        return [resi,resn,atomname,atomnumber,x,y,z,vx,vy,vz]
    else:
        return [resi,resn,atomname,atomnumber,x,y,z]
def incresi(groline,i=1):
    parts = groline2parts(groline)
    parts[0] += i
    return parts2groline(parts)

class MDP:
    def __init__(self):
        self.params = {}
    def write(self,fname):
        f = file(fname,'w')
        for k in self.params:
            try:
                f.write(' '.join([k,'=',]+self.params[k]+['\n',]))
            except TypeError:
                print "Trying to concatenate"
                print [k,'=',]
                print self.params[k]
                print ['\n',]
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
                print "Now parsing",s,line.strip()
            else:
                self.addline(line)
        lines = [i for i in self.parts['atoms'] if (i.strip() and not i.strip().startswith(';'))]
        self.natoms = len(lines)
        assert self.natoms == int(lines[-1].split()[0])
        self.nres = int(lines[-1].split()[2])
        self.incr(0,0) # normalize the lines

    def orderedparts(self):
        keys = self.parts.keys()
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
            print "duplicate for",p
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
                    print "Could not find idx",idx,"in parts",parts
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
        print self._gettxt(includeblanklines=False)
    def printall(self):
        print self._gettxt(includeblanklines=True)
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
                print "Now parsing",s,line.strip()
            else:
                self.addline(line)


        if s == 'moleculetype':
            # We can now expect the rest of the file to be individual molecule definitions.
            for mol in ITP.bymolecules(f):
                self.molecules.append(mol)

    def orderedparts(self):
        keys = self.parts.keys()
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
            print "duplicate for",p
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
        print self._gettxt(includeblanklines=False)
    def printall(self):
        print self._gettxt(includeblanklines=True)
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
                


######
### Using PyMACS
######

try:
    import pymacs
    import numpy
    class PYM:
        @staticmethod
        def getboxes(fname,maxframes=None):
            fp = pymacs.openXTC(fname)
            framenum = 0
            boxes = []
            while 1:
                f = pymacs.readXTCFrame(fp)
                if not f: break
                framenum += 1
                if maxframes and framenum > maxframes: break
                boxes.append(f['box'])
            b = numpy.array(boxes)
            return b[:,0,0],b[:,1,1],b[:,2,2]
            
except ImportError:
    print "No PyMACS!"

