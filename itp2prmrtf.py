#!/usr/bin/env python
from __future__ import division
import os
from pygro import GeneralITP,GromacsError

def convert(itpname,overwrite=True):
    # defaults
    # atomtypes
    rtfname = os.path.splitext(itpname)[0] + '.rtf'
    prmname = os.path.splitext(itpname)[0] + '.prm'
    noename = os.path.splitext(itpname)[0] + '.noe'
    noe = """* NOE for %s
* This mimics GROMACS bond type 6, and is used for pure
* harmonic restraints. You can, of course, do more
* complicated things here if you want.
*
NOE
    RESEt
"""
    itp = GeneralITP(itpname)
    prm = """* PRM for %s
*
"""%os.path.splitext(os.path.split(itpname)[-1])[0]
    rtf = """* RTF for %s
*
   33 1

"""%os.path.splitext(os.path.split(itpname)[-1])[0]
    massidx = 1
    for line in itp.parts['atomtypes']:
        line = line.replace(';','!') # change to CHARMM comments
        if not line.strip() or line.strip().startswith('!'):
            rtf += line
            continue
        name, mass, charge, ptype, c6, c12 = line.split()
        if ptype != 'A' :
            raise GromacsError('Do not know now to deal with ptype S, V, or D (%s)'%ptype)
        if float(c6) != 0:
            raise GromacsError('For CG FFs, I expect c6 0 here (%s)'%c6)
        if float(c12) != 0:
            raise GromacsError('For CG FFs, I expect c12 0 here (%s)'%c12)
        rtf += 'MASS %-5s %-4s %10.6f %-5s\n'%(massidx,name,float(mass),name) # I guess each one is its own particle type/element.
        massidx += 1

    ### nonbond_params
    prm += """!TIM: NONBONDED NBXMOD 2 ATOM CDIEL SHIFT VATOM VDISTANCE VSHIFT -
!TIM: CUTNB 16 CTOFNB 12 EPS 1.0 WMIN 1.5 E14FAC 0.7
NONBONDED NBXMod 2 -
! GROMACS: coulombtype = Shift, rcoulomb_switch=0.0 rcoulomb = 1.2
CUTNB 14 -
CTOFNB 12 -
! ELEC
ATOM CDIEl FSHIft EPS 15 -
! GROMACS: vdw_type = Shift, rvdw_switch=0.9, rvdw=1.2
! VDW
VATOM VSWItch CTONNB 9.0
!V(Lennard-Jones) = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]
!
!epsilon: kcal/mole, Eps,i,j = sqrt(eps,i * eps,j)
!Rmin/2: A, Rmin,i,j = Rmin/2,i + Rmin/2,j
!
!atom  ignored    epsilon      Rmin/2   ignored   eps,1-4       Rmin/2,1-4
"""
    nbfixlines = ''
    for line in itp.parts['nonbond_params']:
        line = line.replace(';','!') # change to CHARMM comments
        if not line.strip() or line.strip().startswith('!'):
            rtf += line
            continue
        a1,a2,funda,c6,c12 = line.split()[:5]
        c6,c12 = float(c6),float(c12)
        # First, get sigma and epsilon from c6,c12
        # We assume combination rule 1
        sigma = (c12/c6)**(1/6)
        epsilon = c6*c6/(4*c12)
        epsilon = epsilon/-4.184
        rminover  =     10*sigma*2**(1/6)
        rminover2 = 0.5*10*sigma*2**(1/6)
        if '!' in line:
            comment = ' '+line[line.index('!'):].strip()
        else:
            comment = ''
        if a1 == a2:
            prm += '%-4s 0.0 %7.5f %7.5f%s\n'%(a1,epsilon,rminover2,comment)
        else:
            nbfixlines += '%-4s %-4s %7.5f %7.5f%s\n'%(a1,a2,epsilon,rminover,comment) #Grr, not rmin/2, unlike NBONd
    if nbfixlines:
        prm += '! You must allow at least %s nbfix lines when compiling CHARMM (demens.fcm).\n'%len(nbfixlines.split('\n'))
        prm += 'NBFIX\n'
        prm += nbfixlines

    ### moleculetype
        
    ## Well isn't this fun. Many parts of the CG model actually do not
    ## have angles. We could take our chances and define dummy 0.00
    ## angle parameters for those, but what if that particular
    ## combination of particle types wants an angle later on? Looks
    ## like we have to explicitly write out the angles. Whee.

    ##rtf += 'AUTOgenerate ANGLe DIHEdral\n\n'
        

    # the sanity checking is explained below, but these dictionaries
    # are used to check all bonds and angles from all molecules
    bonds = {}
    angles = {}
    def removecomments(line):
        if ';' not in line: return line
        return line[:line.index(';')]
    def bycgnr(atomlines):
        currgrp = -1
        currlines = []
        for line in atomlines:
            if not removecomments(line).strip(): continue
            grp = int(line.split()[5])
            if grp != currgrp:
                if currlines:
                    yield currlines
                currgrp = grp
                currlines = [line]
            else:
                currlines.append(line)
        if currlines:
            yield currlines
    bondtxt = 'BONDs\n'
    angletxt = 'ANGLes\n'
    for mol in itp.molecules:
        mollines = [line for line in mol['moleculetype'] if removecomments(line).strip()]
        if len(mollines) != 1:
            print mol
            raise GromacsError("Did not find exactly one line in moleculetype: %s,%s"%(mollines,mol['moleculetype']))
        molname,nrexcl = mollines[0].split()
        print "Processing",molname
        nrexcl = int(nrexcl)
        if nrexcl != 1:
            raise GromacsError("MARTINI should have nrexcl 1. This is why the CHARMM .prm file says NBXMOD 2 (%s)"%(mollines[0]))

        ###  atoms

        # we accumulate the atomtypes map here for use in the bonds
        # section later.
        atomtypemap = {}
        atomnamemap = {}
        chg = sum(float(removecomments(line).split()[6]) for line in mol['atoms'] if removecomments(line).strip())
        "RESI LPPC        0.00 ! deoxylysophosphatidylcholine"
        rtf += 'RESI %-4s        %.2f\n'%(molname,chg)
        for lines in bycgnr(mol['atoms']):
            rtf += 'GROUP\n'
            for line in lines:
                "ATOM C13  CTL5  -0.35 !      H15A-C15-H15C"
                try:
                    atomid,atomtype,resnr,residu,atomname,cgnr,charge = removecomments(line).split()
                except ValueError:
                    raise GromacsError("Trouble unpacking %s\nas [atomid,atomtype,resnr,residu,atomname,cgnr,charge]"%removecomments(line).split())
                if int(resnr) != 1:
                    print 'Multiple residue groups not yet implemented (%s)'%line.strip()
                    print 'Setting residue to 1, perhaps breaking',molname
                atomtypemap[int(atomid)] = atomtype
                atomnamemap[int(atomid)] = atomname
                rtf += 'ATOM %-4s %5s %5.2f\n'%(atomname,atomtype,float(charge))

        ###  bonds

        # one of the differences between CHARMM and GROMACS is that CHARMM
        # defines bonds between atom types in the .prm file, while GROMACS
        # defines bonds on a per-residue basis. So, back to the prm. We'll
        # do lots of sanity checking here. atomtypemap is created as we loop
        # through the atom lines above.
        if 'bonds' in mol:
            for line in mol['bonds']:
                """
                CHARMM:
                !
                !V(bond) = Kb(b - b0)**2
                !
                !Kb: kcal/mole/A**2
                !b0: A
                !
                !atom type Kb          b0

                GROMACS:
                for function type 1 (see section 4.2.1 of the GROMACS manual and table 5.4,
                ; V(bond) = (1/2) Kb(b - b0)**2

                so we will need to multiply by b0 by 10

                and Kb gets divided by 4.184 to convert to Joules, 100 to
                convert from nm**2 to A**2, and 2 to handle the leading
                (1/2). The total conversion factor is then division by 836.8.

                NOTE:
                http://www.ks.uiuc.edu/Research/namd/mailing_list/namd-l/1358.html
                seems to be off .. it says 824.8
                """
                if not removecomments(line).strip(): continue
                i,j,functype,length,forceconst = removecomments(line).split()
                i,j,functype,length,forceconst = int(i),int(j),int(functype),float(length),float(forceconst)
                if functype == 1:
                    length = length*10.
                    forceconst = forceconst / 836.8
                    atoms = [atomtypemap[a] for a in (i,j)]
                    atoms = tuple(sorted(atoms)) # Must be sorted so that we don't
                                                 # double count later, tuple so we
                                                 # can stuff it in a dictionary.
                    if atoms not in bonds:
                        'CTL2  CL    200.0       1.522   ! methyl acetate'
                        bondtxt += '%-5s %-5s %8.2f %8.3f\n'%(atoms[0],atoms[1],forceconst,length)
                        bonds[atoms] = (forceconst,length)
                    else:
                        if bonds[atoms] != (forceconst,length):
                            raise GromacsError("Multiple definitions for bond %s <--> %s : %s vs. %s\nThis line was %s"%(atoms[0],
                                                                                                       atoms[1],
                                                                                                       (forceconst,length),
                                                                                                       bonds[atoms],line))

                    rtf += 'BOND %-4s %-4s\n'%(atomnamemap[i],atomnamemap[j])
                elif functype == 6:
                    # gromacs type 6 is exactly the same as type 1, a
                    # harmonic potential.  CHARMM lets us add this as
                    # a constraint, but note that the definition is
                    # k*(x - x0)**2, rather than 1/2(k) ...
                    #
                    # We cheat and use NOE stuff here, assuming r will
                    # never get close to rlim.
                    #
                    #         /  0.5*KMIN*(RAVE-RMIN)**2    R<RMIN
                    #        /
                    #       /    0.0                        RMIN<R<RMAX
                    #  E(R)=
                    #       \    0.5*KMAX*(RAVE-RMAX)**2    RMAX<RAVE<RLIM
                    #        \
                    #         \  FMAX*(RAVE-(RLIM+RMAX)/2)  RAVE>RLIM
                    #
                    # Where: RLIM = RMAX+FMAX/KMAX (the value of RAVE where the
                    #                               force equals FMAX)
                    #
                    # So, if we want RLIM == 101*RMAX,
                    # (RLIM - RMAX)*KMAX  = FMAX
                    # 100*RMAX*KMAX = FMAX
                    #
                    length = length*10.
                    forceconst = forceconst / 836.8
                    noe += """    ! You must explicitly edit this to insert the correct residue IDs.
    ASSIgn SELE RESID @INSERTRESIDHERE .AND. RESN %(resn)s .AND. TYPE %(a1)s END  SELE RESID @INSERTRESIDHERE .AND. RESN %(resn)s .AND. TYPE %(a2)s END -
           KMIN %(kmin)8.3f RMIN %(rmin)8.3f KMAX %(kmax)8.3f RMAX %(rmax)8.3f  FMAX %(fmax)8.3f
"""%{'kmin':forceconst,'rmin':length,'kmax':forceconst,'rmax':length,'fmax':forceconst*length*100,'a1':atoms[0],'a2':atoms[1],'resn':molname}
                else:
                    raise GromacsError("Bond function type %s not supported: %s"%(functype,line))

        if 'angles' in mol:
            ###  angles #yay
            # Same comments from BONDS section apply in terms of how GROMACS
            # and CHARMM define things in different places, need for sanity
            # checking, etc.
            for line in mol['angles']:
                """
                ANGLES
                !
                !V(angle) = Ktheta(Theta - Theta0)**2
                !
                !V(Urey-Bradley) = Kub(S - S0)**2
                !
                !Ktheta: kcal/mole/rad**2
                !Theta0: degrees
                !Kub: kcal/mole/A**2 (Urey-Bradley)
                !S0: A
                !
                !atom types     Ktheta    Theta0   Kub     S0
                !

                GROMACS supports those as well, but MARTINI uses the
                cosine-based angle potential, defined as

                V(angle) = (1/2)Ktheta(cos(Theta) = cos(Theta0))**2

                This can be found in seciton 4.2.6 and table 5.4.

                Tim Miller put this potential into CHARMM. You access it by
                using a negative angle. E.g. +135 becomes -225. So we'll take
                Theta to Theta-360

                Then we'll divide by 4.184 and by 2 to get CHARMM's Ktheta, a
                total conversion of division by 8.368.

                NOTE: The previously listed mailing list message again seems
                to be wrong.
                """
                if not removecomments(line).strip(): continue
                i,j,k,functype,angle,forceconst = removecomments(line).split()
                i,j,k,functype,angle,forceconst = int(i),int(j),int(k),int(functype),float(angle),float(forceconst)
                if functype != 2:
                    raise GromacsError("Angle function type %s not supported: %s"%(functype,line))
                angle = angle-360
                forceconst = forceconst / 8.368
                outer = [atomtypemap[a] for a in (i,k)]
                outer.sort()
                atoms = outer[0],atomtypemap[j],outer[1]
                rtf += 'ANGLE %-5s %-5s %-5s\n'%(atomnamemap[i],atomnamemap[j],atomnamemap[k])
                if atoms not in angles:
                    'CTL3 CTL2 CL      52.0     108.00   ! alkane'
                    angletxt += '%-5s %-5s %-5s %8.3f %8.3f\n'%(atoms[0],atoms[1],atoms[2],forceconst,angle)
                    angles[atoms] = (forceconst,angle)
                else:
                    if angles[atoms] != (forceconst,angle):
                        raise GromacsError("Multiple definitions for angle %s <--> %s <--> %s : %s vs. %s"%(atoms[0],
                                                                                                            atoms[1],
                                                                                                            (forceconst,length),
                                                                                                            angles[atoms],))
        rtf += 'PATCH FIRST NONE LAST NONE\n\n'
    rtf += 'END\n'
    prm += bondtxt
    prm += angletxt
    prm += 'END\n'
    noe += 'END\n'
            
    def writefile(fn,txt):
        f = file(fn,'w')
        f.write(txt)
        f.close()
    writefile(rtfname,rtf)
    writefile(prmname,prm)
    writefile(noename,noe)
if __name__ == '__main__':
    print "NOTE: This doesn't include connectivity information for anything. Just individual residues."
    import argparse
    parser = argparse.ArgumentParser(description='Converts gromacs itp file to CHARMM prm and rtf files')
    parser.add_argument('-i',dest='itp',default='cglipid.itp',help='Name of ITP file [default %(default)s]')

    args = parser.parse_args()

    convert(args.itp)
