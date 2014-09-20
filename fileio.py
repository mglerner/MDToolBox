#!/usr/bin/env python
"""
Here we convert between different standard formats.

At the moment, everything is expected to be in GROMACS internal units.
"""
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

def crdline2parts(line):
    pass

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

