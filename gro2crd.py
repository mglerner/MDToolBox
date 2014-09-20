#!/usr/bin/env python

from __future__ import division
import pygro as pyg
import argparse,sys

description = """Convert a .gro file into a .crd file. The output file will have the
same name as the input file, but with extension replaced by "crd". You must have
GROMACS (specifically editconf) installed to handle the recentering. Some limited
control is provided for renaming residues and specifying segments.
"""

epilog = """This program also prints out the box size information from
GROMACS. CHARMM doesn't keep box information in the .crd file, so
you'll need this when constructing the box in your CHARMM input files."""

def addtrailingspaces(d):
    """Adds trailing spaces to the values of a dictionary so that
    they are all the same length. Mostly used for remapping dictionaries
    acting on fixed-column files.
    
    Arguments:
     - `d`: the dictionary
    """
    if not d:
        return d
    newd = {}
    maxlen = max(len(v) for v in d.values())
    template = '{{v:{maxlen}s}}'.format(maxlen=maxlen)
    for (k,v) in d.iteritems():
        newd[k] = template.format(v=v)
    return newd

class DictAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        #print('%r %r %r' % (namespace, values, option_string))
        if not divmod(len(values),2)[-1] == 0:
            raise ValueError('You must pass an even number of arguments to {name}. You passed {value}.'.format(name=self.dest,value=values))
        it = iter(values)
        d = dict(zip(it, it))
        setattr(namespace, self.dest, d)


parser = argparse.ArgumentParser(description=description,
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 epilog=epilog,
    )
parser.add_argument("grofile",help="A GROMACS-formatted .gro file")
parser.add_argument('--rename',nargs='+',help='A list of RESNs to rename, e.g. "A AA B BB" would rename A to AA and B to BB',action=DictAction,default={})
parser.add_argument('--segs',nargs='+',help='A list of RESN, SEGI pairs. Any residue of name RESN will get assigned to SEGI. e.g. "DPP L W WAT" would assign DPPs to the L segment and Ws to the W segment. If you are also using --rename, use the un-renamed resis. WARNING: we do not reorder things, so this only works if your residues are listed in chunks corresponding to the segments',action=DictAction,default={})

args = parser.parse_args()

#print args

pyg.Conversion.gro2crd(fn=args.grofile,outfn=None,segmap=args.segs,resnmap=args.rename)
