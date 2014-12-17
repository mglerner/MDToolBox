#!/usr/bin/env python
import sys
from pygro import ITP

if __name__ == '__main__':
    import optparse
    usage = """combineitp.py in_itp1 in_itp2 out_itp

This combines two molecules, one in each of the input itps, into a
single molecule in the output itp file.

Only one molecule may be defined per input file.
"""
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-i','--input',nargs=2,help='List of input files')
    parser.add_option('-m','--moltype',help='Name of output molecule type')
    parser.add_option('-o','--output',help='Name of output ITP file')
    parser.add_option('-B','--skip_blanks',dest='blanks',action='store_false',default=True)
    options,args = parser.parse_args()

    i1 = ITP(options.input[0])
    print i1
    i2 = ITP(options.input[1])
    print i2
    i = i1 + i2
    i.setmoltype(options.moltype)
    print i
    i.write(options.output,options.blanks)
    print options.blanks
    sys.exit()

    if 0:
        i1 = ITP('protein_withele.itp')
        i2 = ITP('dpponly.itp')
        i = i1 + i2
        print i.natoms
        i.moltype = 'PHD'
        i.write('phd.itp')
