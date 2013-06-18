#!/usr/bin/env python
#!/v/bigbox7b/home/mglerner/software/bin/python
#!/v/apps/python-2.5/bin/python


import os,glob,re,sys
import numpy as N, matplotlib as M, pylab as P
M.rcParams['legend.fancybox'] = True
from pylab import plot,xlabel,ylabel,legend,title,clf,subplot,errorbar,savefig,annotate,axis
from numpy import array,average,std
import gzip,bz2
import util as u

def splitseq(seq,size):
    """ Split up seq in pieces of size

    In [34]: u.splitseq(range(30),10)
    Out[34]: 
    [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
    [10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
    [20, 21, 22, 23, 24, 25, 26, 27, 28, 29]]

    In [35]: u.splitseq(range(34),10)
    Out[35]: 
    [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
    [10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
    [20, 21, 22, 23, 24, 25, 26, 27, 28, 29],
    [30, 31, 32, 33]]
    """
    try:
        return [seq[i:i+size] for i in range(0, len(seq), size)]
    except ValueError:
        print "Cannot split this seq",seq
        print "Into this size",size
        raise

def plotdata(x,y=None,yerr=None,numblocks=10,numerrbars=50,colors=('blue','green'),clear=False,
             plottitle=None,xaxislabel=None,yaxislabel=None,description='',grid=False):
    """ Our standard plotting routine.

     - The first argument is the data. We try to extract x and y data from it if possible.
     - Plot 50 error bars if we're given errors.
     - Do block averaging.
     - description gets added to the legend and title.
     """
    txt = ''
    # Figure out x,y,err
    if y is None:
        if len(x.shape) > 1:
            x,y = x[:,0],x[:,1]
        else:
            x,y = array(range(len(x))),x
    if clear: clf()

    if len(x) >= 100000:
        step = divmod(len(x),100000)[0]
    else:
        step = 1
    plot(x[::step],y[::step],color=colors[0],zorder=10,alpha=0.5,label=description+'raw data')
    xmin,xmax,ymin,ymax = x.min(),x.max(),y.min(),y.max()
    #print y
    if yerr is not None:
        error_step = int(len(x)/numerrbars)
        if error_step == 0: error_step = len(x)
        errorbar(x[::error_step],y[::error_step],yerr[::error_step],color=colors[0],zorder=20,label='FLUC>')
        ymin = ymin - abs(yerr.max())
        ymax = ymax + abs(yerr.max())
    else:
        ymin = ymin - abs(.1*ymin)
        ymax = ymax - abs(.1*ymax)
    if numblocks <= len(x):
        blocksize = int(len(x)/numblocks)
    else:
        blocksize = 1
        numblocks = len(x)
        print "Not enough data points, setting blocksize to",blocksize,"and numblocks to",numblocks
    blocksizes = [len(i) for i in splitseq(x,blocksize)] 
    x = [average(i) for i in splitseq(x,blocksize)][:numblocks]
    yerr = array([std(i) for i in splitseq(y,blocksize)])[:numblocks]
    y = array([average(i) for i in splitseq(y,blocksize)])[:numblocks]
    txt += 'Avg: %.4f'%average(y)
    txt += ', std err: %.4f, avg err: %.4f'%(std(y),average(yerr))
    txt += '\nblocks of length %s'%blocksize
    if len(blocksizes) > numblocks:
        txt += ' discarding %s points at the end'%sum(blocksizes[numblocks:])
    txt += '\nslope of best-fit line to block averages:  %s'% u.fitline(x,y)[0]
    #print y
    #errorbar(x,y,yerr,elinewidth=20,color=colors[1],barsabove=True,zorder=30)
    errorbar(x,y,yerr,color=colors[1],barsabove=True,zorder=30,label=description+'block averages')

    # After the last plotting, we need to set the axis explicitly,
    # due to the fact that someone might have previously plotted
    # raw data in plotaver(). In that case, we still want to be
    # zoomed in on the average data.
    axis([xmin,xmax,ymin,ymax])
    if plottitle: title(description+plottitle)
    if xaxislabel: xlabel(xaxislabel)
    if yaxislabel: ylabel(yaxislabel)
    if grid: P.grid('on')
    else: P.grid('off')
    annotation_location = (min(x) + (max(x) - min(x))*0.1,min(y) + (max(y) - min(y))*0.1)
    annotation_location = (xmin + (xmax - xmin)*0.1,ymin + (ymax - ymin)*0.1)
    #print annotation_location,max(y)
    annotate(txt,annotation_location)
    legend(shadow=True)
"""
Display parameters from a dynamics run.
"""
# My own universal file. Just a wrapper around the unzipping classes.
def ufile(fname):
    """ Return uncompressed version of file """
    if fname.endswith('.gz'):
        return gzip.GzipFile(fname)
    elif fname.endswith('.bz2'):
        return bz2.BZ2File(fname)
    return file(fname)

# Extract name/value pairs
def extract_name_value_pairs(line,skip=1):
    # Split up name=value lines that can have multi-word names and
    # non-uniform numbers of spaces. E.g.
    #  FLUC A     =    0.00000 B    =    0.00000 C     =    0.00000
    #  FLUC Alpha =    0.00000 Beta =    0.00000 Gamma =    0.00000
    #  FLUC PIXX =     1645.47 PIYY =    1658.51 PIZZ =    1807.79
    #  FLUC PIXY =      958.14 PIXZ =    1002.77 PIYZ =    1037.63
    #  FLUC Gradient Norm =   19.11575
    line = ' '.join(line.split()[skip:]) # it always starts with DYNA, FLUC, AVER, etc.
    parts = line.split(' = ')
    names = [parts[0],]
    values = []
    for part in parts[1:-1]:
        _part = part.split()
        values.append(float(_part[0]))
        names.append(' '.join(_part[1:]))
    values.append(float(parts[-1]))
    result = dict(zip(names,values))
    return result

class DynProp:
    """This is a generic class for encapsulating the results of a CHARMM
    dynamics run.  You want to initialize it with a directory that
    contains (optionally) a dyn.out and an Out subdir that contains
    other dyn.out files. They can be gzipped.

    You can plot specific properties vs. time with ``plotprop()`` or
    ``plotaver()``. ``plotaver()`` does block averaging.
    
    You can get a list of properties with ``listprops()``.

    You can plot specific subsets of properties with ``plotbox()``,
    ``plotext()``, ``plotstandardff()``, ``plototherff()``, and
    ``plottvp()``.

    """
    def __init__(self, d,numfiles=0,namepat='dyn*.out',skipcom=False):
        """
        Look in directory d for dyn*.out files, initialize
        
        Arguments:
        - `d`:directory to look in
        - `numfiles`:number of files to look at. 0 means everything.
        - `namepat`:must contain a single * (look at definition of frameno below).
                    We will automatically look for thing.gz in addition to thing.
        """

        """
        There are 3 kinds of DYNA lines:

        1. Template lines, which look like DYNA SOMETHING: blah blha blha.
           These tell us what DYNA property lines look like.

        2. Property lines, which look like DYNA SOMETHING> blah blhablh
           These fill in the properties we've defined above.

        3. Crystal lines. These look like DYNA A = 35.  These are like
           property lines, but they follow name = value syntax and
           don't give us templates ahead of time. Note that There are
           words like "Gradient Norm", so they're not just single word
           names.
        """
        self.directory = d
        self.namepat = namepat
        self.skipcom = skipcom
        if not self.skipcom:
            self._comdict = {} # keys are times
        self._dyndict = {}
        self._averdict = {}
        self._flucdict = {}
        self._minidict = {}
        self._dicts = {'DYNA':self._dyndict,
                       'AVER':self._averdict,
                       'FLUC':self._flucdict,
                       'MINI':self._minidict,
                       }
        self._templates = {}
        self._props = {}
        self.time = 0.0

        fnames = glob.glob(os.path.join(d,self.namepat+'*')) + glob.glob(os.path.join(d,'Out',self.namepat+'*'))
        print "fnames",fnames
        print "namepat",namepat
        # dyn1.out .. dyn99.out
        def fnamenumber(fname):
            """A/B/dyn32.out.gz --> 32
                   dyn*.out
                   bh1f-rex-ld-2.out_1
                   bh1f-rex-ld-2.out_*
            """
            try:
                fn = os.path.split(fname)[-1]
                sidx = self.namepat.index('*')
                nums = ''.join([c for c in fn[sidx:] if c in '0123456789'])
                return int(nums)
                return int(os.path.split(fname)[-1].split('.')[0][3:])
            except:
                print "Found no frame number; returning 99999999"
                return 99999999
        fnames = [(fnamenumber(fname),fname) for fname in fnames]
        fnames.sort()
        fnames = [part[1] for part in fnames]
        if 0:
            if os.path.isfile(os.path.join(d,'dyn.out')):
                fnames.append(os.path.join(d,'dyn.out'))
        if numfiles > 0:
            fnames = fnames[:numfiles]

        self._seedtemplates()

        for fname in fnames:
            print "Processing",fname
            self._appenddatafromfile(ufile(fname))
        print "Finalizing data"
        self._finalizedata()
    def _seedtemplates(self):
        # Seed with known templates
        for line in """DYNA DYN: Step         Time      TOTEner        TOTKe       ENERgy  TEMPerature
DYNA PROP:             GRMS      HFCTote        HFCKe       EHFCor        VIRKe
DYNA INTERN:          BONDs       ANGLes       UREY-b    DIHEdrals    IMPRopers
DYNA CROSS:           CMAPs
DYNA EXTERN:        VDWaals         ELEC       HBONds          ASP         USER
DYNA IMAGES:        IMNBvdw       IMELec       IMHBnd       RXNField    EXTElec
DYNA EWALD:          EWKSum       EWSElf       EWEXcl       EWQCor       EWUTil
DYNA CONSTR:       HARMonic    CDIHedral          CIC     RESDistance       NOE
DYNA MMFP:              GEO        MDIP           SSBP        SHEL       DROFfa
DYNA PRESS:            VIRE         VIRI       PRESSE       PRESSI       VOLUme
DYNA XTLE:                       XTLTe         SURFtension  XTLPe        XTLtemp
MINI MIN: Cycle      ENERgy      Delta-E         GRMS    Step-size
MINI MISCEL:         SBOUnd          TSM        PRMS        PANGle        other
MINI QUANTM:        QMELec        QMVDw  """.split('\n'):
            self._addtempl(line)
    
    def _appenddatafromfile(self,f):
        proptypes = 'DYNA AVER FLUC MINI'.split()
        pt = '|'.join(proptypes)
        templpat = re.compile('('+pt+') \w*:')
        proppat = re.compile('('+pt+')[ A-Za-z]*>')
        crystpat = re.compile('('+pt+')[ A-Za-z]*=')
        for line in f:
            if line.strip() == 'DETAILS ABOUT CENTRE OF MASS' and not self.skipcom:
                # Adds the next several lines
                self._addcomprop(f)
            elif line.strip()[:4] in proptypes:
                if templpat.match(line):
                    # It's easier to just seed this with the known
                    # templates.
                    #self._addtempl(line)
                    pass
                elif proppat.match(line):
                    self._addprop(line)
                elif crystpat.match(line.strip()):
                    self._addcryst(line)
    def _addtempl(self,line):
        #DYNA DYN: Step         Time      TOTEner        TOTKe       ENERgy
        parts = line.split()
        linetype = parts[1][:-1] # skip the :
        if linetype not in self._templates:
            self._templates[linetype] = parts[2:]
    def _addprop(self,line):
        """
        These all get converted to arrays later on. To save time, we
        make everything here a list, even if it contains only a single
        value. I.e. COM properties are [x,y,z] and normal properties
        are [x,]. It turns out that the parsing stage (here) is
        relatively fast compared to the data conversion stage later
        on.
        """
        #DYNA>        0      0.00000   1196.49291   1332.28317   -135.79026
        #DYNA PROP>          0.95423   1196.57558   1332.53115      0.08267
        line = line.replace('-',' -')
        line = line.replace('>','> ')
        parts = line.split()
        proptype = parts[0][:4]
        if parts[0][-1] == '>':
            # This one is special because we need to update the time
            # so that everything else is stored correctly. Also, it's
            # missing the linetype.
            if parts[0].startswith('MINI'):
                self.time = float(parts[1])
                parts = [parts[0],'MIN>'] + parts[1:]
            else: # parts[0].startswith('DYN'): # Also handles AVER lines, etc.
                # So here's an annoying issue: sometimes the step and
                # time columns run into eachother for long simulations. E.g.:
                #line DYNA>   5000001000000.00000 -14813.48661   3407.35317  -18220.83978    323.36787
                #line DYNA>        01280000.00000 -14837.90276   3373.91949 -18211.82225    320.19491
                #line AVER>     10001280020.00000 -14814.10356   3405.04348 -18219.14704    323.14867
                # as opposed to
                #line DYNA>        0 220000.00000-235970.65702  54412.32678-290382.98380    322.65757
                
                # in that case, we need to split up that first entry
                # into <step> and <time>. We'll always have a decimal
                # and five numbers after it. So, the length of our
                # time string will be <time>+6. However, we'll need to
                # figure out what the actual time should be first.
                if '.' in parts[1]: #Decimal in step means they're merged
                    try:
                        _dt = self._dt
                    except AttributeError:
                        _times = self._dicts[proptype].keys()
                        _times.sort()
                        _dt = _times[-1] - _times[-2]
                        self._dt = _dt
                    _nexttime = self.time + _dt
                    tslen = len(str(int(_nexttime))) + 6
                    assert len(parts[1]) > tslen #(sanity check)
                    parts = [parts[0], parts[1][:-tslen], parts[1][-tslen:]] + parts[2:]

                    if float(parts[2]) == 0.0:
                        print "Bad time for"
                        print "line",line
                        print "parts",parts
                        print "previous time",self.time
                        print "tslen",tslen
                        print "It shouldn't be zero since the step and time are merged."
                        sys.exit()
                self.time = float(parts[2])
                parts = [parts[0],'DYN>'] + parts[1:]
            if self.time in self._dicts[proptype] and self._dicts[proptype][self.time] != {}:
                # This can show up repeatedly if NTRFRq spits things
                # out for recentering.
                return
                assert self.time not in self._dyndict and self.time not in self._minidict
            self._dicts[proptype][self.time] = {}
        linetype = parts[1][:-1] # strip out the >
        try:
            for (i,prop) in enumerate(self._templates[linetype]):
                try:
                    self._dicts[proptype][self.time][prop] = [float(parts[2+i]),]
                except IndexError:
                    print "line",line
                    print "linetype",linetype
                    print "prop",prop
                    print "i",i,"2+i",2+i
                    print "time",self.time
                    print "proptype",proptype
                    
                    raise
        except KeyError:
            print "Could not parse",line
            print "linetype",linetype
            raise
    def _addcryst(self,line):
        """
        see comments from _addprop
        """
        # It would be better to use regexps to split this up.
        #DYNA A     =   25.16770 B    =   25.16770 C     =   25.16770
        #DYNA Gradient Norm =   67.79419
        #line = ' '.join(line.split()[1:])
        proptype = line.strip().split()[0][:4]
        for (k,v) in extract_name_value_pairs(line).items():
            self._dicts[proptype][self.time][k] = [v,]
        
    def _addcomprop(self,f):
        """
        see comments from _addprop
        """
        line = f.next()
        assert line.strip().startswith('POSITION')
        if self.time in self._comdict:
            print "Found repeated time",self.time
            print "Probably at end of file for COM recentering"
            return
        self._comdict[self.time] = {'Position':[float(i) for i in line.split()[-3:]]}
        line = f.next()
        assert line.strip().startswith('VELOCITY')
        self._comdict[self.time]['Velocity'] = [float(i) for i in line.split()[-3:]]
        line = f.next()
        assert line.strip().startswith('ANGULAR MOMENTUM')
        self._comdict[self.time]['Angular momentum'] = [float(i) for i in line.split()[-3:]]
        line = f.next()
        assert line.strip().startswith('KINETIC ENERGY')
        self._comdict[self.time]['Kinetic energy'] = [float(line.split()[-1]),]
    def _getarrays(self,d):
        # This is a bit of a bottleneck, so I'm breaking my "KISS"
        # strategy and optimizing it a bit.
        #
        # We want to return results[prop] = array(props[prop])
        # Those will be 2d arrays, [step,value].
        
        # 1. Seed props
        props = {}
        #print "Will I die? At this point, I'll seed with times",d
        #
        # Note to self: which step should we seed with? It used to be
        # the last step. That causes trouble if output got truncated
        # in the middle of the last step, so we'll take the 2nd to
        # last instead. I hate clever.
        # 
        for prop in d[d.keys()[-2]]:
            props[prop] = []
        for step in sorted(d.keys()):
            for prop in d[step].keys(): # We can't do for prop in
                                        # props because some may be
                                        # missing. E.g. the first
                                        # frame output will be missing
                                        # the DYNA XTLE line.
                p = d[step][prop]
                #if type(p) in [type([]),type(())]:
                sp = [step,] + p
                #    print "I'm getting a list for",prop,step,p
                #    print p
                #    a  = 1/0
                #else: sp = [step,p]
                try:
                    props[prop].append(sp)
                except:
                    print "Could not add",prop,sp
                    print props.keys()
                    print "step",step
                    print "this",d[step].keys()
                    print "seed",d[d.keys()[-2]].keys()
                    raise
        result = {}
        for prop in props:
            result[prop] = array(props[prop])
        return result
    def _finalizedata(self):
        """ Now that we've collected it, put it into numpy arrays """
        # We have to be a little careful because not everything is
        # output at each step.
        self._props = {}
        if self._dyndict.keys() and self._dyndict[self._dyndict.keys()[0]]: # If there's anything in there.
            print "Converting _dyndict and _comdict to a reasonable format"
            if self.skipcom:
                dirs = (self._dyndict,)
            else:
                dirs = (self._dyndict,self._comdict)
            for d in dirs:
                sys.stdout.write('.')
                sys.stdout.flush()
                arrays = self._getarrays(d)
                for a in arrays:
                    sys.stdout.write(',')
                    sys.stdout.flush()
                    self._props[a] = arrays[a]
        elif self._minidict.keys() and self._minidict[self._minidict.keys()[0]]: # If there's anything in there.
            print "Converting _minidict and _comdict to a reasonable format"
            if skipcom:
                dirs = (self._minidict,)
            else:
                dirs = (self._minidict,self._comdict)
            for d in dirs:
                sys.stdout.write('.')
                sys.stdout.flush()
                arrays = self._getarrays(d)
                for a in arrays:
                    sys.stdout.write(',')
                    sys.stdout.flush()
                    self._props[a] = arrays[a]
        else:
            print "Could not find properties!"
            sys.exit()
            
        print "Converting _averprops to a reasonable format"
        self._averprops = self._getarrays(self._averdict)
        print "Converting _flucprops to a reasonable format"
        self._flucprops = self._getarrays(self._flucdict)
        print self._props.keys()
        #print "Read",len(self._props['time'][:,0]),"Frames",self._props['time'][:,0][-1]/1000.,"ns"
    def listprops(self):
        return self._props.keys()
    def plotprop(self, prop, begin=None,grid=False):
        d = self._props[prop]
        x = d[:,0] 
        # convert ps to ns
        x = x / 1000.
        y = d[:,1]
        if begin != None:
            begin_idx = (x < begin).sum()
            x = x[begin_idx:]
            y = y[begin_idx:]
        if len(x) >= 100000:
            step = divmod(len(x),100000)[0]
        else:
            step = 1
        plot(x[::step],y[::step],)
        xlabel('Time (ns)')
        if self.skipcom and (prop in self._comdict[self._comdict.keys()[0]].keys()):
            title('COM '+prop)
            print 'COM:  ',prop
        else:
            print 'NOCOM:',prop,self._comdict[self._comdict.keys()[0]].keys()
            title(prop)
        #ylabel(prop)
        if grid: P.grid('on')
        else: P.grid('off')
    def plotaver(self, prop,colors=('blue','green'),includeraw=True,clear=False,begin=None,grid=False):
        """
        Plot error bars for 25 points.
        Plot 10 blocks
        """
        # convert ps to ns
        if clear: clf()
        if includeraw:
            d = self._props[prop]
            x = d[:,0] 
            # convert ps to ns
            x = x / 1000.
            y = d[:,1]
            if begin != None:
                begin_idx = (x < begin).sum()
                x = x[begin_idx:]
                y = y[begin_idx:]
            if len(x) >= 100000:
                step = divmod(len(x),100000)[0]
            else:
                step = 1
            plot(x[::step],y[::step],color='black',alpha=0.2,label='raw data')
            
        x = self._averprops[prop][:,0]  / 1000.
        y = self._averprops[prop][:,1]
        yerr = self._flucprops[prop][:,1]
        #print "x starts with",x[:10],begin
        if begin != None:
            begin_idx = (x < begin).sum()
            x = x[begin_idx:]
            y = y[begin_idx:]
            yerr = yerr[begin_idx:]
        #print "x now starts with",x[:10]
        plotdata(x=x,y=y,yerr=yerr,plottitle=prop,xaxislabel='Time (ns)',description = 'AVER ',grid=grid)
    # Plot specific collections of data
    def _pa(self,clear,r,c,props):
        if clear: clf()
        for (i,p) in enumerate(props):
            subplot(r,c,i+1)
            if p in self._averprops:
                self.plotaver(p)
            else:
                self.plotprop(p)
    def _pp(self,clear,r,c,props):
        if clear: clf()
        for (i,p) in enumerate(props):
            subplot(r,c,i+1)
            self.plotprop(p)
        
    def generateFigsAndDatafiles(self,overwrite,outdir=None,begin=None,includeraw=True,grid=False):
        """
        Our plots will start at time=begin (in nanoseconds). If begin is set to None,
        we will just start at the beginning of the data.
        """
        cols = 4
        if outdir is None:
            outdir = os.path.join(self.directory,'props')
        datadir = os.path.join(outdir,'data')
        if not overwrite:
            if os.path.exists(outdir):
                sys.exit('Output directory (%s) exists. Please remove it or set "overwrite" to True.'%outdir)
        os.system('mkdir -p '+datadir)
        propnames = []
        # We use a different definition of this funciton below to get relative links in the HTML file.
        def figname(pn):
            return os.path.join(datadir,pn+'.png')
        def dataname(pn):
            return os.path.join(datadir,pn+'.dat')
        def savedata(p,fname):
            d = self._props[p]
            x = d[:,0] 
            # convert ps to ns
            x = x / 1000.
            y = d[:,1]
            f = file (fname,'w')
            for xy in zip(x,y):
                f.write('%f\t%f\n'%xy)
            f.close()
        for p in self._props.keys():
            pn = ''.join([i for i in p if i in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'])
            print "Doing",p
            clf()
            if p in self._averprops:
                if p.startswith('H') and False: # Hopefully dealt with by not plotting more than 200k points.
                    print "Trying to skip raw for",p
                    self.plotaver(p,begin=begin,includeraw=False,grid=grid)
                else:
                    self.plotaver(p,begin=begin,includeraw=includeraw,grid=grid)
            else:
                self.plotprop(p,begin=begin,grid=grid)
            ### Note: This call to savefig often fails with the
            ### following runtime error:
                
            ### RuntimeError: Agg rendering complexity exceeded. Consider downsampling or decimating your data.

            ### when plotting HFCTote with raw data
            ### shown. Unfortunately, catching the error, calling clf,
            ### and replotting without the raw data produces exactly
            ### the same error. The only thing that seems to help is
            ### just setting includeraw to false the first time through.
            ### Hence the check for p.startswith('H') above. FWIW, I think
            ### I tried it with just checking for HFCTote, but that failed.
            ### I might have gotten the syntax wrong, though.

            ### It looks like I've fixed that completely by just
            ### plotting 200k points, but I've left the comments in
            ### for clarity.
            savefig(figname(pn))
            savedata(p,dataname(p))
            propnames.append(pn)
                    
                
        html = '''<html><head><title>CHARMM Output Properties</title></head><body><h1>CHARMM Output Properties</h1><br><h2>Click on thumbnails for larger pictures</h2><br><table border="0" cellpadding="7" cellspacing="0">'''
        def figname(pn):
            return os.path.join('data',pn+'.png')
        propnames.sort()
        for tr in splitseq(propnames,cols):
            html += '<tr>'
            for td in tr:
                html += '  <td>%s<br><a href="%s"><img src="%s" width="100%%"/></td>\n'%(td,figname(td),figname(td))
                
            html += '</tr>'
        html += '</table></body></html>'
        f = file(os.path.join(outdir,'index.html'),'w')
        f.write(html)
        f.close()


        
def analyzedir(d,overwrite,outdir,begin,namepat,numfiles,grid,skipcom):
    """Reads in the directory, makes images, write out html file with links to images"""
    dp = DynProp(d,namepat=namepat,numfiles=numfiles,skipcom=skipcom)
    dp.generateFigsAndDatafiles(overwrite=overwrite,outdir=outdir,begin=begin,grid=grid)
    return dp

if __name__ == '__main__':
    usage = """This script will look for CHARMM output files in the in the specified
directory as well as a subdirectory called 'Out'. .gz and .bz2 files
are automatically uncompressed.

All properties from the CHARMM output files will be plotted, including
DYNA, AVER, FLUC and MINI lines. LAVE and LFUC lines are ignored, as
the same information is contained in the AVER and FLUC lines.

An html file that displays all of the resulting images will be
generated in props/index.html. Individual files containing the
timeseries data can be found by replacing the image name with .dat
(e.g. TEMPerature.png --> TEMPerature.dat).

KNOWN BUGS:

- for output files containing repeated MINI lines, the reported time
will be step/1000. It would be better to list the x axis as 'step'.

"""
    from optparse import OptionParser
    parser = OptionParser(usage=usage)
    parser.add_option("-d", "--dir", dest="directory",
                      help="Directory to analyze", default='.')
    parser.add_option("-o", "--out-dir", dest="outdirectory",
                      help="Directory to analyze", default='./props')
    parser.add_option("-O", action="store_false", dest="overwrite", default=True,
                      help="Do not overwrite output files")
    parser.add_option('-b','--begin',dest='begin',default=-9999,type='float',
                      help="Time (in nanoseconds) at which to begin the plots (typically set to something > 0 for equilibration purposes). -9999 means begin at whatever the beginning is. [default: %default]")
    parser.add_option('-n','--namepat',dest='namepat',default='dyn*.out',
                      help="Pattern to match for names. If it doesn't contain an asterix, one will be added before the final '.'. If there is no '.', one will be added at the end. [default: dyn*.out]")
    parser.add_option('-N','--numfiles',dest='numfiles',default=0,type='int',
                      help='Number of files to process. 0 means process all files. Mostly for debugging. [default: %default]'
                      )
    parser.add_option('-g','--grid',dest='grid',default=False,action='store_true',
                      help='Show grid lines [default: %default]',
                      )
    parser.add_option('-C','--no-com',dest='skipcom',default=False,action='store_true',
                      help='This means we do not expect COM information in the output file.')

    (options, args) = parser.parse_args()
    if options.begin == -9999: options.begin = None

    if '*' not in options.namepat:
        if '.' in options.namepat:
            i = options.namepat.rfind('.')
            options.namepat = options.namepat[:i]+'*'+options.namepat[i:]
        options.namepat = options.namepat + '*'
    
    analyzedir(d=options.directory,overwrite=options.overwrite,outdir=options.outdirectory,begin=options.begin,namepat=options.namepat,numfiles=options.numfiles,grid=options.grid,skipcom=options.skipcom)
