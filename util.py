#!/usr/bin/env python
from __future__ import division
from __future__ import with_statement
class Placeholder:
    pass
#import psyco
#psyco.full()
DEBUG = False
from contextlib import contextmanager
import sys,os,pprint,time
import numpy as np
from numpy import apply_along_axis,array,vstack,average,std,sqrt
try:
    from numpy import loadtxt
except ImportError:
    if DEBUG: print "NO LoadTxt"
from numpy.linalg import norm

try:
    from pylab import plot,errorbar,title,xlabel,ylabel,annotate,clf
    import scipy
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib import mlab
    from matplotlib.ticker import NullFormatter
    import scipy.special
except ImportError:
    mpl = Placeholder()
    mpl.cm = Placeholder()
    mpl.cm.jet = Placeholder()
    if DEBUG: print "no scipy pylab or matplotlib for you"
try:
    import scipy.optimize
except ImportError:
    if DEBUG: print "no scipy optimize for you"


def run(prog,args,input_txt='',verbose=False):
    """
    wrapper to handle spaces on windows.
    prog is the full path to the program.
    args is a string that we will split up for you.
        or a tuple.  or a list. your call.
    input_txt will get sent to stdin.

    return value is (retval,prog_out)

    e.g.

    (retval,prog_out) = run('/bin/ls','-al /tmp/myusername')

    NOTE: I don't use Windows enough to be sure that this is
    robust. Consider jetisoning it and just using subprocess.
    """
    import subprocess,tempfile

    if type(args) == type(''):
        args = tuple(args.split())
    elif type(args) in (type([]),type(())):
        args = tuple(args)
    args = (prog,) + args
    
    try:
        output_file = tempfile.TemporaryFile(mode="w+")  # <-- shouldn't this point to the temp dir
    except IOError:
        print "Error opening output_file when trying to run the command."

    try:
        input_file = tempfile.TemporaryFile(mode="w+")  # <-- shouldn't this point to the temp dir
        if input_txt and not input_txt.endswith('\n'): input_txt = input_txt + '\n'
        input_file.write(input_txt)
        input_file.seek(0)
    except IOError:
        print "Error opening input_file when trying to run the command."

    if verbose:
        print "Running:\n\tprog=%s\n\targs=%s" % (prog,args)
    retcode = subprocess.call(args,stdout=output_file.fileno(),stderr=subprocess.STDOUT,stdin=input_file.fileno())
    output_file.seek(0)
    #prog_out = output_file.read()
    prog_out = ''.join(output_file.readlines())
    output_file.close() #windows doesn't do this automatically
    input_file.close()
    if DEBUG:
        print "Results were:"
        print "Return value:",retcode
        print "Output:"
        print prog_out
    return (retcode,prog_out)


@contextmanager
def timer(description,s=sys.stdout):
    """A simple timer context. s is the stream to which
    we will log info (defaults to sys.stdout)

    If the need arises, I can obviously extend this to 
    include generic logging capabilities, etc.

    usage:

    In [36]: with timer('hjk'):
       ....:     a = range(1000)
       ....:
       ....:
    hjk
    done with hjk 2.88486480713e-05

    """
    s.write('%s\n'%description)
    s.flush()
    start = time.time()
    yield
    s.write('done with %s. Time: %ss\n'%(description,time.time()-start))
    s.flush()

@contextmanager
def lockfile(fname,delay=1):
    """
    VERY simple lockfile, but works on NFS mounts.

    This is not safe in that it'll leave locks around when you kill
    the process. Sorry!

    *delay* We just sleep until the lockfile can be written. delay
         tells us how long to wait between sleep calls.
    """
    lockname = fname+'.lock'
    while os.path.isfile(lockname):
        print "WAITING FOR LOCK",lockname
        time.sleep(delay)
    f = open(lockname,'w')
    print "LOCK ACQUIRED",lockname
    f.write("under lockdown\n")
    f.close()
    yield
    os.remove(lockname)
    print "LOCK RELEASED",lockname

class Data:
    def __init__(self,*args,**kwargs):
        for a in args:
            for k,v in a.items():
                setattr(self,k,v)
        for k,v in kwargs.items():
            setattr(self,k,v)
    def __str__(self):
        s = 'Data instance\n-------------\n'
        for i in dir(self):
            if not i.startswith('__'):
                s += '%s\t--> %s\n'%(i,getattr(self,i))
        return s
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
# From Python Cookbook 6.18
def attributesFromArguments(d):
    """
    class Bork:
        def __init__(self,thing1,thing2):
            U.attributesFromArguments(locals())

    b = Bork('cow','frog')
    b.thing1 ---> 'cow'
    """
    self=d.pop('self')
    codeObject = self.__init__.im_func.func_code
    argumentNames = codeObject.co_varnames[1:codeObject.co_argcount]
    for n in argumentNames:
        setattr(self,n,d[n])

def fitline(xdata,ydata,uncertainties=1,m0=1,b0=1):
    #print "Fitting",xdata.size,"points"
    # errfunc will break if these are not numpy arrays, because we expect to operate with ufuncs.
    xdata = np.array(xdata)
    ydata = np.array(ydata)
    p0 = [m0,b0]
    def fitfunc(p,x):
        m,b = p
        return m*x + b
    def errfunc(p,x,y,err):
        """
        fit to a line, weighted by positional uncertainties (err).
        """
        return (y - fitfunc(p,x))/err
    # Errors ignore the weighting (for now) TODO: FIXME:
    try:
        if np.std(uncertainties) == 0.:
            #print "WARNING WARNING WARNING no uncertainties"
            uncertainties = 1
        else:
            #print "ALL IS GOOD"
            pass
    except:
        uncertainties = 1
    out = scipy.optimize.leastsq(errfunc,p0[:],
                                 args = (xdata,ydata,uncertainties),
                                 full_output=1)
    pfinal,covar,infodict,mesg,ier = out
    slope,intercept = pfinal
    slopeerror,intercepterror=0,0
    return slope,intercept,slopeerror,intercepterror
def plotbestline(xdata,ydata,uncertainties=1,m0=1,b0=1,start=0,stop=0,numerrbars=30,label=None):
    if stop == 0:
        stop = len(xdata)
    phonyUncertainties=False
    if type(uncertainties) in (type(0),type(0.)):
        phonyUncertainties=True
        uncertainties = np.ones(len(ydata)) * uncertainties
    _x,_y,_u = xdata[start:stop],ydata[start:stop],uncertainties[start:stop]
    (m,b,me,be) = fitline(_x,_y,_u,m0,b0)
    if not phonyUncertainties:
        error_step = int(len(_x)/numerrbars)
        if error_step == 0: error_step = len(_x)
        errorbar(_x[::error_step],_y[::error_step],yerr=_u[::error_step],fmt='k.')
    #plot(_x,m*_x+b,'b-')
    if label is not None:
        plot(xdata,m*xdata+b,'b-',label=label)
    else:
        plot(xdata,m*xdata+b,'b-')
    plot(xdata,ydata,'r.')
    return m,b

    
def colorcycle(colors='bgrcmyk',num=1):
    idx = 0
    yields = range(num)
    while True:
        if idx == len(colors): idx = 0
        for i in yields:
            yield colors[idx]
        idx += 1
color = colorcycle()

def testit():
    pass

def leastsqw(x,y,s):
    """
    Least square fitting to y = a + b*x, with uncertainties s.

    *x,y*
        data

    *s*
        standard deviations of the data

    Given a set of data points x,y with standard deviations s (all
    Numpy arrays), fit them to a straight line $y = ax + b$ by
    minimizing $\chi^2$.

    We return a, b, their probably uncertainties (sigmaa,sigmab), the
    $\chi^2$, and the goodness-of fit probability Q (that the fit
    would have $\chi^2$ this or larger).

    $S \equiv \sum_{i=1}^{N}\frac{1}{\sigma^2}$
    (etc.)
    
    Taken from Numerical Recipes "FIT" subroutine
    """
    overs = 1/s
    overs2 = overs*overs 
    # S is scipy
    S_ = overs2.sum() 
    Sx = (x*overs2).sum() 
    Sy = (y*overs2).sum() 
    #Sxx = (x*x*overs2).sum()
    #Sxy = (x*y*overs2).sum()

    #delta = S_*Sxx - Sx**2
    #a = (Sxx*Sy - Sx*Sxy)/delta
    #b = (S_*Sxy - Sx*Sy)/delta

    t = overs*(x - Sx/S_)
    Stt = (t**2).sum()
    overStt = 1/Stt
    b = overStt*(t*y/s).sum()
    a = (Sy - Sx*b)/S_
    siga2 = (1/S_)*(1 + Sx*Sx*overStt/S_)
    sigb2 = overStt
    siga = np.sqrt(siga2)
    sigb = np.sqrt(sigb2)

    # Now figure out the goodness of fit
    covab = -Sx*overStt/S_
    rab = covab/(siga*sigb)
    chisq = (((y - a - b*x)*overs)**2).sum()
    # note to self, you can do gammainc(4,x) to get an array back)
    # this is GAMMAQ from Numerical Recipes
    q = 1-scipy.special.gammainc(.5*(x.size-2), .5*chisq)
    d = Data(a=a,
             b=b,
             siga=siga,
             sigb=sigb,
             rab=rab,
             chisq=chisq,
             q=q,
             )
    return d

def slowacorr(data):
    return [sum([np.dot(i,j) for (i,j) in zip(data,data[offset:])]) for offset in range(len(data))]

def scaledacorr(x, stepsize, normed=False, detrend=mlab.detrend_none,
                usevlines=False, maxlags=None, **kwargs):
    
    import numpy as np
    x = detrend(np.asarray(x))
    Nx = len(x)
    y = x

    c = np.correlate(x, y, mode=2)

    if normed: c/= np.sqrt(np.dot(x,x) * np.dot(y,y))

    if maxlags is None: maxlags = Nx - 1

    if maxlags >= Nx or maxlags < 1:
        raise ValueError('maglags must be None or strictly '
                         'positive < %d'%Nx)

    lags = np.arange(-maxlags,maxlags+1) * stepsize
    c = c[Nx-1-maxlags:Nx+maxlags]

    if usevlines:
        a = vlines(lags, [0], c, **kwargs)
        b = axhline(**kwargs)
    else:

        kwargs.setdefault('marker', 'o')
        kwargs.setdefault('linestyle', 'None')
        a, = plot(lags, c, **kwargs)
        b = None

    return lags, c, a, b

def normalizeNby3(vector):
    # I can normalize an N,M vector via
    norms = apply_along_axis(norm,1,vector)
    norms = vstack([norms,]*vector.shape[1]).transpose()
    return vector/norms
def normalizeNbyMby3(vector):
    vec2 = vector.copy()
    for i in range(vector.shape[1]):
        vec2[:,i,] = normalizeNby3(vector[:,i,:])
    return vec2

"""
Some remaining ideas:

 - tweak fitline to fit two linear segments
   - maybe
     a----b--c-------d--e
     but maybe that won't converge, so maybe
     a----b-------------c
     or fix the size of b--c and d--e.
 - add the ability to remove outliers.
   - Peirce's criterion
"""

#
# Plotting routines
#

def plotdata(x,y=None,yerr=None,numblocks=10,numerrbars=50,colors=('blue','green'),cumavg=True,clear=False,
             plottitle=None,xaxislabel=None,yaxislabel=None):
    """ Our standard plotting routine.

     - The first argument is the data. We try to extract x and y data from it if possible.
     - Plot 50 error bars if we're given errors.
     - Do block averaging.
     - Plot the cumulative average.

     TODO:
      - Option for block averaging that looks at autocorrelation time to figure out block lengths.
     """
    txt = ''
    # Figure out x,y,err
    if y is None:
        if len(x.shape) > 1:
            x,y = x[:,0],x[:,1]
        else:
            x,y = array(range(len(x))),x
    if clear: clf()

    annotation_location = (min(x) + (max(x) - min(x))*0.1,min(y) + (max(y) - min(y))*0.9)
    #print annotation_location,max(y)
    plot(x,y,color=colors[0],zorder=10,alpha=0.5)

    if cumavg:
        ya = np.cumsum(y)/np.arange(1,len(y)+1)
        plot(x,ya,'k-')

    if yerr is not None:
        error_step = int(len(x)/numerrbars)
        if error_step == 0: error_step = len(x)
        errorbar(x[::error_step],y[::error_step],yerr[::error_step],color=colors[0],zorder=20)
    blocksize = int(len(x)/numblocks)
    blocksizes = [len(i) for i in splitseq(x,blocksize)] 
    x = [average(i) for i in splitseq(x,blocksize)][:numblocks]
    yerr = array([std(i) for i in splitseq(y,blocksize)])[:numblocks]
    y = array([average(i) for i in splitseq(y,blocksize)])[:numblocks]
    txt += 'Avg: %.4f'%average(y)
    txt += ', std err: %.4f, avg err: %.4f'%(std(y)/sqrt(numblocks),average(yerr))
    txt += '\nblocks of length %s'%blocksize
    if len(blocksizes) > numblocks:
        txt += ' discarding %s points at the end'%sum(blocksizes[numblocks:])
    txt += '\nblack: running average. light blue: raw data. green: block averages'
    errorbar(x,y,yerr,elinewidth=20,color=colors[1],barsabove=True,zorder=30)
    if plottitle: title(plottitle)
    if xaxislabel: xlabel(xaxislabel)
    if yaxislabel: ylabel(yaxislabel)
    annotate(txt,annotation_location)

def scatterhist(x,y):
    """ From matplotlib examples """
    nullfmt   = NullFormatter()         # no labels

    # definitions for the axes 
    left, width = min(x) - 0.01*min(x),max(x)+0.01*max(x)
    bottom, height = min(y) - 0.01*min(y),max(y)+0.01*max(y)
    bottom_h = left_h = left+width+0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    ## start with a rectangular Figure
    #plt.figure(1, figsize=(8,8))
    print rect_scatter

    axScatter = plt.gca()
    axScatter = plt.axes(rect_scatter)
    if 0:
        axHistx = plt.axes(rect_histx)
        axHisty = plt.axes(rect_histy)
        
        # no labels
        axHistx.xaxis.set_major_formatter(nullfmt)
        axHisty.yaxis.set_major_formatter(nullfmt)

    # the scatter plot:
    axScatter.scatter(x, y)

    if 0:
        # now determine nice limits by hand:
        binwidth = 0.25
        xymax = np.max( [np.max(np.fabs(x)), np.max(np.fabs(y))] )
        lim = ( int(xymax/binwidth) + 1) * binwidth
        
        axScatter.set_xlim( (-lim, lim) )
        axScatter.set_ylim( (-lim, lim) )
        
        bins = np.arange(-lim, lim + binwidth, binwidth)
        axHistx.hist(x, bins=bins)
        axHisty.hist(y, bins=bins, orientation='horizontal')
    
        axHistx.set_xlim( axScatter.get_xlim() )
        axHisty.set_ylim( axScatter.get_ylim() )

def gradplot(x,y=None,t=None,
             label = 'timeseries colored by time',
             xl="x (Angstroms)",
             yl="y (Angstroms)",
             cmap = mpl.cm.jet,
             linewidth=0.5,
             marker=(4,0,0),
             showlines=True,
             **kwargs):
    """
    Plot x,y colored by t.
    """
    if y is None:
        x,y = x[:,0],x[:,1]
    if t is None:
        t = np.linspace(0,100,len(x))
    #plt.plot(x,y)
    points = zip(x,y)
    segments = zip(points[:-1],points[1:])
    ax = plt.gca()
    if showlines:
        #colors = [cmap(i) for i in mpl.colors.normalize()(t)]
        colors = cmap(mpl.colors.normalize()(t))
        LC = mpl.collections.LineCollection(segments, colors = colors)
        LC.set_linewidth(linewidth)
        ax.add_collection(LC)
    if xl is not None: 
        plt.xlabel(xl)
    if yl is not None:
        plt.ylabel(yl)
    plt.scatter(x,y,c=t,edgecolor='None',cmap=cmap,label=label,marker=marker,**kwargs)
    #plt.scatter(x,y,c=t,cmap=cmap)
    #plt.legend()

def smooth(x,windowlen=10,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        windowlen: the dimension of the smoothing window
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   

    Taken from the scipy cookbook: http://www.scipy.org/Cookbook/SignalSmooth
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < windowlen:
        raise ValueError, "Input vector needs to be bigger than window size."


    if windowlen<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=np.r_[2*x[0]-x[windowlen:1:-1],x,2*x[-1]-x[-1:-windowlen:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(windowlen,'d')
    else:
        w=eval('np.'+window+'(windowlen)')

    y=np.convolve(w/w.sum(),s,mode='same')
    return y[windowlen-1:-windowlen+1]

if 0:
    def smooth_demo():
        t=linspace(-4,4,100)
        x=sin(t)
        xn=x+randn(len(t))*0.1
        y=U.smooth(x)
        ws=31
        subplot(211)
        plot(ones(ws))
        windows=['flat', 'hanning', 'hamming', 'bartlett', 'blackman']
        hold(True)
        for w in windows[1:]:
            eval('plot('+w+'(ws) )')
        axis([0,30,0,1.1])
        legend(windows)
        title("The smoothing windows")
        subplot(212)
        plot(x)
        plot(xn)
        for w in windows:
            plot(U.smooth(xn,10,w))
        l=['original signal', 'signal with noise']
        l.extend(windows)
        legend(l)
        title("Smoothing a noisy signal")
        #show()
def prime_factors(n):
    """ Return the prime factors of the given number. """
    factors = []
    lastresult = n
    # 1 is a special case
    if n == 1:
        return [1]
    while 1:
        if lastresult == 1:
            break
        c = 2
        while 1:
            if lastresult % c == 0:
                break
            c += 1
        factors.append(c)
        lastresult /= c
    return factors
class C(object):
    def __init__(self):
        self._x = None

    def getx(self):
        return self._x
    def setx(self, value):
        self._x = value
    def delx(self):
        del self._x
    x = property(getx, setx, delx, "I'm the 'x' property.")

