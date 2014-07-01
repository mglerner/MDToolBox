#!/usr/bin/env python

# TODO:
# - check cpu load

import os,glob

DEBUG = False

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

def run(prog,args):
    '''
    wrapper to handle spaces on windows.
    prog is the full path to the program.
    args is a string that we will split up for you.
        or a tuple.  or a list. your call.

    return value is (retval,progout)

    e.g.

    (retval,progout) = run("/bin/ls","-al /tmp/myusername")
    '''
    import subprocess,tempfile

    if type(args) == type(''):
        args = tuple(args.split())
    elif type(args) in (type([]),type(())):
        args = tuple(args)
    args = (prog,) + args
    
    try:
        outputfile = tempfile.TemporaryFile(mode="w+")
    except IOError:
        print "Error opening outputfile when trying to run the command."

    if DEBUG: 
        print "Running:\n\tprog=%s\n\targs=%s" % (prog,args)
    retcode = subprocess.call(args,stdout=outputfile.fileno(),stderr=subprocess.STDOUT)
    outputfile.seek(0)
    #progout = outputfile.read()
    progout = ''.join(outputfile.readlines())
    outputfile.close() #windows doesn't do this automatically
    if DEBUG:
        print "Results were:"
        print "Return value:",retcode
        print "Output:"
        print progout
    return (retcode,progout)
def getlastmod(workdir):
    fn = os.path.join(workdir,'dyn.out')
    if not os.path.isfile(fn):
        fn = os.path.join(workdir,'mdcont.trr')
        if not os.path.isfile(fn):
            return '???'
    retval,progout = run('ls',('-lhc',fn))
    dynoutsize = progout.split()[4]
    dynoutmod  = ' '.join(progout.split()[5:8])
    return dynoutmod
def getsimtime(workdir):
    try:
        import charmm
        try:
            inp = charmm.DynamicsInputFile('dyn.inp',workdir)
        except (KeyError,IndexError,TypeError):
            print "Unparsable dyn.inp in",workdir
            return '???'
        ntf = charmm.numtrajfiles(d=workdir,skipfirst=False)
        return '%s ns'%(inp.duration_ns*ntf,)
    except (ImportError, IOError):
        return '???'


def getjobinfo(jobid,load):
    """
    jobid is obviously the job id. 
    if load is true, we'll ssh over to the node and figure out
    the load.
    """
    '''
Job Id: 139273.m1.lobos.nih.gov
    Job_Name = tryitout
    Job_Owner = mglerner@m1.lobos.nih.gov
    resources_used.cput = 00:00:00
    resources_used.mem = 3100kb
    resources_used.vmem = 138456kb
    resources_used.walltime = 00:21:42
    job_state = R
    queue = xeon
    server = m1.lobos.nih.gov
    Checkpoint = u
    ctime = Thu Apr 22 18:07:46 2010
    Error_Path = m1.lobos.nih.gov:/v/bigbox7b/home/mglerner/Gromacs/Diffusion/
	TetheredLipidTrimers/60a/47/tryitout.e139273
    exec_host = h117.lobos.nih.gov/7+h117.lobos.nih.gov/6+h117.lobos.nih.gov/5
	+h117.lobos.nih.gov/4+h117.lobos.nih.gov/3+h117.lobos.nih.gov/2+h117.l
	obos.nih.gov/1+h117.lobos.nih.gov/0
    Hold_Types = n
    Join_Path = oe
    Keep_Files = n
    Mail_Points = a
    Mail_Users = mglerner
    mtime = Thu Apr 22 18:07:52 2010
    Output_Path = m1.lobos.nih.gov:/v/bigbox7b/home/mglerner/Gromacs/Diffusion
	/TetheredLipidTrimers/60a/47/tryitout.o139273
    Priority = 0
    qtime = Thu Apr 22 18:07:46 2010
    Rerunable = False
    Resource_List.nodect = 1
    Resource_List.nodes = 1:ppn=8:htown
    session_id = 4302
    substate = 42
    Variable_List = PBS_O_HOME=/v/bigbox7b/home/mglerner,
	PBS_O_LANG=en_US.UTF-8,PBS_O_LOGNAME=mglerner,
	PBS_O_PATH=/v/apps/gromacs/4.0.5-clovertown-ifort11.openib/bin/:/v/ap
	ps/gromacs/4.0.5-clovertown-ifort11.eth/bin/:/v/apps/intel/11.0/bin/in
	tel64:/v/bigbox7b/home/mglerner/Python:/v/bigbox7b/home/mglerner/bin:/
	v/bigbox7b/home/mglerner/software/bin:/v/apps/bin:/usr/kerberos/bin:/u
	sr/local/bin:/bin:/usr/bin:/v/apps/compilers/x86_64/pathscale/3.2/bin:
	/v/bigbox7b/home/mglerner/work/PyPAT/code/drivers:/v/bigbox12/home/aok
	ur/amber/amber10//exe.htown,PBS_O_MAIL=/var/spool/mail/mglerner,
	PBS_O_SHELL=/bin/bash,PBS_SERVER=m1.lobos.nih.gov,
	PBS_O_HOST=m1.lobos.nih.gov,
	PBS_O_WORKDIR=/v/bigbox7b/home/mglerner/Gromacs/Diffusion/TetheredLip
	idTrimers/60a/47,PBS_O_QUEUE=entry
    etime = Thu Apr 22 18:07:46 2010
    submit_args = ./rungro
    start_time = Thu Apr 22 18:07:52 2010
    start_count = 1
    '''
    (retcode,progout) = run('qstat',('-f',jobid))
    lines = progout.split('\n')
    def iscontinuationline(line):
        return line.startswith('	')
    def isstartline(line):
        return line.startswith('    ')
    def splitnv(line):
        # note that we can't just split on the "=" because of things like Variable_List, Resource_List, etc.
        name = line[:line.find('=')].strip()
        value = line[line.find('=')+1:].strip()
        return name,value
    fulljobid = lines[0].strip()
    jobinfo = {}
    for line in lines[1:]:
        if not line.strip(): continue
        if isstartline(line):
            name,value = splitnv(line)
            jobinfo[name] = value
        elif iscontinuationline(line):
            jobinfo[name] = jobinfo[name] + line.strip()
    # Now add in specific info that we want
    if 'Variable_List' in jobinfo:
        vl = {}
        parts = [i.strip() for i in jobinfo['Variable_List'].split(',')]
        for part in parts:
            vl[part.split('=')[0]] = part.split('=')[1]
        jobinfo['_env_vars'] = vl

    else:
        jobinfo['_env_vars'] = {}
    try:
        jobinfo['workdir'] = jobinfo['_env_vars']['PBS_O_WORKDIR']
    except KeyError:
        try:
            jobinfo['workdir'] = os.path.split(jobinfo['Error_Path'].split(':')[1])[0]
        except:
            jobinfo['workdir'] = 'UNKNOWN'
    # We must preserve order so that the head node is listed first
    nodes = []
    if 'exec_host' in jobinfo:
        nodes = [i.split('.')[0] for i in jobinfo['exec_host'].split('+')]
    _nodes = []
    for n in nodes:
        if n not in _nodes: _nodes.append(n)
    jobinfo['nodes'] = ' '.join(_nodes)
    if load:
        jobinfo['loads'] = {}
        for n in jobinfo['nodes'].split():
            # number of CPUs + 1 for header
            numlines = {'n':9, # nehalems
                        'h':9, # harpertowns
                        'o':5, # opterons
                        'c':9, # clovertowns
                        'b':17,# bigmem
                        }[n[0]]
            retcode,progout = run('/usr/bin/rsh','%s ps -eo pcpu,pid,user,args | sort -r -k1 | head -%s'%(n,numlines))
            jobinfo['loads'][n] = progout.split('\n')
    jobinfo['fulljobid'] = fulljobid
    jobinfo['jobid'] = jobid
    jobinfo['Last Modified'] = getlastmod(jobinfo['workdir'])
    jobinfo['Sim Time'] = getsimtime(jobinfo['workdir'])
            
    return jobinfo

#
# These should return True if they found a job in the directory, False otherwise.
# This is used to figure out which jobs have finished, etc.
#
# Setting "lenient" to true allows it to print as much as possible in case we're
# trying to diagnose a non-running directory.
#
def printcharmm(d,lenient=False,numtrajs=4):
    dynout = os.path.join(d,'dyn.out')
    if os.path.isfile(dynout) or lenient:
        junk = run('tail',('-1',dynout))
        retval,progout = run('ls',('-lhc',dynout))
        dynoutsize = progout.split()[4]
        dynoutmod  = ' '.join(progout.split()[5:8])
        print '      dyn.out size: %9s last modified: %s'%(dynoutsize,dynoutmod)
        trjs = glob.glob(os.path.join(d,'dy*trj')) + glob.glob(os.path.join(d,'*.dcd'))
        if trjs:
            for t in trjs:
                junk = run('tail',('-1',t))
            retval,progout = run('ls',['-sxt']+trjs)
            for line in progout.split('\n')[:numtrajs]:
                if not line.strip(): continue
                fsize,fname = line.strip().split()
                fname = os.path.split(fname)[-1]
                print '        %8s %s'%(fsize,fname)
        return True
    return False

def printamber(d,lenient=False,numtrajs=4):
        outs = glob.glob(os.path.join(d,'prod*out'))
        if not outs:
            return False
        for o in outs:
            junk = run('tail',('-1',o))
        retval,progout = run('ls',['-lhct']+outs)
        progout = progout.split('\n')[0]
        outsize = progout.split()[4]
        outmod  = ' '.join(progout.split()[5:8])
        outname = os.path.split(progout.split()[-1])[-1]
        print '      %11s size: %9s last modified: %s'%(outname,outsize,outmod)

        trjs = glob.glob(os.path.join(d,'prod*crd'))
        for t in trjs:
            junk = run('tail',('-1',t))
        retval,progout = run('ls',['-sxt']+trjs)
        for line in progout.split('\n')[:numtrajs]:
            if not line.strip(): continue
            fsize,fname = line.strip().split()
            fname = os.path.split(fname)[-1]
            print '        %8s %s'%(fsize,fname)
        return True

def printgromacs(d,lenient=False):
    mdcont = os.path.join(d,'mdcont.trr')
    if os.path.isfile(mdcont):
        junk = run('tail',('-1',mdcont))
        retval,progout = run('ls',('-lch',mdcont))
        mdcontsize = progout.split()[4]
        mdcontmod  = ' '.join(progout.split()[5:8])
        print '      mdcont.trr size: %9s last modified: %s'%(mdcontsize,mdcontmod)
    mdcont = os.path.join(d,'mdcont.xtc')
    if os.path.isfile(mdcont):
        junk = run('tail',('-1',mdcont))
        retval,progout = run('ls',('-lch',mdcont))
        mdcontsize = progout.split()[4]
        mdcontmod  = ' '.join(progout.split()[5:8])
        print '      mdcont.xtc size: %9s last modified: %s'%(mdcontsize,mdcontmod)
        return True
    return False

def printjob(jobid,jobinfo,user,lenient=False,load=False,numtrajs=4,verbose=True,fmt=None):
    for res in 'Job_Name Resource_List.nodes job_state resources_used.walltime nodes Resource_List.walltime queue'.split():
        if res not in jobinfo: jobinfo[res] = '???'
    if not verbose:
        # We want the list of nodes to fit nicely in our columns.
        # So we'll split it up into chunks of four

        nodechunks = [' '.join(i) for i in splitseq(jobinfo['nodes'].split(),4)]
        if not nodechunks:
            nodechunks = ["NO ASSIGNED NODE",]

        try:
            if jobinfo['Resource_List.walltime'].endswith(':00:00'):
                wallpart = '%s(%s)'% (jobinfo['resources_used.walltime'][:-3],jobinfo['Resource_List.walltime'][:-6])
            else:
                wallpart = '%s(%s)'% (jobinfo['resources_used.walltime'],jobinfo['Resource_List.walltime'])
            if True:
                statepart = '%s(%s)'%(jobinfo['queue'],jobinfo['job_state'])
            else:
                statepart = jobinfo['job_state']
            if True:
                homedir = os.path.expanduser('~%s/'%user)
                workdirpart = jobinfo['workdir'].replace(homedir,'~%s/'%user)
            else:
                workdirpart = jobinfo['workdir']

            print fmt%(jobid,jobinfo['Job_Name'],jobinfo['Resource_List.nodes'],
                       statepart,
                       wallpart,
                       #jobinfo['Last Modified'],jobinfo['Sim Time'], # MGL this is where you turn mod and sim on
                       nodechunks[0],
                       workdirpart.replace('/mounts/al-salam/software/','...'))
        except IndexError,KeyError:
            print "jobid",jobid
            print "jobinfo['Job_Name']",jobinfo['Job_Name']
            print "jobinfo['Resource_List.nodes']",jobinfo['Resource_List.nodes']
            print "jobinfo['job_state']",jobinfo['job_state']
            print "jobinfo['resources_used.walltime']",jobinfo['resources_used.walltime']
            print "nodechunks[0]",nodechunks[0]
            print "jobinfo['workdir']",jobinfo['workdir']
            print "jobinfo['Last Modified']",jobinfo['Last Modified'],
            print "jobinfo['Sim Time']",jobinfo['Sim Time'],
        except TypeError:
            print fmt
        for nc in nodechunks[1:]:
            print fmt%('','','','','','','',nc,'')
            
        found = True
    else:
        print fmt%(jobid,jobinfo['Job_Name'],jobinfo['Resource_List.nodes'],jobinfo['job_state'],jobinfo['resources_used.walltime'])
        print '    %s'%jobinfo['workdir']
        if jobinfo['nodes']: 
            print '    %s'%jobinfo['nodes']
            found = False
            if printcharmm(jobinfo['workdir'],lenient=lenient,numtrajs=numtrajs): found = True
            if printgromacs(jobinfo['workdir'],lenient=lenient): found = True
            if printamber(jobinfo['workdir'],lenient=lenient,numtrajs=numtrajs): found = True
            if load:
                for n in sorted(jobinfo['loads']):
                    print '    Load on node',n
                    for line in jobinfo['loads'][n]:
                        print '      '+line.strip()
        else:
            # If it's in the queue, we'll say "found" is true, because the
            # job itself might be responsible for some of the setup if
            # it's on step 1, etc.
            found = True
        print fmt%('-'*8,'-'*20,'-'*22,'-'*5,'-'*11,)
    return found

def printjobs(jobinfos,load,user,numtrajs,verbose,fmt,dashline,equalline,starline,titleline):
    """ Prints jobs and returns a dictionary containing directories that are running jobs"""
    print dashline
    print titleline
    print dashline
    # numerical sort
    idinfos = [(i[1],i[2]) for i in sorted([(int(j['jobid']),j['jobid'],j) for j in jobinfos])]
    founddirs = {}
    def sortvalue(jobinfo):
        try:
            return jobinfo['Job_Name']
        except KeyError:
            print "No Job_Name for",jobinfo
            return 1
        #return jobinfo['workdir']
    
    runningidinfos = [(sortvalue(jobinfo),jobid,jobinfo) for (jobid,jobinfo) in idinfos if 'exec_host' in jobinfo]
    runningidinfos.sort()
    runningidinfos = [(jobid,jobinfo) for (jobsort,jobid,jobinfo) in runningidinfos]
    queuedidinfos = [(sortvalue(jobinfo),jobid,jobinfo) for (jobid,jobinfo) in idinfos if 'exec_host' not in jobinfo]
    queuedidinfos.sort()
    queuedidinfos = [(jobid,jobinfo) for (jobsort,jobid,jobinfo) in queuedidinfos]
    if queuedidinfos:
        print equalline
        print 'RUNNING JOBS'
        print equalline
    for (jobid,jobinfo) in runningidinfos:
        if printjob(jobid=jobid,jobinfo=jobinfo,user=user,lenient=False,load=load,numtrajs=numtrajs,verbose=verbose,fmt=fmt):
            founddirs[jobinfo['workdir']] = True
    if queuedidinfos:
        print equalline
        print 'QUEUED JOBS'
        print equalline
    for (jobid,jobinfo) in queuedidinfos:
        if printjob(jobid=jobid,jobinfo=jobinfo,user=user,lenient=False,load=load,numtrajs=numtrajs,verbose=verbose,fmt=fmt):
            founddirs[jobinfo['workdir']] = True
    return founddirs

def printjobstatdirs(founddirs,user,verbose,fmt,dashline,equalline,starline,titleline):
    def normalizedirs(dirs):
        _dirs = {}
        for d in dirs:
            d = os.path.expanduser(d)
            if not d.endswith('/'): d = d + '/'
            _dirs[d] = True
        return _dirs
    founddirs = normalizedirs(founddirs)
                
    jdfile = os.path.expanduser('~%s/.jobstat-dirs'%user)
    if os.path.isfile(jdfile):
        expecteddirs = []
        nondirs = [] #sometimes the wildcard matches files.
        for line in file(jdfile):
            line = line.strip()
            if not line: continue
            if line.startswith('#'): continue
            if line[0] not in ('/','~'):
                line = '~%s/'%user + line
            dirs = glob.glob(os.path.expanduser(line.strip()))
            if not dirs:
                nondirs.append(line.strip())
            for d in dirs:
                if os.path.isdir(d):
                    expecteddirs.append(d)
                else:
                    nondirs.append(d)
        def addslash(x):
            _x = []
            for i in x:
                if i.endswith('/'):
                    _x.append(i)
                else:
                    _x.append(i+'/')
            return _x
        expecteddirs = addslash(expecteddirs)
        founddirs = addslash(founddirs)
        nondirs = addslash(nondirs)
        missingdirs = [d for d in expecteddirs if d not in founddirs]
        if missingdirs:
            print starline
            print "**** No job is running in the following directories."
            print "**** Perhaps the jobs have completed?"
            print starline

            for d in missingdirs:
                jobinfo = {'workdir':d,
                           'jobid':'None',
                           'Last Modified': getlastmod(d),
                           'Sim Time': getsimtime(d),
                           
                           }
                printjob('None',jobinfo,user=user,lenient=True,verbose=verbose,fmt=fmt)
        if nondirs:
            print "NOTE: ~/.jobstat-dirs listed the following files"
            print "that are not directories:"
            for d in nondirs:
                print "    "+d
            print dashline


def printqueuesummary(jobinfos):

    results = {}
    for jobinfo in jobinfos:
        if 'exec_host' not in jobinfo: continue
        # 8:nehalem:ppn=8
        # 16:ppn=8:htown:ipath2
        resources = jobinfo['Resource_List.nodes']
        #nodes = int(resources.split(':')[0])
        nodes = jobinfo['nodes'].split()
        if type(nodes) in [type(()),type([])]:
            nodeset = set(nodes)
        else:
            nodeset = set((nodes,))
        ppn = None
        nodetype = None
        for i in resources.split(':'):
            if 'ppn=' in i:
                ppn = int(i.split('=')[-1])
        
        for nt in 'htown htown:ipath htown:ipath2 nehalem ctown o2200 bar2300'.split():
            if nt in resources:
                nodetype = nt
        if not nodetype:
            nodetype='generic'
        if not ppn:
            ppn = {'ctown':8,'htown:ipath':8,'htown:ipath2':8,'htown':8,
                   'o2200':4,'bar2300':16,'nehalem':8}[nodetype]
        if nodetype not in results:
            results[nodetype] = [nodeset,ppn*len(nodes)]
        else:
            results[nodetype][0] = results[nodetype][0].union(nodeset)
            results[nodetype][1] = results[nodetype][1] + ppn*len(nodes)
    for nt in results:
        print "%-15s: you are using %3s nodes and %4s CPUs"%(nt,len(results[nt][0]),results[nt][1])

    # TODO: FIXME: factor, obviously
    results = {}
    for jobinfo in jobinfos:
        if 'exec_host' in jobinfo: continue
        # 8:nehalem:ppn=8
        # 16:ppn=8:htown:ipath2
        resources = jobinfo['Resource_List.nodes']
        nodes = int(resources.split(':')[0])
        for i in resources.split(':'):
            if 'ppn=' in i:
                ppn = int(i.split('=')[-1])
        
        for nt in 'htown htown:ipath htown:ipath2 nehalem ctown o2200 bar2300'.split():
            if nt in resources:
                nodetype = nt
        if nodetype not in results:
            results[nodetype] = [nodes,ppn*nodes]
        else:
            results[nodetype][0] = results[nodetype][0] + nodes
            results[nodetype][1] = results[nodetype][1] + ppn*nodes
    for nt in results:
        print "%-15s: you've queued %3s nodes and %4s CPUs"%(nt,results[nt][0],results[nt][1])


def qstatchunks(progout):
    lines = progout.split('\n')
    assert lines[4].startswith('-----') # Make sure the first 5 lines are a header
    assert lines[-1].strip() == '' # Last is a blank line we can ignore
    # Now the rest of the lines come in chunks of three
    assert divmod(len(lines) - 5 -1,3)[1] == 0
    for (jobline,nodes,comment) in splitseq(lines[5:-1],3):
        yield jobline, nodes, comment

def jobstat(user,load,numtrajs,includejobstatdirs,clusterid='.lo0.la',verbose=False):
    """Notes from the graveyard: we used to ask for 'clusterid'
    ).m1.lobos for NIH, .as0.al for Earlham, etc). We used it to pick
    out the lines from qstat -ns that actually describe jobs, rather
    than separator lines or descriptions. Instead, I'm now just
    parsing the results all together, assuming there's a common
    format.

    """
    #(retcode,progout) = run('qstat',('-nsu',user))
    (retcode,progout) = run('qstat',('-ns'))
    # first 5 lines are header info
    joblines = [jobline for (jobline,nodes,comment) in qstatchunks(progout)]
    #joblines = [i.split()[0] for i in progout.split('\n')[5:] if i.split() and clusterid in i.split()[0]]
    #print "JOBLINES",joblines
    jobids = [i.split('.')[0] for i in joblines]
    jobinfos = [getjobinfo(jobid=i,load=load) for i in jobids]
    #print "JOBINFOS",jobinfos,"JOBIDS",jobids,"JOBLINES",joblines,"CLUSTERID",clusterid

    if verbose:
        cols = ('Job ID','Job_Name    ','Node Desc.','State','Walltime')
    else:
        #cols = ('Job ID','Job_Name                    ','Node Desc.','Queue(S)','Walltime    ','Modified', 'SimTime ', 'Nodes','Working Directory')
        cols = ('Job ID','Job_Name                    ','Node Desc.  ','Queue(S)','Walltime    ', 'Nodes','Working Directory')
        
    fmt = ' '.join(['%%-%is'%len(i) for i in cols])
    titleline = fmt%cols
    dashline = fmt%tuple(['-'*len(i) for i in cols])
    equalline = fmt%tuple(['='*len(i) for i in cols])
    starline = fmt%tuple(['*'*len(i) for i in cols])

    founddirs = printjobs(jobinfos=jobinfos,load=load,user=user,numtrajs=numtrajs,verbose=verbose,fmt=fmt,dashline=dashline,equalline=equalline,starline=starline,titleline=titleline)
    if includejobstatdirs:
        printjobstatdirs(founddirs,user,verbose=verbose,fmt=fmt,dashline=dashline,equalline=equalline,starline=starline,titleline=titleline)
    printqueuesummary(jobinfos)
                
                
            
if __name__ == '__main__':
    import optparse
    usage = """This expects you to use 

 - Rick's scrips to run CHARMM jobs
 - Michael's scripts to run GROMACS jobs
 - Asim's scripts to run AMBER(pmemd) jobs

The first node reported is the head node.

You may provide a list of directories in ~/.jobstat-dirs. It can
contain wildcards. If any of those directories is *not* running a job
at the, moment an extra message will be printed. Lines in
~/.jobstat-dirs that begin with a # will be ignored.

I find that it's more convenient to list the jobs in alphabetical
order of their working directory. Let me know if you want them listed
some other way (e.g. in order of their jobid).

Running jobs are reported first, followed by queued jobs.
"""
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-u','--user',default=os.environ['USER'],
                      help='User [default %default]',
                      )
    parser.add_option('-l','--load',default=False,action='store_true',
                      help='Also display load information for each node. Could hang if we cannot connect to the node. [default %default]',
                      )
    parser.add_option('-n','--numtrajs',default=4,
                      help='Number of trajectories to show per directory.',
                      type='int',
                      )
    parser.add_option('-J','--ignore-jobstat-dirs',default=True,action='store_false',dest='includejobstatdirs',
                      help='Ignore the .jobstat-dirs file.',
                      )
    parser.add_option('-v','--verbose',default=False,action='store_true',dest='verbose',
                      help='Verbose (multiline) output',
                      )
    parser.add_option('-t','--terse',action='store_false',dest='verbose',
                      help='Terse (single line) output',
                      )
    options,args = parser.parse_args()
    jobstat(user=options.user,load=options.load,numtrajs=options.numtrajs,includejobstatdirs=options.includejobstatdirs,verbose=options.verbose)
