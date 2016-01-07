#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
import sys,os


DEFAULTS = {'annealing':'no',
            'comm_grps':'system',
            'comm_mode':'linear',
            'continuation':'no',
            'cos_acceleration':'0',
            'delta_lambda':'0',
            'disre':'no',
            'disre_weighting':'conservative',
            'disre_fc':'1000',
            'disre_mixed':'no',
            'disre_tau':'0',
            'dt':'0.001',
            'epsilon_rf':'0',
            'epsilon_surface':'0',
            'ewald_rtol':'1e-05',
            'fourier_nx':'0',
            'fourier_ny':'0',
            'fourier_nz':'0',
            'fourierspacing':'0.12',
            'free_energy':'no',
            'gen_seed':'173529',
            'gen_temp':'300',
            'gen_vel':'no',
            'init_lambda':'0',
            'morse':'no',
            'ns_type':'grid',
            'nstcomm':'1',
            'nstdisreout':'100',
            'nstenergy':'100',
            'nsteps':'0',
            'nstlog':'100',
            'nstvout':'100',
            'nstxout':'100',
            'nstxtcout':'0',
            'nstfout':'0',
            'optimize_fft':'no',
            'pbc':'xyz',
            'pcoupl':'no',
            'pcoupltype':'isotropic',
            'periodic_molecules':'no',
            'pme_order':'4',
            'sc_alpha':'0',
            'sc_sigma':'0.3',
            'shake_tol':'0.0001',
            'tau_p':'1',
            'unconstrained_start':'no',
            'userint1':'0',
            'userint2':'0',
            'userint3':'0',
            'userint4':'0',
            'userreal1':'0',
            'userreal2':'0',
            'userreal3':'0',
            'userreal4':'0',
            'xtc_grps':'System',

            # no default provided in manual
            'cpp':'',
            'dispcorr':'no',
            'ref_p':'',
            'ref_t':'',
            'tau_t':'',
            'tc_grps':'',
            'domain_decomposition':'',
            'tinit':'0.0',
            'title':'',
            'energygrps':'',
            'zero_temp_time':'',
            }
def mdpparams(fn):
    f = open(fn)
    result = {}
    for line in f:
        if line.find(';') >= 0:
            line = line[:line.find(';')]
        line = line.strip()
        if not line: continue
        if line.startswith(';'): continue
        try:
            k,v = [i.strip() for i in line.split('=')]
        except ValueError:
            print(fn,line)
            raise
        k = k.lower().replace('-','_')
        v = v.lower()
        v = ' '.join(v.split())
        #if k in DEFAULTS and v == DEFAULTS[k]: v = ''
        if k in DEFAULTS and v == '': 
            v = DEFAULTS[k]
        if v:
            result[k] = v
    return result
def printparm(p,m1,m2):
    pm1,pm2 = None,None
    plen = 0
    if p not in m1: m1[p] = DEFAULTS[p]
    pm1 = m1[p]
    plen = max(plen,len(pm1))
    if p not in m2: m2[p] = DEFAULTS[p]
    pm2 = m2[p]
    plen = max(plen,len(pm2))
    if not pm1:        pm1 = '_'*plen
    if not pm2:        pm2 = '_'*plen
    if pm1 != pm2:
        print('%20s : %20s : %20s'%(p,pm1,pm2))
def compare(fn1,fn2):
    m1 = mdpparams(fn1)
    m2 = mdpparams(fn2)
    params = list(set(m1.keys()).union(set(m2.keys())))
    params.sort()
    print('%20s : %20s : %20s'%(' ',fn1,fn2))
    print('%20s : %20s : %20s'%('-'*20,'-'*20,'-'*20,))
    for p in params:
        printparm(p,m1,m2)
if __name__ == '__main__':
    compare(sys.argv[-2],sys.argv[-1])
