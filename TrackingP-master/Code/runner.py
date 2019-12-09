#!/usr/bin/env python
from subprocess import call
from sys import exit

# number of snapshots
ns=400
# final time
tf=1000

# number of mesh elements

# maximum depth
dMax=20.0

# output directory
dir='SpatialResults'

movieOnly = False

for K0S in ((2.5, 5),):
    name='trackP-Pt-%s-K0-%s' % (0.7, K0S[0])
    snapInterval = K0S[1]
    K0 = K0S[0]
    kappa_u= 0.0004
    kappa_bg = 0.3
    nt = snapInterval*ns
    nx=nt
    # run simulation
    if not movieOnly:
        argPairs = {
        'dir' : dir,
        'dMax' : dMax,
        'K0' : K0,
        'kappa_u': kappa_u,
        'kappa_bg': kappa_bg,
        'tf' : tf,
        'nt' : nt,
        'snap' : snapInterval,
        'nx' : nx,
        'o' : '%s-t' % name
        }
        args = ['./trackP.exe']
        for k in argPairs.keys():
            args.append('--%s=%s' % (k,argPairs[k]))
        print (args)
        rtnCode = call(args)
        if rtnCode != 0:
            print ('Simulation failed!')
            exit(-1)
    # produce movie
    movArgPairs = {
        'dir' : dir,
        'dMax' : dMax,
        'K0' : K0,
        'nt' : ns,
        'tf' : tf,
        'f' : '%s-t' % name,
        'o' : name,
        'title' : '$P_t=%g$,  $K_0=%g$' % (0.7, K0) 
    }
    args = ['./makeMovieSig0.py']
    for k in movArgPairs.keys():
        args.append('--%s=%s' % (k,movArgPairs[k]))
    print (args)
    rtnCode = call(args)
    if rtnCode != 0:
        print ('Animation failed!')
        exit(-1)


        
            
            
        
        
