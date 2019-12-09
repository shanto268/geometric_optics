#!/usr/bin/env python

import argparse
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from matplotlib import cm
import matplotlib.ticker as tkr

parser = argparse.ArgumentParser(description="Animator for Algae-Daphnia")
parser.add_argument('--nt', action='store', default=2000, help='number of timesteps')
parser.add_argument('--o', action='store', default='test', help='output filename')
parser.add_argument('--f', action='store', default='test', help='input filename')
parser.add_argument('--dir', action='store', default='SpaceTimePlotResults', help='data directory')
parser.add_argument('--tf', action='store', default=1000, help='final time')
parser.add_argument('--title', action='store', default='', help='title for plot')
parser.add_argument('--K0', action='store', default=2.0, help='y axis scale')
parser.add_argument('--dMax', action='store', default=20.0, help='x axis scale')

args = parser.parse_args()

nt = int(args.nt)
fname = args.f
dirname = args.dir
outname = args.o
yMax = float(args.K0)
xMax = float(args.dMax)
title = args.title
tf = float(args.tf)

def readfile(dirname, fname):
    
    fn = "%s/%s" % (dirname, fname)
    return np.loadtxt(fn)
tupleData=[]
for K0 in (0.5, ):
    name ='trackP-Pt-%s-K0-%s' % (0.6, K0)
    for t in range(401):
        DATA = readfile(dirname, '%s-t-%d'%(name, t))
        ln = len(DATA[:,0])
        for j in range(ln):
            tupleData.append((t/4.0,DATA[j,0],DATA[j,1]),)
    
    x, y, z = zip(*tupleData)
    z= map(float, z)
    grid_x, grid_y = np.mgrid[min(x):max(x):100j, min(y):max(y):100j]
    grid_z = griddata((x, y), z, (grid_x, grid_y), method='cubic')

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(grid_x, grid_y, grid_z, rstride=1, cstride=1, cmap=cm.jet, linewidth = 0, antialiased= False)
    
    ax.set_title('$P_t=%g$, $K0=%s$'%(0.6, K0))
    ax.set_xlabel('Time (day)', fontsize = 13)
    ax.set_ylabel('Depth (m)', fontsize = 13)
    ax.set_zlabel('Algae density (mg C/L)', fontsize = 13)
    fig.colorbar(surf, shrink=0.1, aspect=2)
    fig.savefig('%s.pdf'%(name))
       
