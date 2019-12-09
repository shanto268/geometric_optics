#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as ani

parser = argparse.ArgumentParser(description="Animator for Algae-Daphnia")
parser.add_argument('--nt', action='store', default=2000, help='number of timesteps')
parser.add_argument('--o', action='store', default='test', help='output filename')
parser.add_argument('--f', action='store', default='test', help='input filename')
parser.add_argument('--dir', action='store', default='SpatialResults', help='data directory')
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


DATA0 = readfile(dirname, "%s-%d" % (fname, 0))


FFMpegWriter = ani.writers['ffmpeg']
metadata = dict(title='Bio', artist='Math')
writer = FFMpegWriter(fps=15, metadata=metadata)

fig, ax1 = plt.subplots()
ax1.set_xlim(0, 20)
ax1.set_ylim(0,1)



ax1.set_title(title)
ax1.set_xlabel('Depth (m)')
ax1.set_ylabel('Population')

ax2 = ax1.twinx()
ax2.set_xlim(0, 20)
ax2.set_ylabel('Light intensity')
ax2.set_ylim(0,2)

A, = ax1.plot([], [], 'r-')
D1, = ax1.plot([], [], 'k-')
Q, = ax1.plot([], [], 'g-')
Pf, = ax1.plot([], [], 'b-')
K, = ax2.plot([], [], 'm-')


plt.legend([A,D1,Q,Pf,K], ['Algae', 'Daphnia','Q', '$P_f$','Light'], loc='upper right')

print ('writing to %s' % outname)
with writer.saving(fig, "%s.mp4" % outname, nt):
    for s in range(nt+1):
        DATA = readfile(dirname, "%s-%d" % (fname, s))
        print ('frame %d, filename=%s/%s-%d' % (s, dirname, fname, s))
        x=DATA[:,0];
        A.set_data(x, DATA[:,1])
        D1.set_data(x, DATA[:,2])
        Q.set_data(x, DATA[:,6])
        Pf.set_data(x, 10*DATA[:,5])
        K.set_data(x, DATA[:,3])
        writer.grab_frame()
