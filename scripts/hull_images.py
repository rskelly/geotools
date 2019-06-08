#!/usr/bin/env python3

'''
Loads data from a csv file with format,

c,r,slope,yint,[...]

where c and r are the column and row in the source image, slope and y-int are the coefficients of the
regression and [...] are the coordinates of a convex hull. Coordinates are delineated by :.

The script draws a hull with the regression line for each column/row.
'''

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import random 

csvfile = sys.argv[1]
outdir = sys.argv[2]

try:
    os.makedirs(outdir)
except: pass

def draw(c, r, slope, yint, coords):
    minx = miny = 99999999.
    maxx = maxy = -99999999.
    for x, y in coords:
        if x < minx: minx = x
        if x > maxx: maxx = x
        if y < miny: miny = y
        if y > maxy: maxy = y
    width = maxx - minx
    height = maxy - miny

    coords = np.array(coords)
    
    ax = plt.subplot()
    ax.plot(coords[...,0], coords[...,1])

    x1 = minx
    y1 = slope * minx + yint
    x2 = maxx
    y2 = slope * maxx + yint
    
    ax.plot((x1, x2), (y1, y2))
    
    #plt.legend()
    plt.savefig(os.path.join(outdir, 'hull_{c}_{r}.png'.format(c=c,r=r)))
    plt.close()

with open(csvfile, 'r') as f:
    queue = []
    head = f.readline()
    line = f.readline()
    while line:
        try:
            c, r, slope, yint, coordss = line.split(',')
            #coords = [list(map(float, x.split(':'))) for x in coords.split(';')]
            coords = []
            try:
                for xy in coordss.split(';'):
                    coords.append(list(map(float, xy.split(':'))))
            except: pass
            queue.append((int(c), int(r), float(slope), float(yint), coords))
        except Exception as e:
            print(e)
            break
        line = f.readline()
    
    random.shuffle(queue)
    for q in queue[:10]:
        draw(*q)
