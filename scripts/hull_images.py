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
import cairo

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
    w = 200.
    h = 200.
    xs = 200. / width
    ys = 200. / height
    hbuf = 20
    path = os.path.join(outdir, 'hull_{c}_{r}.svg'.format(c=c, r=r))
    print(width, height, w, h, xs, ys, hbuf, path)
    with cairo.SVGSurface(path, w + 2 * hbuf, h + 2 * hbuf) as surf:
        ctx = cairo.Context(surf)
        #ctx.scale(width, height)
        ctx.translate(-200., 0)
        
        x, y = coords[0]
        ctx.move_to(x * xs + hbuf, y * ys + hbuf)
        for x, y in coords[1:]:
            ctx.line_to(x * xs + hbuf, y * ys + hbuf)
        ctx.stroke()
        
        ctx.move_to(minx * xs + hbuf, (slope * minx + yint) * ys + hbuf)
        ctx.line_to(maxx * xs + hbuf, (slope * maxx + yint) * ys + hbuf)
        print(minx * xs + hbuf, (slope * minx + yint) * ys + hbuf)
        print(maxx * xs + hbuf, (slope * maxx + yint) * ys + hbuf)
        ctx.stroke()
        
with open(csvfile, 'r') as f:
    head = f.readline()
    line = f.readline()
    i = 0
    while line:
        try:
            c, r, slope, yint, coordss = line.split(',')
            #coords = [list(map(float, x.split(':'))) for x in coords.split(';')]
            coords = []
            try:
                for xy in coordss.split(';'):
                    coords.append(list(map(float, xy.split(':'))))
            except: pass
            draw(int(c), int(r), float(slope), float(yint), coords)
            if i > 10:
                break
            i += 1
        except Exception as e:
            print(e)
            break
        line = f.readline()
        