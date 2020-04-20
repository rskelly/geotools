#!/usr/bin/env python3

'''
Merges outputs of the pc2grid program. Sums the first count layer, then
computes the weighted mean of remaining layers using the counts.

The files should overlap according to the order that results
from a lexicographic sort on the filenames.

If there's no count layer, assume 1 for each weight.

Only works with UTM.
'''

import sys
import gdal
import os
import math
import numpy as np

def get_bounds(infiles):

    bounds = [math.inf, math.inf, -math.inf, -math.inf]
    
    for f in infiles:
        ds = gdal.Open(f)
        drv = ds.GetDriver()
        bands = ds.RasterCount
        trans = ds.GetGeoTransform()
        cols = ds.RasterXSize
        rows = ds.RasterYSize
        proj = ds.GetProjection()
        resx = trans[1]
        resy = trans[5]
        
        x1 = trans[0]
        y1 = trans[3]
        x2 = x1 + cols * resx
        y2 = y1 + rows * resy
        
        bounds[0] = min(bounds[0], min(x1, x2))
        bounds[1] = min(bounds[1], min(y1, y2))
        bounds[2] = max(bounds[2], max(x1, x2))
        bounds[3] = max(bounds[3], max(y1, y2))
    
    cols = int((bounds[2] - bounds[0]) / abs(resx))
    rows = int((bounds[3] - bounds[1]) / abs(resy))
    
    return bounds, resx, resy, cols, rows, drv, bands, proj

def to_x(c, trans):
    return trans[0] + c * trans[1] + 0.5 * trans[1]

def to_y(r, trans):
    return trans[3] + r * trans[5] + 0.5 * trans[5]

def to_c(x, trans):
    return int((x - trans[0]) / trans[1])

def to_r(y, trans):
    return int((y - trans[3]) / trans[5])

def merge(outfile, infiles, has_counts = True):
    
    #infiles.sort(key = lambda a: os.path.basename(a))

    # Get the overall bounds and transform for inputs.
    bounds, resx, resy, cols, rows, drv, bands, proj = get_bounds(infiles)
    trans = [bounds[0] if resx > 0 else bounds[2], resx, 0., bounds[1] if resy > 0 else bounds[3], 0., resy]
    
    # Create the new raster.
    ds = gdal.GetDriverByName("GTiff").Create(outfile, cols, rows, bands if has_counts else bands + 1, gdal.GDT_Float32)
    ds.SetGeoTransform(trans)
    ds.SetProjection(proj)
        
    # 2D arrays for the value sums and counts.
    sums = np.zeros(cols * rows, dtype = np.float32).reshape((rows, cols))
    counts = np.zeros(cols * rows, dtype = np.float32).reshape((rows, cols))
    
    # Output band
    ob = 0

    for b in range(1 if has_counts else 0, bands):
        print('Band', b + 1, 'of', bands)

        sums.fill(0.)
        counts.fill(0)

        for f in infiles:
            
            # Get the properties of the current band for the current file.
            bounds1, resx1, resy1, cols1, rows1, drv1, bands1, proj1 = get_bounds([f])
            ds1 = gdal.Open(f)
            trans1 = ds1.GetGeoTransform()
            band1 = ds1.GetRasterBand(b + 1)
            nd1 = band1.GetNoDataValue()
            
            for r1 in range(rows1):

                # Transform the local coordinate to global.
                x1 = to_x(0, trans1)
                y1 = to_y(r1, trans1)
                c = to_c(x1, trans)
                r = to_r(y1, trans)
                
                # The length of the out grid buffer.
                colst = min(cols - c, cols1)
                vbuf = band1.ReadAsArray(0, r1, colst, 1)
                idx = vbuf[0,...] != nd1

                # Create a counts buffer or a list of 1.
                if has_counts:
                    cbuf = ds1.GetRasterBand(1).ReadAsArray(0, r1, colst, 1)
    
                # Increment the counts array by the count layer value, or by 1.
                counts[r,c:c + colst][idx] += (cbuf[0,idx] if has_counts else 1)
                # Multiply the non-nd values by the count and add to sums.
                sums[r,c:c + colst][idx] += vbuf[0,idx] * (cbuf[0,idx] if has_counts else 1)

        nd = -9999.
        band = ds.GetRasterBand(ob + 1)
        band.SetNoDataValue(nd)
        
        # Write the counts on the first band.
        if ob == 0:
            band.WriteArray(counts, 0, 0)
            ob += 1
            band = ds.GetRasterBand(ob + 1)
    
        # Calculate weighted mean and set nodata and write.
        sums[counts > 0] /= counts[counts > 0]
        sums[counts == 0] = nd
        band.WriteArray(sums, 0, 0)
        ob += 1

if __name__ == '__main__':

    try:    
        has_counts = True
        outfile = None
        infiles = None
        for i in range(1, len(sys.argv)):
            if sys.argv[i] == '-c':
                has_counts = False
            elif not outfile:
                outfile = sys.argv[i]
            else:
                infiles = sys.argv[i:]
                break
        merge(outfile, infiles, has_counts)
    except:
        import traceback
        print(traceback.format_exc())
        print('Usage: pc2grid_merge.py <outfile> <infiles*>')

