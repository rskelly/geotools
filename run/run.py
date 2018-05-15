#!/usr/bin/env python

'''
This script runs the contrem program with arguments supplied by a CSV file. The csv file 
will have the columns, 

    data_file (-d), roi_file (-r), band_map (-b), bm_header (-z), bm_wl_col (-w), bm_idx_col (-i), output_template (-o), low_wl (-l), high_wl (-h), population (-p), driver (-v), extension (-e)

The fields map to inputs for the contrem program,

 -d A GDAL-readable data file containing spectral samples; can contain any number of bands >= 2.
 -r An ENVI ROI text file.
 -b A CSV file containing a mapping from wavelength to (1-based) band index.
 -w An integer giving the (0-based) column index in -b which contains wavelengths.
 -i An integer giving the (0-based) column index in -b which contains the band indices.
 -z If given, indicates the presence of a header in the band map that must be skipped.
 -o An output file template. This is a filename with no extension that will be modified as
    appropriate. Parent directories will be created.
 -l The minimum wavelength to consider.
 -h The maximum wavelength to consider.
 -s The size of the buffer. Default is 256. Larger buffers are possible, but one must
    consider that multiple buffers may be in memory at once.
 -t The number of threads to use. Default 2.
 -p By default, sample statistics are used. This flag forces the use of
    population statistics.
 -v The driver to use for output rasters. Defaults to ENVI, but any GDAL-writable
    format will do.
 -e File extension for raster files. Defaults to .dat for ENVI files.

'''

import csv
import sys
import os
from subprocess import Popen

fields = ['data_file', 'roi_file', 'band_map', 'bm_header', 'bm_wl_col', 'bm_idx_col', 'output_template', 'low_wl', 'high_wl', 'population', 'driver', 'extension']


def run(configfile):
    '''
    Run the configured jobs.
    '''
    
    config = []
    
    with open(configfile, 'rU') as f:
        db = csv.reader(f)
        head = db.next()
        for row in db:
            config.append(dict(zip(head, row)))
    
    for job in config:
        
        print('Running with', job)
        
        params = ['contrem']
        
        if job.get('data_file', False):
            params.extend(('-d', job['data_file']))
        elif job.get('roi_file'):
            params.extend(('-r', job['roi_file']))
        else:
            print('No input file configured for this job.')
            continue
        
        if job.get('output_template', False):
            params.extend(('-o', job['output_template']))
        else:
            print('An output template is required.')
            continue

        if job.get('driver', False) and job.get('extension', False):
            params.extend(('-v', job['driver'], '-e', job['extension']))
        else:
            print('An output driver and extension are required.')
            
        if job.get('band_map', False):
            params.extend(('-b', job['band_map']))
            if job.get('bm_header', False) and job.get('bm_wl_col', False) and job.get('bm_idx_col', False):
                if job['bm_header'] == 't':
                    params.append('-z')
                params.extend(('-w', job['bm_wl_col'], '-i', job['bm_idx_col']))
            else:
                print('If a band map is given, the header, wl and idx col filds must be set.')
                continue
        
        if job.get('low_wl', False) and job.get('high_wl', False):
            params.extend(('-l', job['low_wl'], '-h', job['high_wl']))
            
        if job.get('population', 'f') == 't':
            params.append('-p')
    
        Popen(params, cwd=os.getcwd()).wait()

        
if __name__ == '__main__':
    
    config = sys.argv[1]
    
    run(config)