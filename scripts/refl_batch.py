#!/usr/bin/env python3

import os
import sys
import re

config = [
    {
        'irrad': '/media/rob/robdata/work/hlrg/field/2018_oct_2/sf10_oct_AbsoluteIrradiance_11-36-32-032_conv.txt',
        'fi_dir': '/media/rob/robdata/work/hlrg/field/2018_oct_2/',
        'rad_dir': '/media/rob/robdata/work/hlrg/field/2018_oct_2/rad_MNFinv',
        'imu_gps': '/media/rob/robdata/work/hlrg/field/2018_oct_2/imu_gps.txt',
        'refl_out_dir': '/media/rob/robdata/work/hlrg/field/2018_oct_2/output',
        'imu_offset': 0.,
        'irrad_offset': 0.
    }
]

def expand_config(conf):
    
    fi_dir = conf['fi_dir']
    rad_dir = conf['rad_dir']
    fi = [x for x in os.listdir(fi_dir) if x.startswith('frameIndex')]
    rad = [x for x in os.listdir(rad_dir) if x.startswith('raw_') and x.endswith('.dat')]
    
    ids = {}
    
    pat = re.compile('frameIndex_([0-9]+)\..*')
    for f in fi:
        id = pat.search(f).group(1)
        ids[id] = [os.path.join(fi_dir, f)]
        
    pat = re.compile('raw_([0-9]+)_.*')
    for r in rad:
        id = pat.search(r).group(1)
        try:
            ids[id].append(os.path.join(rad_dir, r))
        except Exception as e: 
            print(e)
        
    pairs = []
    for id, pair in ids.items():
        if len(pair) == 2:
            pair.append(id)
            pairs.append(tuple(pair))
    
    conf['files'] = pairs
    
def run(config):
    
    for conf in config:
        
        expand_config(conf)
        
        if not os.path.exists(conf['refl_out_dir']):
            os.makedirs(conf['refl_out_dir'])
            
        for fi, raw, id in conf['files']:
        
            ofile = os.path.join(conf['refl_out_dir'], 'refl_{id}.tif'.format(id = id))
            
            cmd = 'reflectance -i {i} -io {io} -r {r} -f {f} -c {c} -co {co} -o {o}'.format(
                i = conf['imu_gps'], io = conf['imu_offset'], r = raw, f = fi, c = conf['irrad'], co = conf['irrad_offset'], o = ofile
            )

            os.system(cmd)

if __name__ == '__main__':
    
    run(config)
