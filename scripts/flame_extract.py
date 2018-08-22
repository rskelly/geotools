#!/usr/bin/env python

from datetime import datetime
from calendar import timegm
from getopt import getopt
import math 
import sys 
from uuid import _last_timestamp

# Extracts rows from a spectral file given some user-specified criteria.

def usage():
    print('''Usage: flame_extract.py <options>
        -i <file>        The input file.
        -o <file>        The output file.
        -b               The start time of extraction, in the format yyyy-mm-dd hh:mm:ss.
                         The argument is optional, and the date part is optional if the time is given. 
        -e               The end time of extraction, same format as -s0.
        -d <interval>    The extraction interval in seconds.
        ''')
    
def parse_time(timestr):
    '''
    Take a time/date string in the form YYYY-MM-DD HH:MM:SS.SSS, with or without
    the date part, and parse it into a time stamp -- the number of seconds since the
    epoch.
    ''' 
    try:
        timestamp = datetime.strptime(timestr, '%Y-%m-%d %H:%M:%S.%f')
    except Exception as e:
        print(e)
        try:
            timestamp = datetime.strptime(timestr, '%H:%M:%S.%f')
        except Exception as e:
            print(e)
    return timegm(timestamp.timetuple())
       
 
def extract(infile, outfile, params):
    start_time = 0
    end_time = float('inf')
    interval = 1
    
    #print(params)
    if params.get('-b', False):
        start_time = parse_time(params['-b'])
    if params.get('-e', False):
        end_time = parse_time(params['-e'])
    if params.get('-d', False):
        interval = int(params['-d'])
        
    print(params)
    print(start_time, end_time, interval)
    
    header = False # false if the header has not been finished
    last_timestamp = 0
    
    with open(outfile, 'w') as output:
        with open(infile, 'rU') as input:
            
            line = input.readline()
            while line:
                try:
                    if not header:
                        
                        output.write(line)
                        if line[:3] == '>>>':
                            output.write(input.readline())
                            header = True
                            
                    elif line[0]: # skip the line if there's nothing in it
                         
                        write = False
                        
                        line0 = line.split('\t')
                        timestamp = parse_time(line0[0])
                        timestamp_i = int(timestamp / interval)

                        if timestamp >= start_time and timestamp <= end_time:
                            if timestamp_i != last_timestamp:
                                write = True
                                
                        if write:
                            output.write(line)

                        #print(line0[0], timestamp, timestamp_i, last_timestamp)
                        last_timestamp = timestamp_i
                         
                except Exception as e:
                    print(e)
                    break
                
                line = input.readline()

if __name__ == '__main__':
    
    try:
        optlist, args = getopt(sys.argv[1:], 'b:e:i:o:d:')
        print(optlist)
        params = dict(optlist)
    
        extract(params['-i'], params['-o'], params)        
    except Exception as e:
        print(e)
        usage()