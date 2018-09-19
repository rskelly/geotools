#!/usr/bin/env python3

import sys
import traceback

def sizeof(data_type):
    if data_type == 'float32' or data_type == 'int32' or data_type == 'uint32':
        return 4
    elif data_type == 'float64' or data_type == 'int64' or data_type == 'uint64':
        return 8
    elif data_type == 'int16' or data_type == 'uint16':
        return 2
    elif data_type == 'char' or data_type == 'byte' or data_type == 'int8' or data_type == 'uint8':
        return 1
    else:
        raise Exception('Unknown data type : ' + data_type)

def process(input_file, output_file, input_cols, left_edge, right_edge, data_type):
    input_cols = int(input_cols)
    left_edge = int(left_edge)
    right_edge = int(right_edge)
    data_size = sizeof(data_type)
    with open(output_file, 'wb') as output:
        with open(input_file, 'rb') as input:
            while True:
                row = input.read(input_cols * data_size)
                if not row:
                    break
                output.write(row[left_edge * data_size : (input_cols - right_edge) * data_size])  

def usage():
    print('''
        This program trims the edges off of a raw image whose bytes 
        comprise the sequential rows of an image. The user provides the input file,
        the number of columns, the data type and the width of each edge to be trimmed,
        and a new raw file is produced with the shorter rows.
        
        The data types understood are:
        
        - float32, int32, uint32 - 4 bytes per pixel;
        - float64, int64, uint64 - 8 bytes per pixel;
        - int16, uint16 - 2 bytes per pixel;
        - char, byte, int8, uint8 - 1 byte per pixel.
        
        The program doesn't actually care about how pixel values are representation, only
        how large the representations are.
        
        The left and right edges are the number of pixels to trim from each.

        Call the program as follows:
        
        python edge_replace.py <input file> <output file> <input columns> <left edge> <right edge> <data type>        
        '''
    )
if __name__ == '__main__':
    
    try:
        args = sys.argv[1:]
        process(*args)
    except Exception as e:
        print(traceback.print_tb(e))
        usage()
    