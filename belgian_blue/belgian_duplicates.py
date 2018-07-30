'''Script to determine duplicate sequences from GenBank file.

Example usage:

python3 belgian_blue.py --input-file-folder=stuff_gbs/ 

python3 belgian_blue.py --input-file-folder=stuff_gbs/ --max-files-batch=1000 1>out 2>err

'''
import sys
import argparse
import os

import Bio
from Bio import SeqIO

import pandas as pd

def parse_cmd(args):
    '''Parse command line.

    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-csv',
                        dest='inp_file',
                        help='Path to CSV file with extracted data from GenBank files')
    parser.add_argument('--max-files-batch',
                        dest='n_files_batch',
                        default='100',
                        help='Maximum files to include in batch to read into ' + \
                             'memory; if memory problems, reduce')

    args_data = parser.parse_args()

    if args_data.inp_file is None:
        raise RuntimeError('Missing mandatory input file folder')

    return args_data.inp_file, int(args_data.n_files_batch)

def square_test(df1, df2):
    '''For two non-overlapping chunks of file paths, create the unique pairs

    Parameters
    ----------

    Returns
    -------

    '''
    data_pairs = []

    for index1, row1 in df1.iterrows():
        for index2, row2 in df2.iterrows():

            equal = row1.sequence == row2.sequence

            if equal:
                data_pairs.append([row1.path, row2.path, 
                                   row1.vntname, row2.vntname,
                                   row1.sequence, row2.sequence])

    return data_pairs

def triangle_test(df1, df2):
    '''For a given chunk of file paths, create the unique pairs

    Parameters
    ----------

    Returns
    -------

    '''
    data_pairs = []

    for index1, row1 in df1.iterrows():
        for index2, row2 in df2.iterrows():
            
            if index1 >= index2:
                equal = False

            else:
                equal = row1.sequence == row2.sequence

            if equal:
                data_pairs.append([row1.path, row2.path, 
                                   row1.vntname, row2.vntname,
                                   row1.sequence, row2.sequence])

    return data_pairs

def main(args):

    inp_file, n_batch = parse_cmd(args)

    print ('path1;path2;vntname1;vntname2;sequence1;sequence2', file=sys.stdout)

    df_reader = pd.read_csv(inp_file, sep=';', chunksize=n_batch / 2)

    for k_chunk1, df_1 in enumerate(df_reader):
        df_reader_inner = pd.read_csv(inp_file, sep=';', chunksize=n_batch / 2)

        for k_chunk2, df_2 in enumerate(df_reader_inner):

            print ('Chunk pairs: %s, %s' %(str(k_chunk1), str(k_chunk2)), file=sys.stderr)

            if k_chunk1 == k_chunk2:
                out_data = triangle_test(df_1, df_2)

            elif k_chunk1 < k_chunk2:
                out_data = square_test(df_1, df_2)

            elif k_chunk1 > k_chunk2:
                out_data = []

            for entry in out_data:
                print (';'.join(entry), file=sys.stdout)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
