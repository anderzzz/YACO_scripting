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

def parse_cmd(args):
    '''Parse command line.

    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-file-folder',
                        dest='inp_files',
                        help='Path to folder with all input files')
    parser.add_argument('--max-files-batch',
                        dest='n_files_batch',
                        default='100',
                        help='Maximum files to include in batch to read into ' + \
                             'memory; if memory problems, reduce')

    args_data = parser.parse_args()

    if args_data.inp_files is None:
        raise RuntimeError('Missing mandatory input file folder')

    return args_data.inp_files, int(args_data.n_files_batch)

def get_gb_list(path, suffixes=['gb', 'gbk']):
    '''Determine the list of GenBank files in the give folder

    Parameters
    ----------
    path : str
        Path to folder that contains all GenBank files, but not necessarily
        only Genbank files
    suffixes : list, optional
        List of string suffixes of the GenBank files.

    Returns
    -------
    file_list : list
        List of path and file name to the discovered GenBank files

    '''
    ret_data = []

    in_directory = os.listdir(path)
    for filename in in_directory:
        for suffix in suffixes:
            len_suffix = len(suffix)

            if filename[-1 * (len_suffix + 1):] == '.' + suffix:
                ret_data.append(path + '/' + filename)

    return ret_data

def chunk_(files, n_batch):
    '''Chunk list of files into batches

    Parameters
    ----------
    files : list
        List of string file paths
    n_batch : int
        Maximum batch size to keep in memory simultaneously

    Returns
    -------
    chunks : list
        List of list of file paths

    '''
    ret_chunks = []

    if len(files) <= n_batch:
        ret_chunks.append(files)

    else:
        tmp = []
        for fff in files:
            tmp.append(fff)

            if len(tmp) >= float(n_batch) / 2.0:
                ret_chunks.append(tmp)
                tmp = []

        if len(tmp) > 0:
            ret_chunks.append(tmp)

    return ret_chunks

def extract_vntdata(record):
    '''Extract the salient vntdata from a potentially lengthy comment

    Parameters
    ----------
    record : SeqIO
        Record of data as obtained from SeqIO

    Returns
    -------
    the_name : str
        The VNTNAME field if it exists. If it is missing a warning string is
        returned instead

    '''
    if not 'comment' in record.annotations:
        the_name = 'WARNING:NO VNTNAME FOUND'

    else:
        comment = record.annotations['comment']
        comment_lines = comment.split('\n')
        for line in comment_lines:
            if 'vntname' in line.lower():
                the_name = line.split('|')[1]
                break

        else:
            the_name = 'WARNING:NO VNTNAME FOUND'

    return the_name

def genbank_reader(handle):
    '''Read the sequence record from a given GenBank file

    Parameters
    ----------
    handle : file handle
        File handle to a GenBank file

    Returns
    -------
    record : SeqRecord
        The sequence record obtained from the file

    Raises
    ------
    NotImplementedError
        If more than one sequence record is found in the GenBank file

    '''
    parser_generator = SeqIO.parse(handle, "genbank")

    record = next(parser_generator)

    try:
        dummy = next(parser_generator)
        raise NotImplementedError('File %s contains more than one ' %(handle.name) + \
                                  'sequence record')

    except StopIteration:
        pass

    return record

def cleanse_comments(handle):
    '''Because the raw input file can have lots of irrelevant comments that
    slow down parsing, this function crudely removes these irrelevant comments

    Parameters
    ----------
    handle : file handle
        File handle to a GenBank file

    Returns
    -------
    tmp_filename : str
        Path to the temporary cleansed version of the GenBank file

    '''
    TMP_FILENAME = 'tmp.gb'

    f_tmp = open(TMP_FILENAME, 'w')

    data = handle.read().split('\n')
    for x in data:
        if 'COMMENT' in x:
            if 'VNT' in x:
                print (x, file=f_tmp)

        else:
            print (x, file=f_tmp)

    f_tmp.close()

    return TMP_FILENAME

def read_data(handle):
    '''Wrapper function to read salient data from the Genbank file

    Parameters
    ----------
    handle : file handle
        File handle to a GenBank file

    Returns
    -------
    salient_data : dict
        Dictionary of data extracted from GenBank file, keys: `sequence` for
        the SeqIO record, `vntname` for the VNTNAME in the file, `path` for the
        file path.

    '''
    tmp_filename = cleanse_comments(handle)

    with open(tmp_filename) as handle_tmp:
        record = genbank_reader(handle_tmp)

    vnt_name = extract_vntdata(record)
    seq_data = record.seq
    
    return {'sequence' : seq_data, 'vntname' : vnt_name, 'path' : handle.name}

def _square(chunk1, chunk2):
    '''For two non-overlapping chunks of file paths, create the unique pairs

    Parameters
    ----------
    chunk1 : list
        List of file paths
    chunk2 : list
        List of file paths

    Returns
    -------
    data_pairs : list
        List of tuples of SeqIO records for all unique pairs of files in input

    '''
    data_pairs = []

    for f1 in chunk1:
        for f2 in chunk2:
            
            with open(f1) as handle_1:
                data1 = read_data(handle_1)

            with open(f2) as handle_2:
                data2 = read_data(handle_2)

            data_pairs.append((data1, data2))

    return data_pairs

def _triangle(chunk):
    '''For a given chunk of file paths, create the unique pairs

    Parameters
    ----------
    chunk : list
        List of file paths

    Returns
    -------
    data_pairs : list
        List of tuples of SeqIO records for all unique pairs of files in input

    '''
    data_pairs = []

    for k, f1 in enumerate(chunk):
        for f2 in chunk[k+1:]:

            with open(f1) as handle_1:
                data1 = read_data(handle_1)

            with open(f2) as handle_2:
                data2 = read_data(handle_2)

            data_pairs.append((data1, data2))

    return data_pairs

def lower_triangular(chunks):
    '''Iterator for lower triangular loop of file chunks. The details of how
    the files are pairwise handled is automatically handled and the iterator
    returns all unique pairs of sequence records

    Parameters
    ----------
    chunks : list
        List of list of file paths to process

    Returns
    -------
    d1 : SeqRecord
        SeqRecord for first file
    d2 : SeqRecord
        SeqRecord for second file

    '''
    if len(chunks) == 0:
        raise RuntimeError('No files in the list to consider')

    for k_chunk1, chunk_1 in enumerate(chunks):
        for k_chunk2, chunk_2 in enumerate(chunks):

            if k_chunk1 == k_chunk2:
                data_pairs = _triangle(chunk_1)
                
            elif k_chunk1 < k_chunk2:
                data_pairs = _square(chunk_1, chunk_2)

            elif k_chunk1 > k_chunk2:
                data_pairs = []

            else:
                raise RuntimeError('Strange index relation encountered')

            if len(data_pairs) > 0:
                print ('File Chunk Pair Indeces ' + \
                       '%s %s' %(str(k_chunk1), str(k_chunk2)), file=sys.stderr)

            for d1, d2 in data_pairs:
                yield d1, d2 
        

def main(args):

    inp_file_folder, n_batch = parse_cmd(args)
    genbank_files = get_gb_list(inp_file_folder)
    chunks_of_files = chunk_(genbank_files, n_batch)
    
    n_chunk_pairs = len(chunks_of_files) * (len(chunks_of_files) + 1) / 2
    print ('Total File Chunk Pairs to Process: %s' %(str(n_chunk_pairs)), file=sys.stderr)

    print ('vntname_1,vntname_2,path_1,path_2', file=sys.stdout)
    n_counter = 0
    for data1, data2 in lower_triangular(chunks_of_files):

        if data1['sequence'] == data2['sequence']: 
            n_counter += 1

            out_str = [data1['vntname'], data2['vntname'], 
                       data1['path'], data2['path']]
            print (','.join(out_str), file=sys.stdout)

    print ('Found %s data pair(s) with identical sequences' %(str(n_counter)), file=sys.stderr)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
