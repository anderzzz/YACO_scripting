'''Script to extract most relevant data from GenBank files 

Example usage:

python3 belgian_extract.py --input-file-folder=stuff_gbs/ --output-file=stuff.csv

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
    parser.add_argument('--output-file',
                        dest='out_file',
                        default='out.csv',
                        help='Output file prefix path (default: %default)')

    args_data = parser.parse_args()

    if args_data.inp_files is None:
        raise RuntimeError('Missing mandatory input file folder')

    return args_data.inp_files, args_data.out_file

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

def extract_features(features):
    '''Extract sequence features of interest

    Parameters
    ----------
    features : list
        List of SeqFeature objects

    Returns
    -------

    '''
    c_region = []
    for fet in features:
        if fet.type == 'C_region':
            c_region.append(fet.qualifiers['label'][0])

    return c_region

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
    error = ''

    with open(tmp_filename) as handle_tmp:
        try:
            record = genbank_reader(handle_tmp)
            vnt_name = extract_vntdata(record)
            seq_data = record.seq

            features = extract_features(record.features)

        except ValueError as err:
            vnt_name = 'MALFORMED_GENBANK_FILE'
            seq_data = 'MALFORMED_GENBANK_FILE'
            features = ['MALFORMED_GENBANK_FILE']
            error = str(err)
    
    return {'sequence' : seq_data, 'vntname' : vnt_name, 'path' : handle.name,
            'c_region' : '|'.join(features), 'error' : error}

def main(args):

    inp_file_folder, out_file = parse_cmd(args)
    genbank_files = get_gb_list(inp_file_folder)

    fout = open(out_file, 'w')
    print ('path;vntname;sequence;c_region_features;error_message', file=fout)

    for filepath in genbank_files:
        with open(filepath) as fin:
            file_data = read_data(fin)

        out_data = [file_data['path'], file_data['vntname'],
                    str(file_data['sequence']), file_data['c_region'],
                    file_data['error']]
        print (';'.join(out_data), file=fout)

    fout.close()
    
if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
