#!/usr/bin/env python
#import pdb; pdb.set_trace()
import os, sys
import argparse
#import re

# Only needed for validation of new first line
import Bio
from Bio.GenBank.Scanner import GenBankScanner
from Bio.GenBank import _FeatureConsumer
from Bio.GenBank.utils import FeatureValueCleaner

def parse_cmd(args):
    parser = argparse.ArgumentParser()
    parser.add_argument('--in-dir',
                        dest='inp_dir',
                        help='Path to directory in which the input GenBank ' + \
                        'files exists')
    parser.add_argument('--out-dir',
                        dest='out_dir',
                        help='Path to existing directory to which the altered ' + \
                        'Genbank files are written')
    parser.add_argument('--err-file-path',
                        dest='err_file_path',
                        default='err.out',
                        help='Path to file in which CDS deviants are listed')
    parser.add_argument('--validate-biopython',
                        dest='test_str_bio',
                        default=False,
                        action='store_true',
                        help='Flag to validate against BioPython routines if altered first string passes parser')

    args_data = parser.parse_args()
    print (args_data)

    return args_data.inp_dir, args_data.out_dir, args_data.err_file_path, args_data.test_str_bio


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

def reformat_genbank_first_line(first_line_inp, test_line=True):

    str_out = ''
    error_name = None

    line_parts = first_line_inp.split()
    if not 'LOCUS' in line_parts[0]:
        error_name = 'Massive error: LOCUS not on first line'
        return str_out, error_name

    try:
        bp_index = line_parts.index('bp')
    except ValueError:
        error_name = 'Missing bp'
        return str_out, error_name

    name_slice = line_parts[1:bp_index - 1]
    name = '_'.join(name_slice)

    padding_data = ' '.join([x.lower() for x in line_parts[bp_index -1:]])
    total_len = 12 + len(name) + len(padding_data)
    if total_len >= 80:
        extra_space = ' '
    else:
        extra_space = ' ' * (80 - total_len)

    str_out = 'LOCUS' + '       ' + name + extra_space + padding_data + '\n'
    print (str_out.split())

    if test_line:
        consumer = _FeatureConsumer(use_fuzziness=1, feature_cleaner=FeatureValueCleaner())
        try:
            GenBankScanner(debug=1)._feed_first_line(consumer, str_out)
        except Exception as err:
            error_name = err

    return str_out, error_name

def main(args):
    # Parse command-line
    inp_folder, out_folder, error_file_path, test_str_bio = parse_cmd(args)

    # Determine genbank files available as input
    genbank_files = get_gb_list(inp_folder)

    # Get error file ready
    f_err = open(error_file_path, 'w')

    for i in genbank_files:
        with open(i, 'r', encoding='cp437') as f:
            lines = f.readlines()
        
#        locus_match = re.search(r"LOCUS\s+((\S+\s+\S+)+)\s+\d{3,6}\sbp", lines[0])
#        if locus_match:
#            locus_new = re.sub(r"\s+","_",locus_match.group(1))
#            lines[0] = lines[0].replace(locus_match.group(1),locus_new)
        new_first_line, error_name = reformat_genbank_first_line(lines[0], test_str_bio)

        if error_name is None:
            lines[0] = new_first_line
        else:
            print('{}, {}'.format(i, error_name), file=f_err)

        fout = i.replace(inp_folder, out_folder)
        with open(fout, "w", encoding='utf-8') as f:
            f.writelines(lines)

    f_err.close()

if __name__ == '__main__':
    sys.argv = ["genbank_update_locus.py", "--in-dir=AllconstructGBfiles_old",
           "--out-dir=format_cleaned", "--validate-biopython"]
    sys.exit(main(sys.argv[1:]))