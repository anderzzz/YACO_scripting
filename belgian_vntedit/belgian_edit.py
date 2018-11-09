'''Script to create Genbank copies with VNTNAME field instructive protein name 

Example usage:

python3 belgian_edit.py --spec-file=pairs.csv --inp-source-dir=input_dir/
                        --out-dir=output_dir/ 

'''
import sys
import argparse
import os
import itertools

import pandas as pd

MATCH_DEFAULT = ['HC', 'heavy chain', '\"heavy chain\"', 'heavy\\chain',
                 'LC', 'light chain', '\"light chain\"', 'light\\chain']

def parse_cmd(args):
    '''Parse command line.

    '''
    default_blanks = ','.join(MATCH_DEFAULT)
    parser = argparse.ArgumentParser()
    parser.add_argument('--spec-file',
                        dest='inp_file',
                        help='Path to CSV specification file that provides ' + \
                        'matching between construct name and protein name')
    parser.add_argument('--blank-hc-lc',
                        dest='blank_hc_lc',
                        action='store_true',
                        default=False,
                        help='If provided, switch to note for the HC and LC ' + \
                        'labels in output')
    parser.add_argument('--kill-hc-lc',
                        dest='kill_hc_lc',
                        action='store_true',
                        default=False,
                        help='If provided, remove entirely the HC and LC ' + \
                        'labels in output')
    parser.add_argument('--blank-hc-lc-names',
                        dest='blank_names',
                        default=default_blanks,
                        help='Optional comma-separated string of HC and LC ' + \
                        'labels (case-insensitive) that are to be blanked out')
    parser.add_argument('--inp-source-dir',
                        dest='inp_dir',
                        help='Path to directory in which the input GenBank ' + \
                        'files exists')
    parser.add_argument('--out-dir',
                        dest='out_dir',
                        help='Path to existing directory to which the altered ' + \
                        'Genbank files are written')

    args_data = parser.parse_args()

    blank_names_list = args_data.blank_names.split(',')

    return args_data.inp_file, args_data.inp_dir, args_data.out_dir, \
           args_data.blank_hc_lc, args_data.kill_hc_lc, blank_names_list

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

def main(args):
    '''Main function that is called first

    '''
    # Parse command-line
    inp_file, inp_folder, out_folder, blank_hc_lc, kill_hc_lc, blank_list = parse_cmd(args)

    # Determine genbank files available as input
    genbank_files = get_gb_list(inp_folder)

    # Read mapping data
    mapping_data = pd.read_csv(inp_file)

    for index, row in mapping_data.iterrows():

        # Read construct and protein name
        file_name = row.loc['construct name']
        protein_name = row.loc['protein name']

        print ('Processing file: %s' %(file_name))

        # Get file path to input source Genbank file
        file_path = [path for path in genbank_files if file_name in path].pop()

        with open(file_path) as fin:
            all_data = fin.read().split('\n')

        # Determine if file contains VNTNAME field
        vntname_mask = ['VNTNAME|' in x for x in all_data]
        vntname_present = any(vntname_mask)

        # If VNTNAME field present
        if vntname_present:
            vntname_strings = list(itertools.compress(all_data, vntname_mask))
            if len(vntname_strings) != 1:
                raise RuntimeError('File %s contains more than one VNTNAME' %(file_name))

            # Obtain row in file with VNTNAME and modify
            vntname_index = [n for n, mask in enumerate(vntname_mask) if mask].pop()
            vntname = vntname_strings.pop()
            vntname_pre = vntname.split('|')[0]
            vntname_new = vntname_pre + '|' + protein_name + '|'

            # Remove old VNTNAME field from input data
            all_data.pop(vntname_index)

            print ('... VNTNAME adjustment done.')

        # If VNTNAME field absent
        else:
            # Make a rough guess where among the comments to fit VNTNAME field
            comment_rows = [n for n, row in enumerate(all_data) if 'COMMENT' in row]
            vntname_index = comment_rows[min(6, len(comment_rows))]
            vntname_new = 'COMMENT     VNTNAME|' + protein_name + '|'

            print ('... VNTNAME addition done.')

        # Insert the new VNTNAME field in input data
        all_data.insert(vntname_index, vntname_new)

        # Optionally remove HC LC labels
        if blank_hc_lc or kill_hc_lc:
            count = 0
            all_data_subst = []
            for line in all_data:
                matcher = any([(blank_name.lower() in line.lower()) for blank_name in blank_list])
                if matcher and 'label' in line:
                    count += 1
                    if blank_hc_lc:
                        line_new = line.replace('label', 'note')
                        all_data_subst.append(line_new)
                    elif kill_hc_lc:
                        pass

                else:
                    all_data_subst.append(line)

            all_data = all_data_subst
            print ('... Number of HC-LC label substitutions: %s' %(str(count)))

        # Print data to new file
        total_rows = len(all_data)
        with open(out_folder + '/' + file_name + '.gb', 'w') as fout:
            for n_row, row_data in enumerate(all_data):
                if n_row < total_rows - 1:
                    fout.write(row_data + '\n')
                else:
                    fout.write(row_data)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
