"""
Target Protein Product (TPP) and Vector (VEC) amino acid sequence matching.

Example run:
    
    Data files on disk
    ------------------
        
    python3 aa_check.py --compare-inp=PEB-VEC_test1.txt --compare-inp-sep=$'\t' \
                        --inp-file-path='.' --result-exclude-identical
                        
    This command parses a master input file with columns separated by tab, the
    data files are XML format on disk, and only cases that are not identical
    are included in the output
    
    Data files on web-server
    ------------------------

    python3 aa_check.py --compare-inp=PEB-VEC_test1.txt --compare-inp-sep=$'\t' \
                        --url-root='http://awesome.com:8080' --web-user='foobar'
                        
    This command parses a master input file with columns separated by tab, the
    data files are XML on the Biologics web-server at the root awesome.com and
    the user for said server is foobar. Note that after this command is run, 
    a password is request from the command-line.
    

GD Coding Support by Anders @ YACO Scientific.
August 2019

"""
import csv
import xml.etree.ElementTree as ET

import argparse
import getpass

import pandas as pd
from collections import namedtuple

import urllib.request
import urllib.parse
from pathlib import Path

AAComparison = namedtuple('AAComparison', 'identical, subset, superset')
'''Amino Acid comparison result tuple'''

'''Web defaults and constants'''
URL_ROOT = 'http://puck.bos.us.genedata.com:8080/'
USER = 'admin1'
URL_CONST = 'Biologics/ws/rest/sequence/aa/qid/'

'''IO defaults'''
PATH_OUT = 'result.csv'
SEP_DEFAULT = ';'

'''Output formatting constants'''
IND_KEYS = ['peb_id', 'tpp_id', 'vec_id', 'vec_sub_id']
VALUE_RESULT_NAME = 'Comparison Value'
VALUE_RESULTS = ['Identical', 'Subset', 'Superset', 'NoMatch']

'''Input formatting constants'''
INP_KEYS = ['ID', 'Vector list', 'Target Product Protein ID']

def read_xml(func):
    '''Wrapper function to parse XML data followed by custom selection'''
    
    def wrapper(path, file_not_string=True):
    
        aa_seqs = []
        if file_not_string:
            tree = ET.parse(path)
        else:
            tree = ET.fromstring(path)
            
        for element in tree.iter():
            x = func(element)
            if not x is None:
                aa_seqs.append(x)
            
        return aa_seqs
    
    return wrapper

@read_xml
def read_vec(element):
    '''XML data selection function for VEC file'''
    
    if 'residues' in element.tag:
        return element.text
    
@read_xml
def read_tpp(element):
    '''XML data selection function for TPP file'''
    
    if 'residues' in element.tag:
        return element.text
    
@read_xml
def read_id(element):
    '''XML data id extraction function'''
    
    if 'qualifiedId' == element.tag[-11:]:
        return element.text
    
def comparison(aa_vec, aa_tpp):
    '''Testing amino acid similarities'''
    
    ret = []
    for k_vec, aa_v in enumerate(aa_vec):
        test_identical = False
        test_subset = False
        test_superset = False
        
        for aa_t in aa_tpp:
            
            if aa_v == aa_t:
                test_identical = True
                
            elif (aa_v in aa_t) and aa_v != aa_t:
                test_subset = True
                
            elif aa_t in aa_v:
                test_superset = True
                
        vec_comp = AAComparison(test_identical, test_subset, test_superset)
        ret.append(vec_comp)
        
    return ret

def output_format(result):
    '''Convert result to string output for CSV writing'''

    if result.identical:
        ret_val = {VALUE_RESULT_NAME : VALUE_RESULTS[0]}
        
    else:
        if result.subset:
            ret_val = {VALUE_RESULT_NAME : VALUE_RESULTS[1]}
            
        elif result.superset:
            ret_val = {VALUE_RESULT_NAME : VALUE_RESULTS[2]}
            
        else:
            ret_val = {VALUE_RESULT_NAME : VALUE_RESULTS[3]}
            
    return ret_val

def init_web(root, user, passwd):
    '''Do the one time initializations to access web-server'''
    
    password_mgr = urllib.request.HTTPPasswordMgrWithDefaultRealm()
    password_mgr.add_password(None, root, user, passwd)
    handler = urllib.request.HTTPBasicAuthHandler(password_mgr)
    opener = urllib.request.build_opener(handler)
    urllib.request.install_opener(opener)

def retrieve_web_(id_name, file_type, url_root):
    '''Inner routine to access one given XML on web server'''
    
    url = url_root + '/' + URL_CONST + id_name
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req, timeout=5) as response:
        the_page = response.read()
        
    if file_type == 'vec':
        aas_ = read_vec(the_page, False)
    elif file_type == 'tpp':
        aas_ = read_tpp(the_page, False)
    else:
        raise RuntimeError('Unknown XML file type: {}'.format(file_type))
        
    return aas_

def retrieve_disk_(id_name, file_type, inp_file_path):
    '''Inner routine to access one given XML on disk'''
    
    path = inp_file_path + '/' + id_name + '.xml'
    my_file = Path(path)
    if not my_file.exists():
        raise RuntimeError('Cannot find file {}'.format(path))
        
    if file_type == 'vec':
        aas_ = read_vec(path, True)
    elif file_type == 'tpp':
        aas_ = read_tpp(path, True)
    else:
        raise RuntimeError('Unknown XML file type: {}'.format(file_type))
    
    return aas_

def get_xml_data(file_type, id_list, from_web, url_root, inp_file_path):
    '''Master routine to retrieve data from sequence files'''
    
    ret_vals = []
    for id_name in id_list.split(','):
        id_name_clean = id_name.lower().replace(' ','')
        
        if from_web:
            aa_seq_ = retrieve_web_(id_name_clean, file_type, url_root)
            
        else:
            aa_seq_ = retrieve_disk_(id_name_clean, file_type, inp_file_path)
            
        ret_vals.append(aa_seq_)
        
    return ret_vals
    
def parse_master_iter(path_inp, sep=';', quote='"'):
    '''Iterator for master file parsing'''
    
    with open(path_inp) as fin:
        reader = csv.DictReader(fin, delimiter=sep, quotechar=quote)
        for row in reader:
            yield dict(row)

def analyze_main(web_access, url_root, web_user, web_passwd, 
                 inp_file_path,
                 compare_inp, compare_inp_sep):
    '''Main analysis routine and output data creation'''
    
    df_result = pd.DataFrame()
    
    if web_access:
        init_web(url_root, web_user, web_passwd)
    
    # Loop over rows in master input file
    for peb_row in parse_master_iter(compare_inp, compare_inp_sep):
    
        # Extract all data from given files
        peb_id = peb_row[INP_KEYS[0]]
        vec_xmls = get_xml_data('vec', peb_row[INP_KEYS[1]], 
                                web_access, url_root, inp_file_path)
        tpp_xml = get_xml_data('tpp', peb_row[INP_KEYS[2]], 
                               web_access, url_root, inp_file_path )
    
        vec_ids = [x.strip() for x in peb_row[INP_KEYS[1]].split(',')]  
        tpp_id = peb_row[INP_KEYS[2]]
        
        # Loop over all vector files sequences (outer) and over all tpp 
        # sequences (inner). Most times expect only one sequence per file,
        # but logic can handle multiple sequences per file
        for vec_id, vec_seq in zip(vec_ids, vec_xmls):
            for tpp_seq in tpp_xml:
                ret = comparison(vec_seq, tpp_seq)
            
                # In case multiple sequences in given vec file, disaggregate
                # the comparison results
                for kk, result_row in enumerate(ret):
                    result_formatted = output_format(result_row)
            
                    ind = pd.MultiIndex.from_tuples([(peb_id, tpp_id, vec_id, kk)], 
                                                    names=IND_KEYS)
                    df_result = df_result.append(pd.DataFrame(result_formatted, index=ind))
                
    return df_result

def parse_args():
    '''Parse command-line arguments'''
    
    parser = argparse.ArgumentParser(description='Validate vector sequence files ' + \
                                     'against target protein product files.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    group_web = parser.add_argument_group('Web Arguments', 'Arguments related to ' + \
                                          'web access. Note that password is ' + \
                                          'retrieved separatedly.')
    group_web.add_argument('--url-root', 
                           dest='url_root', 
                           default=URL_ROOT,
                           help='Root URL for web-server to retrieve files from')
    group_web.add_argument('--web-user',
                           dest='web_user',
                           default=USER,
                           help='Web user name for web-server access')
    
    group_io = parser.add_argument_group('IO Arguments', 'Argument related to ' + \
                                         'input and output processing')
    group_io.add_argument('--inp-file-path',
                          dest='inp_file_path',
                          default=None,
                          help='Root directory for files to read and compare. ' + \
                          'If given, this overrides web access options')
    group_io.add_argument('--compare-inp',
                          dest='compare_inp',
                          required=True,
                          help='Input file specifying the vec-tpp tests to perform')
    group_io.add_argument('--compare-inp-sep',
                          dest='compare_inp_sep',
                          default=SEP_DEFAULT,
                          help='Column separator in input file')
    group_io.add_argument('--result-out',
                          dest='result_out',
                          default=PATH_OUT,
                          help='Path to output file with comparison results')
    group_io.add_argument('--result-out-sep',
                          dest='result_out_sep',
                          default=SEP_DEFAULT,
                          help='Column separator in output file')
    group_io.add_argument('--result-exclude-identical',
                          dest='exclude_identical',
                          default=False,
                          action='store_true',
                          help='Flag to exclude identical sequence matches from result file')
    
    args = parser.parse_args()
    
    if args.inp_file_path is None:
        web_access = True
        passwd = getpass.getpass()
        
    else:
        web_access = False
        passwd = None
    
    return web_access, args.url_root, args.web_user, passwd, \
           args.inp_file_path, args.compare_inp, args.compare_inp_sep, \
           args.result_out, args.result_out_sep, \
           args.exclude_identical

def main():
    
    # Parse and interpret command-line options
    web_access, url_root, web_user, web_passwd, \
    inp_file_path, compare_inp, compare_inp_sep, \
    result_out, result_out_sep, \
    exclude_identical = parse_args()
    
    # Perform analysis
    df = analyze_main(web_access, url_root, web_user, web_passwd, 
                      inp_file_path,
                      compare_inp, compare_inp_sep)
    
    # Generate output file
    if exclude_identical:
        df = df.loc[df[VALUE_RESULT_NAME] != VALUE_RESULTS[0]] 
    df.to_csv(result_out, result_out_sep)

if __name__ == '__main__':
    main()