# -*- coding: utf-8 -*-
"""
TPP and VEC amino acid matching.

GD Coding Support by Anders @ YACO Scientific.

How to run
----------


"""
import csv
import xml.etree.ElementTree as ET
import pandas as pd
from collections import namedtuple
import urllib.request
import urllib.parse

AAComparison = namedtuple('AAComparison', 'identical, subset, superset')

'''Path constants'''
PATH_MASTER = '/Users/andersohrn/Documents/ao_dev/GD_coding/vec_tpp_1908/PEB-VEC_test1.txt'
PATH_OUT = '/Users/andersohrn/Documents/ao_dev/GD_coding/vec_tpp_1908/tmp.csv'
PATH_FILES = '/Users/andersohrn/Documents/ao_dev/GD_coding/vec_tpp_1908/'

'''Web constants'''
URL_ROOT = 'http://puck.bos.us.genedata.com:8000/'
URL_CONST = 'Biologics/ws/rest/sequence/aa/qid/'
USER = 'admin1'
PASSWD = ''

'''Result constants'''
IND_KEYS = ['peb_id', 'tpp_id', 'vec_id', 'vec_sub_id']
VALUE_RESULT_NAME = 'Comparison Value'

'''Master input file constants'''
INP_KEYS = ['ID', 'Vector list', 'Target Product Protein ID']

def read_xml(func):
    '''Wrapper function to parse XML data followed by custom selection'''
    
    def wrapper(path):
    
        aa_seqs = []
        tree = ET.parse(path)
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
        ret_val = {VALUE_RESULT_NAME : 'PERFECT'}
        
    else:
        if result.subset:
            ret_val = {VALUE_RESULT_NAME : 'SUBSET'}
            
        elif result.superset:
            ret_val = {VALUE_RESULT_NAME : 'SUPERSET'}
            
        else:
            ret_val = {VALUE_RESULT_NAME : 'POOR'}
            
    return ret_val

def retrieve_web(file_ids):
    '''Get resource from web api'''
    
    file_ids_str = ','.join(file_ids)
    
    url = URL_ROOT + URL_CONST + file_ids_str
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req) as response:
        the_page = response.read()
        
    print (the_page)
    
def init_web():
    '''Bla bla'''
    
    password_mgr = urllib.request.HTTPPasswordMgrWithDefaultRealm()
    password_mgr.add_password(None, URL_ROOT, USER, PASSWD)
    handler = urllib.request.HTTPBasicAuthHandler(password_mgr)
    opener = urllib.request.build_opener(handler)
    urllib.request.install_opener(opener)

def get_xmls_file(cls_str, vec_or_tpp='vec'):
    '''Obtain corresponding vec or tpp files from disk'''
    
    ret_vals = []
    for vec_name in cls_str.split(','):
        path = PATH_FILES + '/' + vec_name.lower().strip() + '.xml'
        
        if vec_or_tpp == 'vec':
            aas_ = read_vec(path)
        elif vec_or_tpp == 'tpp':
            aas_ = read_tpp(path)
        else:
            raise RuntimeError('Unknown XML file type: {}'.format(vec_or_tpp))
            
        ret_vals.append(aas_)
        
    return ret_vals

def get_vec_xmls(csl_str):
    '''Obtain corresponding vec files'''
    
    return get_xmls_file(csl_str, 'vec')

def get_tpp_xml(csl_str):
    '''Obtain corresponding tpp files'''
    
    return get_xmls_file(csl_str, 'tpp')

def parse_master_iter(path_inp, sep=';'):
    '''Iterator for master file parsing'''
    
    with open(path_inp) as fin:
        reader = csv.DictReader(fin, delimiter=sep, quotechar='"')
        for row in reader:
            yield dict(row)

def analyze_main(path_inp, sep_inp):
    '''Main analysis routine and output data creation'''
    
    df_result = pd.DataFrame()
    
    # Loop over rows in master input file
    for peb_row in parse_master_iter(path_inp, sep_inp):
    
        # Extract all data from given files
        peb_id = peb_row[INP_KEYS[0]]
        vec_xmls = get_vec_xmls(peb_row[INP_KEYS[1]])
        tpp_xml = get_tpp_xml(peb_row[INP_KEYS[2]])
    
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

def main():
    
    df = analyze_main(PATH_MASTER, '\t')
    df.to_csv(PATH_OUT, ';')

if __name__ == '__main__':
    main()