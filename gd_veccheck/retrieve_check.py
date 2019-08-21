# -*- coding: utf-8 -*-
"""
TPP and VEC amino acid matching.

GD Coding Support by Anders @ YACO Scientific.

How to run
----------


"""
import csv
import xml.etree.ElementTree as ET
from collections import namedtuple

AAComparison = namedtuple('AAComparison', 'identical, subset, superset')

'''Path constants'''
PATH_MASTER = '/Users/andersohrn/Documents/ao_dev/GD_coding/vec_tpp_1908/example_inp.txt'
PATH_1 = '/Users/andersohrn/Documents/ao_dev/GD_coding/vec_tpp_1908/vec-398.xml'
PATH_2 = '/Users/andersohrn/Documents/ao_dev/GD_coding/vec_tpp_1908/tpp-250.xml'

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

def output_format(result, key={}):
    '''Convert result to string output for CSV writing'''

    if result.identical:
        ret_val = None
        
    else:
        if result.subset:
            ret_val = {'val' : 'SUBSET'}
            
        elif result.superset:
            ret_val = {'val' : 'SUPERSET'}
            
        else:
            ret_val = {'val' : 'POOR'}
            
    for kk, dd in key.items():
        ret_val[kk] = dd
            
    return ret_val

def outputter(results, keys):
    '''Iterator to output result dictionary'''
    
    for result, key in zip(results, keys):
        yield output_format(result, key)

def parse_master_iter(PATH_MASTER, sep=';'):
    '''Iterator for master file parsing'''
    
    with open(PATH_MASTER) as fin:
        reader = csv.DictReader(fin, delimiter=sep)
        for row in reader:
            yield dict(row)

# 
# MAIN
#
for peb_row in parse_master_iter(PATH_MASTER):
    print (peb_row)

aas_vec = read_vec(PATH_1)
ids_vec = read_id(PATH_1)
aas_tpp = read_tpp(PATH_2)
ret = comparison(aas_vec, aas_tpp)
print (ret)
keys = [{}] * len(ret)
for result_row in outputter(ret, keys):
    if result_row is None:
        pass
    else:
        print (result_row)