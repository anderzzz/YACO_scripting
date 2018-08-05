'''Bla bla

'''
import sys
import argparse
import csv

PHRASES_JOIN = ['with', '+', ',', 'w/', ':', ';', '/', ' ', '?', 'and', 'per']
PHRASES_CONC = ['uM', 'mM', \
                'ug/ml', 'mg/l', 'g/l', 'ug per ml', 'mg per l', \
                'ug/mL', 'mg/L', 'g/L', 'ug per mL', 'mg per L']
PHRASES_VOL = ['uL', 'mL', 'ul', 'ml']
PHRASES_AMOUNT = ['ug', 'mg']

def parse_(args):

    parser = argparse.ArgumentParser()
    parser.add_argument('--input-file',
                        dest='inp_file',
                        help='Path to input CSV file with messy data')
    parser.add_argument('--column-index',
                        dest='column_index',
                        default='0',
                        help='The column index in which the data to process is found')
    parser.add_argument('--header-row',
                        dest='header_row',
                        default='0',
                        help='The header row index to ignore')
    parser.add_argument('--sep',
                        dest='sep',
                        default=',',
                        help='Separator for the CSV file')

    inp_args = parser.parse_args(args)

    inp_file = inp_args.inp_file
    col_ind = int(inp_args.column_index)
    header_ind = int(inp_args.header_row)
    sep = inp_args.sep

    return inp_file, col_ind, header_ind, sep

def recursive_splitter(string, seps, sep_index=0, agg_strings=None, ignore_char=['']):

    if agg_strings is None:
        agg_strings = []

    if len(seps) == sep_index:
        agg_strings.append(string)
        return agg_strings

    else:
        sep = seps[sep_index]
        out_strings = string.split(sep)
        sep_index += 1

        for sub_string in out_strings:
            if not sub_string in ignore_char:
                recursive_splitter(sub_string, seps, sep_index, agg_strings)

    return agg_strings

def merge_(str_data, key_phrases, beforeafter, rule_key, padding):

    skip_index = -1
    str_data_merge = []
    for k, s in enumerate(str_data):

        if k == skip_index:
            continue

        if s in key_phrases:

            if beforeafter == 'before':
                other = str_data[k - 1]
                merged = other + padding + s
            elif beforeafter == 'after':
                other = str_data[k + 1]
                merged = s + padding + other
            else:
                raise ValueError('Modifier relation incorrect specified: %s' %(beforeafter))

            if rule_key(other):
                if beforeafter == 'before':
                    str_data_merge.pop(-1)
                elif beforeafter == 'after':
                    skip_index = k + 1

                str_data_merge.append(merged)

            else:
                print ('Warning: Rule Error for %s' %(str_data), file=sys.stderr)
                str_data_merge.append(s)

        else:
            str_data_merge.append(s)

    return str_data_merge

def _be_numeric(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def _be_amount(s):
    part_of = [x in s for x in PHRASES_AMOUNT]
    return any(part_of)

def string_merger(str_data):

    merged_on_unit = merge_(str_data, PHRASES_VOL, 'before', 
                            _be_amount, '/')
    print (merged_on_unit)
    merged_on_vol = merge_(merged_on_unit, PHRASES_VOL, 'before', 
                           _be_numeric, '')
    print (merged_on_vol)
    merged_on_conc = merge_(merged_on_vol, PHRASES_CONC, 'before', 
                            _be_numeric, '')
    print (merged_on_conc)
    

def main(args):

    file_in_name, col_ind, header_ind, sep = parse_(args)

    data = []
    with open(file_in_name, encoding='utf-8') as fin:
        reader = csv.reader(fin, delimiter=sep)

        for k_row, row in enumerate(reader):
            if k_row != header_ind:
                data.append(row[col_ind])

    # LINE 44
    for entry in data:
        print (entry)
        substrings = recursive_splitter(entry, PHRASES_JOIN)
        print (substrings)
        string_merger(substrings)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
