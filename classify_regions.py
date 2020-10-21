import re
import os
import logging
import glob
import sys
import argparse
import json

raw_f_regex = re.compile(
    "([A-z0-9.-]+)\s+-\s+(\w+)\s+(\w+)\s+(\w+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([-+])\s+([-+])\s+(\d+)\s+(\d+[\.\d]*)\s+(\d+[\.\d]*)\s+(\d+[\.\d]*)\s+(.+)\s!\s+ -")

MIN_OVERLAP = 0.95

MIN_MATCH_PROPORTION = 0.01

regions_16S = {
    'V1': [69, 92],
    'V2': [131, 239],
    'V3': [430, 487],
    'V4': [566, 672],
    'V5': [812, 869],
    'V6': [976, 1033],
    'V7': [1107, 1164],
    'V8': [1234, 1285],
    'V9': [1426, 1456]
}

regions_18S = {
    'V1': [69, 105],
    'V2': [135, 190],
    'V3': [474, 545],
    'V4': [627, 727],
    'V5': [1059, 1102],
    'V7': [1366, 1382],
    'V8': [1528, 1608],
    'V9': [1727, 1750]
}

logging.basicConfig(level=logging.INFO)


def calc_overlap(read, reg):
    read_s, read_f = read
    reg_s, reg_f = reg
    overlap = max(min(read_f, reg_f) - max(read_s, reg_s), 0)
    total = max(reg_f - reg_s, 0)
    try:
        return overlap / total
    except ZeroDivisionError:
        return 0


def get_regions(raw_sequence, regions):
    return {region for region, limits in regions.items() if calc_overlap(raw_sequence, limits) >= MIN_OVERLAP}


# Parse, filter empty lines and unpack into 2D array
def load_data(filename):
    with open(filename) as f:
        return [l[0] for l in [raw_f_regex.findall(l) for l in f] if bool(l)]


def normalise_results(region_matches, regions):
    count_matches = len(region_matches)
    results = {}
    for k in regions:
        matches = len([m for m in region_matches if m == k])
        try:
            ratio = matches / count_matches
        except ZeroDivisionError:
            ratio = 0
        logging.debug(k, matches, ratio)

        if ratio > 0:
            results[k] = {'match_proportion': round(ratio, 4), }
    return filter_minimum_match_proportion(results)


def filter_minimum_match_proportion(region_matches):
    return {region: v for region, v in region_matches.items() if v['match_proportion'] >= MIN_MATCH_PROPORTION}


def identify_run(infile_name):
    """
    Args:
        infile_name: The name of the tblout file.
    Return:
        run: Run ID ERR*|SRR*
    """
    run = 'TEMP_NAME'
    return run


def retrieve_regions(tblout_file, outfile, subunit_type):
    regions = regions_16S if subunit_type == '16S' else regions_18S

    data = load_data(tblout_file)
    region_matches = []
    for read in data:
        limits = list(map(int, read[4:6]))
        region_matches.extend(get_regions(limits, regions))
    normalised_matches = normalise_results(region_matches, regions)
    print(tblout_file, normalised_matches)
    with open(outfile, 'w') as f:
        json.dump(normalised_matches, f)
    run_id = identify_run(tblout_file)
    print_to_table(outfile, normalised_matches, run_id)


def print_to_table(json_out, variable_regions, run_id):
    """Prints the variable regions to a tsv file.
    Args:
        json_out: The name of the json outfile.
        variable_regions: The dictionary that contains a list of variable regions for a run and match proportions.
        run_id: Run ID (ERR*|SRR*)
    """
    # generate a file name for the tsv file
    if json_out.endswith('json'):
        tsv_out = re.sub('json$', 'tsv', json_out)
    else:
        tsv_out = '{}.tsv'.format(json_out)
    # determine the variable region to output
    if len(variable_regions.keys()) > 1:
        amplified_region = '{}-{}'.format(min(variable_regions.keys()), max(variable_regions.keys()))
    else:
        amplified_region = next(iter(variable_regions))
    header = 'Run\tAssertionEvidence\tAssertionMethod\tVariable region'
    record = '{}\tECO_0000363\tautomatic assertion\t{}\n'.format(run_id, amplified_region)
    if not (run_id and amplified_region):
        sys.exit('ERROR: Values missing in line: {}'.format(record))
    with open(tsv_out, 'w') as f:
        f.write('\n'.join([header, record]))


def parse_args(argv):
    parser = argparse.ArgumentParser(description='Tool to determine which regions were amplified in 16S data')
    parser.add_argument('file', help='Overlapped tblout file')
    parser.add_argument('subunit_type', choices=['16S', '18S'], help='rRNA subunit type (16S vs 18S)')
    parser.add_argument('-o', '--output_file', default='amplified_regions.json')
    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    retrieve_regions(args.file, args.output_file, args.subunit_type)


if __name__ == '__main__':
    main(sys.argv[1:])
