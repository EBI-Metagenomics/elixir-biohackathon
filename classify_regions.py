import re
import os
import logging
import glob
import sys
import argparse
import json
import time

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
    run = os.path.basename(infile_name).split('_')[0]
    return run


def verify_gene(gene_declared, cm_detected):
    """Check that the declared gene (16S or 18S) matches the cmsearch result.
    This function expects that these specific models (RF00177 and RF01960) are used and ensures
    that 18S is not treated as 16S and vice versa.
    Args:
        gene_declared: 16S or 18S (declared by the user)
        cm_detected: The model that matched the sequence during cmsearch.
    Return:
        check: True or false
    """
    if gene_declared == '16S' and cm_detected == 'RF00177':
        check = True
    elif gene_declared == '18S' and cm_detected == 'RF01960':
        check = True
    else:
        check = False
    return check


def print_to_table(tsv_out, results):
    """Prints the variable regions to a tsv file.
    Args:
        tsv_out: The name of the tsv outfile.
        results: The dictionary that contains a list of variable regions for a run and their match proportions.
    """
    test_table = tsv_out + '.test'
    t = open(test_table, 'w')
    f = open(tsv_out, 'w')
    # print the table header to file
    f.write('Run\tAssertionEvidence\tAssertionMethod\tVariable region\n')
    for run in results.keys():
        # determine the variable region to output
        if len(results[run].keys()) > 1:
            amplified_region = '{}-{}'.format(min(results[run].keys()), max(results[run].keys()))
        elif len(results[run].keys()) == 1:
            amplified_region = next(iter(results[run]))
        else:
            continue
        record = '{}\tECO_0000363\tautomatic assertion\t{}\n'.format(run, amplified_region)
        if not (run and amplified_region):
            sys.exit('ERROR: Values missing in line: {}'.format(record))
        f.write(record)
        t.write('\n{}\t'.format(run))
        for key in results[run].keys():
            t.write('{}\t{}\t'.format(key, results[run][key]['match_proportion']))
    f.close()
    t.close()


def retrieve_regions(tblout_file_list, outfile, subunit_type):
    regions = regions_16S if subunit_type == '16S' else regions_18S
    file_counter = 0  # count how many files were analyzed
    sequence_counter_total = 0  # count how many sequences in total were analyzed
    sequence_counter_useful = 0  # count how many sequences an output was generated for
    normalised_matches = dict()  # dictionary that will contain results for all runs
    for tblout_file in tblout_file_list:
        data = load_data(tblout_file)
        run_id = identify_run(tblout_file)
        region_matches = []
        data_type_unsupported = 0
        for read in data:
            sequence_counter_total += 1
            if not verify_gene(subunit_type, read[2]):
                print('ERROR: provided subunit type ({}) does not match cmsearch result:\n{}'.format(subunit_type,
                                                                                                        read))
                data_type_unsupported = 1
                break
            limits = list(map(int, read[4:6]))
            region_matches.extend(get_regions(limits, regions))
            sequence_counter_useful += 1
        if data_type_unsupported:
            continue
        else:
            normalised_matches[run_id] = normalise_results(region_matches, regions)
            #print(tblout_file, normalised_matches[run_id])
            file_counter += 1

    json_outfile = '{}.json'.format(outfile)
    tsv_outfile = '{}.tsv'.format(outfile)
    with open(json_outfile, 'w') as f:
        json.dump(normalised_matches, f)
    print_to_table(tsv_outfile, normalised_matches)
    print('Analyzed {} files and {} sequences. Output generated for {} sequences'.format(file_counter,
                                                                                         sequence_counter_total,
                                                                                         sequence_counter_useful))


def parse_args(argv):
    parser = argparse.ArgumentParser(description='Tool to determine which regions were amplified in 16S data')
    parser.add_argument('files', nargs='+', help='A list of overlapped tblout files')
    parser.add_argument('subunit_type', choices=['16S', '18S'], help='rRNA subunit type (16S vs 18S)')
    parser.add_argument('-o', '--output_file', default='amplified_regions', help='Prefix for the outfile name')
    return parser.parse_args(argv)


def main(argv):
    t_start = time.perf_counter()  # time the run
    args = parse_args(argv)
    retrieve_regions(args.files, args.output_file, args.subunit_type)
    t_stop = time.perf_counter()
    t_fact = t_stop - t_start
    print('Elapsed time:', '{0:.2f}'.format(t_fact), 'seconds')


if __name__ == '__main__':
    main(sys.argv[1:])

