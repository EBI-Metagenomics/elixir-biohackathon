import re
import os
import logging
import glob

raw_f_regex = re.compile(
    "([A-z0-9.-]+)\s+-\s+(\w+)\s+(\w+)\s+(\w+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([-+])\s+([-+])\s+(\d+)\s+(\d+[\.\d]*)\s+(\d+[\.\d]*)\s+(\d+[\.\d]*)\s+(.+)\s!\s+ -")
RAW_DATA_FILE = 'SRR4341428_MERGED_FASTQ.cmsearch.all.tblout.deoverlapped_V1_V2'

MIN_OVERLAP = 0.95
regions = {
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

logging.basicConfig(level=logging.INFO)


def calc_overlap(read, reg):
    read_s, read_f = read
    reg_s, reg_f = reg
    if read_s>read_f:
        print(read_s, read_f)
        input()
    overlap = max(min(read_f, reg_f) - max(read_s, reg_s), 0)
    total = max(reg_f - reg_s, 0)
    try:
        return overlap / total
    except ZeroDivisionError:
        return 0


def get_regions(raw_sequence):
    return {region for region, limits in regions.items() if calc_overlap(raw_sequence, limits) >= MIN_OVERLAP}


# Parse, filter empty lines and unpack into 2D array
def load_data(filename=RAW_DATA_FILE):
    with open(filename) as f:
        return [l[0] for l in [raw_f_regex.findall(l) for l in f] if bool(l)]


def print_matches(region_matches):
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
            results[k] = round(ratio, 4)
    results = sorted(results.items(), key=lambda x: x[1])
    results.reverse()
    return results


def main():
    for f in glob.glob(os.path.join(os.getcwd(), 'examples', '*')):
        data = load_data(f)
        region_matches = []
        for read in data:
            limits = list(map(int, read[4:6]))
            region_matches.extend(get_regions(limits))
        logging.info(f + " " + str(print_matches(region_matches)))


if __name__ == '__main__':
    main()
