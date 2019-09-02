#!/usr/bin/env python3
import argparse
import tempfile
import shutil
import sys

def parse_args():
    '''Argument parsin'''
    description = """
    parsing cap3 assembly aln output
    """

    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        '-g',
        '--gff_file',
        default=None,
        required=True,
        help="input gff3 file for appending coverage information",
        type=str,
        action='store')
    parser.add_argument(
        '-p',
        '--profile',
        default=None,
        required=True,
        help="output file for coverage profile",
        type=str,
        action="store")
    return parser.parse_args()

def read_coverage(profile):
    with open(profile) as p:
        d = {}
        for name, prof in zip(p, p):
            d[name[1:].strip()] = [int(i) for i in prof.split()]
    print(d, file=sys.stderr)
    return d


def main():
    args = parse_args()
    coverage_hash = read_coverage(args.profile)
    gff_tmp = tempfile.NamedTemporaryFile()
    with open(args.gff_file) as f, open(gff_tmp.name, 'w') as out:
        for line in f:
            if line[0] == "#":
                out.write(line)
            else:
                line_parts = line.split()
                start = int(line_parts[3])
                end = int(line_parts[4])
                coverage = round( sum(coverage_hash[line_parts[0]][(
                    start - 1):end]) / (end - start + 1), 3)
                new_line = "{};Coverage={}\n".format(line.strip(), coverage)
                out.write(new_line)

    shutil.copyfile(gff_tmp.name, args.gff_file)


if __name__ == "__main__":

    main()
