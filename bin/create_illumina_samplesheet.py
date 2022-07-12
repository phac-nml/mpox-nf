#!/usr/bin/env python3

"""
This small script takes a directory of paired illumina fastq files and outputs
a samplesheet that has the sample name and paired files in it
"""

import argparse
import csv
import glob
import logging
import os
import re

from collections import defaultdict
from pathlib import Path

# Make file regexes
FULL_ILLUMINA_REGEX = re.compile(r'([^/]*)(_S\d+)(_L\d{3})?(_R[12]|\.pair[12])(_\d{3})?\.fastq(\.gz)?$')
SIMPLE_ILLUMINA_REGEX = re.compile(r'([^/]*)(_R[12]|\.pair[12])(_\d{3})?\.fastq(\.gz)?$')

def init_parser() -> argparse.ArgumentParser:
    '''
    Purpose
    -------
    Parse CL inputs to be used in script
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-d',
        '--directory',
        required = True,
        help = 'Directory of sent data paired illumina fastq data'
    )
    parser.add_argument(
        '-o',
        '--out_csv',
        required = False,
        default = 'samplesheet.csv',
        help = 'Output csv file name'
    )

    return parser

def _add_to_dict(dict, full_path, name, read, gzip) -> None:
    '''
    Helper function to add found values to read dict
    '''
    # Add to key
    if name in dict:
        if (read == '_R1') or (read == '.pair1'):
            dict[name]['forward'] = full_path
        elif (read == '_R2') or (read == '.pair2'):
            dict[name]['reverse'] = full_path
        else:
            return
    # If not in dict, create key
    else:
        if (read == '_R1') or (read == '.pair1'):
            dict[name]['forward'] = full_path
        elif (read == '_R2') or (read == '.pair2'):
            dict[name]['reverse'] = full_path
        else:
            return
        
        # Gzipped?
        if gzip == '.gz':
            dict[name]['gzipped'] = True
        else:
            dict[name]['gzipped'] = False

def parse_indir_to_dict(indir) -> dict:
    '''
    Purpose
    -------
    Parse input 

    Parameters
    ----------
    indir : PosixPath
        Absolute path to input fastq directory

    Returns
    -------
    read_dict : dict
        Sample name dict containing forward and reverse read pairs along with gzipped status
            Ex. {'sampleName': {'forward': '/Path/to/sampleName_R1.fastq.gz', 'reverse': '/Path/to/sampleName_R2.fastq.gz', 'gzipped;: True}}
    '''
    paired_reads_dict = defaultdict(lambda: defaultdict(str))
    for file in glob.glob('{}/*.fastq*'.format(indir)):
        full_path = os.path.join(indir, file)
        
        # Match and add to output
        full_match = re.search(FULL_ILLUMINA_REGEX, file)
        if full_match:
            _add_to_dict(paired_reads_dict, full_path, full_match.group(1), full_match.group(4), full_match.group(6))
            continue

        simple_match = re.search(SIMPLE_ILLUMINA_REGEX, file)
        if simple_match:
            _add_to_dict(paired_reads_dict, full_path, simple_match.group(1), simple_match.group(2), simple_match.group(4))
            continue

        # Not recognized extension
        logging.warning('File: {} does not currently match any paired file patterns'.format(file))

    return paired_reads_dict

def output_samplesheet(read_dict, outfile) -> None:
    '''
    Purpose
    -------
    Output samples with found pairs to csv file

    Parameters
    ----------
    read_dict : dict
        Sample name dict containing forward and reverse read pairs along with gzipped status
            Ex. {'sampleName': {'forward': '/Path/to/sampleName_R1.fastq.gz', 'reverse': '/Path/to/sampleName_R2.fastq.gz', 'gzipped;: True}}

    outfile : str
        Outfile name or path
    '''
    logging.info('Writing data to outfile {}'.format(outfile))
    with open(outfile, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['sample', 'read1', 'read2', 'gzipped']) # header
        # Write rows if they have a forward and reverse file found
        for sample, value in read_dict.items():
            if ('forward' in value.keys()) and ('reverse' in value.keys()):
                writer.writerow([sample, value['forward'], value['reverse'], value['gzipped']])
            else:
                logging.warning('Sample {} is unpaired'.format(sample))

def main():
    ## Init Parser and Set Arguments ##
    parser = init_parser()
    args = parser.parse_args()

    ## Logging and Tracking Initialize ##
    log_name = 'setup.log'
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s: %(message)s',
        handlers=[logging.FileHandler(log_name), logging.StreamHandler()]
        )

    # Get full path
    indir = Path(args.directory).resolve()
    if (not indir.exists()) or (not indir.is_dir()):
        logging.error('Input dir {} is not a valid directory'.format(indir.name))
        exit(1)
    logging.info('Input dir: "{}" has absolute path of "{}"'.format(indir.name, indir))
    
    # Create dict of files found
    read_dict = parse_indir_to_dict(indir)

    # Output samplesheet csv
    output_samplesheet(read_dict, args.out_csv)

if __name__ == "__main__":
    main()
