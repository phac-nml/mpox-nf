#!/usr/bin/env python3

import argparse
import os
import pathlib
import re
import subprocess
import pandas as pd

from collections import defaultdict

def init_parser() -> argparse.ArgumentParser:
    '''
    Purpose
    -------
    Parse CL inputs to be used in script
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--directory',
        '-d',
        required=False,
        default=False,
        help='Path to directory containing dehosted paired fastq files or consensus sequences to rename and add an upload SampleList.csv file to'
    )
    parser.add_argument(
        '--samplesheet',
        '-s',
        required=True,
        help='File path to comma separated samplesheet with the following columns: ["sample", "project_id", "sequencing_date"]'
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--fastq', action='store_true', help="Create SampleList for paired illumina fastq files")
    group.add_argument('--fasta', action='store_true', help="Create SampleList for consensus fasta files")
    return parser

def _add_to_paired_dict(dict: dict, full_path: str, name: str, read: str) -> None:
    '''
    Purpose
    -------
        Helper function to add found values to read dict

    Parameters
    ----------
        dict: dict
            Dictionary format paired read values to
        full_path: str
            Path to input file
        name: str
            Name of input file
        read: str
            Read number. Either _R1 or _R2
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

def create_sample_file_df(samplesheet: str, directory: str, file_type='') -> pd.DataFrame:
    '''
    Purpose
    -------
        Take input samplesheet, directory, and file type to create a df of values needed for IRIDA Uploads along
        with renaming the files to include the sequencing date for better tracking

    Parameters
    ----------
        samplesheet: path
            CSV samplesheet containing the columns ["sample", "project_id", "sequencing_date"]
        directory: path
            Path to directory containing the files wanting to upload
        file_type: str
            Eithe fasta or fastq for what type of files we want to upload

    Returns
    -------
        df_out: pd.DataFrame
            Dataframe with files formatted to be used by IRIDA Uploader
    '''
    # Read in input TSV file
    df = pd.read_csv(samplesheet, sep=',')

    # Setup Dict of files -  Ex. {'sampleName': {'forward': '/Path/to/sampleName_R1.fastq', 'reverse': '/Path/to/sampleName_R2.fastq'}}
    file_dict = defaultdict(lambda: defaultdict(str))
    
    # Set file regex based on data type. Regex statements are based on pipeline output names so if those change we must adjust here
    if file_type == 'fastq':
        FILE_WANTED_REGEX = re.compile(r'([^/]*)_dehosted.*(_R[12])\.fastq$')
    else:
        FILE_WANTED_REGEX = re.compile(r'(^.+)\.consensus\.(fa|fasta)$')

    # Populate dict with files
    for file_found in os.listdir(directory):
        match = re.search(FILE_WANTED_REGEX, file_found)
        if match and file_type =='fastq':
            _add_to_paired_dict(file_dict, file_found, match.group(1), match.group(2))
        elif match and file_type =='fasta':
            file_dict[match.group(1)]['forward'] = file_found
            file_dict[match.group(1)]['reverse'] = ''
        else:
            print('WARNING: No matches for file {}'.format(file_found))

    # Generate Out DF from values in input table that have matching files found
    row_list = []
    for sample, proj_id, seq_date in zip(df['sample'], df['project_id'], df['sequencing_date']):
        # Check if sample string a is a substring of any files and then add it to df_out if so
        if sample in file_dict.keys():
            ext = ''.join(pathlib.Path(file_dict[sample]['forward']).suffixes)

            # Rename for better tracking
            if file_type == 'fastq':
                full_path_f = pathlib.Path(directory, file_dict[sample]['forward'])
                new_file_name_f = '{}_{}_R1{}'.format(sample, seq_date, ext)
                new_path_f = pathlib.Path(directory, new_file_name_f)
                
                full_path_r = pathlib.Path(directory, file_dict[sample]['reverse'])
                new_file_name_r = '{}_{}_R2{}'.format(sample, seq_date, ext)
                new_path_r = pathlib.Path(directory, new_file_name_r)

                subprocess.run(['mv', full_path_f, new_path_f])
                subprocess.run(['mv', full_path_r, new_path_r])
            else:
                full_path_f = pathlib.Path(directory, file_dict[sample]['forward'])
                new_file_name_f = '{}_{}{}'.format(sample, seq_date, ext)
                new_path_f = pathlib.Path(directory, new_file_name_f)

                new_file_name_r = ''

                subprocess.run(['mv', full_path_f, new_path_f])
            row_dict = {
                    'Sample_Name': sample,
                    'Project_ID': proj_id, 
                    'File_Forward': new_file_name_f, 
                    'File_Reverse': new_file_name_r
                }
            row_list.append(row_dict)
    df_out = pd.DataFrame(row_list, columns=['Sample_Name', 'Project_ID', 'File_Forward', 'File_Reverse'])
    return df_out

def main() -> None:
    """
    Purpose
    -------
        Main process. Checks and uploads metadata to IRIDA
    """
    # Init Parser and set arguments
    parser = init_parser()
    args = parser.parse_args()

    # Create df_out using directory files and the sample TSV file based on the input file type wanted
    if args.fastq:
        df_out = create_sample_file_df(args.samplesheet, args.directory, file_type='fastq')
    else:
        df_out = create_sample_file_df(args.samplesheet, args.directory, file_type='fasta')

    # Output
    with open('{}/SampleList.csv'.format(args.directory), 'w') as handle:
        handle.write('[Data]\n')
    
    df_out.to_csv('{}/SampleList.csv'.format(args.directory), mode='a', header=True, index=False)

if __name__ == "__main__":
    main()
