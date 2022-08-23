#!/usr/bin/env python3

import argparse
import configparser
import pandas as pd

from iridauploader import api, model

def init_parser() -> argparse.ArgumentParser:
    '''
    Purpose
    -------
    Parse CL Parameters to be used in script
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-c',
        '--config',
        required=True,
        help='Irida Uploader config file to parse'
    )
    parser.add_argument(
        '-m',
        '--metadata_csv',
        required=True,
        help='Metadata qc csv for upload containing matching ids to SampleList.csv in the FIRST column and a project id in the SECOND'
    )
    parser.add_argument(
        '--no_sample_creation',
        required=False,
        default=False,
        action='store_true',
        help='Turns off sample creation for sample names not found in IRIDA'
    )
    return parser

def generate_api_instance(config_f) -> api.ApiCalls:
    '''
    Purpose
    -------
        Generate api instance for IRIDA platform based on input keys found in csv file

    Parameters
    ----------
        config: path
            IRIDA Uploader config file to be used to connect to IRIDA

    Returns
    -------
        irida_api: ApiCalls object
            Connection to IRIDA API based on keys given in config
    '''
    config = configparser.ConfigParser()
    config.read(config_f)
    settings = config['Settings']
    irida_api = api.ApiCalls(settings['client_id'], settings['client_secret'], settings['base_url'], settings['username'], settings['password'])
    return irida_api

def _create_sample_status_dict(sample: str, project_id: int, uploaded: bool, status: str) -> dict:
    """
    Purpose
    -------
        Create and return individual dictionaries for tracking status of uploaded metadata samples

    Parameters
    ----------
        sample: str
            Sample name
        project_id: int
            IRIDA project sample was sent to
        uploaded: bool
            Was data uploaded
        status: str
            Info on upload bool

    Returns
    -------
        d: dict
            Dictionary containing above combined information
    """
    d = {
        'sample': sample,
        'project_id': project_id,
        'uploaded': uploaded,
        'status': status
    }
    return d

def send_metadata(api_instance: api.ApiCalls, metadata_csv: str, no_sample_creation: bool) -> list:
    '''
    Purpose
    -------
        Send metadata for each sample in the input CSV file given a valid project ID

    Parameters
    ----------
        irida_api: ApiCalls object
            Connection to IRIDA API based on keys given in config
        metadata_csv: path
            CSV file that contains all of the metadata along with the sample name and project id
        no_sample_creation: bool
            If true, samples not already in IRIDA are skipped instead of created

    Returns
    -------
        tracking_dict_list: list(dicts)
            List of dictionarys containing the status of each samples upload
    '''
    df_in_dict = pd.read_csv(metadata_csv).fillna('NA').to_dict(orient='records')
    tracking_dict_list = [] # List to append dicts that track the status of samples for data uploads

    # Process and Upload each entry in input csv
    for metadata_dict in df_in_dict:
        sample_name = str(metadata_dict.pop('sample'))
        project_id = metadata_dict.pop('project_id')

        # Check if we have a project ID that is an integer for IRIDA
        try:
            project_id_int = int(project_id)
        except ValueError:
            print('Unknown Project ID for sample {}, moving to next sample'.format(project_id, sample_name))
            tracking_dict_list.append(_create_sample_status_dict(sample_name, project_id, False, 'Unknown Project ID'.format()))
            continue

        # Check that sample exists and if the new values are better than previous
        sample_id = api_instance.get_sample_id(sample_name=sample_name, project_id=project_id_int)
        if sample_id:
            # We have a sample in IRIDA already. Can add whatever checks we want to restrict uploads here
            pass
        # If given, skip samples that are not in IRIDA already
        elif no_sample_creation:
            print('Skipped sample {} metadata as it is not in IRIDA and --no_sample_creation arg passed'.format(sample_name))
            tracking_dict_list.append(_create_sample_status_dict(sample_name, project_id_int, False, 'Not in IRIDA and --no_sample_creation arg passed'))
            continue
        # If sample does not already exist, create it and grab its sample_id
        else:
            irida_sample = model.Sample(sample_name=sample_name)
            response = api_instance.send_sample(sample=irida_sample, project_id=project_id_int)
            sample_id = response['resource']['identifier']

        # Push data to IRIDA
        upload_metadata = model.Metadata(metadata=metadata_dict, project_id=project_id_int, sample_name=sample_name)
        status = api_instance.send_metadata(upload_metadata, sample_id)
        print('Uploaded {} metadata to {}'.format(sample_name, project_id_int))
        tracking_dict_list.append(_create_sample_status_dict(sample_name, project_id_int, True, ''))

    return tracking_dict_list

def main() -> None:
    """
    Purpose
    -------
        Main process. Checks and uploads metadata to IRIDA
    """
    # Init Parser and set arguments
    parser = init_parser()
    args = parser.parse_args()

    # Connect and upload metadata
    irida_api = generate_api_instance(args.config)
    tracking_dict_list = send_metadata(irida_api, args.metadata_csv, args.no_sample_creation)

    # Output tracking df
    df = pd.DataFrame.from_dict(tracking_dict_list)
    df.to_csv('metadata_upload_status.csv', index=False)

if __name__ == "__main__":
    main()
