import os
import re

import pandas as pd

from Pima.pima_data import PimaData
from Pima.utils.settings import Settings

from Pima.utils.utils import (
    print_and_log,
    print_and_run,
    validate_utility,
    validate_file_and_size_or_error,
    make_start_file,
    make_finish_file,
    std_files,
    find_checkpoint,
)

def validate_contamination_check(pima_data: PimaData, settings: Settings):

    if not pima_data.contam_check:
        return

    print_and_log(
        pima_data,
        'Validating contamination check', 
        pima_data.main_process_verbosity, 
        pima_data.main_process_color,
    )
    
    if not pima_data.will_have_ont_fastq and pima_data.illumina_fastq is None:
        pima_data.errors.append('--contamination requires a set of FASTQ reads')
    
    if validate_utility(pima_data, 'kraken2', 'kraken2 is not on the PATH (required by --contamination-check).'):
        command = 'kraken2 --version'
        pima_data.versions['kraken2'] = re.search(r'[0-9]+\.[0-9.]+', print_and_run(pima_data, command)[0]).group(0)

    if os.path.isdir(settings.kraken_database_default): 
        pima_data.kraken_database = settings.kraken_database_default
    elif os.path.isdir(settings.DockerPathKraken):
        pima_data.kraken_database = settings.DockerPathKraken
    else:
        pima_data.errors.append("No kraken2 database detected, try and run pima with --download. Exiting now.")
    
    pima_data.analysis.append(['fastq_contamination', pima_data, settings])


def fastq_contamination(pima_data: PimaData, settings: Settings):

    print_and_log(
        pima_data,
        'Running Kraken2 to check for contamination', 
        pima_data.main_process_verbosity, 
        pima_data.main_process_color,
    )

    pima_data.kraken_dir = os.path.join(pima_data.output_dir, 'contamination')

    if find_checkpoint(pima_data, pima_data.kraken_dir):
        print_and_log(
            pima_data,
            'Using existing kraken2 report',
            pima_data.sub_process_verbosity,
            pima_data.sub_process_color,
        )
        pima_data.did_kraken_fastq = True
        if os.path.isdir(os.path.join(pima_data.kraken_dir, "ont")):
            pima_data.kraken_fracs['ONT'] = read_kraken_report(os.path.join(pima_data.kraken_dir, "ont", "kraken.report"))
        
        if os.path.isdir(os.path.join(pima_data.kraken_dir, "illumina")):
            pima_data.kraken_fracs['Illumina'] = read_kraken_report(os.path.join(pima_data.kraken_dir, "illumina", "kraken.report"))
        return

    os.makedirs(pima_data.kraken_dir)
    make_start_file(pima_data, pima_data.kraken_dir)

    if not (pima_data.ont_fastq is None):
        print_and_log(
            pima_data,
            'Running Kraken2 on ONT data', 
            pima_data.sub_process_verbosity, 
            pima_data.sub_process_color,
        )
        ont_kraken_dir = os.path.join(pima_data.kraken_dir, 'ont')
        pima_data.kraken_fracs['ONT'] = kraken_fastq(pima_data, settings, pima_data.ont_fastq, ont_kraken_dir)

    if not (pima_data.illumina_fastq is None):
        print_and_log(
            pima_data,
            'Running Kraken2 on Illumina data', 
            pima_data.sub_process_verbosity, 
            pima_data.sub_process_color,
        )
        illumina_kraken_dir = os.path.join(pima_data.kraken_dir, 'illumina')
        pima_data.kraken_fracs['Illumina'] = kraken_fastq(pima_data, settings, pima_data.illumina_fastq, illumina_kraken_dir)

    pima_data.did_kraken_fastq = True
    make_finish_file(pima_data, pima_data.kraken_dir)

            
def kraken_fastq(pima_data: PimaData, settings: Settings, fastq, fastq_dir: str):

    os.makedirs(fastq_dir)
    
    kraken_files = [os.path.join(fastq_dir, 'kraken.' + i) for i in ['report', 'out', 'class', 'unclass']]
    kraken_report, kraken_out, kraken_class, kraken_unclass = kraken_files
    kraken_stdout, kraken_stderr = std_files(os.path.join(fastq_dir, 'kraken'))

    fastq_arg = fastq
    if isinstance(fastq, list):
        fastq_arg = ' '.join(fastq)

    command = " ".join(
        [
            'kraken2',
            '--threads', str(pima_data.threads),
            '--report', kraken_report,
            '--out', kraken_out,
            '--class', kraken_class,
            '--unclass', kraken_unclass,
            '--db', pima_data.kraken_database,
            fastq_arg,
            '1>', kraken_stdout, '2>', kraken_stderr,
        ]
    )
    print_and_run(pima_data, command)

    [validate_file_and_size_or_error(pima_data, i, i + ' missing after Kraken2', i + ' file is size 0 after Kraken2', )
        for i in kraken_files]

    # Read in the Kraken fractions and pull out the useful parts
    kraken_fracs = read_kraken_report(kraken_report)
    pima_data.files_to_clean.append([kraken_class, kraken_unclass, kraken_out])
    return(kraken_fracs)

def read_kraken_report(kraken_report: str):
    kraken_fracs = pd.read_csv(kraken_report, delimiter = '\t', header = None)
    kraken_fracs.index = kraken_fracs.iloc[:, 4].values
    kraken_fracs = kraken_fracs.loc[kraken_fracs.iloc[:, 3].str.match('[UG]1?'), :]
    kraken_fracs = kraken_fracs.loc[(kraken_fracs.iloc[:, 0] >= 1) | (kraken_fracs.iloc[:, 3] == 'U'), :]
    kraken_fracs = kraken_fracs.iloc[:, [0, 1, 3, 5]]
    kraken_fracs.columns = ['Fraction', 'Reads', 'Level', 'Taxa']
    kraken_fracs['Fraction'] = (kraken_fracs['Fraction'] / 100).round(4)
    kraken_fracs.sort_values(by = 'Fraction', inplace = True, ascending = False)
    kraken_fracs['Taxa'] = kraken_fracs['Taxa'].str.lstrip()
    return kraken_fracs