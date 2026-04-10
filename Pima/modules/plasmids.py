import os
import re
import shutil

import pandas as pd

from Pima.pima_data import PimaData
from Pima.utils.settings import Settings

from Pima.utils.utils import (
    print_and_log,
    print_and_run,
    validate_utility,
    validate_file_and_size,
    validate_file_and_size_or_error,
    make_start_file,
    make_finish_file,
    find_checkpoint,
    std_files,
)

from .annotations import make_blast_database

def validate_plasmids(pima_data: PimaData, settings: Settings):

    if not pima_data.plasmids:
        return

    print_and_log(
        pima_data,
        'Validating plasmid database and utilities', 
        pima_data.main_process_verbosity, 
        pima_data.main_process_color,
    )
    
    for utility in ['minimap2']:
        if validate_utility(pima_data, utility, utility + ' is not on the PATH.'):
            command = utility + ' --version'
            pima_data.versions[utility] = re.search(r'[0-9]+\.[0-9.]*', print_and_run(pima_data, command)[0]).group(0)

    for utility in ['Rscript', 'R']:
        if validate_utility(pima_data, utility, utility + ' is not on the PATH.'):
            command = utility + ' --version 2>&1'
            pima_data.versions[utility] = re.search(r'[0-9]+\.[0-9.]*', print_and_run(pima_data, command)[0]).group(0)
        
    if not validate_file_and_size(pima_data, the_file = pima_data.plasmid_database, min_size = 2000):
        pima_data.errors.append(f"Can't find plasmid database {pima_data.plasmid_database} or is empty. Try --download?")
        
    if not pima_data.will_have_genome_fasta:
        pima_data.errors.append("Can't call plasmids without a genome or an assembly")

    pima_data.analysis.append(['call_plasmids', pima_data, settings])


def call_plasmids(pima_data: PimaData, settings: Settings):

    print_and_log(
        pima_data,
        'Calling plasmids', 
        pima_data.main_process_verbosity, 
        pima_data.main_process_color,
    )

    # Make a directory for plasmid stuff
    pima_data.plasmid_dir = os.path.join(pima_data.output_dir, 'plasmids')
    
    #Checkpoint
    if find_checkpoint(pima_data, pima_data.plasmid_dir):
        pima_data.did_call_plasmids = True
        pima_data.plasmid_tsv = os.path.join(pima_data.plasmid_dir, 'plasmids.tsv')
        try:
            pima_data.plasmids = pd.read_csv(filepath_or_buffer = pima_data.plasmid_tsv, sep = '\t', header = 0)
        except:
            pima_data.plasmids = None  
        
        #Retrieve notes
        smaller_contigs_fasta = os.path.join(pima_data.plasmid_dir, 'small_contigs.fasta')
        small_contigs = pima_data.load_fasta(smaller_contigs_fasta)
        align_stats = pima_data.query_alignment_stats.copy(deep=True)
        align_stats['Size (bp)'] = pd.to_numeric(align_stats['Size (bp)'])
        smaller_contigs = align_stats[(align_stats['Size (bp)'] < 500000.) & (align_stats['Perc Align'] < 75.)]['Contig'].to_list()
        if len(smaller_contigs) == 0 and pima_data.reference_fasta:
            message = 'No contigs smaller than 500kb that do not align to the reference provided found, skipping plasmid search'
            pima_data.plasmid_notes = pd.concat([pima_data.plasmid_notes, pd.Series(message, dtype='object')])        
        elif len(small_contigs) == 0:
            message = 'No contigs smaller than 500kb found, skipping plasmid search'
            pima_data.plasmid_notes = pd.concat([pima_data.plasmid_notes, pd.Series(message, dtype='object')])
        elif pima_data.plasmids is None:
            message = 'No small contigs matched plasmids within the database after excluding pX01 and pX02'
            pima_data.plasmid_notes = pd.concat([pima_data.plasmid_notes, pd.Series(message, dtype='object')])
        return
    
    os.makedirs(pima_data.plasmid_dir)
    make_start_file(pima_data, pima_data.plasmid_dir)

    # Some of the plasmids in the database have regions that align to the Ba reference genome at the SSU rRNA region
    # these can show contigs that align to "ecoli" plasmids, but they actually align better to the Ba genome
    # The plasmid annotation step was reporting all small contigs, but lets actually remove ones that align well to the reference genome

    if pima_data.reference_fasta:
        print_and_log(
            pima_data,
            'Finding contigs that have at least 25% of their length not aligning to the reference genome and are < 500000 bp',
            pima_data.sub_process_verbosity,
            pima_data.sub_process_color,
        )
        smaller_contigs_fasta = os.path.join(pima_data.plasmid_dir, 'small_contigs.fasta')
        align_stats = pima_data.query_alignment_stats.copy(deep=True)
        align_stats['Size (bp)'] = pd.to_numeric(align_stats['Size (bp)'])
        smaller_contigs = align_stats[(align_stats['Size (bp)'] < 500000.) & (align_stats['Perc Align'] < 75.)]['Contig'].to_list()
        smaller_contigs_fasta = os.path.join(pima_data.plasmid_dir, 'small_contigs.fasta')
        if len(smaller_contigs) == 0:
            message = 'No contigs smaller than 500kb that do not align to the reference provided found, skipping plasmid search'
            print_and_log(
                pima_data,
                message, 
                pima_data.sub_process_verbosity, 
                pima_data.sub_process_color,
            )
            pima_data.did_call_plasmids = True
            pima_data.plasmids = None
            make_finish_file(pima_data, pima_data.plasmid_dir)
            pima_data.plasmid_notes = pd.concat([pima_data.plasmid_notes, pd.Series(message, dtype='object')])
            return
        else:
            command = " ".join(
                [
                    'faidx', pima_data.genome_fasta, " ".join(smaller_contigs), '>', smaller_contigs_fasta,
                ]
            )
            print_and_run(pima_data, command)
    ## If there isn't a reference provided check all the short contigs
    # Take very large things out of the assembly.  They aren't plasmids and take a long time to run
    else:
        print_and_log(
            pima_data,
            'Finding contigs < 500000 bp', 
            pima_data.sub_process_verbosity, 
            pima_data.sub_process_color,
        )
        smaller_contigs_fasta = os.path.join(pima_data.plasmid_dir, 'small_contigs.fasta')
        command = " ".join(
            [
                'faidx -i chromsizes', 
                pima_data.genome_fasta,
                '| awk \'($2 <= 500000){print $1}\'',
                '| parallel -n1 -n1 faidx', pima_data.genome_fasta, '>', smaller_contigs_fasta,
            ]
        )
        print_and_run(pima_data, command)

        small_contigs = pima_data.load_fasta(smaller_contigs_fasta)
        if len(small_contigs) == 0:
            message = 'No contigs smaller than 500kb found, skipping plasmid search'
            print_and_log(
                pima_data,
                message, 
                pima_data.sub_process_verbosity, 
                pima_data.sub_process_color,
            )
            pima_data.did_call_plasmids = True
            pima_data.plasmids = None
            make_finish_file(pima_data, pima_data.plasmid_dir)
            pima_data.plasmid_notes = pd.concat([pima_data.plasmid_notes, pd.Series(message, dtype='object')])
            return
                                
    # Query plasmid sequences against the assembly using minimap2
    print_and_log(
        pima_data,
        'Running minimap2 against the plasmid database', 
        pima_data.sub_process_verbosity, 
        pima_data.sub_process_color,
    )
    plasmid_sam = os.path.join(pima_data.plasmid_dir, 'plasmid_hits.sam')
    _, minimap_stderr = std_files(os.path.join(pima_data.plasmid_dir, 'minimap'))
    command = " ".join(
        [
            'minimap2',
            '-k 20 -p .2 -a',
            '-t', str(pima_data.threads),
            smaller_contigs_fasta,
            pima_data.plasmid_database,
            '1>', plasmid_sam,
            '2>', minimap_stderr,
        ]
    )
    print_and_run(pima_data, command)
    validate_file_and_size_or_error(pima_data, plasmid_sam, 'Plasmid v. contig SAM', 'cannot be found', 'is empty')

    pima_data.files_to_clean.append(plasmid_sam)
    
    # Turn the SAM file in to a PSL file using the modified sam2psl script
    print_and_log(
        pima_data,
        'Converting the SAM file to a PSL file', 
        pima_data.sub_process_verbosity, 
        pima_data.sub_process_color,
    )
    plasmid_psl = os.path.join(pima_data.plasmid_dir, 'plasmid_hits.psl')
    if not validate_file_and_size(pima_data, plasmid_psl):
        message = 'No small contigs matched plasmids within the database after excluding pX01 and pX02'
        print_and_log(
            pima_data,
            message, 
            pima_data.sub_process_verbosity, 
            pima_data.sub_process_color,
        )
        pima_data.did_call_plasmids = True
        pima_data.plasmids = None
        make_finish_file(pima_data, pima_data.plasmid_dir)
        pima_data.plasmid_notes = pd.concat([pima_data.plasmid_notes, pd.Series(message, dtype='object')])
        return

    sam2psl_stdout, sam2psl_stderr = std_files(os.path.join(pima_data.plasmid_dir, 'sam2psl'))
    path2sam2psl = os.path.join(settings.pima_path, "accessory_scripts", "sam2psl.py")
    command = " ".join(
        [
            'python3',
            path2sam2psl,
            '-i', plasmid_sam,
            '-o', plasmid_psl,
            '1>', sam2psl_stdout, '2>', sam2psl_stderr,
        ]
    )
    print_and_run(pima_data, command)
    validate_file_and_size_or_error(pima_data, plasmid_sam, 'Plasmid v. contig PSL', 'cannot be found', 'is empty')
    
    # Make a BLAST database of the plasmid sequences
    make_blast_database(pima_data, pima_data.plasmid_database)
    
    # Pass the data onto pChunks
    print_and_log(
        pima_data,
        'Running pChunks', 
        pima_data.sub_process_verbosity, 
        pima_data.sub_process_color,
    )
    pima_data.pchunks_dir = os.path.join(pima_data.plasmid_dir, 'pChunks')
    os.makedirs(pima_data.pchunks_dir)
    
    pima_data.plasmid_tsv = os.path.join(pima_data.pchunks_dir, 'plasmids.tsv')
    stdout_file, stderr_file = std_files(os.path.join(pima_data.plasmid_dir, "pchunks"))
    path2pChunks = os.path.join(settings.pima_path, "accessory_scripts", "pChunks.R")
    command = " ".join(
        [
            'Rscript',
            path2pChunks, '--plasmid-psl', plasmid_psl,
            '--output', pima_data.pchunks_dir,
            '--no-amr', '--no-inc',
            '--plasmid-database', pima_data.plasmid_database,
            '--threads', str(pima_data.threads),
            '1>', stdout_file, '2>', stderr_file,
        ]
    )
    print_and_run(pima_data, command)
    pima_data.plasmid_tsv = os.readlink(os.path.join(pima_data.pchunks_dir, 'plasmids.tsv'))
    if not validate_file_and_size(pima_data, pima_data.plasmid_tsv):
        message = 'No small contigs matched plasmids within the database after excluding pX01 and pX02'
        print_and_log(
            pima_data,
            message, 
            pima_data.sub_process_verbosity, 
            pima_data.sub_process_color,
        )
        pima_data.did_call_plasmids = True
        pima_data.plasmids = None
        make_finish_file(pima_data, pima_data.plasmid_dir)
        shutil.rmtree(pima_data.pchunks_dir, ignore_errors=True)
        pima_data.plasmid_notes = pd.concat([pima_data.plasmid_notes, pd.Series(message, dtype='object')])
        return

    # delete orphan amrFinal.bed file that pChunks creates even though we turn off amr calling
    os.remove(os.path.join(pima_data.pchunks_dir, 'amrFinal.bed'))
    
    # The final file is in pChunks
    new_plasmid_tsv = os.path.join(pima_data.plasmid_dir, 'plasmids.tsv')
    shutil.copy2(pima_data.plasmid_tsv, new_plasmid_tsv)
    pima_data.plasmid_tsv = new_plasmid_tsv

    try:
        pima_data.plasmids = pd.read_csv(filepath_or_buffer = pima_data.plasmid_tsv, sep = '\t', header = 0)
    except:
        pima_data.plasmids = None

    pima_data.did_call_plasmids = True
    make_finish_file(pima_data, pima_data.plasmid_dir)