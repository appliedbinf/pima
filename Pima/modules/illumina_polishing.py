import os
import re

import subprocess
import pandas as pd
import pathlib

from Pima.pima_data import PimaData
from Pima.utils.utils import (
    print_and_log,
    print_and_run,
    add_warning,
    validate_utility,
    validate_file_and_size_or_error,
    find_checkpoint,
    make_start_file,
    make_finish_file,
    make_report_info_file,
    std_files,
)
from Pima.utils.mapping import (
    bwa_short_illumina_fastq_and_sort,
    minimap_and_sort,
    bwa_mem_all_aln_illumina,
)

def validate_illumina_polish(pima_data: PimaData):
    if not pima_data.will_have_ont_assembly and not pima_data.will_have_genome_fasta:
        return

    if pima_data.illumina_polisher == "pilon":
        validate_pilon(pima_data)
    elif pima_data.illumina_polisher == "polypolish":
        validate_polypolish(pima_data)
    else:
        return


def validate_pilon(pima_data: PimaData):

    # skip conditions
    if pima_data.no_assembly:
        return
    if not (pima_data.will_have_genome_fasta and pima_data.illumina_fastq):
        return
    if not pima_data.illumina_polisher == "pilon":
        return
    
    print_and_log(
        pima_data,
        'Validating Pilon and memory arguments', 
        pima_data.main_process_verbosity, 
        pima_data.main_process_color,
    )

    for utility in ['minimap2', 'pilon']:
        if validate_utility(pima_data, utility, f"{utility} is not on the PATH (required by --illumina-fastq)."):
            command = f"{utility} --version"
            pima_data.versions[utility] = re.search(r'[0-9]+\.[0-9.]+', print_and_run(pima_data,command)[0]).group(0)

    pima_data.analysis.append(["pilon_assembly", pima_data])


def validate_polypolish(pima_data: PimaData):
    # skip conditions
    if pima_data.no_assembly:
        return
    if not (pima_data.will_have_genome_fasta and pima_data.illumina_fastq):
        return
    if not pima_data.illumina_polisher == "polypolish":
        return
    
    print_and_log(
        pima_data,
        'Validating polypolish', 
        pima_data.main_process_verbosity, 
        pima_data.main_process_color,
    )

    if validate_utility(pima_data, "polypolish", "polypolish is not on the PATH (required by --illumina-polisher polypolish)."):
        command = "polypolish --version"
        pima_data.versions['polypolish'] = re.search(r'[0-9]+\.[0-9.]+', print_and_run(pima_data,command)[0]).group(0)

    if validate_utility(pima_data, "bwa", "bwa is not on the path (required by --illumina-polisher polypolish)."):
        command = "bwa"
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        version_element = [i for i,x in enumerate(result.stderr.split()) if re.search("Version", x)]
        pima_data.versions['bwa'] = re.search(r'[0-9]+\.[0-9.]+', result.stderr.split()[version_element[0] + 1]).group(0)
    
    pima_data.analysis.append(["polypolish_assembly", pima_data])
    

def polypolish_assembly(pima_data: PimaData):

    print_and_log(
        pima_data,
        'Running polypolish on genome assembly', 
        pima_data.main_process_verbosity, 
        pima_data.main_process_color,
    )

    pima_data.illumina_polish_dir = os.path.join(pima_data.output_dir, 'illumina_polish')

    ## Check if pilon has been run before to completion
    if find_checkpoint(pima_data, pima_data.illumina_polish_dir):
        pima_data.genome_fasta = os.path.join(pima_data.illumina_polish_dir, 'assembly.fasta')
        pima_data.load_genome()

        print_and_log(
            pima_data,
            'Polypolish had previously been run and finished successfully', 
            pima_data.sub_process_verbosity, 
            pima_data.sub_process_color,
        )
        pima_data.did_polypolish_ont_assembly = True
        pima_data.files_to_clean.extend(
            list(pathlib.Path(pima_data.illumina_polish_dir).glob("*.bam"))
        )
        pima_data.files_to_clean.append(os.path.join(pima_data.illumina_polish_dir, 'polypolish.fasta'))
        return

    os.makedirs(pima_data.illumina_polish_dir)
    make_start_file(pima_data, pima_data.illumina_polish_dir)

    # Map illumina reads onto the assembly 
    #We need to do this again because polypolish needs separate bam files 
    #and all alignments from bwa
    print_and_log(
        pima_data,
        'Mapping Illumina reads to assembly', 
        pima_data.sub_process_verbosity, 
        pima_data.sub_process_color,
    )

    #This gets turned into 2 files by the bwa function if there are R1/R2 files provided
    polypolish_bam = os.path.join(pima_data.illumina_polish_dir, 'mapping.bam')

    bams = bwa_mem_all_aln_illumina(pima_data, pima_data.genome_fasta, pima_data.illumina_fastq, polypolish_bam)
    polypolish_cmd = "polypolish polish"
    # If mean coverage is low then give some warning & use correct form of polypolish
    if pima_data.mean_coverage['Illumina'] < 25:
        warning = f"Mean illumina coverage ({pima_data.mean_coverage['Illumina']}X) is below 25X, using careful mode)."
        pima_data.add_warning(warning)
        pima_data.assembly_notes = pd.concat([pima_data.assembly_notes, pd.Series(warning, dtype = 'object')])
        polypolish_cmd = "polypolish polish --careful"

    # Actually run polypolish
    print_and_log(
        pima_data,
        'Running Polypolish', 
        pima_data.sub_process_verbosity, 
        pima_data.sub_process_color,
    )

    pfilt_stdout, pfilt_stderr = std_files(os.path.join(pima_data.illumina_polish_dir, 'poly_filter'))
    polypolish_prefix = os.path.join(pima_data.illumina_polish_dir, 'polypolish')
    polished_fasta = polypolish_prefix + '.fasta'

    if len(bams) == 2:
        filt_bams = [re.sub(r"\.bam", r"_filt.bam", x) for x in bams]
        #you can't filter unpaired bams
        command = " ".join(
            [
                "polypolish filter",
                "--in1", bams[0],
                "--in2", bams[1],
                "--out1", filt_bams[0],
                "--out2", filt_bams[1],
                "1>", pfilt_stdout, "2>", pfilt_stderr,
            ]
        )
        print_and_run(pima_data, command)
        validate_file_and_size_or_error(pima_data, filt_bams[0], 'Polypolish filtered bam', 'cannot be found after filtering', 'is empty')
        _, polish_stderr = std_files(os.path.join(pima_data.illumina_polish_dir, 'polypolish'))
        command = " ".join(
            [
                polypolish_cmd,
                pima_data.genome_fasta,
                " ".join(filt_bams),
                "1>", polished_fasta,
                "2>", polish_stderr,
            ]
        )
        print_and_run(pima_data, command)
        pima_data.files_to_clean.extend(bams + filt_bams)

    else:
        #We have only a single fastq file so we can't filter the alignments
        _, polish_stderr = std_files(os.path.join(pima_data.illumina_polish_dir, 'polypolish'))
        command = " ".join(
            [
                polypolish_cmd,
                pima_data.genome_fasta,
                " ".join(bams),
                "1>", polished_fasta,
                "2>", polish_stderr,
            ]
        )
        print_and_run(pima_data, command)
        pima_data.files_to_clean.extend(bams)
    validate_file_and_size_or_error(pima_data, polished_fasta, 'Polypolished assembly', 'cannot be found after running polypolish', 'is empty')

    print_and_log(
        pima_data,
        'Repairing contig names after polypolish', 
        pima_data.sub_process_verbosity, 
        pima_data.sub_process_color,
    )
    pima_data.genome_fasta = os.path.join(pima_data.illumina_polish_dir, 'assembly.fasta')
    command = " ".join(
        [
            'cat', 
            polished_fasta,
            '| awk \'{if($0 ~ /^>/){gsub("_pilon", "", $0)}print}\'',
            '>', pima_data.genome_fasta,
        ]
    )
    print_and_run(pima_data, command)
    validate_file_and_size_or_error(pima_data, pima_data.genome_fasta, 'Genome assembly', 
                                            'cannot be found after fixing names', 'is empty')
    pima_data.files_to_clean.append(polished_fasta)
    make_finish_file(pima_data, pima_data.illumina_polish_dir)
    pima_data.did_polypolyish_ont_assembly = True #track for the report

def pilon_assembly(pima_data: PimaData):
    print_and_log(
        pima_data,
        'Running Pilon on genome assembly', 
        pima_data.main_process_verbosity, 
        pima_data.main_process_color,
    )
    pima_data.illumina_polish_dir = os.path.join(pima_data.output_dir, 'illumina_polish')

    ## Check if pilon has been run before to completion
    if find_checkpoint(pima_data, pima_data.illumina_polish_dir):
        pima_data.genome_fasta = os.path.join(pima_data.illumina_polish_dir, 'assembly.fasta')
        pilon_bam = os.path.join(pima_data.illumina_polish_dir, 'mapping.bam')
        pima_data.files_to_clean.append(pilon_bam)
        pima_data.load_genome()

        print_and_log(
            pima_data,
            'Pilon had previously been run and finished successfully', 
            pima_data.main_process_verbosity, 
            pima_data.main_process_color,
        )
        pima_data.did_pilon_ont_assembly = True
        return

    os.makedirs(pima_data.illumina_polish_dir)
    make_start_file(pima_data, pima_data.illumina_polish_dir)

    # Map illumina reads onto the assembly
    print_and_log(
        pima_data,
        'Mapping Illumina reads to assembly', 
        pima_data.sub_process_verbosity, 
        pima_data.sub_process_color,
    )
    pilon_bam = os.path.join(pima_data.illumina_polish_dir, 'mapping.bam')
    pima_data.files_to_clean.append(pilon_bam)

    # See what mapping method to use - bwa aln or minimap 2
    if pima_data.illumina_length_mean <= 50:
        bwa_short_illumina_fastq_and_sort(pima_data, pima_data.genome_fasta, pima_data.illumina_fastq, pilon_bam)
    else: # We have longer short reads
        minimap_and_sort(pima_data, pima_data.genome_fasta, pilon_bam, pima_data.illumina_fastq, ont=False)


    # Figure out the depth here.  If it's too low, give some sort of warning?
    coverage_tsv = os.path.join(pima_data.illumina_polish_dir, 'coverage.tsv')
    command = " ".join(
        [
            'samtools depth -a', 
            pilon_bam,
            '| awk \'{s += $3; c++}END{printf "%s\\t%i\\t%.0f\\n", i, c, (s / c)}\'',
            '>', coverage_tsv,
        ]
    )
    print_and_run(pima_data, command)
    validate_file_and_size_or_error(pima_data, coverage_tsv, 'Coverage TSV', 'cannot be found after samtools', 'is empty')
    pilon_coverage = pd.read_csv(coverage_tsv, header = None, index_col = None, sep = '\t').iloc[0,2]

    # If mean coverage is low then give some warning
    if pilon_coverage < pima_data.pilon_coverage_min:
        warning = f"Illumina coverage for Pilon {pilon_coverage}X is below the recommended minimum ({pima_data.pilon_coverage_min}X), we recommend using polypolish instead."
        add_warning(pima_data, warning)
        
    # Actually run pilon
    print_and_log(
        pima_data,
        'Running Pilon', 
        pima_data.sub_process_verbosity, 
        pima_data.sub_process_color,
    )
    pilon_stdout, pilon_stderr = std_files(os.path.join(pima_data.illumina_polish_dir, 'pilon'))
    pilon_prefix = os.path.join(pima_data.illumina_polish_dir, 'pilon')
    polished_fasta = pilon_prefix + '.fasta'
    bam_option = '--frags'
    if len(pima_data.illumina_fastq) == 1:
        bam_option = '--unpaired'

    command = " ".join(
        [
            'pilon',
            '-Xms4g -Xmx4g',
            '--genome', pima_data.genome_fasta,
            bam_option, pilon_bam,
            '--output', pilon_prefix,
            '1>', pilon_stdout, '2>', pilon_stderr,
        ]
    )
    print_and_run(pima_data, command)
    validate_file_and_size_or_error(pima_data, polished_fasta, 'Pilon FASTA', 'cannot be found after pilon', 'is empty')

    print_and_log(
        pima_data,
        'Repairing contig names after Pilon', 
        pima_data.sub_process_verbosity, 
        pima_data.sub_process_color,
    )
    pima_data.genome_fasta = os.path.join(pima_data.illumina_polish_dir, 'assembly.fasta')
    command = " ".join(
        [
            'cat', 
            polished_fasta,
            '| awk \'{if($0 ~ /^>/){gsub("_pilon", "", $0)}print}\'',
            '>', pima_data.genome_fasta,
        ]
    )
    print_and_run(pima_data, command)
    validate_file_and_size_or_error(pima_data, pima_data.genome_fasta, 'Genome assembly', 
                                            'cannot be found after fixing names', 'is empty')
    pima_data.load_genome()
    pima_data.did_pilon_ont_assembly = True
    make_finish_file(pima_data, pima_data.illumina_polish_dir)
