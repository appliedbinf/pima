import re
import os
import subprocess

import pandas as pd

from Pima.pima_data import PimaData
from Pima.utils.utils import (
    add_warning,
    print_and_log,
    validate_utility,
    validate_file_and_size_or_error,
    print_and_run,
    find_checkpoint,
    make_start_file,
    make_finish_file,
    std_files,
)

def validate_medaka(pima_data: PimaData):
    # Skip conditions
    if pima_data.no_assembly:
        return

    if pima_data.no_medaka:
        return

    if not pima_data.will_have_genome_fasta:
        return

    if not pima_data.will_have_ont_fastq:
        return

    print_and_log(
        pima_data,
        "Validating medaka",
        pima_data.main_process_verbosity,
        pima_data.main_process_color,
    )

    if validate_utility(
        pima_data,
        "medaka_consensus",
        "medaka_consensus is not on the PATH (required by --medaka)",
    ):
        command = "medaka --version"
        pima_data.versions["medaka"] = re.search(
            r"[0-9]+\.[0-9.]+", print_and_run(pima_data, command)[0]
        ).group(0)

    if pima_data.ont_model == 'auto':
        print_and_log(
            pima_data,
            "Trying to determine the ont basecalling model for medaka",
            pima_data.sub_process_verbosity,
            pima_data.sub_process_color,
        )

        determine_ont_model(pima_data)
        return
    
    pima_data.analysis.append(["medaka_ont_assembly", pima_data])

    # Assume this means we have an ONT assembly
    pima_data.will_have_ont_assembly = True

def medaka_ont_assembly(pima_data: PimaData):

    print_and_log(
        pima_data,
        "Running Medaka on ONT assembly",
        pima_data.main_process_verbosity,
        pima_data.main_process_color,
    )


    medaka_dir = os.path.join(pima_data.output_dir, "medaka")
    if find_checkpoint(pima_data, medaka_dir):
        print_and_log(
            pima_data,
            "Medaka had previously been run and finished successfully",
            pima_data.main_process_verbosity,
            pima_data.main_process_color,
        )
        pima_data.genome_fasta = os.path.join(medaka_dir, "assembly.fasta")
        pima_data.did_medaka_ont_assembly = True
        return

    os.makedirs(medaka_dir)
    make_start_file(pima_data, medaka_dir)

    # Actually run Medaka
    print_and_log(
        pima_data,
        f"Starting medaka using model: '{pima_data.ont_model}'",
        pima_data.sub_process_verbosity,
        pima_data.sub_process_color,
    )
    medaka_threads = min(pima_data.threads, 2)
    medaka_fasta = os.path.join(medaka_dir, "consensus.fasta")
    medaka_stdout, medaka_stderr = std_files(
        os.path.join(medaka_dir, "medaka")
    )
    command = " ".join(
        [
            "medaka_consensus",
            "-i",
            pima_data.ont_fastq,
            "-d",
            pima_data.genome_fasta,
            "-o",
            medaka_dir,
            "-t",  
            str(medaka_threads), # Medaka throttles anything more than 2 due to poor scaling
            "-b 50",  # MUCH more efficient on scicomp than the default [-b 100] (10x faster, 1/10th the RAM)
            "-m",
            pima_data.ont_model,
            "1>", medaka_stdout,
            "2>", medaka_stderr,
        ]
    )
    print_and_run(pima_data, command)
    validate_file_and_size_or_error(
        pima_data,
        medaka_fasta,
        "Medaka FASTA",
        "cannot be found after Medaka",
        "is empty",
    )

    medaka_bam = os.path.join(medaka_dir, "calls_to_draft.bam")

    print_and_log(
        pima_data,
        "Repairing contig names after Medaka",
        pima_data.sub_process_verbosity,
        pima_data.sub_process_color,
    )
    pima_data.genome_fasta = os.path.join(medaka_dir, "assembly.fasta")
    command = " ".join(
        [
            'awk \'{if($0 ~ /^>/){gsub(":.*", "", $0);gsub("_segment", "_", $0)}print}\'',
            medaka_fasta,
            ">",
            pima_data.genome_fasta,
        ]
    )
    print_and_run(pima_data, command)
    validate_file_and_size_or_error(
        pima_data,
        pima_data.genome_fasta,
        "Genome assembly",
        "cannot be found after fixing names",
        "is empty",
    )

    pima_data.load_genome()
    make_finish_file(pima_data, medaka_dir)
    pima_data.did_medaka_ont_assembly = True
    pima_data.files_to_clean.extend([medaka_bam, medaka_bam + ".bai"])
    pima_data.files_to_clean.append(os.path.join(medaka_dir, "consensus_probs.hdf"))

def determine_ont_model(pima_data: PimaData):
    #skip if this ran during the info_fastq assessment step ()
    if pima_data.ont_model != "auto":
        return
    
    command = " ".join(
        [
            "medaka tools resolve_model --auto_model consensus",
            pima_data.ont_fastq,
        ]
    )
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if result.returncode == 0:
        pima_data.ont_model = result.stdout.strip() #tell medaka to automatically determine the model
        pima_data.analysis.append(["medaka_ont_assembly", pima_data])
        print_and_log(
            pima_data,
            f"Identified basecalling model: {pima_data.ont_model}",
            pima_data.sub_process_verbosity,
            pima_data.sub_process_color,
        )
        
    elif result.returncode != 0:
        warn = "Medaka could not determine the basecalling model used to generate the fastq files and '--ont-model' was not provided, continuing PiMA without medaka polishing."
        add_warning(pima_data, warn)
        pima_data.assembly_notes = pd.concat([pima_data.assembly_notes, pd.Series(warn, dtype='object')])