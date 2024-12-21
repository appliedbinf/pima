#!/usr/bin/env python
from __future__ import annotations  # for "|" in docstrings python 3.7-3.9
import os
import sys
import matplotlib as mpl

from Pima.pima_data import PimaData
from Pima.utils.utils import print_and_log, stop_logging, clean_up
from Pima.utils.settings import Settings
from Pima.utils.cli import parse_args
from Pima import modules

mpl.use("Agg")

def run_prechecks(pima_data: PimaData, settings: Settings, pima_cmdtext: str):
    print_and_log(
        pima_data,
        "STARTING VALIDATION STEPS",
        pima_data.main_process_verbosity,
        pima_data.main_process_color,
    )
    modules.validate_download(pima_data, settings)
    modules.validate_output_dir(pima_data, settings, pima_cmdtext)
    modules.validate_organism(pima_data)
    modules.validate_multiplex_fastq(pima_data)

    if len(pima_data.errors) > 0:
        print_and_log(
            pima_data,
            "Errors were found during validation.",
            pima_data.fail_verbosity,
            pima_data.error_color,
        )
        for error in pima_data.errors:
            print_and_log(
                pima_data,
                error,
                pima_data.fail_verbosity,
                pima_data.error_color,
            )
        print_and_log(
            pima_data,
            "Aborting.",
            pima_data.fail_verbosity,
            pima_data.error_color,
        )
        sys.exit(1)

def run_validation(pima_data: PimaData, settings: Settings):
    modules.validate_ont_fastq(pima_data, settings)
    modules.validate_illumina_fastq(pima_data)
    modules.validate_genome_fasta(pima_data)
    modules.validate_genome_assembly_size(pima_data)
    modules.validate_contamination_check(pima_data, settings)
    modules.validate_assembler(pima_data)
    modules.validate_assembly_info(pima_data)
    modules.validate_medaka(pima_data)
    modules.validate_illumina_polish(pima_data)
    modules.validate_evaluate_assembly(pima_data)
    modules.validate_plasmids(pima_data, settings)
    modules.validate_features(pima_data, settings)
    modules.validate_blast(pima_data)
    modules.validate_reference_fasta(pima_data)
    modules.validate_quast(pima_data)
    modules.validate_mutations(pima_data)
    modules.validate_draw_features(pima_data)
    modules.validate_draw_amr_matrix(pima_data)
    modules.validate_draw_circos(pima_data, settings)
    modules.validate_make_report(pima_data, settings)

    if len(pima_data.errors) > 0:
        print_and_log(
            pima_data,
            "Errors were found during validation.",
            pima_data.fail_verbosity,
            pima_data.error_color,
        )
        for error in pima_data.errors:
            print_and_log(
                pima_data,
                error,
                pima_data.fail_verbosity,
                pima_data.error_color,
            )
        print_and_log(
            pima_data,
            "Aborting.",
            pima_data.fail_verbosity,
            pima_data.error_color,
        )
        sys.exit(1)

    #log all the versions
    version_log = "Utility versions:"
    for utility, version in pima_data.versions.items():
        version_log = version_log + "\n" + "{:<25} {:<10}".format(utility, version)
    print_and_log(
        pima_data,
        version_log,
        pima_data.main_process_verbosity,
        pima_data.main_process_color,
    )

def define_workflow(pima_data: PimaData):
    """
    Step through modules and run pima

    The run_validation steps parse the cmdline args, check all necessary tools / files are availabe for the
    requested steps in the pipeline, and queues up the modules by adding the steps to 
    pima_data.analysis object
    """
    print_and_log(
        pima_data,
        "STARTING PiMA Analysis",
        pima_data.main_process_verbosity,
        pima_data.main_process_color,
    )
    while (True):
        step = pima_data.analysis[0]
        pima_data.analysis = pima_data.analysis[1:]
        
        ## See if we have arguments to pass to our function
        if type(step) is list:
            arguments = []
            if len(step) > 1:
                arguments = step[1:]
            step = step[0]
            function = getattr(modules, step)
            function(*arguments)
        else:
            function = getattr(modules, step)
            function()

        if (len(pima_data.analysis) == 0):
            break

def run_workflow(pima_data: PimaData, settings: Settings):
    run_validation(pima_data, settings)
    define_workflow(pima_data)
    clean_up(pima_data)
    stop_logging(pima_data, "PiMA completed successfully")

def main():
    """ """

    settings = Settings()
    opts, unknown_args = parse_args(settings)

    # Start the analysis
    pima_data = PimaData(opts, unknown_args)

    # Capture commandline options used
    pima_data.run_command = ' '.join(sys.argv)
    run_prechecks(pima_data, settings, [f"PiMA command used: {' '.join(sys.argv)}"])

    ##Initialize serial multiplex analysis
    if pima_data.multiplexed:
        modules.initialize_multiplex_analysis(pima_data, settings)

    else:
        run_workflow(pima_data, settings)


if __name__ == "__main__":
    main()
