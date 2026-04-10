from __future__ import annotations
import os
import re
import copy
import datetime
import glob
import shutil
import csv

from collections import defaultdict
from dataclasses import dataclass

from Pima.pima_data import PimaData
from Pima.utils.settings import Settings
from Pima import modules
import Pima.pima

from Pima.utils.utils import (
    print_and_log,
    print_and_run,
    error_out,
    std_files,
    stop_logging,
    validate_utility,
    validate_file_and_size_or_error,
)

@dataclass
class barcode_data:
    barcode_id: str
    from_sample_sheet: bool = False
    #sample_sheet_specific vals
    illumina_r1: str | None = None
    illumina_r2: str | None = None
    genome: str | None = None
    genome_size: str | None = None
    reference_organism: str | None = None
    reference_genome: str | None = None
    reference_mutation_bed_file: str | None = None
    self_circos: str | None = None
    
    #multiplex decoding vals
    ont_fastq: str | None = None
    barcode_root_path: str | list | None = None
    barcode_fastq_list: list = None
    barcode_fastq_paths: list = None
    barcode_size_bytes: int = None
    barcode_size_bp: int | None = None

    def update_barcode(self, barcode_root_path: str, barcode_fastq_list: list, barcode_fastq_paths: list, barcode_size_bytes: int):
        self.barcode_root_path = [self.barcode_root_path, barcode_root_path]
        self.barcode_fastq_list = [*self.barcode_fastq_list, *barcode_fastq_list]
        self.barcode_fastq_paths = [*self.barcode_fastq_paths, *barcode_fastq_paths]
        self.barcode_size_bytes = self.barcode_size_bytes + barcode_size_bytes

    def create_concat_fastq(self, pima_data, fastq_path:str = None):
        if self.from_sample_sheet:
            return

        if len(self.barcode_fastq_list) == 1:
            if fastq_path:
                #for nf multiplexing
                if re.search(r'\.(gz|gzip)$', self.barcode_fastq_list[0]):
                    fastq_path = fastq_path + ".gz"
                os.symlink(self.barcode_fastq_paths[0], fastq_path)
                self.ont_fastq = fastq_path

            else:
                self.ont_fastq = self.barcode_fastq_paths[0]
            return
        
        print_and_log(
                pima_data,
                "Concatenating barcode fastq files",
                pima_data.sub_process_verbosity,
                pima_data.sub_process_color,
        )

        if re.search(r'\.(gz|gzip)$', self.barcode_fastq_list[0]):
            if fastq_path:
                self.ont_fastq = fastq_path + ".gz"
                pima_data.ont_fastq = self.ont_fastq
            else:
                self.ont_fastq = pima_data.ont_fastq + ".gz"
                pima_data.ont_fastq = self.ont_fastq
        command = " ".join(
            [
                "cat",
                " ".join(self.barcode_fastq_paths),
                f"> {self.ont_fastq}",
            ]
        )
        print_and_run(pima_data, command)

    def report_sample(self):
        message = f"Running PiMA on {self.barcode_id}"
        return message


def validate_multiplex_fastq(pima_data: PimaData):
    if not (pima_data.multiplexed or pima_data.sample_sheet):
        return
    
    if pima_data.nextflow and not (pima_data.multiplexed or pima_data.sample_sheet):
        error = "Nextflow is only useful for multiplexed data"
        error_out(
            pima_data,
            error,
        )

    if pima_data.sample_sheet and any([pima_data.ont_fastq, pima_data.illumina_fastq]):
        error_out(
            pima_data,
            "No other data inputs allowed when providing a sample_sheet. See an example by running 'pima --example-sample-sheet'",
        )

    if pima_data.resume and pima_data.nextflow:
        error_out(
            pima_data,
            "--resume does not currently work with nextflow multiplexing. If the assemblies were completed in the previous attempt, you can resume the multiplex run in serial by removing '--nextflow' or resume each sample independently without '--multiplexed', otherwise just use '--overwrite'",
        )
    if not any([pima_data.ont_fastq, pima_data.sample_sheet]):
        error_out(
            pima_data,
            "--multiplexed requires that a directory of FASTQ files or directories of FASTQ files be given",
        )

    if pima_data.sample_sheet:
        validate_file_and_size_or_error(pima_data, pima_data.sample_sheet)
        pima_data.multiplexed = True
        return
         
    # Conditions for when not using a sample_sheet
    if pima_data.illumina_fastq:
        error_out(
            pima_data,
            "--multiplexing with illumina data only works by providing a --sample-sheet as input. Exiting",
        )

    if not os.path.exists(pima_data.ont_fastq):
        error_out(
            pima_data,
            f"The provided ont-fastq path: {pima_data.ont_fastq} does not exist. Please check your input. Exiting."
        )


    if pima_data.barcode_min_fraction >= 1:
        error_out(
            pima_data,
            f"--barcode_min_fraction is greater than 1, did you mean to use {pima_data.barcode_min_fraction / 100}?",
        ) 

    if pima_data.genome_assembly_size is not None and pima_data.genome_assembly_size != "estimate":
        print_and_log(
            pima_data,
            f"Using the same --genome_size {pima_data.genome_assembly_size} for every sample in the multiplex run. If you do not expect all samples to have the same genome size (+/-10%), please cancel (ctrl+c) and re-run using '--genome-size estimate' or leave it blank (prevents downsampling)'",
            pima_data.warning_verbosity,
            pima_data.warning_color,
        )

    if os.path.isfile(pima_data.ont_fastq) and pima_data.barcode_kit:
        #fastq data has not been demultiplexed - try to demulitplex it ?
        ## currently will just exit if user gives a single fastq file
        
        message = ("You provided a single fastq file and indicated it is multiplexed. "
                    "PiMA currently doesn't demultiplex a fastq file since this file type is not common. "
                    "If you need to demultiplex, we recommend using dorado. "
                    "Please let us know if this is a feature you'd like to see added."
                    "Exiting now"
        )
        error_out(
            pima_data,
            message,
        )

    pima_data.ont_fastq = os.path.realpath(pima_data.ont_fastq)
    pima_data.output_dir = os.path.realpath(pima_data.output_dir)


def identify_multiplexed_fastq_files(pima_data: PimaData):

    print_and_log(
        pima_data,
        "Starting Multiplex Analysis",
        pima_data.main_process_verbosity,
        pima_data.main_process_color,
    )
    
    multiplexed_dirs = defaultdict()
    total_dir_size = 0

    if not any(os.path.isdir(os.path.join(pima_data.ont_fastq, item)) for item in os.listdir(pima_data.ont_fastq)) and any([re.search(r"fastq", item) for item in os.listdir(pima_data.ont_fastq)]):
        #user is providing a directory containing fastq files for each sample
        for fastq in os.listdir(pima_data.ont_fastq):
            name = os.path.splitext(os.path.basename(fastq))[0]

            if re.search(r"\.(gz|gzip)$", fastq):
                name = os.path.splitext(os.path.splitext(fastq)[0])[0]

            fastq_size = os.path.getsize(os.path.join(pima_data.ont_fastq, fastq))
            multiplexed_dirs[name] = barcode_data(
                        barcode_id = name,
                        barcode_root_path = pima_data.ont_fastq,
                        barcode_fastq_list = [fastq],
                        barcode_fastq_paths = [os.path.join(pima_data.ont_fastq, fastq)],
                        barcode_size_bytes = fastq_size,
                        ont_fastq = os.path.join(os.path.join(pima_data.output_dir, name), f"{name}.fastq"))
            total_dir_size = total_dir_size + fastq_size

    else: 
        for root, dirs, files in os.walk(pima_data.ont_fastq):
            if any([re.search(r"fail", dir) for dir in dirs]):
                print_and_log(
                    pima_data,
                    "There are directories with the string 'fail' in the name. We will use both passing and failing reads for each sample. If you wish to use only passing reads, please cancel (ctrl+c) and re-run giving just the fastq_pass directory as input",
                    pima_data.warning_verbosity,
                    pima_data.warning_color,
                )

            if len(files) > 0:
                if re.search(r"fastq",files[0]) and not re.search(r"unclassified", root):
                    if os.path.basename(root) in multiplexed_dirs:
                        dir_size = sum(os.path.getsize(os.path.join(root, f)) for f in files)
                        multiplexed_dirs[os.path.basename(root)].update_barcode(root, files, [os.path.join(root, f) for f in files], dir_size)
                        total_dir_size = total_dir_size + dir_size

                    else:
                        dir_size = sum(os.path.getsize(os.path.join(root, f)) for f in files)
                        multiplexed_dirs[os.path.basename(root)] = barcode_data(
                                                    barcode_id = os.path.basename(root),
                                                    barcode_root_path = root,
                                                    barcode_fastq_list = files,
                                                    barcode_fastq_paths = [os.path.join(root, f) for f in files],
                                                    barcode_size_bytes = dir_size,
                                                    ont_fastq = os.path.join(pima_data.output_dir, f"{os.path.basename(root)}.fastq"))
                        total_dir_size = total_dir_size + dir_size

    ignored_barcodes = dict()
    for barcode in multiplexed_dirs.copy().values():
        perc_data = barcode.barcode_size_bytes / total_dir_size
        if perc_data < pima_data.barcode_min_fraction:
            ignored_barcodes[barcode.barcode_id] = perc_data
            del multiplexed_dirs[barcode.barcode_id]

    print_and_log(
            pima_data,
            f"Running PiMA on barcodes: {', '.join(multiplexed_dirs.keys())}",
            pima_data.main_process_verbosity,
            pima_data.main_process_color,
    )

    if len(ignored_barcodes) > 0:
        message = (
            "The following barcodes were found in the input directory but were NOT analyzed "
            f"because they contained less than {pima_data.barcode_min_fraction*100}% (default=0.025 [2.5%]) of the fastq data:\n"
            "If you need to change the min_fraction, please re-run pima with the following flag '--barcode_min_fraction <fractional value>'\n"
        )
        for k, v in ignored_barcodes.items():
            message = message + "{:<15} {:>.1%}".format(k, v) + "\n"
        print_and_log(
            pima_data,
            message,
            pima_data.warning_verbosity,
            pima_data.warning_color,
        )
    pima_data.barcodes = multiplexed_dirs


def parse_input_sample_sheet(pima_data: PimaData):
    """
    Probably need to use the Barcode class similar to identify_multiplexed_fastq_files fx
    Loop over each row in the sample sheet & set variables
    Return pima_data.barcodes [ dict of barcode classes ]

    May need to update barcode class to allow for individual parameters (e.g. illumina reads, reference, circos type...)
    These need to be optionally handed over to nf_template OR to the serial pima command
    """
    validate_file_and_size_or_error(pima_data, pima_data.sample_sheet, min_size = 0)

    #all allowed columns
    mandatory_cols = ["sample_name"] # This gets saved in the data class as "barcode_id"
    optional_cols = [
        "ont_fastq", 
        "illumina_r1", 
        "illumina_r2", 
        "genome", 
        "genome_size", 
        "reference_organism", 
        "reference_genome", 
        "reference_mutation_bed_file",
        "self_circos",
        ]

    samples = defaultdict()
    with open(pima_data.sample_sheet, 'r') as f:
        reader = csv.DictReader(f, delimiter = '\t')
        header_columns = reader.fieldnames

        if not set(mandatory_cols).issubset(set(header_columns)):
            missing_cols = set(mandatory_cols) - set(header_columns)
            message = (f"Your sample sheet is missing the required column {missing_cols}")
            error_out(
                pima_data,
                message
            )

        if not any(data_col in header_columns for data_col in ["ont_fastq", "illumina_r1", "genome"]):
            message = (f"PiMA needs data to run, please include at least an ont_fastq, illumina_r1, or assembled genome as input")
            error_out(
                pima_data,
                message
            )

        if all(val in header_columns for val in ["reference_organism", "reference_genome"]) or all(val in header_columns for val in ["reference_organism", "reference_mutation_bed_file"]):
            message = (f"Please provide either a reference organism ('e.g. Bacillus_anthracis') OR a reference_genome and/or reference_mutation_bed_file. The sample sheet cannot have both")
            error_out(
                pima_data,
                message
            )       
    
        for sample in reader:
            sample_data = {}
            for col in mandatory_cols:
                sample_data['barcode_id'] = sample[col]
            for col in optional_cols:
                sample_data[col] = sample.get(col, None)
            
            try:
                cur_samp = barcode_data(**sample_data)
                cur_samp.from_sample_sheet = True

            except TypeError as e:
                message = (
                    f"""Failed to parse input sample sheet, unexpected column name: {str(e).split("'")[1]}\n"""
                    "Please correct the sample sheet and try again. "
                    f"The allowed column names are:\n{', '.join([x for x in mandatory_cols and optional_cols])}"
                )
                error_out(
                    pima_data,
                    message
                )        
            samples[cur_samp.barcode_id] = cur_samp

        pima_data.barcodes = samples


def initialize_multiplex_analysis(pima_data: PimaData, settings: Settings):
    
    if pima_data.sample_sheet:
        parse_input_sample_sheet(pima_data)

    else:
        identify_multiplexed_fastq_files(pima_data)

    if pima_data.nextflow or isinstance(pima_data.nextflow, str):
        validate_nextflow(pima_data)

        #This is pretty hacky, but on rosalind I can't get it to pass the modules / environment through correctly
        ## Here we basically use the pima run environment to feed the nextflow.config script user-specific values
        
        #generate nextflow config with the users conda environment that pima ran in
        try:
            conda_env = os.environ['CONDA_PREFIX']
        except KeyError:
            conda_env = "None"

        #use the activate script that was used in the parent environment
        try:
            activate_sh = os.environ['CONDA_EXE'].replace(r"/conda", "/activate")
        except KeyError:
            activate_sh = "None"

        nextflow_dir = os.path.join(settings.pima_path, "nextflow_parallelization")
        nextflow_config_template = os.path.join(nextflow_dir, "nextflow.config.template")
        user_nextflow_config = os.path.join(nextflow_dir, "nextflow.config")
        find_replace = {
            "conda = None": f"conda = '{conda_env}'",
            "beforeScript = None": f"beforeScript = 'source {activate_sh}'"
        }
        with open(nextflow_config_template, "rt") as fin:
            with open(user_nextflow_config, "wt") as fout:
                for line in fin:
                    for key in find_replace:
                        if key in line:
                            line = line.replace(key, find_replace[key])
                    fout.write(line)

        if isinstance(pima_data.nextflow, str):
            nextflow_args = pima_data.nextflow.replace("'","").replace('"','')
        else:
            nextflow_args = ""
        
        print_and_log(
            pima_data,
            "Handing off multiplexing to Nextflow",
            pima_data.main_process_verbosity,
            pima_data.main_process_color,
        )
        
        stop_logging(pima_data, "Sample specific logs are found in their respective directories, closing multiplex log now.")

        nf_file = os.path.join(pima_data.output_dir, "nf_singleplex_inputs.csv")
        header = [
            "sample_name",
            "output_directory",
            "ont_fastq",
            "illumina_fastq",
            "genome",
            "genome_size",
            "reference_organism",
            "reference_genome",
            "reference_mutation_bed_file",
            "self_circos",
            'original_command',
        ]
        with open(nf_file, "w") as nf_handle:
            nf_handle.write(",".join(header) + '\n')
            for barcode in pima_data.barcodes.keys():
                pima_data.barcodes[barcode].create_concat_fastq(pima_data, pima_data.barcodes[barcode].ont_fastq)
                updated_cmd = strip_pima_cmd(pima_data, pima_data.run_command)
                illumina_list = [fastq for fastq in [pima_data.barcodes[barcode].illumina_r1, pima_data.barcodes[barcode].illumina_r2] if fastq]
                
                line = (
                    f"{barcode},"
                    f"{os.path.join(pima_data.output_dir, barcode)},"
                    f"{pima_data.barcodes[barcode].ont_fastq or ''},"
                    f"{' '.join(illumina_list) if len(illumina_list) > 0 else ''},"
                    f"{pima_data.barcodes[barcode].genome or ''},"
                    f"{pima_data.barcodes[barcode].genome_size or pima_data.genome_assembly_size},"
                    f"{pima_data.barcodes[barcode].reference_organism or ''},"
                    f"{pima_data.barcodes[barcode].reference_genome or ''},"
                    f"{pima_data.barcodes[barcode].reference_mutation_bed_file or ''}," 
                    f"{pima_data.barcodes[barcode].self_circos.lower() if pima_data.barcodes[barcode].self_circos else ''},"                     
                    f"{updated_cmd}\n"
                )
                nf_handle.write(line)

        #run nextflow
        #specify nextflow workdir to be within the pima output dir so we can delete everything after a successful run
        nextflow_stdout, nextflow_stderr = std_files(f"{pima_data.output_dir}/nextflow")
        command = " ".join(
            [
                "nextflow run",
                os.path.join(settings.pima_path, "nextflow_parallelization/main.nf"),
                "--sample_sheet",
                nf_file,
                "--output",
                pima_data.output_dir,
                "-w",
                os.path.join(pima_data.output_dir, "work"),
                nextflow_args,
                "1>",
                nextflow_stdout,
                "2>",
                nextflow_stderr,
            ]
        )
        print_and_run(pima_data, command,  change_exe_dir = pima_data.output_dir)
        cleanup_nextflow(pima_data)

    #not using nextflow, running pima in serial
    else:
        stop_logging(pima_data, "Sample specific logs are found in their respective directories, closing multiplex log now.")
        for barcode in pima_data.barcodes.keys():
            barcode_pima_data = copy.deepcopy(pima_data)
            barcode_pima_settings = copy.deepcopy(settings)
            barcode_pima_data.multiplexed = None
            barcode_pima_data.sample_sheet = None

            barcode_pima_data.analysis_name = barcode
            barcode_pima_data.output_dir = os.path.join(barcode_pima_data.output_dir, barcode)
            barcode_pima_data.ont_fastq = pima_data.barcodes[barcode].ont_fastq
            barcode_pima_data.ont_raw_fastq = barcode_pima_data.ont_fastq
            barcode_pima_data.illumina_fastq = [fastq for fastq in [pima_data.barcodes[barcode].illumina_r1, pima_data.barcodes[barcode].illumina_r2] if fastq is not None]
            barcode_pima_data.genome_fasta = pima_data.barcodes[barcode].genome
            barcode_pima_data.genome_assembly_size = pima_data.barcodes[barcode].genome_size or pima_data.genome_assembly_size
            barcode_pima_data.organism = pima_data.barcodes[barcode].reference_organism or pima_data.organism
            barcode_pima_data.reference_fasta = pima_data.barcodes[barcode].reference_genome or pima_data.reference_fasta
            barcode_pima_data.mutation_region_bed = pima_data.barcodes[barcode].reference_mutation_bed_file or pima_data.mutation_region_bed
            barcode_pima_data.logging_file = os.path.join(barcode_pima_data.output_dir, "pima.log")

            if pima_data.barcodes[barcode].self_circos is not None and pima_data.barcodes[barcode].self_circos.lower() == "yes":
                barcode_pima_data.self_circos = True

            #Tell user what options are passed to the singleplex command
            updated_cmd = strip_pima_cmd(pima_data, pima_data.run_command)
            settings_from_multiplex = ' '.join([arg for arg in updated_cmd.split(' ') if not re.search("pima",arg)])
            pima_options = (
                    f"output_directory: {barcode_pima_data.output_dir}\n"
                    f"ont-fastq: {barcode_pima_data.ont_fastq or ''}\n"
                    f"illumina-fastq: {barcode_pima_data.illumina_fastq if len(barcode_pima_data.illumina_fastq) > 0 else ''}\n"
                    f"genome: {barcode_pima_data.genome_fasta or ''}\n"
                    f"genome-size: {barcode_pima_data.genome_assembly_size}\n"
                    f"organism: {barcode_pima_data.organism or ''}\n"
                    f"reference-genome: {barcode_pima_data.reference_fasta or ''}\n"
                    f"mutation-regions: {barcode_pima_data.mutation_region_bed or ''}\n" 
                    f"self-circos: {barcode_pima_data.self_circos or 'No'}\n"                     
                    f"flags from the original multiplex command: {settings_from_multiplex}"
                )

            log_messages = [
                ("main", f'[{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}]', f"{barcode_pima_data.barcodes[barcode].report_sample()}"),
                ("main", f'[{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}]', f"Sample inputs to pima:\n{pima_options}"),
                ]
            
            modules.validate_output_dir(barcode_pima_data, settings, log_messages)
            barcode_pima_data.barcodes[barcode].create_concat_fastq(barcode_pima_data)
            Pima.pima.run_workflow(barcode_pima_data, settings)


def strip_pima_cmd(pima_data, system_args: list):
    """
    Regenerate the pima command except singleplex with the specific ONT reads, updated output directory, removed multiplexing/nextflow statements
    """

    params_to_change = ['--output', '--ont-fastq', '--threads', '--sample-sheet']
    params_to_remove = ['--multiplexed', '--nextflow']
    if isinstance(pima_data.nextflow, str):
        system_args = re.sub(pima_data.nextflow, "", system_args)
    params_to_fix_path = ['--reference-genome', '--mutation-regions']
    param_iter = iter(system_args.split(" "))
    new_cmd = []
    for param in param_iter:
        if not param in [*params_to_change, *params_to_remove, *params_to_fix_path]:
            new_cmd.append(param)
        elif param in params_to_change:
            next(param_iter)
        elif param in params_to_fix_path:
            new_cmd.append(param)
            new_cmd.append(os.path.realpath(next(param_iter)))
        elif param in params_to_remove:
            continue
    return " ".join(i for i in new_cmd)


def validate_nextflow(pima_data):
    if not pima_data.nextflow:
        return
    if pima_data.no_assembly:
        return
      
    if validate_utility(
        pima_data, "nextflow", "nextflow is not on the PATH (required by --multiplexed --nextflow)."
    ):
        command = "nextflow -v"
        pima_data.versions["nextflow"] = re.search(
            r"[0-9]+\.[0-9.]+", print_and_run(pima_data, command)[0]
        ).group(0)


def cleanup_nextflow(pima_data):
    nf_temp_files = glob.glob(os.path.join(pima_data.output_dir, ".nextflow*"))
    nf_work_dir = os.path.join(pima_data.output_dir, "work")
    for file in nf_temp_files:
        try: 
            shutil.rmtree(file)
        except NotADirectoryError:
            os.remove(file)

    shutil.rmtree(nf_work_dir)