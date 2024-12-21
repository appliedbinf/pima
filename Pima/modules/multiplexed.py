from __future__ import annotations
import os
import re
import copy
import datetime
import glob
import shutil

from collections import defaultdict

from Pima.pima_data import PimaData
from Pima.utils.settings import Settings
from Pima import modules
import Pima.pima

from Pima.utils.utils import (
    print_and_log,
    print_and_run,
    std_files,
    stop_logging,
    validate_utility,
)

class barcode_data:
    def __init__(self, barcode_id: str, barcode_root_path: str | list, barcode_fastq_list: list, barcode_fastq_paths: list, barcode_size_bytes: int):
        self.barcode_id = barcode_id
        self.barcode_root_path = barcode_root_path
        self.barcode_fastq_list = barcode_fastq_list
        self.barcode_fastq_paths = barcode_fastq_paths
        self.barcode_size_bytes = barcode_size_bytes
        self.barcode_size_bp = None

    def update_barcode(self, barcode_root_path: str, barcode_fastq_list: list, barcode_fastq_paths: list, barcode_size_bytes: int):
        self.barcode_root_path = [self.barcode_root_path, barcode_root_path]
        self.barcode_fastq_list = [*self.barcode_fastq_list, *barcode_fastq_list]
        self.barcode_fastq_paths = [*self.barcode_fastq_paths, *barcode_fastq_paths]
        self.barcode_size_bytes = self.barcode_size_bytes + barcode_size_bytes

    def create_concat_fastq(self, pima_data, fastq_path:str = None):
        if len(self.barcode_fastq_list) == 1:
            if fastq_path:
                #for nf multiplexing
                if re.search(r'\.(gz|gzip)$', self.barcode_fastq_list[0]):
                    fastq_path = fastq_path + ".gz"
                os.symlink(self.barcode_fastq_paths[0], fastq_path)
                pima_data.ont_fastq = fastq_path

            else:
                pima_data.ont_fastq = self.barcode_fastq_paths[0]
            return
        
        print_and_log(
                pima_data,
                "Concatenating barcode fastq files",
                pima_data.sub_process_verbosity,
                pima_data.sub_process_color,
        )

        if re.search(r'\.(gz|gzip)$', self.barcode_fastq_list[0]):
            pima_data.ont_fastq = pima_data.ont_fastq + ".gz"

        command = " ".join(
            [
                "cat",
                " ".join(self.barcode_fastq_paths),
                f"> {pima_data.ont_fastq}",
            ]
        )
        print_and_run(pima_data, command)

    def report_multiplex_sample(self):
        message = f"Running PiMA on {self.barcode_id}"
        return message


def validate_multiplex_fastq(pima_data: PimaData):
    if not pima_data.multiplexed:
        return
    
    if pima_data.resume and pima_data.nextflow:
        pima_data.errors.append("--resume does not currently work with nextflow multiplexing. If the assemblies were completed in the previous attempt, you can resume the multiplex run in serial by removing '--nextflow' or resume each sample independently without '--multiplexed', otherwise just use '--overwrite'")
        
    if not pima_data.ont_fastq:
        pima_data.errors.append("--multiplexed requires that a directory of FASTQ files or directories of FASTQ files be given")

    if pima_data.illumina_fastq:
        pima_data.errors.append("--multiplexing does not currently work with illumina data. Exiting")

    if pima_data.barcode_min_fraction >= 1:
        pima_data.errors.append(f"--barcode_min_fraction is greater than 1, did you mean to use {pima_data.barcode_min_fraction / 100}?") 

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
        pima_data.errors.append(message)

    pima_data.will_have_ont_fastq = True
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
            #user is providing a directory contains fastq files for each sample (hopefully)
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
                            barcode_size_bytes = fastq_size)
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
                                                        barcode_size_bytes = dir_size)
                            total_dir_size = total_dir_size + dir_size

        #ignored_barcodes = []
        ignored_barcodes = dict()
        for barcode in multiplexed_dirs.copy().values():
            perc_data = barcode.barcode_size_bytes / total_dir_size
            if perc_data < pima_data.barcode_min_fraction:
                #ignored_barcodes.append(barcode.barcode_id)
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


def initialize_multiplex_analysis(pima_data: PimaData, settings: Settings):
    
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

        nf_file = os.path.join(pima_data.output_dir, "nf_singplex_inputs.csv")
        
        with open(nf_file, "w") as nf_handle:
            for barcode in pima_data.barcodes.keys():
                barcode_pima_data = copy.deepcopy(pima_data)
                barcode_pima_data.output_dir = os.path.join(pima_data.output_dir, barcode)
                barcode_pima_data.ont_fastq = os.path.join(pima_data.output_dir, f"{barcode}.fastq")
                barcode_pima_data.barcodes[barcode].create_concat_fastq(barcode_pima_data, barcode_pima_data.ont_fastq)
                updated_cmd = strip_pima_cmd(pima_data, pima_data.run_command)
                line = (
                    f"{barcode},"
                    f"{barcode_pima_data.output_dir},"
                    f"{barcode_pima_data.ont_fastq},"
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
            barcode_pima_data.output_dir = os.path.join(barcode_pima_data.output_dir, barcode)
            barcode_pima_data.logging_file = os.path.join(barcode_pima_data.output_dir, "pima.log")
            barcode_pima_data.ont_fastq = os.path.join(barcode_pima_data.output_dir, f"{barcode}.fastq")
            barcode_pima_data.multiplexed = None
            log_message = [("main", f'[{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}]', f"{barcode_pima_data.barcodes[barcode].report_multiplex_sample()}")]
            modules.validate_output_dir(barcode_pima_data, settings, log_message)
            barcode_pima_data.barcodes[barcode].create_concat_fastq(barcode_pima_data)
            Pima.pima.run_workflow(barcode_pima_data, settings)


def strip_pima_cmd(pima_data, system_args: list):
    """
    Regenerate the pima command except singleplex with the specific ONT reads, updated output directory, removed multiplexing/nextflow statements
    """
    #TODO: How should Illumina data be handled for a multiplexed run??? 
    #      - Just need to add a sample_sheet mode

    params_to_change = ['--output', '--ont-fastq', '--threads']
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
    if not pima_data.will_have_ont_fastq:
        error = "Nextflow is only setup for using ont-fastq data, not assemblies"
        pima_data.errors.append(error)
        print_and_log(
            pima_data,
            error,
            pima_data.fail_verbosity,
            pima_data.error_color,
        )
        return
    
    if pima_data.nextflow and not pima_data.multiplexed:
        error = "Nextflow is only useful for multiplexed data"
        pima_data.errors.append(error)
        print_and_log(
            pima_data,
            error,
            pima_data.fail_verbosity,
            pima_data.error_color,
        )
    
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