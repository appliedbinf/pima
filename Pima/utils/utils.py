import sys
import os
import re
import datetime
import shutil
import subprocess
from pathlib import Path

import numpy as np

from Pima.pima_data import PimaData
from Pima.pima_colors import Colors


def nicenumber(x: float, round: int):
    exp = np.floor(np.log10(x))
    f = x / 10**exp

    if round:
        if f < 1.5:
            nf = 1.0
        elif f < 3.0:
            nf = 2.0
        elif f < 7.0:
            nf = 5.0
        else:
            nf = 10.0
    else:
        if f <= 1.0:
            nf = 1.0
        elif f <= 2.0:
            nf = 2.0
        elif f <= 5.0:
            nf = 5.0
        else:
            nf = 10.0

    return nf * 10.0**exp


def pretty(low, high, n):
    range = nicenumber(high - low, False)
    d = nicenumber(range / (n - 1), True)
    miny = np.floor(low / d) * d
    maxy = np.ceil(high / d) * d
    return np.arange(miny, maxy + 0.5 * d, d)


def format_kmg(number: float, decimals: int = 0):
    if number == 0:
        return "0"

    magnitude_powers = [10**9, 10**6, 10**3, 1]
    magnitude_units = ["G", "M", "K", ""]
    for i in range(len(magnitude_units)):
        if number >= magnitude_powers[i]:
            magnitude_power = magnitude_powers[i]
            magnitude_unit = magnitude_units[i]
            return f"{round(number/magnitude_power,decimals)}{magnitude_unit}"


def print_and_log(
    pima_data: PimaData, text: str, verbosity: int, color: str = Colors.ENDC, time_string:str = None
    ):
    if not time_string:
        time_string = f'[{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}]'
    if verbosity <= pima_data.verbosity:
        sys.stderr.write(f"{time_string} {color}{text}{Colors.ENDC}\n")

    if pima_data.logging_handle:
        pima_data.logging_handle.write(f"{time_string} {text}\n")
        #otherwise log doesn't get created until it is closed
        ## I don't think this is too expensive of a process since we are not logging that frequently
        pima_data.logging_handle.flush()


def start_logging(pima_data: PimaData):

    if pima_data.resume:
        # if user has specified a new logfile for us
        if pima_data.logging_file and not os.path.isfile(pima_data.logging_file):
            pima_data.logging_handle = open(pima_data.logging_file, 'w')

        # if user has resumed a run using the same log path
        elif pima_data.logging_file and os.path.isfile(pima_data.logging_file):
            backup_log_file(pima_data.logging_file)
            pima_data.logging_handle = open(pima_data.logging_file, 'w')

        # if user is using the default pima.log file
        else:
            backup_log_file(os.path.join(pima_data.output_dir, "pima.log"))
            pima_data.logging_file = os.path.join(pima_data.output_dir, "pima.log")
            pima_data.logging_handle = open(pima_data.logging_file, 'w')

    else:
        if pima_data.logging_file:
            pima_data.logging_handle = open(pima_data.logging_file, 'w')

        else:
            pima_data.logging_file = os.path.join(pima_data.output_dir, "pima.log")
            pima_data.logging_handle = open(pima_data.logging_file, 'w')


def backup_log_file(existing_log_path: str):
    log_path = Path(existing_log_path)

    if os.path.isfile(log_path):
        rename_log_path = os.path.join(log_path.parent, f"previous_{log_path.name}")
        
        if os.path.isfile(rename_log_path):
            with open(rename_log_path, 'a') as prev_log_file:
                prev_log_file.write("\n\n")
                with open(existing_log_path) as cur_log_file:
                    shutil.copyfileobj(cur_log_file, prev_log_file)
            os.remove(existing_log_path)
            return

        os.rename(existing_log_path, rename_log_path)


def stop_logging(pima_data: PimaData, message: str = None):

    if message:
        print_and_log(
            pima_data,
            message,
            pima_data.main_process_verbosity,
            pima_data.main_process_color,
        )
        
    pima_data.logging_handle.close()
    pima_data.logging_handle = None


def validate_file(the_file: str):
    return os.path.isfile(the_file)


def validate_file_size(pima_data: PimaData, the_file: str, min_size: int = 0):
    if not pima_data.fake_run:
        return os.stat(the_file).st_size >= min_size
    else:
        return True


def validate_file_and_size(pima_data: PimaData, the_file: str, min_size: int = 0):
    return validate_file(the_file) and validate_file_size(pima_data, the_file, min_size)


def validate_file_and_size_or_error(
    pima_data,
    the_file,
    error_prefix="The file",
    presence_suffix="doesn't exist",
    size_suffix="is size 0",
    min_size=0,
):
    if not validate_file(the_file) and not pima_data.fake_run:
        error_out(pima_data, " ".join([error_prefix, the_file, presence_suffix]))

    if not validate_file_size(pima_data, the_file, min_size) and not pima_data.fake_run:
        error_out(pima_data, " ".join([error_prefix, the_file, size_suffix]))


def validate_utility(pima_data: PimaData, utility: str, error: str):
    if not shutil.which(utility):
        pima_data.errors.append(error)
        print_and_log(
            pima_data,
            error,
            pima_data.fail_verbosity,
            pima_data.error_color,
        )
        return False
    else:
        return True


def print_and_run(pima_data: PimaData, command: str, change_exe_dir: str = None):
    print_and_log(pima_data, command, pima_data.command_verbosity)
    return run_command(pima_data, command, change_exe_dir)


def run_command(pima_data: PimaData, command: str, change_exe_dir: str = None):
    if not pima_data.fake_run:
        if change_exe_dir:
            result = subprocess.run(command, shell=True, capture_output=True, text=True, cwd = change_exe_dir)
        else:
            result = subprocess.run(command, shell=True, capture_output=True, text=True)

        if result.returncode == 0:
            return result.stdout.split("\n")
        elif re.search(r"\.stderr$", command):
            sterr_f = [x for x in command.split(" ") if re.search(r"\.stderr$", x)][0]
            message = f"Command {command} failed with the following error. exiting\n{Path(sterr_f).read_text()}"
            error_out(pima_data, message)
        else:
            message = f"Command {command} failed with the following error. exiting\n{result.stderr}"
            error_out(pima_data, message)


def error_out(pima_data: PimaData, message: str):
    print_and_log(
        pima_data,
        message,
        pima_data.fail_verbosity,
        pima_data.error_color,
    )
    sys.exit(1)


def print_warning(pima_data: PimaData, warning: str):
    print_and_log(
        pima_data, warning, pima_data.warning_verbosity, pima_data.warning_color
    )


def add_warning(pima_data: PimaData, warning: str):
    print_warning(pima_data, warning)
    pima_data.warnings += [warning]


def find_checkpoint(pima_data: PimaData, dir: str):
    """Searches for the .finish file generated after each analysis step is completed

    Args:
        dir (path): Path to analysis directory

    Returns:
        True if 'resume' flag provided & ".finish" found
        False if ".finish" is not found, existing directory deleted
    """
    if not pima_data.resume:
        return False

    if os.path.exists(os.path.join(dir, ".finish")):
        return True

    else:
        if os.path.exists(dir):
            shutil.rmtree(dir)
        return False


def std_files(prefix: str):
    return [prefix + ".std" + i for i in ["out", "err"]]


def touch_file(pima_data: PimaData, a_file: str):
    command = " ".join(["touch", a_file])
    print_and_run(pima_data, command)


def make_start_file(pima_data: PimaData, a_dir: str):
    start_file = os.path.join(a_dir, ".start")
    touch_file(pima_data, start_file)


def make_finish_file(pima_data: PimaData, a_dir: str):
    finish_file = os.path.join(a_dir, ".finish")
    touch_file(pima_data, finish_file)


def make_report_info_file(pima_data: PimaData, a_dir: str):
    report_info_file = os.path.join(a_dir, ".report_info")
    touch_file(pima_data, report_info_file)
    return report_info_file


def clean_up(pima_data: PimaData):

    print_and_log(
        pima_data,
        "Cleaning up",
        pima_data.main_process_verbosity,
        pima_data.main_process_color,
    )

    if (pima_data.genome_fasta):
        final_fasta = os.path.join(pima_data.output_dir, 'assembly.fasta')
        command = ' '.join(['cp', pima_data.genome_fasta, final_fasta])
        print_and_run(pima_data, command)

    #create a shortcut to the report in the primary output_dir
    if os.path.isfile(os.path.join(pima_data.output_dir, "report.pdf")):
        os.remove(os.path.join(pima_data.output_dir, "report.pdf"))
    
    os.symlink(os.path.join(pima_data.output_dir, "report", "report.pdf"), os.path.join(pima_data.output_dir, "report.pdf"))

    if not pima_data.keep_intermediates:
        if len(pima_data.files_to_clean) > 0:
            for file in pima_data.files_to_clean :
                if os.path.isfile(file):
                    try:
                        os.remove(file)
                    except OSError:
                        continue