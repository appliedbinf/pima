import os
import shutil
import datetime
from pathlib import Path

from Pima.pima_data import PimaData
from Pima.utils.settings import Settings
from Pima.utils.utils import print_and_log, start_logging


def validate_output_dir(pima_data: PimaData, settings: Settings, log_messages:list=[]):

    log_messages.append(("main", f'[{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}]', f"PiMA version: {settings.pima_version}"))
    log_messages.append(("main", f'[{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}]', "Validating output dir"))

    if not pima_data.output_dir:
        pima_data.errors += ["No output directory given (--output)"]
    elif pima_data.overwrite and pima_data.resume:
        pima_data.errors += ["--overwrite and --resume are mutually exclusive"]
    elif os.path.exists(pima_data.output_dir) and not (
        pima_data.overwrite or pima_data.resume
    ):
        pima_data.errors += [
            "Output directory "
            + pima_data.output_dir
            + " already exists.  Add --overwrite OR --resume to ignore"
        ]
    
    else:
        pima_data.output_dir = os.path.realpath(pima_data.output_dir)
        make_outdir(pima_data, log_messages)


def make_outdir(pima_data: PimaData, log_messages: list):
    if pima_data.resume and os.path.isdir(pima_data.output_dir):
        start_logging(pima_data)
        log_messages.append(("main", f'[{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}]', f"Resuming from previous run, previous log has been renamed to 'previous_<logname>.log.'"))
        report_logs(pima_data, log_messages)
        return

    if os.path.isdir(pima_data.output_dir):
        log_messages.append(("warn", f'[{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}]', f"Output directory {pima_data.output_dir} already exists. It will be removed."))
        shutil.rmtree(pima_data.output_dir)
    elif os.path.isfile(pima_data.output_dir):
        log_messages.append(("warn", f'[{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}]', f"Output directory {pima_data.output_dir} already exists. It will be removed."))
        os.remove(pima_data.output_dir)

    os.makedirs(pima_data.output_dir)
    start_logging(pima_data)
    report_logs(pima_data, log_messages)

def report_logs(pima_data, log_messages):
    for message in log_messages:
        if isinstance(message, str):
            print_and_log(
                pima_data,
                message,
                pima_data.main_process_verbosity,
                pima_data.main_process_color,
            )
        else:
            if message[0] == "main":
                print_and_log(
                    pima_data,
                    message[2],
                    pima_data.main_process_verbosity,
                    pima_data.main_process_color,
                    message[1],
                )
            elif message[0] == "warn":
                print_and_log(
                    pima_data,
                    message[2],
                    pima_data.warning_verbosity,
                    pima_data.warning_color,
                    message[1],
                )
