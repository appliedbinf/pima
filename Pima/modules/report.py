import os
import shutil
import re
import string

from Pima.pima_data import PimaData
from Pima.utils.settings import Settings
from Pima.utils.utils import (
    print_and_log,
    print_and_run,
    validate_utility,
    find_checkpoint,
    validate_file_and_size_or_error,
    std_files,
)
from Pima.accessory_scripts.MarkdownReport import PimaReport

def validate_make_report(pima_data: PimaData, settings: Settings):
    if pima_data.no_report:
        return

    print_and_log(
        pima_data,
        "Validating reporting utilities",
        pima_data.main_process_verbosity,
        pima_data.main_process_color,
    )

    if pima_data.bundle:
        if not os.path.isdir(pima_data.bundle):
            pima_data.errors.append(f"Can't find pandoc bundle {pima_data.bundle}")

    validate_utility(
        pima_data, "pandoc", "pandoc is not on the PATH (required for reporting)."
    )
    
    pima_data.analysis.append(["make_report", pima_data, settings])

def make_report(pima_data: PimaData, settings: Settings):

    print_and_log(
        pima_data,
        "Making report", 
        pima_data.main_process_verbosity, 
        pima_data.main_process_color,
    )
    pima_data.report_dir = os.path.join(pima_data.output_dir, 'report')

    if find_checkpoint(pima_data, pima_data.report_dir):
        ##Always regenerate report in case downstream steps have changed the results
        shutil.rmtree(pima_data.report_dir) 
    os.mkdir(pima_data.report_dir)

    pima_data.report_prefix = os.path.join(pima_data.report_dir, 'report')
    pima_data.report_md = pima_data.report_prefix + '.md'

    pima_data.markdown_report = PimaReport(pima_data, settings)
    pima_data.markdown_report.make_report()

    ## Add appendices to the pima report
    if len(pima_data.markdown_report.appendices) > 0:
        #assign each new AMR class its own appendix ID (A-Z)
        for i, appendix in zip(string.ascii_uppercase, pima_data.markdown_report.appendices):
            png = re.sub(r"\.md", ".png", appendix)
            mod_appendix = os.path.join(pima_data.report_dir, os.path.basename(appendix))
            shutil.copyfile(png, os.path.join(pima_data.report_dir, os.path.basename(png)))
            with open(appendix, "rt") as fin:
                with open(mod_appendix, "wt") as fout:
                    for line in fin:
                        fout.write(line.replace("LETTER", i))

            # translate the appendix into the markdown report
            with open(pima_data.report_md, "a") as file:
                with open(mod_appendix, "r") as temp_file:
                    file.write(temp_file.read())
            #pima_data.files_to_clean.append(mod_appendix)
            #pima_data.files_to_clean.append(os.path.join(pima_data.report_dir, os.path.basename(png)))

    pima_data.report_pdf = pima_data.report_prefix + '.pdf'
    validate_file_and_size_or_error(pima_data, pima_data.report_md, 'Report MD', 'cannot be found', 'is empty')
    
    tectonic_stdout, tectonic_stderr = std_files(os.path.join(pima_data.report_dir, 'markdown2pdf'))
    command = ' '.join(
        [
            'pandoc -f gfm',
            pima_data.report_md,
            '-o', pima_data.report_pdf,
            '--pdf-engine=weasyprint',
            '--css ' + settings.pima_css,
            '1>', tectonic_stdout, 
            '2>', tectonic_stderr,
        ]
    )
    print_and_run(pima_data, command, change_exe_dir=pima_data.report_dir)
    validate_file_and_size_or_error(pima_data, pima_data.report_pdf, 'Report PDF', 'cannot be found', 'is empty')