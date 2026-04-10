import re
import os

import pandas as pd

from Pima.pima_data import PimaData
from Pima.utils.utils import (
    print_and_log,
    validate_utility,
    validate_file_and_size,
    validate_file_and_size_or_error,
    add_warning,
    print_and_run,
    find_checkpoint,
    std_files,
    make_start_file,
    make_finish_file,
)
from Pima.utils.mapping import (
    minimap_and_sort,
    filter_bam,
)

def validate_genome_fasta(pima_data: PimaData):
    if not pima_data.genome_fasta:
        return

    if pima_data.no_assembly:
        warning = f"Providing a genome file and using the flag '--no-assembly' is redundant."
        add_warning(pima_data, warning)
        pima_data.no_assembly = None
        
    print_and_log(
        pima_data,
        'Validating genome FASTA', 
        pima_data.main_process_verbosity, 
        pima_data.main_process_color,
    )

    if not validate_file_and_size(pima_data, pima_data.genome_fasta, min_size = 1000) :
        pima_data.errors.append('Input genome FASTA ' + pima_data.genome_fasta + ' cannot be found')
    else :
        pima_data.load_genome()
        
    pima_data.will_have_genome_fasta = True
    

def validate_genome_assembly_size(pima_data: PimaData):
    # skip conditions
    if pima_data.no_assembly:
        return

    if pima_data.genome_assembly_size is None:
        return

    print_and_log(
        pima_data,
        "Validating assembly size",
        pima_data.main_process_verbosity,
        pima_data.main_process_color,
    )

    if not re.match(r"^[0-9]+(\.[0-9]+)?[mkgMKG]?$", pima_data.genome_assembly_size):
        if not pima_data.genome_assembly_size == 'estimate':
            error = f"--genome-size needs to be a floating point number, accepting k/m/g, got {pima_data.genome_assembly_size}"
            pima_data.errors.append(error)
            print_and_log(
                pima_data,
                error,
                pima_data.fail_verbosity,
                pima_data.error_color,
            )
            return

    if pima_data.genome_assembly_size == 'estimate':
        #pima_data.genome_assembly_raw size is set from the fastq.py script
        return
    
    power = 0
    if re.match(r"^[0-9]+(\.[0-9]+)?[kK]$", pima_data.genome_assembly_size):
        power = 3
    elif re.match(r"^[0-9]+(\.[0-9]+)?[mM]$", pima_data.genome_assembly_size):
        power = 6
    elif re.match(r"^[0-9]+(\.[0-9]+)?[gG]$", pima_data.genome_assembly_size):
        power = 9

    pima_data.genome_assembly_raw_size = (
        float(re.sub("[mkgMKG]", "", pima_data.genome_assembly_size)) * 10**power
    )


def validate_assembler(pima_data: PimaData):
    if pima_data.no_assembly:
        return

    if pima_data.genome_fasta is not None:
        return

    if not (pima_data.will_have_ont_fastq or pima_data.illumina_fastq):
        return

    if pima_data.ont_fastq and pima_data.assembler in ["flye", "flye_sup"]:
        if pima_data.assembler == "flye_sup":
            pima_data.assembler == "flye"
            pima_data.flye_sup == True
        validate_flye(pima_data)

    elif pima_data.ont_fastq and pima_data.assembler == "raven":
        validate_raven(pima_data)

    elif pima_data.illumina_fastq and not pima_data.will_have_genome_fasta:
        #only run spades if we know we don't have a better assembly
        validate_spades(pima_data)

    else:
        return
    

def validate_flye(pima_data: PimaData):

    print_and_log(
        pima_data,
        "Validating flye utilities",
        pima_data.main_process_verbosity,
        pima_data.main_process_color,
    )

    if validate_utility(
        pima_data, "flye", "flye is not on the PATH (required by --assembler flye)."
    ):
        command = "flye --version"
        pima_data.versions["flye"] = re.search(
            r"[0-9]+\.[0-9.]+", print_and_run(pima_data, command)[0]
        ).group(0)

    pima_data.will_have_ont_assembly = True
    pima_data.will_have_genome_fasta = True
    pima_data.analysis.append(["flye_ont_fastq", pima_data])


def validate_raven(pima_data: PimaData):

    print_and_log(
        pima_data,
        "Validating Raven assembler utilities",
        pima_data.main_process_verbosity,
        pima_data.main_process_color,
    )
    
    if validate_utility(
        pima_data, "raven", "raven is not on the PATH (required by --assembler raven)."
    ):
        command = "raven --version"
        pima_data.versions["raven"] = re.search(
            r"[0-9]+\.[0-9.]+", print_and_run(pima_data, command)[0]
        ).group(0)

    pima_data.will_have_ont_assembly = True
    pima_data.will_have_genome_fasta = True
    pima_data.analysis.append(["raven_ont_fastq", pima_data])


def validate_spades(pima_data: PimaData):

    print_and_log(
        pima_data,
        "Validating SPAdes assembler utilities",
        pima_data.main_process_verbosity,
        pima_data.main_process_color,
    )
    
    if validate_utility(pima_data, 'spades.py', 'spades.py is not on the PATH (required by --illumina-fastq)') :
        command = 'spades.py --version'
        pima_data.versions['spades'] = re.search(r'[0-9]+\.[0-9.]+', print_and_run(pima_data, command)[0]).group(0)
            
    pima_data.analysis.append(['spades_illumina_fastq', pima_data])
    pima_data.will_have_genome_fasta = True


def validate_assembly_info(pima_data):
    if not (pima_data.will_have_ont_fastq or pima_data.illumina_fastq):
        return

    if not pima_data.will_have_genome_fasta:
        return

    pima_data.contig_info = pd.Series(dtype=object)

    validate_utility(pima_data, "samtools", "is not on the PATH")
    pima_data.analysis.append(["check_assembly_coverages", pima_data])


def flye_ont_fastq(pima_data: PimaData):
    print_and_log(
        pima_data,
        "Assembling ONT reads using flye",
        pima_data.main_process_verbosity,
        pima_data.main_process_color,
    )

    # Make the directory for new assembly files
    pima_data.ont_assembly_dir = os.path.join(pima_data.output_dir, "ont_assembly")
    if find_checkpoint(pima_data, pima_data.ont_assembly_dir):
        pima_data.genome_fasta = pima_data.ont_assembly_dir + "/assembly.fasta"
        pima_data.load_genome()
        pima_data.did_flye_ont_fastq = True

        print_and_log(
            pima_data,
            "Flye previously finished successfully",
            pima_data.main_process_verbosity,
            pima_data.main_process_color,
        )
        return

    os.makedirs(pima_data.ont_assembly_dir)
    make_start_file(pima_data, pima_data.ont_assembly_dir)

    # Assemble with Flye
    print_and_log(
        pima_data,
        "Running flye",
        pima_data.sub_process_verbosity,
        pima_data.sub_process_color,
    )

    flye_output_dir = pima_data.ont_assembly_dir
    flye_stdout, flye_stderr = std_files(
        os.path.join(pima_data.ont_assembly_dir, "flye")
    )
    flye_fasta = os.path.join(flye_output_dir, "assembly.fasta")

    # Actually run Flye // need to update this for R9 vs R10 data
    raw_or_corrected = "--nano-raw"
    if pima_data.flye_sup or re.search(r"sup", os.path.basename(pima_data.ont_model)):
        raw_or_corrected = "--nano-hq"

    if pima_data.genome_assembly_size == "estimate":
        genome_size = f"--genome-size {str(round(pima_data.genome_assembly_raw_size))}"
    elif pima_data.genome_assembly_size is None:
        genome_size = ""
    else:
        genome_size = f"--genome-size {pima_data.genome_assembly_size}"
    command = " ".join(
        [
            "flye",
            raw_or_corrected,
            pima_data.ont_fastq,
            genome_size,
            "--out-dir",
            flye_output_dir,
            "--threads",
            str(pima_data.threads),
            "--no-alt-contigs",
            "1>",
            flye_stdout,
            "2>",
            flye_stderr,
        ]
    )
    #########
    #########   TEMP TO AVOID ASSEMBLING DURING DEV
    #########
    #command = " ".join(["cp", "-r", os.path.join(pima_data.backup_dir, "ont_assembly/*"), flye_output_dir])
    print_and_run(pima_data, command)
    validate_file_and_size_or_error(
        pima_data, flye_fasta, "Flye fasta", "cannot be found after flye", "is empty"
    )

    pima_data.genome_fasta = pima_data.ont_assembly_dir + "/assembly.fasta"
    validate_file_and_size_or_error(
        pima_data,
        pima_data.genome_fasta,
        "Genome fasta",
        "cannot be found after copying Flye output",
        "is empty",
    )
    pima_data.load_genome()

    # Pull in the assembly summary and look at the coverage
    assembly_info_txt = os.path.join(pima_data.ont_assembly_dir, "assembly_info.txt")
    assembly_info = pd.read_csv(assembly_info_txt, header=0, index_col=0, sep="\t")

    # Look for non-circular contigs
    open_contigs = assembly_info.loc[assembly_info["circ."] == "N", :]
    if open_contigs.shape[0] > 0:
        open_contig_ids = open_contigs.index.values
        warning = f"Flye reported {open_contigs.shape[0]} open contigs ({', '.join(open_contig_ids)}); assembly may be incomplete."
        add_warning(pima_data, warning)
        pima_data.assembly_notes = pd.concat([pima_data.assembly_notes, pd.Series(warning, dtype='object')])

    make_finish_file(pima_data, pima_data.ont_assembly_dir)
    pima_data.did_flye_ont_fastq = True


def raven_ont_fastq(pima_data: PimaData):
    print_and_log(
        pima_data,
        "Assembling ONT reads using raven",
        pima_data.main_process_verbosity,
        pima_data.main_process_color,
    )

    # Make the directory for new assembly files
    pima_data.ont_assembly_dir = os.path.join(pima_data.output_dir, "ont_assembly")
    if find_checkpoint(pima_data, pima_data.ont_assembly_dir):
        pima_data.genome_fasta = pima_data.ont_assembly_dir + "/assembly.fasta"
        pima_data.load_genome()
        pima_data.did_raven_ont_fastq = True

        print_and_log(
            pima_data,
            "Raven previously finished successfully",
            pima_data.main_process_verbosity,
            pima_data.main_process_color,
        )
        return

    os.makedirs(pima_data.ont_assembly_dir)
    make_start_file(pima_data, pima_data.ont_assembly_dir)

    # Assemble with Flye
    print_and_log(
        pima_data,
        "Running raven",
        pima_data.sub_process_verbosity,
        pima_data.sub_process_color,
    )

    raven_output_dir = pima_data.ont_assembly_dir
    _, raven_stderr = std_files(
        os.path.join(pima_data.ont_assembly_dir, "raven")
    )
    raven_fasta = os.path.join(raven_output_dir, "assembly.fasta")
    raven_gfa = os.path.join(raven_output_dir, "assembly.gfa")
    command = " ".join(
        [
            "raven",
            "--threads",
            str(pima_data.threads),
            "--min-unitig-size",
            "4999g",
            "--graphical-fragment-assembly",
            raven_gfa,
            pima_data.ont_fastq,
            "--disable-checkpoints",
            "1>",
            raven_fasta,
            "2>",
            raven_stderr,
        ]
    )

    print_and_run(pima_data, command)
    validate_file_and_size_or_error(
        pima_data, raven_fasta, "raven fasta", "cannot be found after raven", "is empty"
    )
    
    print_and_log(
        pima_data,
        'Adjusting contig names after Raven', 
        pima_data.sub_process_verbosity, 
        pima_data.sub_process_color,
    )
    pima_data.genome_fasta = pima_data.ont_assembly_dir + "/assembly.fasta"
    
    command = " ".join(
        [
            'awk \'BEGIN {i=0} {if($0 ~ /^>/){print ">contig_" i; i++}else{print} }\'',
            pima_data.genome_fasta,
            '>', os.path.join(pima_data.ont_assembly_dir, "temp.fasta"),
            '&& mv', os.path.join(pima_data.ont_assembly_dir, "temp.fasta"), pima_data.genome_fasta,
        ]
    )
    print_and_run(pima_data, command)
    validate_file_and_size_or_error(pima_data, pima_data.genome_fasta, 'Genome assembly', 
                                            'cannot be found after fixing names', 'is empty')
    pima_data.load_genome()
    make_finish_file(pima_data, pima_data.ont_assembly_dir)
    pima_data.did_raven_ont_fastq = True


def spades_illumina_fastq(pima_data: PimaData):

    print_and_log(pima_data,
        'Assembling Illumina FASTQ data with SPAdes', 
        pima_data.main_process_verbosity, 
        pima_data.main_process_color,
    )
    
    # Figure out where the assembly is going
    pima_data.illumina_asm_dir = os.path.join(pima_data.output_dir, 'illumina_assembly')
    
    if find_checkpoint(pima_data, pima_data.illumina_asm_dir):
        assembly_fasta = os.path.join(pima_data.illumina_asm_dir, 'assembly.fasta')
        pima_data.genome_fasta = assembly_fasta
        pima_data.did_spades_illumina_fastq = True
        pima_data.load_genome()

        print_and_log(
            pima_data,
            "Spades previously finished successfully",
            pima_data.main_process_verbosity,
            pima_data.main_process_color,
        )
        return
    
    os.makedirs(pima_data.illumina_asm_dir)
    make_start_file(pima_data, pima_data.illumina_asm_dir)
    
    # Figure out if we have single-end or paired-end reads
    fastq_input = '-s ' + pima_data.illumina_fastq[0]
    if len(pima_data.illumina_fastq) == 2:
        fastq_input = ' '.join(['-1', pima_data.illumina_fastq[0],
                                '-2', pima_data.illumina_fastq[1]])

    spades_stderr, spades_stdout = std_files(os.path.join(pima_data.illumina_asm_dir, 'spades'))
    command = ' '.join(
        [
            'spades.py',
            fastq_input,
            '-o', pima_data.illumina_asm_dir,
            '-t', str(pima_data.threads),
            '-k auto --isolate',
            '1>', spades_stdout, '2>', spades_stderr, #flipped because spades writes the error to stdout, but we need to report this if there is an issue
        ]
    )
    print_and_run(pima_data, command)

    pima_data.genome_fasta = os.path.join(pima_data.illumina_asm_dir, 'scaffolds.fasta')
    if len(pima_data.illumina_fastq) == 1:
        pima_data.genome_fasta = os.path.join(pima_data.illumina_asm_dir, 'contigs.fasta')

    # Fix the silly contig names
    print_and_log(
        pima_data,
        'Repairing contig names after SPAdes', 
        pima_data.sub_process_verbosity, 
        pima_data.sub_process_color,
    )
    fixed_names_fasta = os.path.join(pima_data.illumina_asm_dir, 'fixed_names.fasta')
    command = ' '.join(
        [
            'awk \'{if($0 ~ /^>/){gsub("NODE", "contig", $0);gsub("_length.*", "", $0)}print}\'', 
            pima_data.genome_fasta,
            '>', fixed_names_fasta,
        ]
    )
    print_and_run(pima_data, command)
    validate_file_and_size_or_error(pima_data, fixed_names_fasta, 'Genome assembly', \
                                            'cannot be found after fixing names', 'is empty')
        
    # Filter small contigs
    assembly_fasta = os.path.join(pima_data.illumina_asm_dir, 'assembly.fasta')
    command = ' '.join(
        [
            'faidx -i chromsizes',
            fixed_names_fasta,
            '| awk \'($2 > 1000){print $1}\'',
            '| parallel faidx',
            fixed_names_fasta,
            '>', assembly_fasta,
        ]
    )
    print_and_run(pima_data, command)
    pima_data.genome_fasta = assembly_fasta
    validate_file_and_size_or_error(pima_data, pima_data.genome_fasta, 'Genome assembly', 
                                            'cannot be found after filtering short contigs', 'is empty')
    pima_data.load_genome()
    pima_data.did_spades_illumina_fastq = True
    make_finish_file(pima_data, pima_data.illumina_asm_dir)


def check_assembly_coverages(pima_data: PimaData):
    print_and_log(
        pima_data,
        "Getting assembly description/coverage",
        pima_data.main_process_verbosity,
        pima_data.main_process_color,
    )

    pima_data.info_dir = os.path.join(pima_data.output_dir, "info")

    if find_checkpoint(pima_data, pima_data.info_dir):
        if pima_data.ont_fastq:
            parse_assembly_info_fastq(pima_data, read_type="ONT", output_prefix="ont_coverage")
        if pima_data.illumina_fastq:
            parse_assembly_info_fastq(
                pima_data, read_type="Illumina", output_prefix="illumina_coverage"
            )
        return

    os.makedirs(pima_data.info_dir)
    make_start_file(pima_data, pima_data.info_dir)

    if pima_data.ont_fastq:
        assembly_info_fastq(
            pima_data,
            fastq=pima_data.ont_fastq,
            read_type="ONT",
            output_prefix="ont_coverage",
        )
    if pima_data.illumina_fastq:
        assembly_info_fastq(
            pima_data,
            fastq=pima_data.illumina_fastq,
            read_type="Illumina",
            output_prefix="illumina_coverage",
        )


def parse_assembly_info_fastq(pima_data: PimaData, read_type: str, output_prefix: str):
    """Read the assembly coverages and populate the information into the Assembly class for the report generation"""
    coverage_tsv = os.path.join(pima_data.info_dir, output_prefix + ".tsv")
    validate_file_and_size_or_error(
        pima_data,
        coverage_tsv,
        "Coverage TSV",
        "cannot be found after samtools",
        "is empty",
    )

    pima_data.contig_info[read_type] = pd.read_csv(
        coverage_tsv, header=None, index_col=None, sep="\t"
    ).sort_values(1, axis=0, ascending=False)
    pima_data.contig_info[read_type].columns = ["contig", "size", "coverage"]

    mean_coverage = (
        pima_data.contig_info[read_type].iloc[:, 1]
        * pima_data.contig_info[read_type].iloc[:, 2]
    ).sum() / pima_data.contig_info[read_type].iloc[:, 1].sum()
    pima_data.mean_coverage[read_type] = mean_coverage
    
    # TODO make this use ONT/Illumina specific coverage
    if mean_coverage <= pima_data.ont_coverage_min:
        warning = (
            read_type
            + " mean coverage ({:.0f}X) is less than the recommended minimum ({:.0f}X).".format(
                mean_coverage, pima_data.ont_coverage_min
            )
        )
        print_and_log(
            pima_data, warning, pima_data.warning_verbosity, pima_data.warning_color
        )

    # If some contigs have low coverage, report that
    # TODO make this use ONT/Illumina specific coverage
    low_coverage = pima_data.contig_info[read_type].loc[
        pima_data.contig_info[read_type]["coverage"] < pima_data.ont_coverage_min, :
    ]
    if low_coverage.shape[0] >= 0:
        for contig_i in range(low_coverage.shape[0]):
            warning = (
                read_type
                + " coverage of {:s} ({:.0f}X) is less than the recommended minimum ({:.0f}X).".format(
                    low_coverage.iloc[contig_i, 0],
                    low_coverage.iloc[contig_i, 2],
                    pima_data.ont_coverage_min,
                )
            )
            print_and_log(
                pima_data, warning, pima_data.warning_verbosity, pima_data.warning_color
            )

    # See if some contigs have anolously low coverage
    fold_coverage = pima_data.contig_info[read_type]["coverage"] / mean_coverage
    low_coverage = pima_data.contig_info[read_type].loc[fold_coverage < 1 / 5, :]
    if low_coverage.shape[0] >= 0:
        for contig_i in range(low_coverage.shape[0]):
            warning = (
                read_type
                + " coverage of {:s} ({:.0f}X) is less than 1/5 the mean coverage ({:.0f}X).".format(
                    low_coverage.iloc[contig_i, 0],
                    low_coverage.iloc[contig_i, 2],
                    mean_coverage,
                )
            )
            print_and_log(
                pima_data, warning, pima_data.warning_verbosity, pima_data.warning_color
            )


def assembly_info_fastq(pima_data, fastq, read_type, output_prefix):
    """Generate the bam files & calculate assembly coverages"""
    coverage_tsv = os.path.join(pima_data.info_dir, output_prefix + ".tsv")
    # Map reads to the assembly
    unfilt_coverage_bam = os.path.join(pima_data.info_dir, output_prefix + "_unfilt.bam")
    coverage_bam = os.path.join(pima_data.info_dir, output_prefix + ".bam")
    pima_data.files_to_clean.extend([coverage_bam, unfilt_coverage_bam])

    if read_type == "ONT":
        minimap_and_sort(
            pima_data, pima_data.genome_fasta, coverage_bam, fastq, ont=True
        )
    else:
        minimap_and_sort(
            pima_data,
            pima_data.genome_fasta,
            unfilt_coverage_bam,
            fastq,
            ont=False,
        )
        filter_bam(pima_data, inbam=unfilt_coverage_bam, outbam=coverage_bam, F="0x0100", q="30")
    # Generate a table of contig/coverage
    command = " ".join(
        [
            "samtools depth -a",
            coverage_bam,
            "| awk '{s[$1] += $3; c[$1]++}",
            'END{for(i in s){printf "%s\\t%i\\t%.0f\\n", i, c[i], (s[i] / c[i])}}\'',
            ">",
            coverage_tsv,
        ]
    )
    print_and_run(pima_data, command)
    parse_assembly_info_fastq(pima_data, read_type, output_prefix)
    make_finish_file(pima_data, pima_data.info_dir)
