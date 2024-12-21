import re
import os

import pandas as pd

from Pima.pima_data import PimaData
from Pima.utils.settings import Settings

from Pima.utils.utils import (
    validate_file_and_size,
    validate_file_and_size_or_error,
    print_and_log,
    print_and_run,
    error_out,
    add_warning,
    find_checkpoint,
    make_start_file,
    make_finish_file,
    format_kmg,
)
from Pima.utils.mapping import (
    minimap_and_sort,
    filter_bam,
)


def validate_ont_fastq(pima_data: PimaData, settings: Settings):
    # skip condition
    if not pima_data.ont_fastq:
        return

    if os.path.isdir(pima_data.ont_fastq):
        pima_data.errors.append(f"The provided '--ont-fastq' is a directory, did you mean to include the '--multiplexed' flag?\nOffending input: {pima_data.ont_fastq}")
        return
    
    print_and_log(
        pima_data,
        "Validating ONT FASTQ",
        pima_data.main_process_verbosity,
        pima_data.main_process_color,
    )

    if not validate_file_and_size(
        pima_data, the_file=pima_data.ont_fastq, min_size=1000
    ):
        error = f"Input ONT FASTQ file {pima_data.ont_fastq} cannot be found"
        pima_data.errors.append(error)

    pima_data.ont_fastq = os.path.realpath(pima_data.ont_fastq)
    pima_data.will_have_ont_fastq = True
    pima_data.analysis.append(["info_given_ont_fastq", pima_data, settings])


def validate_illumina_fastq(pima_data: PimaData):
    if not pima_data.illumina_fastq:
        return
    
    print_and_log(
        pima_data,
        'Validating Illumina data', 
        pima_data.main_process_verbosity, 
        pima_data.main_process_color,
    )
    
    if not type(pima_data.illumina_fastq) is list:
        pima_data.illumina_fastq = [pima_data.illumina_fastq]

    for r_fastq in pima_data.illumina_fastq:
        if not validate_file_and_size(pima_data, the_file = r_fastq, min_size = 1000):
            pima_data.errors.append('Illumina FASTQ file ' + r_fastq + ' cannot be found or is size 0')

    pima_data.analysis.append(['info_illumina_fastq', pima_data])


def validate_genome_estimate(pima_data: PimaData, settings: Settings):
    # skip conditions
    if pima_data.no_assembly:
        return

    if pima_data.genome_fasta is not None:
        return

    #skip this if we already have a genome assembly, let the assembly module re-initiate assembly info
    if pima_data.resume:
        if find_checkpoint(pima_data, os.path.join(pima_data.output_dir, "ont_assembly")):
            return

    if not pima_data.will_have_ont_fastq:
        return

    if not pima_data.genome_assembly_size == "estimate":
        return
    
    print_and_log(
        pima_data,
        "Estimating expected genome size using median single copy gene coverages",
        pima_data.main_process_verbosity,
        pima_data.main_process_color,
    )
    median_cov = estimate_genome_size(pima_data, settings)
    result = (
        f"Estimated genome size: {round(pima_data.genome_assembly_raw_size)}\t"   
        f"Median coverage: {round(median_cov,1)}"
    )
    print_and_log(
        pima_data,
        result,
        pima_data.main_process_verbosity,
        pima_data.main_process_color,
    )


def info_ont_fastq(pima_data: PimaData, settings: Settings):
    print_and_log(
        pima_data,
        "Getting ONT FASTQ info",
        pima_data.main_process_verbosity,
        pima_data.main_process_color,
    )

    opener = "cat"
    if re.search(r"\.(gz|gzip)$", pima_data.ont_fastq):
        opener = "gunzip -c"

    command = " ".join(
        [
            opener,
            pima_data.ont_fastq,
            "| awk '{getline;print length($0);s += length($1);getline;getline;} END{print \"+\"s}'",
            "| sort -gr",
            "| awk 'BEGIN{bp = 0;f = 0}",
            '{if(NR == 1){sub(/+/,"", $1);s=$1} else{bp += $1;if(bp > s / 2 && f == 0){n50 = $1;f = 1}}}',
            'END{printf "%d\\t%d\\t%d\\n", n50, (NR - 1), s;exit}\'',
        ]
    )
    result = list(re.split(r"\t", print_and_run(pima_data, command)[0]))
    if result[1] == "0":
        error_out(pima_data, "No ONT reads found")
    (pima_data.ont_n50, pima_data.ont_read_count, pima_data.ont_raw_bases) = [
        int(i) for i in result
    ]

    command = " ".join(
        [
            opener,
            pima_data.ont_fastq,
            "| awk '{getline;print length($0);getline;getline;}'",
        ]
    )
    result = print_and_run(pima_data, command)
    result = list(filter(lambda x: x != "", result))
    pima_data.ont_read_lengths = [int(i) for i in result]


def info_given_ont_fastq(pima_data: PimaData, settings: Settings):
    info_ont_fastq(pima_data, settings)
    validate_genome_estimate(pima_data, settings)

    if pima_data.ont_n50 <= pima_data.ont_n50_min:
        warning = (
            "ONT N50 ("
            + str(pima_data.ont_n50)
            + ") is less than the recommended minimum ("
            + str(pima_data.ont_n50_min)
            + ")."
        )
        add_warning(pima_data, warning)
        pima_data.assembly_notes = pd.concat([pima_data.assembly_notes, pd.Series(warning, dtype='object')])
    
    # See if we should downsample the ONT FASTQ file
    if not (pima_data.genome_assembly_size):
        print_and_log(
            pima_data,
            "Cannot downsample since --genome-size was not provided",
            pima_data.sub_process_verbosity,
            pima_data.sub_process_color,
        )
        return
    
    if pima_data.genome_assembly_raw_size is not None:
        present_coverage = pima_data.ont_raw_bases / pima_data.genome_assembly_raw_size
        if present_coverage >= pima_data.assembly_coverage + 1:
            downsample_ont_fastq(pima_data)
            return

    print_and_log(
        pima_data,
        "No downsampling of reads performed",
        pima_data.sub_process_verbosity,
        pima_data.sub_process_color,
    )


def downsample_ont_fastq(pima_data: PimaData):
    print_and_log(
        pima_data,
        f"Downsampling ONT data to {pima_data.assembly_coverage}X coverage",
        pima_data.sub_process_verbosity,
        pima_data.sub_process_color,
    )

    opener = "cat"
    if re.search(r"\.(gz|gzip)$", pima_data.ont_fastq):
        opener = "gunzip -c"

    pima_data.downsample_dir = os.path.join(pima_data.output_dir, "downsample")
    
    if find_checkpoint(pima_data, pima_data.downsample_dir):
        downsampled_fastq = os.path.join(pima_data.downsample_dir, "downsampled.fastq")
        pima_data.ont_fastq = downsampled_fastq
        return

    os.makedirs(pima_data.downsample_dir)
    make_start_file(pima_data, pima_data.downsample_dir)

    desired_bases = int(
        pima_data.genome_assembly_raw_size * pima_data.assembly_coverage
    )

    read_length_tsv = os.path.join(pima_data.downsample_dir, "read_sizes.tsv")
    command = " ".join(
        [
            opener,
            pima_data.ont_fastq,
            "| awk '{printf $1\"\t\";getline;print length($1);getline;getline}'",
            "| sort -k 2gr",
            "| awk 'BEGIN{bp = 0}{bp += $2;print; if (bp >= ",
            str(desired_bases),
            "){exit}}'",
            ">",
            read_length_tsv,
        ]
    )
    print_and_run(pima_data, command)
    validate_file_and_size_or_error(
        pima_data, read_length_tsv, "Read length index", "doesn't exist", "is empty"
    )

    downsampled_fastq = os.path.join(pima_data.downsample_dir, "downsampled.fastq")
    command = " ".join(
        [
            opener,
            pima_data.ont_fastq,
            "| awk 'BEGIN{while(getline < \"" + read_length_tsv + '"){l[$1] = 1}}',
            "{if($1 in l){print;getline;print;getline;print;getline;print}}'",
            ">",
            downsampled_fastq,
        ]
    )
    print_and_run(pima_data, command)
    validate_file_and_size_or_error(
        pima_data, downsampled_fastq, "Downsampled FASTQ", "doesn't exist", "is empty"
    )

    # Keep this new FASTQ
    pima_data.ont_fastq = downsampled_fastq
    make_finish_file(pima_data, pima_data.downsample_dir)


def info_illumina_fastq(pima_data: PimaData) :

    print_and_log(
        pima_data,
        'Getting Illumina FASTQ info', 
        pima_data.main_process_verbosity, 
        pima_data.main_process_color,
    )

    pima_data.illumina_length_mean, pima_data.illumina_read_count, pima_data.illumina_bases = 0,0,0
    opener = 'cat'
    if re.search(r'\.(gz|gzip)$', pima_data.illumina_fastq[0]):
        opener = 'gunzip -c'

    for r in range(len(pima_data.illumina_fastq)):
        r_fastq = pima_data.illumina_fastq[r]
        command = ' '.join(
            [
                opener,
                r_fastq,
                '| awk \'{getline;s += length($1);getline;getline;}END{print s/(NR/4)"\t"(NR/4)"\t"s}\'',
            ]
        )
        
        values = [float(i) for i in re.split(r'\t', print_and_run(pima_data, command)[0])]
        pima_data.illumina_length_mean += values[0]
        pima_data.illumina_read_count += int(values[1])
        pima_data.illumina_bases += int(values[2])
        
    pima_data.illumina_length_mean /= len(pima_data.illumina_fastq)
    pima_data.illumina_bases = format_kmg(pima_data.illumina_bases, decimals = 1)


def estimate_genome_size(pima_data: PimaData, settings: Settings):
    """Estimate expected genome size using median coverage of a set of 20 universal single copy genes (pulled from Bacillus anthracis)
    
    The set was pulled from this reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10631056/
    """

    reference_ffn = os.path.join(settings.data_dir, 'ames_single_copy_genes.fna')

    temp_bam = os.path.join(pima_data.output_dir, 'single_copy_genes.bam')
    minimap_and_sort(
        pima_data,
        reference_ffn,
        temp_bam,
        pima_data.ont_fastq,
        ont=True,
    )
    filter_bam(
        pima_data,
        inbam = temp_bam,
        F = "4",
        q = "0",
    )

    command = " ".join(
        [
            "samtools coverage",
            temp_bam,
            "| grep -v '^#'",
            "| cut -f 7 | sort -n",
            "| awk '{ a[i++]=$1; } END { print a[int(i/2)]; }'",
        ]
    )
    median_cov = float(print_and_run(pima_data, command)[0])
    pima_data.genome_assembly_raw_size = pima_data.ont_raw_bases / median_cov
    [os.remove(file) for file in [temp_bam, temp_bam + '.bai', os.path.splitext(temp_bam)[0] + '.stderr']]
    return median_cov