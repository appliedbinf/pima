import re
import os

from Pima.pima_data import PimaData
from Pima.utils.utils import (
    print_and_log,
    print_and_run,
    validate_file_and_size,
    validate_file_and_size_or_error,
    std_files,
)


def minimap_and_sort(
    pima_data: PimaData,
    genome: str,
    bam: str,
    fastq: list,
    ont: bool = True):

    if ont:
        map_mode = "map-ont"
    else:
        map_mode = "sr"

    if not type(fastq) is list:
        fastq = [fastq]
    
    std_prefix = re.sub(r"\.bam$", "", bam)
    _, minimap_stderr = std_files(std_prefix)
    command = " ".join(
        [
            "minimap2 -a",
            "-t",
            str(pima_data.threads),
            "-x",
            map_mode,
            genome,
            " ".join(fastq), #handle 1 or 2 fastqs
            "2>", minimap_stderr,
            "| samtools sort",
            "-@",
            str(pima_data.threads),
            "-o",
            bam,
            "-T reads.tmp -",
            "1>/dev/null 2>/dev/null",
        ]
    )

    print_and_run(pima_data, command)
    validate_file_and_size_or_error(pima_data, the_file=bam, error_prefix="The file", presence_suffix="doesn't exist", size_suffix="is below min expected size", min_size=1000)
    index_bam(pima_data, bam)

def filter_bam(pima_data: PimaData, inbam: str, outbam: str = None, F: str = None, q: str = None):
    """Filter the bam file, if outbam not provided, we filter in-place"""
    if not outbam:
        outbam = inbam
    command = " ".join(
        [
            "samtools view -h",
            "-F",
            F,
            "-q",
            q,
            inbam,
            "| samtools sort",
            "-@",
            str(pima_data.threads),
            "-o",
            outbam,
            "-T reads.tmp -",
            "1>/dev/null 2>/dev/null",
        ]
    )
    print_and_run(pima_data, command)
    validate_file_and_size_or_error(pima_data, the_file=outbam, min_size=1000)
    index_bam(pima_data, outbam)

def index_bam(pima_data: PimaData, bam: str):
    command = " ".join(["samtools index", bam, "1>/dev/null 2>/dev/null"])
    print_and_run(pima_data, command)
    index_bai = bam + ".bai"
    validate_file_and_size_or_error(pima_data, the_file=index_bai, min_size=1000)

def mpileup_bam(pima_data: PimaData, reference_genome: str, bam: str, mpileup: str, output_dir: str):

    print_and_log(
        pima_data,
        "Making mpileup from BAM", 
        pima_data.sub_process_verbosity, 
        pima_data.sub_process_color,
    )
    
    mpileup_stdout, mpileup_stderr = std_files(os.path.join(output_dir, 'mpileup'))
    command = " ".join(
        [
            'samtools mpileup',
            '-B',
            '-a',
            '-f', reference_genome,
            '-o' + mpileup,
            bam,
            '1>', mpileup_stdout, 
            '2>', mpileup_stderr,
        ]
    )
    print_and_run(pima_data, command)
    validate_file_and_size_or_error(pima_data,
                                    mpileup, 
                                    'Region MPILEUP file', 
                                    'cannot be found', 
                                    'is empty',
    )

def bwa_index_fasta(pima_data: PimaData, fasta: str):
    
    print_and_log(
        pima_data,
        'Indexing FASTA with bwa index', 
        pima_data.sub_process_verbosity, 
        pima_data.sub_process_color,
    )
        
    # Check for an index already there
    bwa_index = f"{fasta}.bwt"
    if validate_file_and_size(pima_data, bwa_index):
        return
    
    # Make the bwa index
    std_prefix = re.sub(r'\.f(na|asta)$', '', fasta)
    bwa_index_stdout, bwa_index_stderr = std_files(std_prefix + '_index')
    command = " ".join(
        [
            'bwa',
            'index',
            fasta,
            '1>', bwa_index_stdout, '2>', bwa_index_stderr,
        ]
    )
    print_and_run(pima_data, command)

    # Check that the index was built
    validate_file_and_size_or_error(pima_data, bwa_index, 'BWA index', 'doesn\'t exist', 'is empty')
                 
def bwa_short_illumina_fastq_and_sort(pima_data: PimaData, genome: str, fastq: str, bam: str):

    std_prefix = re.sub(r'\.bam$', '', bam)
        
    bwa_index_fasta(pima_data, genome)

    # Align the reads 
    sai = []
    for i in range(len(fastq)):
        bwa_stdout, bwa_stderr = std_files(std_prefix + '_aln')
        this_sai = std_prefix + '_aln_' + str(i) + '.sai'
        command = " ".join(
            [
                'bwa aln',
                '-t', str(pima_data.threads),
                genome,
                fastq[i],
                '1>', this_sai,
                '2>', bwa_stderr,
            ]
        )
        sai.append(this_sai)
        print_and_run(pima_data, command)
        validate_file_and_size_or_error(pima_data, this_sai)
        
    # And turn the SAI into a proper SAM file
    read_type = 'samse'
    if len(fastq) > 1:
        read_type = 'sampe'
    bwa_stdout, bwa_stderr = std_files(std_prefix + '_sam')
    tmp_file = std_prefix + '.tmp'
    command = " ".join(
        [
            'bwa',
            read_type,
            genome,
            ' '.join(sai),
            ' '.join(fastq),
            '2>', bwa_stderr,
            '| samtools',
            'sort',
            '-T', tmp_file,
            '-o', bam,
            '-',
            '1>/dev/null 2>/dev/null',
        ]
    )
    print_and_run(pima_data, command)
    validate_file_and_size_or_error(pima_data, the_file = bam, min_size = 100)
    index_bam(pima_data, bam)

def bwa_mem_all_aln_illumina(pima_data: PimaData, genome: str, fastq: list, bam: str):
    std_prefix = re.sub(r'\.bam$', '', bam)
        
    bwa_index_fasta(pima_data, genome)

    bams = []
    # Align the reads 
    for i in range(len(fastq)):
        read = f"_R{i+1}"
        _, bwa_stderr = std_files(std_prefix + read + '_aln')
        this_bam = std_prefix + read + '.bam'
        #Polypolish warns not to sort the bam files
        command = " ".join(
            [
                'bwa mem',
                '-a',
                '-t', str(pima_data.threads),
                genome,
                fastq[i],
                '1>', this_bam,
                '2>', bwa_stderr,
            ]
        )
        bams.append(this_bam)
        print_and_run(pima_data, command)
        validate_file_and_size_or_error(pima_data, this_bam)
    return bams