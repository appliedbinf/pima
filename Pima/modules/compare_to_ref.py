import os
import re
import shutil
from typing import Literal

import pandas as pd
import numpy as np

from Pima.pima_data import PimaData
from Pima.utils.utils import (
    print_and_log,
    print_and_run,
    validate_utility,
    validate_file_and_size,
    validate_file_and_size_or_error,
    add_warning,
    find_checkpoint,
    make_start_file,
    make_finish_file,
    std_files,
    error_out,
)
from Pima.utils.mapping import (
    minimap_and_sort,
    filter_bam,
    mpileup_bam,
)

def validate_reference_fasta(pima_data: PimaData):
    # skip conditions
    if pima_data.only_assemble:
        return

    if not pima_data.reference_fasta:
        return

    pima_data.reference_fasta = os.path.realpath(pima_data.reference_fasta)

    print_and_log(
        pima_data,
        "Validating reference FASTA",
        pima_data.main_process_verbosity,
        pima_data.main_process_color,
    )

    if not validate_file_and_size(pima_data, pima_data.reference_fasta, min_size=1000):
        pima_data.errors.append(
            f"Reference FASTA {pima_data.reference_fasta} can't be found or is size 0."
        )
        # prevents nonsensical error when checking if the mutation regions are present in a nonexistant file
        pima_data.will_have_reference_fasta = False

    else:
        pima_data.load_reference()

        # Check for utils
        for utility in ['nucmer', 'mummer']:
            validate_utility(pima_data, utility, f"{utility} is not on the PATH.")  
        if validate_utility(pima_data, 'dnadiff', 'dnadiff is not on the PATH.'):
            command = 'dnadiff -version 2>&1'
            pima_data.versions['dnadiff'] = re.search(r'[0-9]+\.[0-9.]*', print_and_run(pima_data, command)[1]).group(0)

    pima_data.analysis.append(['call_insertions', pima_data])

def validate_quast(pima_data: PimaData):
    # skip conditions
    if pima_data.only_assemble:
        return

    if not pima_data.reference_fasta:
        return
    print_and_log(
        pima_data,
        "Validating running QUAST using reference FASTA",
        pima_data.main_process_verbosity,
        pima_data.main_process_color,
    )

    if validate_utility(pima_data, 'quast.py', 'quast is not on the PATH.'):
        command = 'quast.py --version 2>&1'
        pima_data.versions['quast'] = re.search(r'[0-9]+\.[0-9.]*', print_and_run(pima_data, command)[0]).group(0)

    pima_data.analysis.append(['quast_genome', pima_data])

def validate_mutations(pima_data: PimaData):
    # skip conditions
    if pima_data.only_assemble:
        return

    if not (pima_data.organism or pima_data.mutation_region_bed):
        return

    if pima_data.mutation_region_bed and not pima_data.reference_fasta:
        pima_data.errors.append("--mutation-regions requires --reference-genome")

    if not (pima_data.ont_fastq or pima_data.illumina_fastq):
        return

    print_and_log(
        pima_data,
        "Validating mapping and variant utilities",
        pima_data.main_process_verbosity,
        pima_data.main_process_color,
    )

    if not pima_data.mutation_region_bed:
        return

    # Either from the built in set or given as an argument, make sure the file is there
    pima_data.mutation_region_bed = os.path.realpath(pima_data.mutation_region_bed)
    if validate_file_and_size(pima_data, pima_data.mutation_region_bed, min_size=100):
        #correct_bed_column_names = ['#contig', 'start', 'stop', 'name', 'type', 'drug', 'note']
        correct_bed_column_names = ['#contig', 'start', 'stop', 'name', 'type', 'drug', 'priority', 'note'] #updated for ardap style reporting
        pima_data.mutation_regions = pd.read_csv(
            pima_data.mutation_region_bed, header=0, sep="\t", index_col=False
        )

        if pima_data.mutation_regions.shape[1] != 8:
            message = (
                f"Mutation region file not formatted correctly. \nIt should be an eight column, tab delimited, file with the following headers:\n"
                f"{', '.join(correct_bed_column_names)}"
            )
            pima_data.errors.append(message)

        elif pima_data.mutation_regions.shape[0] == 0:
            pima_data.errors.append("No rows in mutation regions file.")

        elif pima_data.mutation_regions.columns.tolist() != correct_bed_column_names:
            pima_data.errors.append(
                "Mutation regions bed file does not contain the expected column names.\n"
                f"bed file column names: {', '.join(pima_data.mutation_regions.columns.tolist())}\n"
                f"Expected column names: {', '.join(correct_bed_column_names)}\n"                
            )

        elif pima_data.will_have_reference_fasta:
            # Make sure that the positions in the BED file fall within the chromosomes provided in the reference sequence
            for mutation_region in range(pima_data.mutation_regions.shape[0]):
                mutation_region = pima_data.mutation_regions.iloc[mutation_region,:]
                if mutation_region['#contig'] not in pima_data.reference:
                    pima_data.errors.append(
                        f"Mutation region: {' '.join(mutation_region.drop(['type', 'drug', 'priority', 'note']).astype(str))}"
                        + " not found in reference genome."
                    )
                    continue

                if not isinstance(mutation_region["start"], np.int64):
                    pima_data.errors.append(
                        f"Non-integer found in mutation region start (column 2): {mutation_region[1]}"
                    )
                    break
                elif not isinstance(mutation_region["stop"], np.int64):
                    pima_data.errors.append(
                        f"Non-integer found in mutation region stop (column 3): {mutation_region[2]}"
                    )
                    break

                if mutation_region["start"] <= 0 or mutation_region["stop"] <= 0:
                    pima_data.errors.append(
                        f"Mutation region: {' '.join(mutation_region.drop(['type', 'drug', 'priority', 'note']).astype(str))}"
                        + " starts before the reference sequence."
                    )
                if mutation_region["start"] > len(
                    pima_data.reference[mutation_region["#contig"]].seq
                ) or mutation_region["stop"] > len(
                    pima_data.reference[mutation_region["#contig"]].seq
                ):
                    pima_data.errors.append(
                        f"Mutation region: {' '.join(mutation_region.drop(['type', 'drug', 'priority', 'note']).astype(str))}"
                        + " ends after the reference sequence."
                    )
    else:
        pima_data.errors.append(
            f"Mutation region BED {pima_data.mutation_region_bed} can't be found or is size 0."
        )

    # We need reads to make calls
    if not any(
        [
            pima_data.ont_fastq,
            pima_data.illumina_fastq,
            pima_data.will_have_genome_fasta,
        ]
    ):
        pima_data.errors.append(
            "Can't call mutations without a FASTQ dataset or assembly."
        )

    for utility in ["minimap2", "samtools"]:
        if validate_utility(
            pima_data,
            utility,
            f"{utility} is not on the PATH (required for AMR mutations)",
        ):
            command = f"{utility} --version"
            pima_data.versions[utility] = re.search(
                r"[0-9]+\.[0-9.]*", print_and_run(pima_data, command)[0]
            ).group(0)

    if validate_utility(
        pima_data, "varscan", "varscan is not on the PATH (required for AMR mutations)"
    ):
        command = "varscan 2>&1 | head -1"
        quickcheck = print_and_run(pima_data, command)[0]
        pima_data.versions["varscan"] = re.search(r"[0-9]+\.[0-9.]*", quickcheck).group(
            0
        )
    pima_data.amr_region_names = pima_data.mutation_regions.get('name').unique().tolist()
    pima_data.analysis.append(["call_amr_mutations", pima_data])
  
def dnadiff_fasta(pima_data: PimaData, output_prefix: str):

    dnadiff_stdout, dnadiff_stderr = std_files(output_prefix)
    command = ' '.join(['dnadiff',
                        '-p', output_prefix,
                        pima_data.reference_fasta,
                        pima_data.genome_fasta,
                        '1>', dnadiff_stdout, '2>', dnadiff_stderr])
    output_1coords = output_prefix + '.1coords'
    print_and_run(pima_data, command)
    validate_file_and_size_or_error(
        pima_data,
        the_file = output_1coords, 
        min_size = 5,
    )

#switch to using sniffles if we have ONT data for INDEL/SV hunting
def call_insertions(pima_data: PimaData):
    """Use dnadiff to identify insertions / deletions between the assembled genome and the reference"""

    if not pima_data.genome_fasta:
        return
           
    print_and_log(
        pima_data,
        'Calling insertions', 
        pima_data.main_process_verbosity, 
        pima_data.main_process_color,
    )

    pima_data.insertions_dir = os.path.join(pima_data.output_dir, "insertions")

    #Enable checkpointing
    if find_checkpoint(pima_data, pima_data.insertions_dir):
        print_and_log(
            pima_data,
            'Reusing previously identified insertions', 
            pima_data.main_process_verbosity, 
            pima_data.main_process_color,
        )
        pima_data.dnadiff_prefix = os.path.join(pima_data.insertions_dir, "vs_reference")
        parse_insertion_results(pima_data)
        return
    
    os.makedirs(pima_data.insertions_dir)
    make_start_file(pima_data, pima_data.insertions_dir)

    # Align the assembly against the reference sequence
    print_and_log(
        pima_data,
        "Running dnadiff against the reference", 
        pima_data.sub_process_verbosity, 
        pima_data.sub_process_color,
    )
    pima_data.dnadiff_prefix = os.path.join(pima_data.insertions_dir, "vs_reference")
    dnadiff_fasta(pima_data, pima_data.dnadiff_prefix)

    parse_insertion_results(pima_data)
    make_finish_file(pima_data, pima_data.insertions_dir)

#TODO: Handle large indels that intersect with the mutations file correctly
def parse_insertion_results(pima_data: PimaData):
    """Read the results from dnadiff and populate the information into the Analysis class & final report"""
    # Figure out the %-identity with the reference sequence
    dnadiff_report = pima_data.dnadiff_prefix + '.report'
    command = ' '.join(['grep AvgIdentity', dnadiff_report,
                        '| head -1',
                        '| awk \'{print $2}\''])
    pima_data.reference_identity = float(
        print_and_run(pima_data, command)[0]
    )

    command = ' '.join(['grep AlignedBases', dnadiff_report,
                        '| head -1',
                        '| awk \'{sub("\\\\(.*", "", $2);print $2}\''])
    pima_data.reference_aligned_bases = int(
        print_and_run(pima_data, command)[0]
    ) * 100

    # Figure out the % of the query 
    command = ' '.join(['grep AlignedBases', dnadiff_report,
                '| head -1',
                '| awk \'{print $3}\'',
                '| sed \'s/(.*//\''])
    pima_data.query_aligned_bases = int(
        print_and_run(pima_data, command)[0]
    ) * 100
    
    # Figure out the aligned fraction
    pima_data.reference_aligned_fraction = round(pima_data.reference_aligned_bases / pima_data.reference_size,2)
    pima_data.query_aligned_fraction = round(pima_data.query_aligned_bases / pima_data.genome_size,2)
    # Give a warning if the distance is to high
    if pima_data.reference_identity <= pima_data.reference_identity_min:
        warning = f"Identity with the reference ({round(pima_data.reference_identity,2)}%) is less than {round(pima_data.reference_identity_min,2)}%"
        add_warning(pima_data, warning)
        pima_data.alignment_notes = pd.concat([pima_data.alignment_notes, pd.Series(warning, dtype='object')])
    
    if pima_data.reference_aligned_fraction <= pima_data.reference_alignment_min:
        warning = f"Fraction of the reference alignment ({pima_data.reference_aligned_fraction}%) is less than {pima_data.reference_alignment_min}%."
        add_warning(pima_data, warning)
        pima_data.alignment_notes = pd.concat([pima_data.alignment_notes, pd.Series(warning, dtype='object')])

    if pima_data.query_aligned_fraction <= pima_data.query_alignment_min:
        warning = f"Fraction of the query alignment ({pima_data.query_aligned_fraction}%) is less than {pima_data.query_alignment_min}%."
        add_warning(pima_data, warning)
        pima_data.alignment_notes = pd.concat([pima_data.alignment_notes, pd.Series(warning, dtype='object')])

    # Pull out the aligned regions of the two genomes
    print_and_log(
        pima_data,
        'Finding reference and query specific insertions', 
        pima_data.sub_process_verbosity, 
        pima_data.sub_process_color,
    )
    pima_data.one_coords = pima_data.dnadiff_prefix + '.1coords'
    reference_aligned_bed =  os.path.join(pima_data.insertions_dir, 'reference_aligned.bed')
    genome_aligned_bed = os.path.join(pima_data.insertions_dir, 'genome_aligned.bed')
    command = " ".join(
        [
            'awk \'{OFS = "\\t"; if ($2 < $1){t = $2; $2 = $1; $1 = t} print $12,$1,$2}\'', 
            pima_data.one_coords,
            '| sort -k 1,1 -k 2,2n',
            '>', 
            reference_aligned_bed,
        ]
    )
    print_and_run(pima_data, command)
    command = " ".join(
        [
            'awk \'{OFS = "\\t"; if ($4 < $3){t = $4; $4 = $3; $3 = t} print $13,$3,$4}\'', 
            pima_data.one_coords,
            '| sort -k 1,1 -k 2,2n',
            '>', 
            genome_aligned_bed,
        ]
    )
    print_and_run(pima_data, command)

    # Find the unaligned regions.  These are our insertions and deletions
    pima_data.reference_sizes = os.path.join(pima_data.insertions_dir, 'reference.sizes')
    reference_insertions_bed = os.path.join(pima_data.insertions_dir, 'reference_insertions.bed')
    genome_sizes = os.path.join(pima_data.insertions_dir, 'genome.sizes')
    genome_insertions_bed = os.path.join(pima_data.insertions_dir, 'genome_insertions.bed')

    command = " ".join(
        [
            'faidx -i chromsizes', 
            pima_data.reference_fasta, 
            ' | sort -k 1,1 -k 2,2n >', 
            pima_data.reference_sizes,
        ]
    )
    print_and_run(pima_data, command)
    command = " ".join(
        [
            'faidx -i chromsizes', 
            pima_data.genome_fasta, 
            ' | sort -k1,1 -k2,2n >', 
            genome_sizes,
        ]
    )
    print_and_run(pima_data, command)

    command = " ".join(
        [
            'bedtools complement',
            '-i', reference_aligned_bed,
            '-g', pima_data.reference_sizes,
            '| awk \'($3 - $2 >= 25){OFS = "\\t";print $1,$2,$3,($3 - $2)}\'',
            '>', 
            reference_insertions_bed,
        ]
    )
    print_and_run(pima_data, command)
    pima_data.reference_contig_order = [contig for contig in pd.read_csv(pima_data.reference_sizes, sep="\t", header = None).iloc[:,0]]

    # Check how much each contig in the query aligns with the reference
    query_alignment_stats = []
    with open(genome_sizes, 'r', encoding="UTF-8") as fp:
        for line in fp:
            contig, size = line.strip().split()
            command = " ".join(
                [
                    "grep", 
                    contig, 
                    genome_aligned_bed,
                    "| awk \'BEGIN {sum = 0} {sum+=($3-$2)} END {print sum}\'",
                ]
            )
            contig_bases_align = print_and_run(pima_data, command)[0]
            contig_perc_align = round(int(contig_bases_align) / int(size) * 100,1)
            query_alignment_stats.append(
                {
                    'Contig': contig,
                    'Size (bp)': size,
                    'Bases Aligned to Ref': contig_bases_align,
                    'Perc Align': contig_perc_align
                }
            )
    pima_data.query_alignment_stats = pd.DataFrame(query_alignment_stats)

    # Check how much each reference element aligns with the query
    reference_alignment_stats = []
    with open(pima_data.reference_sizes, 'r', encoding="UTF-8") as fp:
        for line in fp:
            contig, size = line.strip().split()
            command = " ".join(
                [
                    "grep", 
                    contig, 
                    reference_aligned_bed, 
                    "| awk \'BEGIN {sum = 0} {sum+=($3-$2)} END {print sum}\'",
                ]
            )
            contig_bases_align = print_and_run(pima_data, command)[0]
            contig_perc_align = round(int(contig_bases_align) / int(size) * 100,1)
            reference_alignment_stats.append(
                {
                    'Contig': contig,
                    'Size (bp)': size,
                    'Bases Aligned to Query': contig_bases_align,
                    'Perc Align': contig_perc_align
                }
            )
    pima_data.reference_alignment_stats = pd.DataFrame(reference_alignment_stats)

    # There may or may not be any insertions seen
    try:
        pima_data.reference_insertions = pd.read_csv(filepath_or_buffer = reference_insertions_bed, sep = '\t', header = None,
                                                     names=["contig_name", "start", "stop", "length"])
    except Exception:
        pima_data.reference_insertions = pd.DataFrame()

    command = " ".join(
        [
            'bedtools complement',
            '-i', genome_aligned_bed,
            '-g', genome_sizes,
            '| awk \'($3 - $2 >= 25){OFS = "\\t";print $1,$2,$3,($3 - $2)}\'',
            '>', 
            genome_insertions_bed,
        ]
    )
    print_and_run(pima_data, command)

    # There may or may not be any insertions seen
    try:
        pima_data.genome_insertions = pd.read_csv(filepath_or_buffer = genome_insertions_bed, sep = '\t', header = None,
                                                  names=["contig_name", "start", "stop", "length"])
    except FileNotFoundError:
        pima_data.genome_insertions = pd.DataFrame()
    except pd.errors.EmptyDataError:
        pima_data.genome_insertions = pd.DataFrame()
    except Exception as e:
        error_out(
            pima_data, f"Unexpected exception when processing insertions relative to the reference genome: {e}"
        )
    # Report the large indels
    ## About 40 get reported / page currently -> will restrict this table to 1 pg if there are a lot of large INDELs
    if pima_data.genome_insertions.shape[0] > 40:
        genome_insertions_csv = os.path.join(pima_data.insertions_dir, "query_insertions.csv")
        warning = (f"There are {pima_data.genome_insertions.shape[0]} large (>25 bp) query-based insertions detected. "
                    f"40 longest insertions reported. A complete list can be accessed at: \n"
                    f"{genome_insertions_csv}"
        )
        add_warning(pima_data, warning)
        pima_data.genome_insertions.to_csv(
            genome_insertions_csv, 
            header = ["Query Contig", "Start", "End", "Size (bp)"], 
            index = False,
        )
        pima_data.large_indel_notes = pd.concat([pima_data.large_indel_notes, pd.Series(warning, dtype='object')])
        pima_data.genome_insertions = pima_data.genome_insertions.sort_values('length', ascending = False).head(40).sort_index()

    if pima_data.reference_insertions.shape[0] > 40:
        reference_insertions_csv = os.path.join(pima_data.insertions_dir, "reference_insertions.csv")
        warning = (f"There are {pima_data.reference_insertions.shape[0]} large (>25 bp) reference insertions detected. "
                    f"40 longest insertions reported. A complete list can be accessed at: \n"
                    f"{reference_insertions_csv}"
        )
        pima_data.reference_insertions.to_csv(
            reference_insertions_csv, 
            header = ["Reference Contig", "Start", "End", "Size (bp)"], 
            index = False,
        )
        add_warning(pima_data, warning)
        pima_data.large_indel_notes = pd.concat([pima_data.large_indel_notes, pd.Series(warning, dtype='object')])
        pima_data.reference_insertions = pima_data.reference_insertions.sort_values('length', ascending = False).head(40).sort_index()

    pima_data.large_indels['Reference insertions'] = pima_data.reference_insertions
    pima_data.large_indels['Query insertions'] = pima_data.genome_insertions

    # Also pull in the number of SNPs and small indels
    genome_snps = pima_data.dnadiff_prefix + '.snps'
    try:
        snps = pd.read_csv(filepath_or_buffer = genome_snps, sep = '\t', header = None)
        pima_data.small_indels = snps.loc[(snps.iloc[:, 1] == '.') | (snps.iloc[:, 2] == '.'), :]
        pima_data.snps = snps.loc[(snps.iloc[:, 1] != '.') & (snps.iloc[:, 2] != '.'), :]
    except FileNotFoundError:
        snps = pd.DataFrame()
    except pd.errors.EmptyDataError:
        snps = pd.DataFrame()
    except Exception as e:
        error_out(
            pima_data, f"Unexpected exception when processing SNPs and small indels from the reference to assembly genome alignment: {e}"
        )

    #TODO: rewrite this to handle the large INDELs & report them correctly
    ## here the large indels are >25 bp, should I filter the varscan variants so that only <25 indels are reported there?
    ### handling this with the call_amr_mutations below
    #### We can't handle large insertions right now!!

    # See if any mutation regions intersect with deletions
    ## Isn't this superfluous? Read mapping will identify INDELs
    # if pima_data.mutation_region_bed:
    #     amr_deletion_bed = os.path.join(pima_data.insertions_dir, 'amr_deletions.bed')
    #     command = " ".join(
    #         [
    #             'bedtools intersect',
    #             '-nonamecheck',
    #             '-a', pima_data.mutation_region_bed,
    #             '-b', reference_insertions_bed,
    #             '>', 
    #             amr_deletion_bed,
    #         ]
    #     )
    #     print_and_run(pima_data, command)

    #     # There may be no deletions, let's see
    #     try:
    #         pima_data.amr_deletions = pd.read_csv(filepath_or_buffer = amr_deletion_bed, sep = '\t', header = None)
    #     except FileNotFoundError:
    #         pima_data.amr_deletions = pd.DataFrame()
    #     except pd.errors.EmptyDataError:
    #         pima_data.amr_deletions = pd.DataFrame()
    #     except Exception as e:
    #         error_out(
    #             pima_data, f"Unexpected exception when processing deletions within locations listed in the mutations_regions.bed file: {e}"
    #         )  

    #     if pima_data.amr_deletions.shape[0] > 0:
    #         pima_data.amr_deletions.columns = ['contig', 'start', 'stop', 'name', 'type', 'drug', 'note']
    #         #only report regions within the mutation_bed file that are 'any' or 'large-deletion'
    #         pima_data.amr_deletions = pima_data.amr_deletions.loc[pima_data.amr_deletions['type'].isin(['large-deletion', 'any']), :]
    
    pima_data.did_call_large_indels = True

def quast_genome(pima_data: PimaData):
    
    if not pima_data.genome_fasta:
        return
    

    def parse_quast_output(pima_data: PimaData):
        # Look at the quast TSV file
        quast_report_tsv = os.path.join(pima_data.quast_dir, 'report.tsv')
        quast_report = pd.read_csv(filepath_or_buffer=quast_report_tsv, header = 0, index_col = 0, sep = '\t')
        
        pima_data.quast_mismatches = int(float(quast_report.loc['# mismatches per 100 kbp', :].iloc[0]) * \
                                    (float(quast_report.loc['Total length (>= 0 bp)', :].iloc[0]) / 100000.))
        pima_data.quast_indels = int(float(quast_report.loc['# indels per 100 kbp', :].iloc[0]) * \
                                (float(quast_report.loc['Total length (>= 0 bp)', :].iloc[0]) / 100000.))
    
    print_and_log(
        pima_data,
        'Running Quast of the genome vs. the reference sequence',
        pima_data.main_process_verbosity, 
        pima_data.main_process_color,
    ) 

    pima_data.quast_dir = os.path.join(pima_data.output_dir, 'quast')
    if find_checkpoint(pima_data, pima_data.quast_dir):
        print_and_log(
            pima_data,
            'Reading previously completed quast report',
            pima_data.main_process_verbosity, 
            pima_data.main_process_color,
        ) 
        parse_quast_output(pima_data)
        return

    # Quast output directory
    os.makedirs(pima_data.quast_dir)
    make_start_file(pima_data, pima_data.quast_dir)
    
    # Run Quast
    quast_stdout, quast_stderr = std_files(os.path.join(pima_data.quast_dir, 'quast'))
    command = " ".join(
        [
            'quast.py',
            '-R', pima_data.reference_fasta,
            '-o', pima_data.quast_dir,
            pima_data.genome_fasta,
            '1>', quast_stdout, 
            '2>', quast_stderr,
        ]
    )
    print_and_run(pima_data, command)

    # Collect results (refactored as function to enable checkpointing)
    parse_quast_output(pima_data)
    make_finish_file(pima_data, pima_data.quast_dir)

def call_amr_mutations(pima_data: PimaData):
    #TODO: Handle multiple reads as input (e.g. still use ONT even if Illumina provided)
    
    if pima_data.no_assembly and not pima_data.genome_fasta:
        warning = (
            f"Since an assembly was not generated or provided, we do not check and see how similar the input organism is to the reference. "
            f"If it is < {pima_data.reference_identity_min}% then amr mutation identification can be misleading and will include many false positive reports. "
        )
        add_warning(
            pima_data,
            warning,
        )
    
    print_and_log(
        pima_data,
        'Calling AMR mutations', 
        pima_data.main_process_verbosity, 
        pima_data.main_process_color,
    )

    # Make the directory for small mutation related work
    pima_data.mutations_dir = os.path.join(pima_data.output_dir, 'mutations')

    # Checkpoint; parse results and continue
    if find_checkpoint(pima_data, pima_data.mutations_dir):

        varscan_prefix = os.path.join(pima_data.mutations_dir, 'varscan')
        print_and_log(
            pima_data,
            'Reading previously identified AMR mutations', 
            pima_data.main_process_verbosity, 
            pima_data.main_process_color,
        )
        #need to clean-up after the last run if this one completes successfully
        if pima_data.illumina_fastq:
            kind_of_reads = "Illumina"
            reference_mapping_unfilt_illumina_bam = os.path.join(pima_data.mutations_dir, 'reference_mapping_illumina_unfilt.bam')
            pima_data.reference_mapping_illumina_bam = os.path.join(pima_data.mutations_dir, 'reference_mapping_illumina.bam')
            reference_mapping_mpileup = os.path.join(pima_data.mutations_dir, 'reference_mapping_illumina.mpileup')
            pima_data.files_to_clean.extend([reference_mapping_unfilt_illumina_bam, pima_data.reference_mapping_illumina_bam, reference_mapping_mpileup])
        else:
            kind_of_reads = "ONT"
            pima_data.reference_mapping_ont_bam = os.path.join(pima_data.mutations_dir, 'reference_mapping_ont.bam')
            reference_mapping_mpileup = os.path.join(pima_data.mutations_dir, 'reference_mapping_ont.mpileup')
            pima_data.files_to_clean.extend([reference_mapping_mpileup, pima_data.reference_mapping_ont_bam])

        merged_hits_tsv = os.path.join(pima_data.mutations_dir, "intersect_amr-regions_variants_deletions.tsv")
        pima_data.mutations_read_type = kind_of_reads
        merged_hits = pd.read_csv(
            merged_hits_tsv,
            sep="\t",
            header=0,
        )
        
        pima_data.amr_mutations = merged_hits
        pima_data.did_call_mutations = True
        return
    
    os.makedirs(pima_data.mutations_dir)
    make_start_file(pima_data, pima_data.mutations_dir)

    # Map the reads to the reference sequence
    ## Do we want to map both ONT & Illumina when they are provided - allows later drawing of coverage profiles
    if pima_data.illumina_fastq:
        reference_mapping_unfilt_illumina_bam = os.path.join(pima_data.mutations_dir, 'reference_mapping_illumina_unfilt.bam')
        pima_data.reference_mapping_illumina_bam = os.path.join(pima_data.mutations_dir, 'reference_mapping_illumina.bam')
        minimap_and_sort(
            pima_data,
            pima_data.reference_fasta,
            reference_mapping_unfilt_illumina_bam,
            pima_data.illumina_fastq,
            ont = False,
        )
        filter_bam(
            pima_data,
            inbam = reference_mapping_unfilt_illumina_bam,
            outbam = pima_data.reference_mapping_illumina_bam,
            F = "0x0100",
            q = "30",
        )
        pima_data.files_to_clean.extend([reference_mapping_unfilt_illumina_bam, pima_data.reference_mapping_illumina_bam])
        kind_of_reads = 'Illumina'
    
        # Make an mpileup file out of the reads
        reference_mapping_mpileup = os.path.join(pima_data.mutations_dir, 'reference_mapping_illumina.mpileup')
        pima_data.files_to_clean.append(reference_mapping_mpileup)
        mpileup_bam(
            pima_data,
            reference_genome = pima_data.reference_fasta,
            bam = pima_data.reference_mapping_illumina_bam, 
            mpileup = reference_mapping_mpileup, 
            output_dir = pima_data.mutations_dir,
        )
    else: #We'll generate the ONT mapping information during the circos step to build the coverage, but not to call mutations
        pima_data.reference_mapping_ont_bam = os.path.join(pima_data.mutations_dir, 'reference_mapping_ont.bam')
        pima_data.files_to_clean.append(pima_data.reference_mapping_ont_bam)
        minimap_and_sort(
            pima_data,
            pima_data.reference_fasta, 
            pima_data.reference_mapping_ont_bam,
            pima_data.ont_fastq, 
            ont = True,
        )
        kind_of_reads = 'ONT'

        # Make an mpileup file out of the reads
        reference_mapping_mpileup = os.path.join(pima_data.mutations_dir, 'reference_mapping_ont.mpileup')
        pima_data.files_to_clean.append(reference_mapping_mpileup)
        mpileup_bam(
            pima_data,
            pima_data.reference_fasta, 
            pima_data.reference_mapping_ont_bam, 
            reference_mapping_mpileup, 
            pima_data.mutations_dir
        )
    pima_data.mutations_read_type = kind_of_reads

    ##rplV calling has issues due to repetitive regions where the INDELs occur 
    ##     - true variants are often at a lower percentage of population (< 30%)
    ## ONT data is actually better at calling these long-non-homopolymer mutations than illumina in these cases
    ## For now, will enable indel calling with ONT-only & filter by INDEL type
    ### I don't know how this is reported when both ONT & Illumina data are provided? I suspect only Illumina data are used!!

    # Now call and filter variants with Varscan and filter
    varscan_raw_prefix = os.path.join(pima_data.mutations_dir, 'varscan_raw')
    varscan_mpileup(
        pima_data,
        reference_mapping_mpileup, 
        varscan_raw_prefix, 
        'snp', 
        0.8,
    )

    varscan_mpileup(
        pima_data,
        reference_mapping_mpileup, 
        varscan_raw_prefix, 
        'indel', 
        1./4., #changed from 1/3 to be more sensitive to long repetitive indels
    ) 
        
    varscan_prefix = os.path.join(pima_data.mutations_dir, 'varscan')
    filter_varscan(pima_data,
                   varscan_raw_prefix, 
                   varscan_prefix, 
                   'snp', 
                   kind_of_reads,
    )

    filter_varscan(pima_data,
                   varscan_raw_prefix, 
                   varscan_prefix, 
                   'indel', 
                   kind_of_reads,
    )
        
    # And dump as a VCF
    pima_data.varscan_vcf = vcf_varscan(pima_data, varscan_prefix)
    verified_hits = intersect_variants_with_mutations(pima_data)
    verified_large_dels = intersect_dnadiff_large_dels(pima_data)

    #combine the two datasets
    ## Bug where variants in the "any" var_type get collapsed since they have the same start/stop, added 'loc' for the specific, but might not work? 
    merged_hits = pd.concat([verified_hits, verified_large_dels]).drop_duplicates(subset=['GE', 'start', 'stop', 'loc'], keep='first')
    merged_hits = pd.concat([verified_hits, verified_large_dels])
    merged_hits_tsv = os.path.join(pima_data.mutations_dir, "intersect_amr-regions_variants_deletions.tsv")
    merged_hits.to_csv(
        path_or_buf = merged_hits_tsv,
        sep="\t",
        header=True,
        index=False,
    )
    pima_data.amr_mutations = merged_hits
    pima_data.did_call_mutations = True
    make_finish_file(pima_data, pima_data.mutations_dir)

def intersect_variants_with_mutations(pima_data: PimaData):
    command = " ".join(
        [
            "bedtools intersect",
            "-a",
            pima_data.mutation_region_bed,
            "-b",
            pima_data.varscan_vcf,
            "-loj",
            "-nonamecheck",
        ]
    )
    intersect_bed = print_and_run(pima_data, command)
    intersect_df = pd.DataFrame(
        [x.split("\t") for x in intersect_bed],
        )
    #drop columns we don't care about
    intersect_df.drop([8,10,13,14,15,16,17], axis=1, inplace=True)
    df_column_names=["GE", "start", "stop", "region_name", "var_type", "amr_class", "classification_type", "note", "loc","ref", "var", "AF"]
    intersect_df = intersect_df.set_axis(df_column_names, axis=1) 

    #drop reference_regions that had no detected variants
    hits = intersect_df.query('loc != ["-1", None]')

    #make sure the intersect matches the correct kinds of mutations with the confirmed mutations
    cond1 = hits.query("var_type == 'large-indel' & var.str.len() > 25") #this maybe needs to be handled by the dnadiff codeblock
    cond2 = hits.query("var_type == 'snp' & var.str.len() == 1")
    cond3 = hits.query("var_type == 'indel' & 1 < var.str.len() <= 25")
    cond4 = hits.query("var_type == 'any'")
    verified_hits = pd.concat([cond1, cond2, cond3, cond4])
    #mutations file is organized by confirmed variants -> potential variants in key regions
    ## here we want to remove a variant that intersects both locations and only report the confirmed location
    ### note: we don't actually check to see if it is the same exact variant, just the same location?
    verified_hits.drop_duplicates(subset=['loc', 'ref', 'var'], keep='first', inplace=True)
    return verified_hits

def intersect_dnadiff_large_dels(pima_data: PimaData):
    ##check if dnadiff-based large deletions intersect our mutations bed file
    ### We can't pull out the insertions due to how the genome alignments are saved, would need to rework that section
    ### our >25bp deletions
    reference_insertions_bed = os.path.join(pima_data.insertions_dir, 'reference_insertions.bed')
    command = " ".join(
        [
            'bedtools intersect',
            '-nonamecheck',
            '-loj',
            '-a', pima_data.mutation_region_bed,
            '-b', reference_insertions_bed,
        ]
    )
    intersect_large_del = print_and_run(pima_data, command)
    large_del_df = pd.DataFrame(
        [x.split("\t") for x in intersect_large_del],
        )
    large_del_df.drop([8,10], axis=1, inplace=True)
    #Use same column names as the variant caller so we can merge & reconcile them
    df_column_names=["GE", "start", "stop", "region_name", "var_type", "amr_class", "classification_type", "note", "loc", "var"]
    large_del_df = large_del_df.set_axis(df_column_names, axis=1)

    #drop reference_regions that had no detected variants
    large_del_hits = large_del_df.query('loc != ["-1", None] & (var_type == "large-indel" | classification_type == "potential_confer_amr")')
    large_del_hits = large_del_hits.drop_duplicates(subset=['start', 'stop'], keep='first')
    large_del_hits['ref'] = "large-deletion"
    large_del_hits['AF'] = "detected-by-assembly"
    return large_del_hits

def varscan_mpileup(pima_data: PimaData, mpileup: str, prefix: str, var_type: Literal['snp', 'indel'], min_var_freq: float):

    print_and_log(
        pima_data,
        f"Calling {var_type} with varscan",
        pima_data.sub_process_verbosity, 
        pima_data.sub_process_color,
    )
    varscan_var = f"{prefix}.{var_type}"
    varscan_stderr = ''.join([prefix, '_', var_type, '.stderr'])

    command = " ".join(
        [
            'varscan',
            'pileup2' + var_type,
            mpileup,
            '-p-value 0.01',
            '--min-coverage 15',
            '--min-var-freq', str(min_var_freq),
            '--variants',
            '1>', varscan_var, 
            '2>', varscan_stderr,
        ]
    )
    print_and_run(pima_data, command)
    validate_file_and_size_or_error(
        pima_data,
        varscan_var, 
        " ".join(['Varscan', var_type, 'file']), 
        'cannot be found', 
        'is empty',
    )               

def filter_varscan(pima_data: PimaData, in_prefix: str, out_prefix: str, var_type: Literal['snp', 'indel'], read_type: Literal['Illumina', 'ONT']):
    
    print_and_log(
        pima_data,
        f"Filtering Varscan {var_type} calls",
        pima_data.sub_process_verbosity, 
        pima_data.sub_process_color,
    )
    varscan_raw_var = ''.join([in_prefix, '.',var_type])
    varscan_var = ''.join([out_prefix, '.', var_type])
    
    ## Add in that ONT detection for homopolymers is allowed if R10 data is used
    ### Changing the logic, lets call all mutations the same way, and change how we report them. Knowledgeable users can inspect these intermediate files
    ## Use this statement for all Illumina sets, and the ONT SNP mode, and R10_sup basecalled data
    #if (read_type == "Illumina") or (var_type == "snp") or (re.search(r"r1041.*sup.*", pima_data.ont_model)):
    command = " ".join(
        [
            'awk \'(NR > 1 && $9 == 2 && $5 + $6 >= 15)',
            '{OFS = "\\t";f = $6 / ($5 + $6); gsub(/.*\\//, "", $4);s = $4;gsub(/[+\\-]/, "", s);$7 = sprintf("%.2f%%", f * 100);',
            'min = 1 / log(length(s) + 2) / log(10) + 2/10;if(f > min){print}}\'', 
            varscan_raw_var,
            '1>' + varscan_var,
        ]
    )
        ##That complicated 'min' formula scales the minimum percent frequency required by the length of the variant
        ## Variants of length 1 (SNPs) require at least 59% of reads
        ## Variants of length 10 (INDELs) require 37% of reads
        ## Formula is: 1 / (ln(length+2) / ln(10)) + 0.2
        ## This formula converges on 20% (x -> inf, y -> 0.2); drops from 0.6 -> 0.3 by length=10

    # only use long-non homopolymer INDELs if ONT is used
    # elif (read_type == "ONT") and (var_type == "indel"):
    #     ## add the statement that the sequence cannot be a repetition of the same letter (avoid homopolymers)
    #     ## also required the INDEL to be at least 6 bp long (length(s) > 5)
    #     command = " ".join(
    #         [
    #             'awk \'(NR > 1 && $9 == 2 && $5 + $6 >= 15)',
    #             '{OFS = "\\t";f = $6 / ($5 + $6); gsub(/.*\\//, "", $4);s = $4;gsub(/[+\\-]/, "", s);$7 = sprintf("%.2f%%", f * 100);',
    #             'min = 1 / log(length(s) + 2) / log(10) + 2/10;if(f > min && length(s) > 5 && s~"[^"substr(s,1,1)"]"){print}}\'',
    #             varscan_raw_var,
    #             '1>' + varscan_var, 
    #             '2>' + varscan_var,
    #         ]
    #     )
    print_and_run(pima_data, command)
    validate_file_and_size_or_error(
        pima_data, varscan_var, 
        f"Filtered Varscan {var_type} file", "cannot be found", "is empty",
        )
    
def vcf_varscan(pima_data: PimaData, varscan_prefix):

    #Need to include allele frequency in the data structures
    print_and_log(
        pima_data,
        'Making VCF from varscan', 
        pima_data.sub_process_verbosity, 
        pima_data.sub_process_color,
    )

    varscan_snp, varscan_indel = [varscan_prefix + '.' + i for i in ['snp', 'indel']]
    snp_vcf, indel_vcf = [varscan_prefix + '_' + i + '.vcf' for i in ['snp', 'indel']]
    found_vcf = []
    varscan_vcf = varscan_prefix + '.vcf'

    if os.path.isfile(varscan_snp):
        print_and_log(
            pima_data,
            'Making SNP VCF', 
            pima_data.sub_process_verbosity, 
            pima_data.sub_process_color,
        )
        command = " ".join(
            [
                'awk \'{OFS = "\\t"; print $1,$2,".",$3,$4,-log($14),"PASS",".","GT","1|1",$6"/"$5+$6}\'',
                varscan_snp,
                '1>' + snp_vcf,
            ]
        )
        print_and_run(pima_data, command)
        found_vcf += [snp_vcf]

    if os.path.isfile(varscan_indel):
        print_and_log(pima_data,
                      'Making indel VCF', 
                      pima_data.sub_process_verbosity, 
                      pima_data.sub_process_color,
                    )
        command = " ".join(
            [
                'awk \'{OFS = "\\t"; print $1,$2,".",$3,$4,-log($14),"PASS",".","GT","1|1",$6"/"$5+$6}\'',
                varscan_indel,
                '1>' + indel_vcf,
            ]
        )
        print_and_run(pima_data, command)
        found_vcf += [indel_vcf]

    command = " ".join(
        [
            'cat', 
            " ".join(found_vcf),
            '| sort -k 1,1 -k 2n,2n',
            '| awk \'BEGIN{OFS = "\\t";print "##fileformat=VCFv4.2";',
            'print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\tAF/DP"}{print}\'',
            '1>' + varscan_vcf,
        ]
    )
    print_and_run(pima_data, command)
    
    return(varscan_vcf)
