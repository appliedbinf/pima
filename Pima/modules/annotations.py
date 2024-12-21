import csv
import re
import os
import subprocess

import pandas as pd

from Pima.pima_data import PimaData
from Pima.utils.settings import Settings
from Pima.utils.utils import (
    print_and_log,
    print_and_run,
    validate_utility,
    validate_file_and_size,
    find_checkpoint,
    make_start_file,
    make_finish_file,
    std_files,
    error_out,
)


def validate_features(pima_data: PimaData, settings: Settings):
    # skip conditions
    if pima_data.only_assemble:
        return

    if pima_data.feature_fastas is None:
        pima_data.feature_fastas = []

    if not pima_data.no_amr:
        pima_data.feature_fastas.append(pima_data.amr_database)
        pima_data.feature_colors.append(settings.amr_default_color)
        if validate_file_and_size(pima_data, settings.amr_gene_drug_tsv):
            pima_data.amr_gene_drug = pd.read_csv(
                settings.amr_gene_drug_tsv,
                index_col=None,
                sep="\t",
                quoting=csv.QUOTE_NONE,
                header=None,
            )
            pima_data.drug_categories = pima_data.amr_gene_drug.iloc[:, 1].unique()

    if not pima_data.no_inc:
        pima_data.feature_fastas.append(pima_data.inc_database)
        pima_data.feature_colors.append(settings.inc_default_color)

    if len(pima_data.feature_fastas) == 0:
        return

    if not pima_data.will_have_genome_fasta:
        return

    print_and_log(
        pima_data,
        "Validating feature sets",
        pima_data.main_process_verbosity,
        pima_data.main_process_color,
    )

    for feature_fasta in pima_data.feature_fastas:
        if not validate_file_and_size(pima_data, feature_fasta):
            # See if the missing database can be downloaded
            if feature_fasta in settings.included_databases:
                if not pima_data.download:
                    pima_data.errors.append(
                        f"Can't find feature database {feature_fasta} or it is empty.  Try --download?"
                    )
            else:
                pima_data.errors.append(f"Can't find feature database {feature_fasta}")


def validate_blast(pima_data: PimaData):
    # skip conditions
    if pima_data.only_assemble:
        return

    if len(pima_data.feature_fastas) == 0:
        return

    if not pima_data.will_have_genome_fasta:
        return

    print_and_log(
        pima_data,
        "Validating blast utilities",
        pima_data.main_process_verbosity,
        pima_data.main_process_color,
    )

    for utility in ["makeblastdb", "blastn", "bedtools"]:
        if validate_utility(pima_data, utility, f"{utility} isn't on the PATH."):
            command = utility + " -version"
            pima_data.versions[utility] = re.search(
                r"[0-9]+\.[0-9.]+", print_and_run(pima_data, command)[0]
            ).group(0)
    pima_data.analysis.append(["blast_feature_sets", pima_data])

    
def blast_feature_sets(pima_data: PimaData):
    """Find genes within both 'amr' and 'inc' databases within the assembly
    
    Generates a dictionary of dataframes, 1 dataframe for amr and 1 for inc
    """

    print_and_log(
        pima_data,
        "BLASTing feature sets",
        pima_data.main_process_verbosity,
        pima_data.main_process_color,
    )

    # Keep track of feature hits for reporting
    pima_data.features_dir = os.path.join(pima_data.output_dir, "features")

    # Check if results already exist
    if find_checkpoint(pima_data, pima_data.features_dir):
        print_and_log(
            pima_data,
            "BLASTing features had previously been run and finished successfully",
            pima_data.main_process_verbosity,
            pima_data.main_process_color,
        )
        pima_data.did_blast_feature_sets = True
        found_feature_dirs = [
            feature_dir.path
            for feature_dir in os.scandir(pima_data.features_dir)
            if feature_dir.is_dir()
        ]
        for feature_dir in found_feature_dirs:
            feature_name = os.path.basename(feature_dir)
            best_bed = os.path.join(feature_dir, "best.bed")
            parse_blast_features(pima_data, best_bed, feature_name)
        return

    os.makedirs(pima_data.features_dir)
    make_start_file(pima_data, pima_data.features_dir)

    # Make a blast database of the genome
    make_blast_database(pima_data, pima_data.genome_fasta)

    for feature_number in range(len(pima_data.feature_fastas)):
        feature_fasta = pima_data.feature_fastas[feature_number]
        feature_name = re.sub(r"\.f.*", "", os.path.basename(feature_fasta))
        feature_dir = os.path.join(pima_data.features_dir, feature_name)
        blast_features(pima_data, feature_fasta, feature_dir, feature_name)
        pima_data.feature_dirs += [feature_dir]
        pima_data.feature_names += [feature_name]

    pima_data.did_blast_feature_sets = True
    make_finish_file(pima_data, pima_data.features_dir)


def make_blast_database(pima_data: PimaData, database_fasta: str):

    if os.path.isfile(f"{database_fasta}.nin"):
        command = " ".join(
            [
                'blastdbcmd -info -db',
                database_fasta,
            ]
        )
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        if result.returncode == 0:
            return

    print_and_log(
        pima_data,
        "Making a BLAST database for " + database_fasta,
        pima_data.sub_process_verbosity,
        pima_data.sub_process_color,
    )
    std_prefix = re.sub(r"\.[^.]*$", "", database_fasta)
    stdout_file, stderr_file = std_files(std_prefix)
    command = " ".join(
        [
            "makeblastdb -in",
            database_fasta,
            "-dbtype nucl -parse_seqids",
            "1>", stdout_file,
            "2>", stderr_file,
        ]
    )
    print_and_run(pima_data, command)


def blast_features(
    pima_data: PimaData, feature_fasta: str, feature_dir: str, feature_name: str
):
    # Make a directory for the new features
    os.makedirs(feature_dir)

    # BLASTn the feature set
    blast_output = os.path.join(feature_dir, "blast_output.tsv")
    print_and_log(
        pima_data,
        "BLASTing features against the assembly",
        pima_data.sub_process_verbosity,
        pima_data.sub_process_color,
    )
    blastn_stdout, blastn_stderr = std_files(os.path.join(feature_dir, "blastn"))
    command = " ".join(
        [
            "blastn -db",
            pima_data.genome_fasta,
            "-query",
            feature_fasta,
            "-perc_identity 95.0",
            '-outfmt "6',
            'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident qlen"',
            "-evalue 1e-10 -out ",
            blast_output,
            "1>", blastn_stdout,
            "2>", blastn_stderr,
        ]
    )
    print_and_run(pima_data, command)

    # Clean up the results into a handy BED file
    print_and_log(
        pima_data,
        "Converting feature hits to BED",
        pima_data.sub_process_verbosity,
        pima_data.sub_process_color,
    )
    all_bed = os.path.join(feature_dir, "all.bed")
    command = " ".join(
        [
            "cat",
            blast_output,
            "| awk -F '\\t' '($3 >= 95) && ($4 / $14 >= .90){OFS = \"\\t\";"
            + 'print $2,($9 < $10 ? $9 : $10),($9 < $10 ? $10 : $9),$1,$3/100,($9 < $10 ? "+" : "-")}\'',
            "| sort -k 1,1 -k 2,2n >",
            all_bed,
        ]
    )
    print_and_run(pima_data, command)

    # Make clusters of hits
    print_and_log(
        pima_data,
        "Clustering feature hits",
        pima_data.sub_process_verbosity,
        pima_data.sub_process_color,
    )
    merge_bed = os.path.join(feature_dir, "merge.bed")
    _, merge_stderr = std_files(os.path.join(feature_dir, "bedtools_merge"))
    command = " ".join(
        ["bedtools merge -d -30 -i", all_bed, "1>", merge_bed, "2>", merge_stderr]
    )
    print_and_run(pima_data, command)

    # Pick the best hit for each cluster
    print_and_log(
        pima_data,
        "Finding the best hit for each feature cluster",
        pima_data.sub_process_verbosity,
        pima_data.sub_process_color,
    )
    best_bed = os.path.join(feature_dir, "best.bed")
    command = " ".join(
        [
            "bedtools intersect",
            "-a",
            all_bed,
            "-b",
            merge_bed,
            "-f .9 -F .9 -wao",
            "| awk '$7 != \".\"'",
            '| awk \'{OFS="\\t";locus=$7"\\t"$8"\\t"$9; if($5 > s[locus]){s[locus]=$5;id = sprintf("%.3f", $5); b[locus] = $1"\\t"$2"\\t"$3"\\t"$4"\\t"id"\\t"$6}}',
            "END{for(i in b){print b[i]}}'",
            "| sort -k 1,1 -k2,2n",
            ">" + best_bed,
        ]
    )
    print_and_run(pima_data, command)
    parse_blast_features(pima_data, best_bed, feature_name)


def parse_blast_features(pima_data: PimaData, best_bed: str, feature_name: str):
    # Keep the feature hits for later drawing.  It may be empty, i.e., no feature hits
    try:
        best = pd.read_csv(filepath_or_buffer=best_bed, sep="\t", header=None)
    except FileNotFoundError:
        best = pd.DataFrame()
    except pd.errors.EmptyDataError:
        best = pd.DataFrame()
    except Exception as e:
        error_out(
            pima_data, f"Unexpected exception when processing BLAST features: {e}"
        )

    pima_data.feature_hits[feature_name] = best
