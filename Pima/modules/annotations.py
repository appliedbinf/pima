import csv
import re
import os
import subprocess
import shutil
import pandas as pd
import pybedtools as bd
from pyfaidx import Fasta

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

def validate_virulence(pima_data:PimaData, settings: Settings):
    if not pima_data.will_have_genome_fasta:
        return
    if pima_data.only_assemble:
        return

    virulence_fasta = None
    if pima_data.organism_dir:
        candidate = os.path.join(pima_data.organism_dir, "virulence_genes.fasta")
        if os.path.isfile(candidate):
            virulence_fasta = candidate

    if virulence_fasta is None and pima_data.organism == "Bacillus_anthracis":
        if not validate_file_and_size(pima_data, settings.ba_virulence_genes):
            pima_data.errors.append("Can't find Bacillus anthracis virulence gene database")
            return
        virulence_fasta = settings.ba_virulence_genes

    if virulence_fasta is None:
        return

    pima_data.ba_virulence_genes = virulence_fasta


def validate_inc(pima_data:PimaData, settings: Settings):
    if not pima_data.will_have_genome_fasta:
        return
    
    if pima_data.only_assemble:
        return
    
    if not pima_data.no_inc:
        pima_data.inc_database = settings.inc_database

    if not validate_file_and_size(pima_data, pima_data.inc_database):
        pima_data.errors.append(f"Can't find feature database {pima_data.inc_database}, this indicates a problem with the pima install as this dataset should be packaged")
    

def validate_resfinder(pima_data: PimaData, settings: Settings):
    # skip conditions
    if pima_data.only_assemble:
        return

    if not pima_data.will_have_genome_fasta:
        return
    
    if not pima_data.no_amr:
        pima_data.amr_database = settings.amr_database
        if validate_file_and_size(pima_data, settings.amr_gene_drug_tsv):
            pima_data.amr_gene_drug = pd.read_csv(
                settings.amr_gene_drug_tsv,
                index_col=None,
                sep="\t",
                quoting=csv.QUOTE_NONE,
                header=None,
            )
            pima_data.drug_categories = pima_data.amr_gene_drug.iloc[:, 1].unique()
        else:
            pima_data.errors.append(f"Can't find resource {pima_data.amr_gene_drug}, this indicates a problem with the pima install as this dataset should be packaged")


        if not validate_file_and_size(pima_data, pima_data.amr_database):
            pima_data.errors.append(f"Can't find feature database {pima_data.amr_database}, this indicates a problem with the pima install as this dataset should be packaged")


def validate_blast(pima_data: PimaData, settings: Settings):
    """
    Check for all tools necessary to run the resfinder, plasmidfinder, virulence screens
    """
    # skip conditions
    if pima_data.only_assemble:
        return

    if not pima_data.will_have_genome_fasta:
        return

    if not any([pima_data.amr_database, pima_data.inc_database, pima_data.ba_virulence_genes]):
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
    pima_data.analysis.append(["blast_feature_sets", pima_data, settings])


def validate_amrfinder(pima_data: PimaData, settings: Settings):
    """Check if AMRFinder and it's database are available"""

    if not pima_data.amrfinder:
        return

    if pima_data.no_amr:
        print_and_log(
            pima_data,
            "Running AMRFinder even though you requested '--no-amr'",
            pima_data.sub_process_verbosity,
            pima_data.sub_process_color,
        )

    print_and_log(
        pima_data,
        "Validating AMRFinder",
        pima_data.main_process_verbosity,
        pima_data.main_process_color,
    )    

    if validate_utility(pima_data, "amrfinder", "amrfinder isn't on the PATH."):
        #check tool
        command = "amrfinder --version"
        pima_data.versions['amrfinder'] = print_and_run(pima_data, command)[0]
        pima_data.analysis.append(["run_amrfinder", pima_data, settings])

        command = f"amrfinder -d {settings.amrfinder_database}/latest --database_version"
        #database: settings.amrfinder_database

        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        if result.returncode == 0:
            result = [x for x in result.stdout.split("\n") if re.search("Database version", x)]
            version = result[0].split(": ")[1]
            pima_data.versions['amrfinder database'] = version

        else:
            message = "AMRFinder was unable to find it's database. Try running 'pima --download' to fix"
            error_out(pima_data, message)


def blast_feature_sets(pima_data: PimaData, settings: Settings):
    """Find genes within 'amr', 'inc', 'virulence' databases within the assembly
    
    Generates a separate dataframes saved within pima_data
    """

    #split into separate results: resfinder amr, plasmidfinder inc, virulence
    print_and_log(
        pima_data,
        "BLASTing feature sets",
        pima_data.main_process_verbosity,
        pima_data.main_process_color,
    )

    # Keep track of feature hits for reporting
    features_dir = os.path.join(pima_data.output_dir, "features")
    
    #AMRfinder can make this dir
    os.makedirs(features_dir, exist_ok=True)
    make_start_file(pima_data, features_dir)

    # Make a blast database of the genome
    make_blast_database(pima_data, pima_data.genome_fasta)

    #Remove loop & run each independently
    #Resfinder Blast
    if pima_data.amr_database: #set during validation, otherwise None
        amr_dir = os.path.join(features_dir, 'amr')
        #delete if re-running
        if os.path.exists(amr_dir):
            shutil.rmtree(amr_dir)
        blast_features(pima_data, pima_data.amr_database, amr_dir, 'amr', 90.)

    #Inc/plasmidfinder Blast
    if pima_data.inc_database:
        inc_dir = os.path.join(features_dir, 'inc')
        if os.path.exists(inc_dir):
            shutil.rmtree(inc_dir)        
        blast_features(pima_data, pima_data.inc_database, inc_dir, 'inc', 95.)

    #Virulence Blast
    if pima_data.ba_virulence_genes:
        vir_dir = os.path.join(features_dir, 'ba_virulence_genes')
        if os.path.exists(vir_dir):
            shutil.rmtree(vir_dir)  
        blast_features(pima_data, pima_data.ba_virulence_genes, vir_dir, 'ba_virulence_genes', 90.)

    overlap_features(pima_data)
    pima_data.did_blast_feature_sets = True
    make_finish_file(pima_data, features_dir)


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
    pima_data: PimaData, feature_fasta: str, feature_dir: str, feature_name: str, percent_identity: float,
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
            f"| awk -F '\\t' '($3 >= {str(percent_identity)}) && ($4 / $14 >= .90){{OFS = \"\\t\";"
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
    # INC groups might get dropped here if they are on opposite strands on overlapping genomic positions
    # Throw a warning if the merged bed file has data and the best bed file is blank
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

    #TODO:Need to rework this - only triggers if every single region gets dropped, would not correctly report some regions geting dropped
    if os.path.getsize(best_bed) == 0 and os.path.getsize(merge_bed) > 0:
        note =  (
                    "PIMA may not be reporting true features due to larger than expected overlap. "
                    f"Check the intermediate files: {blast_output} and {all_bed}."
                )
        print_and_log(
            pima_data,
            note,
            pima_data.warning_verbosity,
            pima_data.warning_color,
        )
        pima_data.annotation_notes = pd.concat([pima_data.annotation_notes, pd.Series(note, dtype='object')])

    parse_blast_features(pima_data, best_bed, feature_name)


def parse_blast_features(pima_data: PimaData, best_bed: str, feature_name: str):
    # Keep the feature hits for later drawing.  It may be empty, i.e., no feature hits
    try:
        best = pd.read_csv(filepath_or_buffer=best_bed, sep='\t', header=None, skip_blank_lines=True)
        best[0] = best[0].astype(str)
    except (FileNotFoundError,  pd.errors.EmptyDataError):
        best = pd.DataFrame()
    except Exception as e:
        error_out(
            pima_data, f"Unexpected exception when processing BLAST features: {e}"
        )

    if feature_name == "ba_virulence_genes":
        pima_data.ba_virulence_hits = best

    elif feature_name == "amr":
        best = pd.merge(
            best,
            pima_data.amr_gene_drug, 
            how="left", left_on=3, right_on=0, suffixes=("","_y")
        ).set_axis([0,1,2,3,4,5,7,6], axis=1)[[0,1,2,3,4,5,6]] if not best.empty else best
        pima_data.resfinder_hits = best

    elif feature_name == 'inc':
        if not best.empty:
            best[6] = "inc"
        pima_data.inc_hits = best

    else:
        error_out(pima_data, "Should not see this, something wrong with the annotations")


def run_amrfinder(pima_data: PimaData, settings: Settings):
    """Run AMRFinder to search for whole genes potentially conferring resistance
    
    Generates a pandas dataframe
    """

    print_and_log(
        pima_data,
        "Running AMRFinder",
        pima_data.main_process_verbosity,
        pima_data.main_process_color,
    )

    # Keep track of feature hits for reporting
    features_dir = os.path.join(pima_data.output_dir, "features", 'amrfinder')   
    
    # Check if results already exist
    if find_checkpoint(pima_data, features_dir):
        expected_result = os.path.join(features_dir, 'amrfinder.tsv')
        if validate_file_and_size(pima_data, expected_result):
            print_and_log(
                pima_data,
                "Using previously identified amrfinder results",
                pima_data.sub_process_verbosity,
                pima_data.sub_process_color,
            )
            parse_amrfinder(pima_data, expected_result), 
            return
    
    os.makedirs(features_dir)
    make_start_file(pima_data, features_dir)
    amrfinder_results = os.path.join(features_dir, "amrfinder.tsv")
    amrfinder_stdout, amrfinder_stderr = std_files(os.path.join(features_dir, "amrfinder"))
    command = " ".join(
        [
            "amrfinder --nucleotide",
            pima_data.genome_fasta,
            "--database",
            os.path.join(settings.amrfinder_database, "latest"),
            "--output",
            amrfinder_results,
            "--threads",
            str(pima_data.threads), 
            "1>", amrfinder_stdout,
            "2>", amrfinder_stderr,
        ]
    )
    print_and_run(pima_data, command)
    make_finish_file(pima_data, features_dir)
    parse_amrfinder(pima_data, amrfinder_results)


def parse_amrfinder(pima_data: PimaData, path_to_raw_amrfinder_results: str):
    
    with open(path_to_raw_amrfinder_results, 'r') as f:
        lines = f.readlines()
    header_idx = 0
    for i, line in enumerate(lines):
        if "Protein identified" in line or "Gene symbol" in line:
            header_idx = 1
            break
    df = pd.read_csv(
        path_to_raw_amrfinder_results,
        sep="\t",
        skipinitialspace=True,
        skip_blank_lines=True,
        skiprows=header_idx,   
    )
    df.columns = [col.lstrip('#').strip() for col in df.columns]
    subset = df[['Contig id', 'Start', 'Stop', 'Strand', '% Coverage of reference', '% Identity to reference', 'Element symbol', 'Class']]
    pima_data.amrfinder_results = subset


def overlap_features(pima_data: PimaData):
    """
    Takes as input the 3 possibile annotation results and 
    generate a pandas dataframe of non-overlapping features for the visualizations
    """
    
    merged_df = pd.DataFrame(columns = ['contig', 'start', 'stop', 'strand', 'gene', 'class', 'source'])

    #Ugly - concatenates the results into 1 dataframe, and handles the case where there might be no annotations or the analysis wasn't run
    merged_df = pd.concat(
        [
        merged_df,
        (
            pd.DataFrame(pima_data.amrfinder_results[['Contig id', 'Start', 'Stop', 'Strand', 'Element symbol', 'Class']].assign(source='amr').values, columns=merged_df.columns) 
            if isinstance(pima_data.amrfinder_results, pd.DataFrame) and not pima_data.amrfinder_results.empty
            else pd.DataFrame(columns = merged_df.columns)
        ),
        (
            pd.DataFrame(pima_data.resfinder_hits[[0,1,2,5,3,6]].assign(source='amr').values, columns=merged_df.columns) 
            if isinstance(pima_data.resfinder_hits, pd.DataFrame) and not pima_data.resfinder_hits.empty
            else pd.DataFrame(columns = merged_df.columns)
        ),
        (
            pd.DataFrame(pima_data.inc_hits[[0,1,2,5,3,6]].assign(source='inc').values, columns=merged_df.columns) 
            if isinstance(pima_data.inc_hits, pd.DataFrame) and not pima_data.inc_hits.empty
            else pd.DataFrame(columns = merged_df.columns)
        )], ignore_index=True)       
    
    merged_df['class'] = merged_df['class'].str.lower()
    merged_df['gene'] = merged_df['gene'].str.replace(r'_[^_]+$', '', regex=True)
    if merged_df.empty:
        pima_data.unique_hits = merged_df
        return
    
    bed = bd.BedTool.from_dataframe(merged_df)
    merged_bed = bed.sort().merge(c = [4,5,6,7], o='collapse')
    merged_df = merged_bed.to_dataframe()
    merged_df['name'] = merged_df['name'].str.split(",").str[0]
    merged_df['score'] = merged_df['score'].str.split(",").str[0]
    merged_df['strand'] = merged_df['strand'].str.split(",").str[0]
    merged_df['thickStart'] = merged_df['thickStart'].str.split(",").str[0]
    merged_df.rename(columns = {'name': 'strand', 'score': 'gene', 'strand': 'class', 'thickStart': 'source'}, inplace=True)
    pima_data.unique_hits = merged_df
