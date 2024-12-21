import os
import sys
import shutil
import glob

from Pima.pima_data import PimaData
from Pima.utils.settings import Settings

from Pima.utils.utils import (
    print_and_log,
    print_and_run,
    validate_file_and_size,
)

def validate_download(pima_data: PimaData, settings: Settings):

    if not pima_data.download:
        return

    pima_data.errors = []

    pima_data.verbosity = 3
    
    download_databases(pima_data, settings)
    # We don't need this as part of the analysis pipeline, if user tries to donwload, run the download function and then end
    #pima_data.analysis.append(['download_databases', pima_data, settings])

    print_and_log(
        pima_data,
        "Finished validating databases. Re-run analysis without '--download' argument",
        pima_data.main_process_verbosity,
        pima_data.main_process_color,
    )
    sys.exit(0)


def validate_organism(pima_data: PimaData):
    if pima_data.only_assemble:
        return
    
    if not pima_data.organism and not pima_data.list_organisms:
        return

    list_of_org = [
        "Bacillus_anthracis",
    ]

    if pima_data.list_organisms:
        print_and_log(
            pima_data,
            f"List of available reference organisms:\n{' '.join(list_of_org)}",
            pima_data.main_process_verbosity,
            pima_data.main_process_color,
        )
        sys.exit(0)

    print_and_log(
        pima_data, 
        'Validating organism', 
        pima_data.main_process_verbosity, 
        pima_data.main_process_color,
    )
    
    if pima_data.organism and pima_data.reference_fasta:
        pima_data.errors.append("--organism and --reference-genome are mutually exclusive")

    if pima_data.organism and pima_data.mutation_region_bed:
        pima_data.errors.append("--organism and --mutation-regions are mutually exclusive")

    if not pima_data.organism in list_of_org:
        pima_data.errors.append(
            f"--organism {pima_data.organism} is not available, please specify a specifc --reference-genome and --mutations-regions"
            f" or run PiMA without a reference"
        )
        return
    
    if not os.path.isdir(pima_data.reference_dir):
        os.mkdir(pima_data.reference_dir)

    pima_data.organism_dir = os.path.join(pima_data.reference_dir, pima_data.organism)
    if not os.path.isdir(pima_data.organism_dir):
        os.mkdir(pima_data.organism_dir)
    
    pima_data.reference_fasta = os.path.join(pima_data.organism_dir, "genome.fasta")
    pima_data.mutation_region_bed = os.path.join(pima_data.organism_dir, "confirmed_amr_mutations.bed")
    #pima_data.mutation_regions = os.path.join(pima_data.organism_dir, "mutation_regions.bed")
    pima_data.organism_amr_appendices = glob.glob(os.path.join(pima_data.organism_dir, "amr_appendices","*md"))
    
    if not validate_file_and_size(pima_data, pima_data.reference_fasta):
        print_and_log(
            pima_data,
            f"Downloading reference genome for {pima_data.organism}",
            pima_data.sub_process_verbosity, 
            pima_data.sub_process_color,
        )
        download_organism(pima_data, pima_data.organism)
        pima_data.load_reference()

    pima_data.will_have_reference_fasta = True

def download_organism(pima_data: PimaData, organism: str):

    print_and_log(
        pima_data, 
        f"Downloading references specific for {organism}", 
        pima_data.main_process_verbosity, 
        pima_data.main_process_color,
    )

    if organism == "Bacillus_anthracis":
        genome_temp = os.path.join(pima_data.organism_dir, "genome_temp.fasta")
        genome = os.path.join(pima_data.organism_dir, "genome.fasta")
        command = " ".join(
            [
                "wget -O",
                f"{genome_temp}.gz",
                "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/008/445/GCA_000008445.1_ASM844v1/GCA_000008445.1_ASM844v1_genomic.fna.gz",
                "1> /dev/null 2> /dev/null;",
                "gunzip", f"{genome_temp}.gz",
            ]
        )
        print_and_run(pima_data, command)

        replace_dict = {
            ">AE017334.2 Bacillus anthracis str. 'Ames Ancestor', complete genome": ">chromosome",
            ">AE017336.2 Bacillus anthracis str. 'Ames Ancestor' plasmid pXO1, complete sequence": ">pX01",
            ">AE017335.3 Bacillus anthracis str. 'Ames Ancestor' plasmid pXO2, complete sequence": ">pX02",
        }
        with open(genome_temp, "r") as f:
            with open(genome, "w") as w:
                for line in f:
                    for key in replace_dict:
                        if key in line:
                            line = line.replace(key, replace_dict[key])
                    w.write(line)

        os.remove(genome_temp)

def download_databases(pima_data: PimaData, settings: Settings):

    print_and_log(
        pima_data,
        "Checking for missing databases", 
        pima_data.main_process_verbosity, 
        pima_data.main_process_color,
    )

    database_fasta = settings.plasmid_database_default_fasta
    if not validate_file_and_size(pima_data, database_fasta) and validate_file_and_size(pima_data, settings.DockerPathPlasmid):
        pima_data.plasmid_database = settings.DockerPathPlasmid

    elif not validate_file_and_size(pima_data, database_fasta) and not validate_file_and_size(pima_data, settings.DockerPathPlasmid):
        print_and_log(
            pima_data,
            'Downloading plasmid database', 
            pima_data.sub_process_verbosity, 
            pima_data.sub_process_color,
        )
        command = " ".join(
            [
                'wget',
                '-O', database_fasta,
                'http://pima.appliedbinf.com/data/plasmids_and_vectors.fasta',
            ]
        )
        print_and_run(pima_data, command)
    else:
        print_and_log(
            pima_data,
            'Plasmid database present', 
            pima_data.sub_process_verbosity, 
            pima_data.sub_process_color,
        )
    if not os.path.isdir(settings.kraken_database_default) and validate_file_and_size(pima_data, settings.DockerPathKraken):
        pima_data.kraken_database = settings.DockerPathKraken
    elif not validate_file_and_size(pima_data, os.path.join(settings.kraken_database_default, "hash.k2d")) and not validate_file_and_size(pima_data, settings.DockerPathKraken):
        print_and_log(
            pima_data,
            'Downloading and the prebuilt 8gb kraken2 database, 20230605, (may take some time)', 
            pima_data.sub_process_verbosity, 
            pima_data.sub_process_color,
        )
        #if download was corrupted and the files (hash, opts, taxo) are not present, we need to delete and try again
        if os.path.isdir(settings.kraken_database_default):
            shutil.rmtree(settings.kraken_database_default)

        os.makedirs(settings.kraken_database_default)
        command = " ".join(
            [
                'wget -O', f"{settings.kraken_database_default}.tar.gz",
                'https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20230605.tar.gz;',
                'tar xvf', 
                f"{settings.kraken_database_default}.tar.gz", 
                '-C', 
                settings.kraken_database_default,
                '1> /dev/null 2> /dev/null',
            ]
        )
        print_and_run(pima_data, command)
        command = " ".join(
            [
                'rm', 
                f"{settings.kraken_database_default}.tar.gz",
                f"{settings.kraken_database_default}/*kmer_distrib",
                f"{settings.kraken_database_default}/inspect.txt",
                f"{settings.kraken_database_default}/ktaxonomy.tsv",
                f"{settings.kraken_database_default}/seqid2taxid.map"
            ]
        )
        print_and_run(pima_data, command)
    else:
        print_and_log(
            pima_data,
            'Kraken2 database found', 
            pima_data.sub_process_verbosity, 
            pima_data.sub_process_color,
        )