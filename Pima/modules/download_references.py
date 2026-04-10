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
    calc_md5,
)

def validate_download(pima_data: PimaData, settings: Settings):

    if not pima_data.download:
        return

    pima_data.errors = []

    pima_data.verbosity = 3
    
    download_databases(pima_data, settings)

    print_and_log(
        pima_data,
        "Finished validating databases. Re-run analysis without '--download' argument",
        pima_data.main_process_verbosity,
        pima_data.main_process_color,
    )
    sys.exit(0)


def validate_organism(pima_data: PimaData, settings: Settings):
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
    
    if not os.path.isdir(settings.reference_dir):
        os.mkdir(settings.reference_dir)

    pima_data.organism_dir = os.path.join(settings.reference_dir, pima_data.organism)
    if not os.path.isdir(pima_data.organism_dir):
        os.mkdir(pima_data.organism_dir)
    
    pima_data.reference_fasta = os.path.join(pima_data.organism_dir, "genome.fasta")
    pima_data.mutation_region_bed = os.path.join(pima_data.organism_dir, "confirmed_amr_mutations.bed")
    pima_data.organism_amr_appendices = glob.glob(os.path.join(pima_data.organism_dir, "amr_appendices","*md"))

    #READ version info & grab md5sums for the key files
    with open(os.path.join(pima_data.organism_dir, "version.txt"), 'r') as f:
        version = f.readline().strip().split(",")[1]
        bed_md5 = f.readline().strip().split(",")[2]
        genome_md5 = f.readline().strip().split(",")[2]
    pima_data.organism_version = {"version": version, "genome": genome_md5, "bed": bed_md5}

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
    #Verify checksums
    genome_version = calc_md5(pima_data, os.path.join(pima_data.organism_dir, "genome.fasta"))
    bed_version = calc_md5(pima_data, os.path.join(pima_data.organism_dir, "confirmed_amr_mutations.bed"))

    if not genome_version == pima_data.organism_version['genome']:
        pima_data.errors.append(
            f"The {pima_data.reference_fasta} checksum does not match the expected version. "
            f"If this was expected, please update the {os.path.join(pima_data.organism_dir, 'version.txt')} file."
        )
        return
    if not bed_version == pima_data.organism_version['bed']:
        pima_data.errors.append(
            f"The {pima_data.mutation_region_bed} checksum does not match the expected version. "
            f"If this was expected, please update the {os.path.join(pima_data.organism_dir, 'version.txt')} file."
        )
        return

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
        dl_genome_md5 = calc_md5(pima_data, genome)
        if not pima_data.organism_version['genome'] == dl_genome_md5:
            pima_data.errors.append(
                f"The newly downloaded Banthracis genome has a md5sum that does not match expectations. "
                f"This should not have happened. If you are confident the genome is identical you can update "
                f"the {os.path.join(pima_data.organism_dir, 'version.txt')} file with the new genome's md5sum."
                f"If the positions are not identical the mutations that PiMA searches for will not be valid."
            )
            return

def download_databases(pima_data: PimaData, settings: Settings):

    print_and_log(
        pima_data,
        "Checking for missing databases", 
        pima_data.main_process_verbosity, 
        pima_data.main_process_color,
    )

    database_fasta = settings.plasmid_database_default_fasta

    if not validate_file_and_size(pima_data, database_fasta):
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

    if not validate_file_and_size(pima_data, os.path.join(settings.kraken_database_default, "hash.k2d")):
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