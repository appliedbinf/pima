# PIMA: Plasmid, Integrations, Mutations, and Antibiotic resistance annotation pipeline

  * [PIMA: Plasmid, Integrations, Mutations, and Antibiotic resistance annotation pipeline](#pima-plasmid-integrations-mutations-and-antibiotic-resistance-annotation-pipeline)
   * [Installation](#installation)
      * [Conda environments](#conda-environments)
      * [Dependencies](#dependencies)
         * [Software](#software)
      * [Cloning PiMA pipeline](#Cloning PiMA pipeline)
         
         
      * [Manual install (not recommended)](#manual-install-not-recommended)
         * [Software](#software-1)
   * [Quickstart / demo](#quickstart--demo)
   * [Usage](#usage)
   * [Using your own data](#using-your-own-data)


PIMA (Plasmid, Integrations, Mutations, and Antibiotic resistance annotation pipeline) is a high-throughput sequence analysis pipeline.  PIMA is an end-to-end solution encompassing all the steps required to transform raw ONT output into high-quality annotated assemblies.  PIMA supports providing optional Illumina paired end reads for polishing and error correcting assemblies.

The PIMA pipeline performs the following actions [Analysis program options in brackets, default first]:
1. Performs basecalling from FAST5 files [guppy or albacore]
2. Demultiplexes basecalled ONT data (FASTQ files) [qcat or porechop]
3. Performs basic quality control [Porechop]
4. Optional base quality correction [LorDEC]
5. Assemble reads [Flye, miniasm, or wtdbg2]
6. Polish and improve assembly [Optional, racon, Medaka, nanopolish, Pilon]
7. Annotate plasmids
8. Annotate plasmid incompatiblity (INC) groups
9. Annotate antimicrobial resistance (AMR) genes
10. Annotation of other user-supplied genomic features (from FASTA)
11. Comparison to reference genome, including variant calling and annotation

Codebase stage: development   
Developers and maintainers: Andrew Conley, Lavanya Rishishwar   
Testers: Lavanya Rishishwar, Andrew Conley, Vasanta Chivukula, Nolan English

# Installation
We are working on improving the user experience with PiMA. This section will be updated with easier installation procedure for PiMA and it's dependencies.

## Conda environments

We **highly** recommend using the supplied Conda environments for installing and managing all the dependencies for PIMA.  

```bash
# Install the PIMA base environment (Python 3.6)
conda env create -f pima.yml

```
## Dependencies

### Software
* BCFtools & Samtools: http://www.htslib.org/
* bedtools: https://bedtools.readthedocs.io/en/latest/
* NCBI BLAST+ (v 2.8+): https://blast.ncbi.nlm.nih.gov/Blast.cgi
* bwa: https://github.com/lh3/bwa
* flye: https://github.com/fenderglass/Flye
* Medaka: https://github.com/nanoporetech/medaka
* miniasm: https://github.com/lh3/miniasm
* minimap2: https://github.com/lh3/minimap2
* nanopolish: https://github.com/jts/nanopolish
* parallel: https://www.gnu.org/software/parallel/
* Pilon: https://github.com/broadinstitute/pilon/wiki
* Porechop: https://github.com/rrwick/Porechop
* qcat: https://github.com/nanoporetech/qcat
* SAMtools: http://htslib.org/
* SPAdes: http://cab.spbu.ru/software/spades/
* wtdbg2: https://github.com/ruanjue/wtdbg2

```bash
# You can also install the dependencies while creating the PiMA base environment.
conda create -n pima medaka=1.4.3 varscan r flye blast circos minimap2 bwa samtools bedtools pandas pathos joblib pylatex tectonic mummer qcat -y 

# Some more dependencies (si_prefix, dna_features_viewer, and quast) installed using pip installer
pip install si_prefix dna_features_viewer
pip install git+https://github.com/ablab/quast
```

## Cloning PiMA pipeline

```bash
# Activate the conda environement using
conda activa pima

# clone the PiMA pipeline using
git clone https://github.com/abconley/pima.git

# If you want to download the reference file, kraken database, and other required folders, run
pima.py --download

```
## Quickstart guide
## Basic usage
### General program structure overview
PiMA is a multi-step process that can be run as a single command - base calling, checking for contamination, assembly, generating a quast report, identifying any mutations based on the reference file provided.

```bash
usage: pima.py [--help] 

pima.py - Plasmid, Integrations, Mutations, and Antibiotic resistance annotation pipeline

optional arguments:
  --help, -h, --h
  
pima.py command: pima.py --output <your output folder name> --ont-fast5 <where the input fast5 files are located> --genome-size <5M for Ba; genome size in Mb> --threads <number of parallel threads; 12 is the maximum>

```
### Input file format

The input file can be in fast5 format (use the --ont-fast5 flag for fast5 files) or fastq format (the flag would be --ont-fastq). 

### A typical run

Here is an example of a basic run using fast5 files as input.

```bash
# basic command
./pima.py --output output_folder --ont-fast5 fast5_files --genome-size 5M --threads 12

## to run kraken, plasmids, quast, and mutations, the following command should be used
pima.py --out out --ont-fast5 barcodes_folder --threads 16 --overwrite --genome-size 5m --verb 3 --reference-genome ref_genome.fasta --mutation-regions mutation_regions.bed --plasmids

```
Consider an example scenario where you want to assemble Bacillus anthracis ont reads. If the reference files is named ref.fasta 
and the query fast5 files are in the folder named barcodes_folder, the mutation regions bed file is names mutation_regions.bed 
and the output folder you named is ont_output then your pima command will look as follows:

```bash
pima.py --out ont_output --ont-fast5 barcodes_folder --threads 16 --overwrite --genome-size 5m --verb 3 --reference-genome ref.fasta --mutation-regions mutation_regions.bed

```
I prefer to run PiMA using the verbose flag and my hardware can support upto 20 threads. 
Feel free to turn off the verbose flag if you so desire, and change the number of threads in accordance with your hardware 
capability. This command will create the following outputs listed in no particular order:

1. ont_assembly: Folder with the assembly metrics and the fasta file
2. ont_fastq: Folder with the raw_fastq file that is used to generate the assembly and the guppy run.
3. quast: Folder with the quality metrics and the quast report.
4. downsample: folder with a fastq with a subset of a high-coverage data set to a given average coverage.
5. medaka: folder with the consensus sequence.
6. drawing: provides the contig plot.
7. insertions: provides information on the alignment with the reference, alignment coordinates, any snps and report generated with this information.
8. circos: provides svg and png files with visuals of the chromosome and plasmids (if present).
9. mutations: includes files with information on mapping with the reference and the snps identified.
10. report: A final report in a pdf format describing the assemblies and other information pertaining to the assembled genome (such as mutations, amr genes identified).
11. info: provides the coverage information based on mapping to a reference.
12. features: provides information on the presence of amr genes and any inserions, if present in the assembled genome.


## Documentation for the full run
The full description of each commandline option is provided below.

```bash

usage: pima.py [--help] [--version] [--ont-watch <ONT_DIR>] [--ont-watch-min-reads <INT>] [--ont-watch-max-time <HOURS>]
               [--ont-watch-between-time <MINUTES>] [--ont-watch-min-coverage <X>] [--ont-fast5 <ONT_DIR>]
               [--basecaller {guppy}] [--ont-fastq <FASTQ|GZ>] [--multiplexed] [--error-correct] [--contamination]
               [--only-basecall] [--illumina-fastq <R1> [R2] [<R1> [R2] ...]] [--genome <GENOME_FASTA>]
               [--output <OUTPUT_DIR>] [--overwrite] [--assembler {wtdbg2,flye}] [--genome-size <GENOME_SIZE>]
               [--assembly-coverage <X>] [--racon] [--racon-rounds <NUM_ROUNDS>] [--no-medaka] [--only-assemble]
               [--no-assembly] [--download] [--plasmids] [--plasmid-database <PLASMID_FASTA>]
               [--amr-database <AMR_FASTA>] [--no-amr] [--inc-database <INC_FASTA>] [--no-inc]
               [--feature <FEATURE_FASTA>] [--no-drawing] [--reference-dir <REFERNCE_DIR>] [--organism <ORGANISM>]
               [--list-organisms] [--auto-reference] [--reference-genome <GENOME_FASTA>]
               [--mutation-regions <REGION_BED>] [--name <NAME>] [--threads <NUM_THREADS>] [--verbosity <INT>]
               [--bundle <PATH>] [--fake-run]

P.I.M.A. bacterial genome analysis pipeline

Help and version:
  --help                                      Print this help and exit.
  --version                                   Print the software version.

Input and basecalilng options:
  --ont-watch <ONT_DIR>                       Directory from an ONT run.
  --ont-watch-min-reads <INT>                 Minimum number of ONT reads to watch for.
  --ont-watch-max-time <HOURS>                Maxmium time to watch for new reads.
  --ont-watch-between-time <MINUTES>          Maximum to to wait between new FAST5 files
  --ont-watch-min-coverage <X>                Minimum genome coverage achieved before starting analysis.
  --ont-fast5 <ONT_DIR>                       Directory containing ONT FAST5 files
  --basecaller {guppy}                        The basecaller for ONT FAST5 data (default : guppy)
  --ont-fastq <FASTQ|GZ>                      File containing basecalled ONT reads
  --multiplexed                               The ONT data are multiplexed (default : False)
  --error-correct                             Use LORMA to error-correct ONT reads (default : False)
  --contamination                             Use Kraken2 to look for contamination in the input read (default : False)
  --only-basecall                             Just basecall and demultiplex/trim/correct. No downstream analysis.
  --illumina-fastq <R1> [R2] [<R1> [R2] ...]  Files containing R1 & R2 Illumina reads
  --genome <GENOME_FASTA>                     A genome FASTA file to be used in place of assembly

Output options:
  --output <OUTPUT_DIR>                       Output directory for the analysis
  --overwrite                                 Overwrite an existing output directory (default : False)

Assembly options:
  --assembler {wtdbg2,flye}                   Assembler to use (default : flye)
  --genome-size <GENOME_SIZE>                 Genome size estimate for the assembly & downsampling (default : None))
  --assembly-coverage <X>                     Downsample the provided reads to this coverage (default : 200X)
  --racon                                     Force the generation a racon consenus (default : False)
  --racon-rounds <NUM_ROUNDS>                 Number of RACON rounds used to generate a consensus (default : 4)
  --no-medaka                                 Skip Medaka polising of the ONT assembly (faster) (default : False)
  --only-assemble                             Only carry out assembly steps with the given data; no downstream analysis.
  --no-assembly                               Don't attempt to assembly/polish a given genome/set of reads (default :
                                              False)

Database downloading arguments:
  --download                                  Attempt to download Kraken/Plasmid databases if not found locally.Use
                                              witout other options.

Plasmid and vector search options:
  --plasmids                                  Do a plasmid search (default : False)
  --plasmid-database <PLASMID_FASTA>          Path to a FASTA file with reference plasmid sequences

AMR gene search options:
  --amr-database <AMR_FASTA>                  Path to a FASTA file with AMR gene sequences (default :
                                              /data/home/vchivukula7/pima/data/amr.fasta)
  --no-amr                                    Skip AMR search (default : False)

Incompatibility group search options:
  --inc-database <INC_FASTA>                  Path to a FASTA file with incompatibility group sequences (default :
                                              /data/home/vchivukula7/pima/data/inc.fasta)
  --no-inc                                    Skip incompatibility group search (default : False)

Other feature search options:
  --feature <FEATURE_FASTA>                   Path to a FASTA file with feature sequences

Drawing options:
  --no-drawing                                Skip drawing of contigs & Features (default : False)

Reference options:
  --reference-dir <REFERNCE_DIR>              Directory containing refrence organisms (default :
                                              /data/home/vchivukula7/pima/data/reference_sequences)
  --organism <ORGANISM>                       Reference organism to compare against
  --list-organisms                            List the reference organisms available to this pipeline
  --auto-reference                            Automatically choose an appropriate reference if available.
  --reference-genome <GENOME_FASTA>           Reference genome to compare against (default : None)
  --mutation-regions <REGION_BED>             Regions in the reference genome to screen for mutations (default : None)

Other options:
  --name <NAME>                               Name of this analysis for reporting.
  --threads <NUM_THREADS>                     Number of worker threads to use (default : 1)
  --verbosity <INT>                           How much information to print as PiMA runs (default : 1)
  --bundle <PATH>                             Local Tectonic bundle (default : None)
  --fake-run                                  Don't actually run the pipeline, just pretend to (default : False)
```




### Python modules and programs


## Literature Citations
- McLaughlin HP, Bugrysheva JV, Conley AB, Gulvik CA, Kolton CB, Marston C, Swaney E, Lonsway DR, Cherney B, Gargis AS, Kongphet-Tran T, Lascols C, Michel P, Villanueva J, Hoffmaster ER, Gee JE, Sue D. 2020. When minutes matter: rapid nanopore whole genome sequencing for anthrax emergency preparedness. Emerging Infectious Diseases
- [BCFtools and SAMtools](https://www.ncbi.nlm.nih.gov/pubmed/19505943)
- [BEDTools](https://www.ncbi.nlm.nih.gov/pubmed/25199790)
- [BLAST+](https://www.ncbi.nlm.nih.gov/pubmed/20003500)
- [GNU parallel](https://www.usenix.org/publications/login/february-2011-volume-36-number-1/gnu-parallel-command-line-power-tool)
- [Flye](https://www.ncbi.nlm.nih.gov/pubmed/30936562)
- [miniasm and minimap2](https://www.ncbi.nlm.nih.gov/pubmed/27153593)
- [MUMmer](https://www.ncbi.nlm.nih.gov/pubmed/14759262)
- [racon](https://www.ncbi.nlm.nih.gov/pubmed/28100585)
- [SPAdes](https://www.ncbi.nlm.nih.gov/pubmed/22506599)
- [Wtdbg2](https://www.nature.com/articles/s41592-019-0669-3)
