# PIMA: Plasmid, Integrations, Mutations, and Antibiotic resistance annotation pipeline

- [PIMA: Plasmid, Integrations, Mutations, and Antibiotic resistance annotation pipeline](#pima-plasmid-integrations-mutations-and-antibiotic-resistance-annotation-pipeline)
- [Installation](#installation)
  - [Install using prebuilt conda package](#install-using-prebuilt-conda-package)
- [Quickstart Guide](#quickstart-guide)
  - [A typical run](#a-typical-run)
  - [Using fastq files generated by MinKNOW](#using-fastq-files-generated-by-minknow)
    - [General program structure overview](#general-program-structure-overview)
    - [Input file format](#input-file-format)
  - [All Available Arguments](#all-available-arguments)
  - [Additional Notes on Optional Arguments](#additional-notes-on-optional-arguments)
  - [Literature Citations](#literature-citations)
  - [FAQ](#faq)
    - [Error: "(wkhtmltopdf:\<\>): Gtk-WARNING \*\*: \<\>: cannot open display: ".](#error-wkhtmltopdf-gtk-warning---cannot-open-display-)
  - [Software Dependencies](#software-dependencies)
  
  
PIMA (Plasmid, Integrations, Mutations, and Antibiotic resistance annotation pipeline) is a high-throughput sequence analysis pipeline.  PIMA is an end-to-end solution encompassing all the steps required to transform raw ONT output into high-quality annotated assemblies.  PIMA supports providing optional Illumina paired end reads for polishing and error correcting assemblies.
  
The PIMA pipeline performs the following actions [Analysis program options in brackets, default first]:
1. Assess FASTQ files and performs basic quality control and downsampling
2. Check for contamination or multiple species present [Optional, kraken2]
2. Assemble reads [Flye, Raven]
4. Polish and improve assembly [Optional, Medaka, Pilon, Polypolish]
5. Identify and Annotate plasmids [Optional, pChunks]
6. Annotate plasmid incompatiblity (INC) groups
7. Annotate antimicrobial resistance (AMR) genes [resFinder database]
8. Annotation of user-supplied genomic features (from input FASTA file)
9. Comparison to reference genome provided:
     - including variant calling relative to reference [Varscan]
     - annotation of mutations within regions provided by a bed file [custom Bacillus anthracis file provided]
10. Generate visualizations of the location of detected AMR or user-supplied features [dna_features_viewer]
11. Generate visualization of genome (based on provided reference OR the denovo genome generated) [py_circos]
12. Generate PDF report summarizing all results
  
Codebase stage: development  
Developers and maintainers: [Will Overholt](https://github.com/waoverholt), [Alan Collins](https://github.com/Alan-Collins), [Andrew Conley](https://github.com/abconley), [Lavanya Rishishwar](https://github.com/lavanyarishishwar), [Nolan English](https://github.com/Kyxsune), [Vasanta Chivukula](https://github.com/vchivukula7)

### Disclaimer
PiMA has not been cleared or approved by the FDA. The results are intended for public health purposes only and must NOT be communicated to the patient, their care provider, or placed in the patient’s medical record. These results should NOT be used for diagnosis, treatment, or assessment of patient health or management.

# Main differences in version 2

PiMA was rewritten and restructured to make it easier to maintain and add new features. With some notable exceptions, the code base is similar and has just been split into mutiple modules to separate steps into different scripts. Upgrading from pima1 to pima2 will not drastically change the results.

## Significant changes

### Fast5 / pod5 files

PiMA2 does not handle raw ONT data anymore. This change was motivated by ONT consolidating the basecalling, demultiplexing, model selection, and QAQC steps into their sequencing software MinKNOW. As MinKNOW transitioned from using guppy to (dorado)[https://github.com/nanoporetech/dorado] to basecall, we found it more difficult to keep PiMA insync with minKNOW. CDC-users preferred working with MinKNOW and including all the options to basecall with dorado in pima made the software more complicated to run. It also made choosing default settings in pima consequential.
Currently, PiMA 2 only accepts fastq (or fastq.gz) files as input. Since we almost exclusively utilize minKNOW produced fastq files, pima can accept a sequencing folder as input and will generate the sample-specific fastq files for you (see [Using fastq files generated by MinKNOW](#using-fastq-files-generated-by-minknow) below).

### ONT Model

We strongly recommend using SUP models for basecalling in MinKNOW (or separately using dorado). You can ask PiMA to try and determine which model was used with the flag `--ont-model auto`. On the backend, PiMA2 runs `medaka tools resolve_model --auto_model consensus <fastq_file>` on the input sequencing file(s). This command will fail if run information is not encoded in the fastq file heads. We noticed that following the recommended basecalling procedures in (dorado)[https://github.com/nanoporetech/dorado] and emitting a fastq_file (`--emit-fastq`) does not include the run specific information in the sequence headers and will not correctly report to medaka which basecalling model was used. Using `samtools split` instead to create the bams with long sequence IDs and then converting those to fastq (`samtools fastq <input_bam_path> > <output_fastq_path>`) works (as of dorado 0.7.3).

It is always an option to tell PiMA which model was used to basecall:  
`--ont-model r1041_e82_400bps_sup_variant_v4.3.0`  

### ONT watch commands

PiMA2 will not monitor a sequencing folder and initiate an analysis after sufficient reads, as pima1 could. Due to yield differences, we found this application challenging to parameterize. Users preferred either monitoring the data output on MinKNOW during a sequencing run and then manually kickstarting pima if time was pressing, or more likely just running PiMA after sequencing had completed. We are currently evaluating if we should re-introduce this method.

### Estimating the input genome size

PiMA1 required the users to esimate the final size of the expected genome assembly. This information was used for Flye (used to be required) and to guide downsampling. PiMA2 has a method in place to estimate a genome size from the input fastq files. This method is not very accurate, but has been successful in our tests. It works by mapping the ONT reads to a set of [20 universal single copy genes](https://pmc.ncbi.nlm.nih.gov/articles/PMC10631056/) (extracted from *Bacillus anthracis* Ames ancestor). 
The median coverage across all positions from these genes is extracted and the estimated genome size is calculated using the formula (genome_size = number of bases sequenced / median coverage).

If possible we recommend setting this based on your target organism (Bacillus anthracs ~ 5.5m).

### Parallelizing PiMA on HPCs

PiMA2 has Nextflow included to handle parallelizing the analysis of multiplex runs. Since PiMA is a stand-alone python program, we rather crudely wrapped the multiplex handling functions to serve the entire sample-run to nextflow (rather than the more typical parallelizing steps in nextflow-based pipelines). We have exclusively tested this functionality on a QSUB (Oracle Grid / Sun Grid job scheduler). PiMA2 contains a nextflow.config template in the `nextflow_parallelization` folder that will need to be updated for different systems. This works on conda-installed versions of PiMA2 on our system, AFTER modifying the nextflow.config script to point to the location of the conda environment (e.g. `conda = /full/path/to/conda/envs/pima`). You may need to specify a `beforeScript` block to use conda on your HPC system (e.g. we need to have a line in the process block that points to the system conda activate script similar to: `beforeScript = 'source <path/to/conda/bin/activate'`)

### Parallelizing PiMA on local machines

Similar to PiMA1, PiMA2 will analyze a multiple sequencing run from a single command (provide the `--multiplexed` flag). In this mode, each sample will be analyzed sequentially and the `--threads` flag will tell each process how many threads to use (each sample will be analyzed with the same number of threads, one after another until all are complete).
This section of pima was entirely rewritten so please let us known if (when) you identify bugs. 

### Calling AMR-conferring variants

PiMA has always been targetted at the analysis of *Bacillus anthracis* and PiMA1 provided a `mutation-regions.bed` file containing the AMR-conferring locations in the standard reference (*Bacillus anthracis* Ames ancestor)["https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/008/445/GCA_000008445.1_ASM844v1/GCA_000008445.1_ASM844v1_genomic.fna.gz"]. We have re-worked these functions to differentiate between specific SNPs/INDELs that have been experimentally shown to confer mutations and mutations that have not been observed before but are located in genes where other AMR-conferring mutations have been identified. 
The new report format separates these experimentally confirmed mutations (will be bolded) from potentially-conferring mutations.
There has also been a lot of processing done on reporting identified variants if your input organism is too distant from your provided reference (or the PiMA2 Ba_Ames_ancestor reference).

#### Calling AMR-conferring INDELS

Some of the experimentally confirmed AMR-conferring mutations are due to 1bp deletions within homopolymer regions within some genes (e.g. rsiP). ONT has traditionally struggled to detect these INDELs and PiMA1 would only report these mutations if Illumina data was also provided. PiMA2 will report these mutations IF a SUP (super accuracy) model was used to generate the fastq files. If pima cannot determine which model was used OR a fast or hac model was used, PiMA2 will suppress these from the report as they have a high chance of being a false positive. The intermediate files within the mutations directory in the output folder can be scrutinized, and the `reference_mapping_ont.bam` can be imported into IGV along with the Ames ancestor reference genome to inspect the raw sequence mapping results.

### Appendices

If AMR-conferring mutations are detected, PiMA will include an Appendix specific to that mutation and the class of antibiotics it confers resistance towards.

### Changes to the report

The report for PiMA2 has been reworked and new tools for visualizing the results been incorporated. These changes are aimed at improving readability, reproducbility, and overall clarity. Genome visualizations will be generated when the reference *Bacillus anthracs* Ames ancestor is used, as in PiMA1. However, PiMA2 now includes a visualization of the depth of coverage (ONT and/or Illumina) and will create these plots for the de novo assembly, if a reference is not provided.

### Resuming a run
PiMA2 will register a "checkpoint" after each successful completion of a step in the pipeline. Users can re-run a sample that failed part way by including the `--resume` flag. However, if you change how you want the run to proceed, you can trick PiMA as previously completed steps are based off the original command and PiMA does not check to make sure that a re-run command would not change these older results.  

We find this function most useful if the de novo genome assembly finishes (which takes the most time) but you realize you wanted to include a plasmids analysis or a kraken2 analysis.  

We also find it useful at times to delete completed folders within the output to force PiMA to regenerate those results (e.g. changing the reference or mutations file and regenerating the mutations / insertions / circos / quast results).

### Singularity container
While a singularity container for PiMA2 has been generated, we are still working out a few kinks with how nextflow-based multiplexing is performed. Standard multiplexing and singleplex analysis works (most of the time), but we have not tested these containers as well as we have tested the conda environment based installs. The steps to build a docker container and then create a singularity container are outlined in the `dockerbuild/Dockerfile` file. 

# Installation
## Install using prebuilt conda package
A conda package of PiMA is available on our conda channel (appliedbinf)[https://anaconda.org/appliedbinf/dashboard]. Using this package, you can install PiMA into a new conda environment using mamba (preferred) or conda by:
```{bash}
mamba create -n pima -c appliedbinf -c conda-forge -c bioconda pima
mamba activate pima
pima -h

# If you want to download the Bacillus anthracis reference files, kraken database, and other required folders, run:
pima --download
```

## Install using a singularity container
A singularity/docker container is not currently pubicly available, although all the steps to generate one are present in the codebase (`dockerbuild/Dockerfile`). Please note that the nextflow-based multiplexing DOES NOT WORK in a container.

# Quickstart Guide

## A typical run

Here is an example of a basic run using fastq files as input.

```
# basic command
pima --output output_folder --ont-fastq fastq_file --genome-size 5M --threads 4 --organism Bacillus_anthracis

```

We prefer to run PiMA using the verbose flag (--verbosity) and my hardware can support upto 20 threads. 
Feel free to turn off the verbose flag if you so desire, and change the number of threads in accordance with your hardware 
capability. This command will create the following outputs listed in no particular order:
  
1. ont_assembly: Folder with the assembly FASTA file (assembly.fasta) and assembly metrics.
2. ont_fastq: Folder with the raw FASTQ file that is used to generate the assembly and the guppy run.
3. quast: Folder with the quality metrics and the quast report.
4. downsample: folder with a FASTQ with a subset of a high-coverage data set to a given average coverage.
5. medaka: folder with the consensus sequence (consensus.fasta).
6. drawing: provides the contig plot.
7. insertions: provides information on the alignment with the reference, alignment coordinates, any snps and report generated with this information.
8. circos: provides svg and png files with visuals of the chromosome and plasmids (if present).
9. mutations: includes files with information on mapping with the reference and the snps identified.
10. report: A final report in a pdf format describing the assemblies and other information pertaining to the assembled genome (such as mutations, amr genes identified).
11. info: provides the coverage information based on mapping to a reference (ont_coverage.tsv).
12. features: provides information on the presence of amr genes and any inserions, if present in the assembled genome. If there are no insertions or amr genes then you will get an empty file.

## Using fastq files generated by MinKNOW
Pima can also take directories as input and search for the fastq files. We usually work with ONT data generated by minknow, and therefore give pima the fastq_pass folder. PiMA checks each of the folders within fastq_pass (barcodes if multiplexed). If a folder contains fastq files that contain at least 2.5% (tunable using --barcode-min-fraction) of all data in the provided folder (e.g. that barcode contains at least 2.5% of the sequence data) it will be included in the run.

```
pima --output output_folder --ont-fastq  ../20241030_1637_MD-100828_FBA24753_235158ef/fastq_pass --threads 4 --genome-size estimate --organism Bacillus_anthracis --ont-model r1041_e82_400bps_sup_variant_v4.3.0 --verbosity 3 --multiplexed
```

When providing PiMA an input directory, you need to allow provide the "--multiplexed" flag. If the data is all contained within a single barcode pima will identify which folder has the sequencing data (ignoring the unclassified folder).

You can also just provide the entire run folder, in which case pima will merge the reads in the fastq_fail folders with the same names from the fastq_pass folders. Again, the final detected set will need to be at least 2.5% of the total amount of fastq data generated. 

If there are multiple folders containing reads (e.g. a multiplexed run), PiMA will analyze each sample separately, in serial. This works great for sequencing computers or individual computers. If you are on an high performance computing system, you can use Nextflow to handle parallelizing the sample analyses. We have used this with the CDC HPC system and have example nextflow.config files that will need to be modified (Pima/nextflow_parallelization/nextflow.config).

```
pima --output output_folder --ont-fastq  ../20241030_1637_MD-100828_FBA24753_235158ef/fastq_pass --threads 4 --genome-size estimate --organism Bacillus_anthracis --ont-model r1041_e82_400bps_sup_variant_v4.3.0 --verbosity 3 --multiplexed --nextflow
```

## Screening for potential contaminants
PiMA can be run with a kraken2 step to check what organisms are mostly likely present in the sequence dataset by enabling the "--contamination" flag.

```
pima --output output_folder --ont-fastq fastq_file --genome-size 5M --threads 4 --organism Bacillus_anthracis --contamination
```

## Screening for unknown plasmids
PiMA can also search for non-pX01/px02 plasmids that might be present in a sample. With this enabled PiMA will screen each contig to see if it is similar to a plasmid contained the [plasmids_and_vectors.fasta](http://pima.appliedbinf.com/data/plasmids_and_vectors.fasta) file. If multiple contigs belong to the same plasmid, pima will report then together (using pChunks).

```
pima --output output_folder --ont-fastq fastq_file --genome-size 5M --threads 4 --organism Bacillus_anthracis --plasmids
```

## Overwriting a previous run
```
pima --output output_folder --ont-fastq fastq_file --genome-size 5M --threads 4 --organism Bacillus_anthracis --overwrite
```
This will tell pima to delete the <output_folder> and restart each step.

## Resuming a previous run that did not complete
```
pima --output output_folder --ont-fastq fastq_file --genome-size 5M --threads 4 --organism Bacillus_anthracis --resume
```
This will tell pima to continue from the last completed step within the <output_folder>.

Users can control this by manually deleting each folder they want to regenerate.
```
rm -r ./ouput_folder/insertions
rm -r ./output_folder/mutations
rm -r ./output_folder/circos

pima --output output_folder --ont-fastq fastq_file --genome-size 5M --threads 4 --organism Bacillus_anthracis --resume
```
The above steps will re-run all the steps necessary to generate the insertions, mutations, and circos results.

# General program structure overview
PiMA is a multi-step process that can be run as a single command - assessing input fastq files, checking for contamination, assembly, generating a quast report, identifying any mutations based on the reference file provided, and summarizing the results in a report PDF.

## Input file format

The input file can be a single fastq file (or fastq.gz file) representing all the generated sequences from a sample. The input can also be a directory or directories of fastq.gz files (see below). This example directory structure is generated by the rapid barcoding kit, and only some of the barcode folders have sequence data while most have erroneously assigned sequences (due to errors in the barcodes). Running PiMA in a multiplex mode will search through the directory structure, pulling out folders with fastq(.gz) files.  

If you provide the fastq_pass folder, pima will analyze each barcode with at least 2.5% (--barcode-min-fraction 0.025) of the total fastq data generated.  

If you provide the Example_Run folder, pima will merge the fastq_fail and fastq_pass folders with the same name (e.g. barcode01 will comprise reads from both fastq_pass and fastq_fail folders). The combined set of reads will need to be above the 2.5% threshold to be included.

```
20240918_MinION_Example_Run
└── no_sample_id
    └── 20240918_XXXX_XX-XXXXXX_XXXXXXX_XXXXXXX
        ├── fastq_fail
        │   ├── barcode01
        │   ├── barcode02
        │   ├── barcode04
        │   ├── barcode05
        │   ├── barcode06
        │   ├── barcode07
        │   ├── barcode08
        │   ├── barcode09
        │   ├── barcode10
        │   ├── barcode12
        │   ├── barcode13
        │   ├── barcode14
        │   ├── barcode16
        │   ├── barcode17
        │   ├── barcode18
        │   ├── barcode19
        │   ├── barcode20
        │   ├── barcode22
        │   ├── barcode24
        │   └── unclassified
        ├── fastq_pass
        │   ├── barcode01
        │   ├── barcode02
        │   ├── barcode03
        │   ├── barcode04
        │   ├── barcode05
        │   ├── barcode06
        │   ├── barcode07
        │   ├── barcode08
        │   ├── barcode09
        │   ├── barcode10
        │   ├── barcode11
        │   ├── barcode12
        │   ├── barcode13
        │   ├── barcode16
        │   ├── barcode17
        │   ├── barcode19
        │   ├── barcode20
        │   ├── barcode21
        │   ├── barcode22
        │   ├── barcode23
        │   ├── barcode24
        │   └── unclassified
        ├── other_reports
        └── pod5
```
  
## All Available Arguments
The full description of each commandline option is provided below.

```

usage: pima.py [-h] [-v] [--ont-watch <ONT_DIR>] [--ont-watch-min-reads <INT>] [--ont-watch-max-time <HOURS>]
               [--ont-watch-between-time <MINUTES>] [--ont-watch-min-coverage <X>] [--ont-model <ONT BASECALLING MODEL]
               [--ont-fastq <FASTQ|GZ>, fastq_pass] [--multiplexed] [--nextflow [nextflow configs]]
               [--barcode-min-fraction] [--contamination] [--illumina-fastq <R1> [R2] [<R1> [R2] ...]]
               [--genome <GENOME_FASTA>] [--output <OUTPUT_DIR>] [--overwrite] [--keep-intermediates]
               [--assembler {flye,flye_sup,raven}] [--genome-size <GENOME_SIZE>] [--assembly-coverage <X>] [--no-medaka]
               [--illumina-polisher {pilon,polypolish,skip}] [--only-assemble] [--no-assembly] [--download] [--plasmids]
               [--plasmid-database <PLASMID_FASTA>] [--amr-database <AMR_FASTA>] [--no-amr] [--inc-database <INC_FASTA>]
               [--no-inc] [--feature <FEATURE_FASTA>] [--no-drawing] [--self-circos] [--reference-dir <REFERNCE_DIR>]
               [--organism <Genus_species>] [--list-organisms] [--reference-genome [<GENOME_FASTA>]]
               [--mutation-regions [<REGION_BED>]] [--name <NAME>] [--threads <NUM_THREADS>] [--verbosity <INT>]
               [--logfile <FILENAME>] [--resume] [--bundle <PATH>] [--fake-run]

P.I.M.A. bacterial genome analysis pipeline

Help and version:
  -h, --help                                   Print this help and exit.
  -v, --version                                Print the software version.

Input and basecalilng options:
  --ont-watch <ONT_DIR>                        Directory from an ONT run.
  --ont-watch-min-reads <INT>                  Minimum number of ONT reads to watch for.
  --ont-watch-max-time <HOURS>                 Maxmium time to watch for new reads.
  --ont-watch-between-time <MINUTES>           Maximum to to wait between new FAST5 files
  --ont-watch-min-coverage <X>                 Minimum genome coverage achieved before starting analysis.
  --ont-model <ONT BASECALLING MODEL           ONT model used for base calling. 'auto' will try and determine which
                                               basecalling model was used and if 'sup' is detected will tell flye to use
                                               '--nano-hq'. Available models can be listed using 'medaka tools
                                               list_models' as described by medaka
                                               (https://github.com/nanoporetech/medaka) (default : auto)
  --ont-fastq <FASTQ|GZ>, fastq_pass           File containing basecalled ONT reads (fastq|gz), or folder containing
                                               directories of demultiplexed basecalled reads (fastq_pass)
  --multiplexed                                The ONT data are multiplexed (default : False)
  --nextflow [nextflow configs]                Use nextflow to parallelize PiMA (otherwise samples will be run in
                                               serial). Can also provide additional commands to the nextflow run wrapped
                                               in single quotes ('), e.g. --nextflow '-c <my_hpc_nextflow.config>
                                               -resume' (default : False)
  --barcode-min-fraction                       The minimum fraction of data necessary to include a barcode in the
                                               analysis (default : 0.025)
  --contamination                              Use Kraken2 to look for contamination in the input read (default : False)
  --illumina-fastq <R1> [R2] [<R1> [R2] ...]   Files containing R1 & R2 Illumina reads
  --genome <GENOME_FASTA>                      A genome FASTA file to be used in place of assembly

Output options:
  --output <OUTPUT_DIR>                        Output directory for the analysis
  --overwrite                                  Overwrite an existing output directory (default : False)
  --keep-intermediates                         Keep all intermediate files (mostly retains large bam/mpileup files)
                                               (default : False)

Assembly options:
  --assembler {flye,flye_sup,raven}            Assembler to use. Choose 'flye_sup' if you are providing R10 fastq data
                                               base called with a super_accuracy model and wish to use flye (default :
                                               flye)
  --genome-size <GENOME_SIZE>                  Genome size estimate for the assembly & downsampling, can use M,m,K,k
                                               shorthand (e.g. 5.5m). If unknown, you can let pima estimate the genome
                                               size using 'estimate' (default : estimate)
  --assembly-coverage <X>                      Downsample the provided reads to this coverage (default : 200X)
  --no-medaka                                  Skip Medaka polising of the ONT assembly (faster) (default : False)
  --illumina-polisher {pilon,polypolish,skip}  Polish the genome assembly using short-reads, to skip use '--illumina-
                                               polisher skip' (default : pilon)
  --only-assemble                              Only carry out assembly steps with the given data; no downstream
                                               analysis.
  --no-assembly                                Don't attempt to assembly/polish a given genome/set of reads (default :
                                               False)

Database downloading arguments:
  --download                                   Attempt to download Kraken/Plasmid databases if not found locally.Use
                                               witout other options.

Plasmid and vector search options:
  --plasmids                                   Do a plasmid search (default : False)
  --plasmid-database <PLASMID_FASTA>           Path to a FASTA file with reference plasmid sequences

AMR gene search options:
  --amr-database <AMR_FASTA>                   Path to a FASTA file with AMR gene sequences (default : /scicomp/home-
                                               pure/tsz0/.conda/envs/pima2/lib/python3.10/site-
                                               packages/Pima/data/amr.fasta)
  --no-amr                                     Skip AMR search (default : False)

Incompatibility group search options:
  --inc-database <INC_FASTA>                   Path to a FASTA file with incompatibility group sequences (default :
                                               /scicomp/home-pure/tsz0/.conda/envs/pima2/lib/python3.10/site-
                                               packages/Pima/data/inc.fasta)
  --no-inc                                     Skip incompatibility group search (default : False)

Other feature search options:
  --feature <FEATURE_FASTA>                    Path to a FASTA file with feature sequences

Drawing options:
  --no-drawing                                 Skip drawing of contigs & Features (default : False)
  --self-circos                                Use the assembled genome as the reference to draw circos images. NOT
                                               recommended for Illumina assemblies (default : False)

Reference options:
  --reference-dir <REFERNCE_DIR>               Directory containing refrence organisms (default : /scicomp/home-
                                               pure/tsz0/.conda/envs/pima2/lib/python3.10/site-
                                               packages/Pima/data/reference_sequences)
  --organism <Genus_species>                   Reference organism to compare against
  --list-organisms                             List the reference organisms available to this pipeline
  --reference-genome [<GENOME_FASTA>]          Reference genome to compare against (default : False)
  --mutation-regions [<REGION_BED>]            Regions in the reference genome to screen for mutations (default : False)

Other options:
  --name <NAME>                                Name of this analysis for reporting.
  --threads <NUM_THREADS>                      Number of worker threads to use (default : 1)
  --verbosity <INT>                            How much information to print as PiMA runs (default : 1)
  --logfile <FILENAME>                         Name of logfile written in output directory (default : pipeline.log)
  --resume                                     Restart pipeline from last completed step (default : False)
  --bundle <PATH>                              Local Tectonic bundle (default : None)
  --fake-run                                   Don't actually run the pipeline, just pretend to (default : False)

```

## Literature Citations
If you are using PIMA, please cite the following literature:
McLaughlin HP, Bugrysheva JV, Conley AB, Gulvik CA, Kolton CB, Marston C, Swaney E, Lonsway DR, Cherney B, Gargis AS, Kongphet-Tran T, Lascols C, Michel P, Villanueva J, Hoffmaster ER, Gee JE, Sue D. 2020. When minutes matter: rapid nanopore whole genome sequencing for anthrax emergency preparedness. Emerging Infectious Diseases  
