# PIMA: Plasmid, Integrations, Mutations, and Antibiotic resistance annotation pipeline

  * [PIMA: Plasmid, Integrations, Mutations, and Antibiotic resistance annotation pipeline](#pima-plasmid-integrations-mutations-and-antibiotic-resistance-annotation-pipeline)
   * [Installation](#installation)
      * [Conda environments](#conda-environments)
      * [Dependencies](#dependencies)
         * [Software](#software)
         * [Python modules and programs](#python-modules-and-programs)
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



# Installation

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
* LoRDEC: http://www.atgc-montpellier.fr/lordec/
* Medaka: https://github.com/nanoporetech/medaka
* miniasm: https://github.com/lh3/miniasm
* minimap2: https://github.com/lh3/minimap2
* nanopolish: https://github.com/jts/nanopolish
* parallel: https://www.gnu.org/software/parallel/
* Pilon: https://github.com/broadinstitute/pilon/wiki
* Porechop: https://github.com/rrwick/Porechop
* qcat: https://github.com/nanoporetech/qcat
* racon: https://github.com/isovic/racon
* SAMtools: http://htslib.org/
* SPAdes: http://cab.spbu.ru/software/spades/
* wtdbg2: https://github.com/ruanjue/wtdbg2

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