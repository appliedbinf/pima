import sys
import os
from argparse import ArgumentParser, HelpFormatter

from Pima.pima_colors import Colors
from Pima.utils.settings import Settings


def parse_args(settings: Settings):
    # with open(os.path.join(settings.pima_path, "VERSION"), "r") as version_fp:
    #     VERSION = version_fp.read().strip()

    parser = ArgumentParser(
        allow_abbrev=False,
        prog="pima.py",
        add_help=False,
        description="""
                            P.I.M.A. bacterial genome analysis pipeline
                            """,
        formatter_class=lambda prog: HelpFormatter(
            prog, width=120, max_help_position=120
        ),
    )

    parser._optionals.title = "Help and version"
    parser.add_argument(
        "-h", "--help", action="store_true", help="Print this help and exit."
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        help="Print the software version.",
        version=f"PIMA microbial genome analysis pipeline (version {settings.pima_version})",
    )

    # Input arguments
    input_group = parser.add_argument_group("Input and basecalilng options")

    input_group.add_argument(
        "--ont-model",
        required=False,
        default="auto",
        metavar="<ONT BASECALLING MODEL",
        help="ONT model used for base calling. 'auto' will try and determine which basecalling model was used and if 'sup' is detected will tell flye to use '--nano-hq'. Available models can be listed using 'medaka tools list_models' as described by medaka (https://github.com/nanoporetech/medaka) (default : %(default)s)",
    )
    input_group.add_argument(
        "--ont-fastq",
        required=False,
        default=None,
        metavar="<FASTQ|GZ>, fastq_pass",
        help="File containing basecalled ONT reads (fastq|gz), or folder containing directories of demultiplexed basecalled reads (fastq_pass)",
    )
    input_group.add_argument(
        "--multiplexed",
        required=False,
        default=False,
        action="store_true",
        help="The ONT data are multiplexed (default : %(default)s)",
    )
    input_group.add_argument(
        "--nextflow",
        required=False,
        default=False,
        const=True,
        nargs='?',
        metavar='nextflow configs',
        help="Use nextflow to parallelize PiMA (otherwise samples will be run in serial). Can also provide additional commands to the nextflow run wrapped in single quotes ('), e.g. --nextflow '-c <my_hpc_nextflow.config> -resume' (default : %(default)s)",
    )
    input_group.add_argument(
        "--barcode-min-fraction",
        required=False,
        default=0.025,
        type=float,
        metavar="",
        help="The minimum fraction of data necessary to include a barcode in the analysis (default : %(default)s)",
    )
    input_group.add_argument(
        "--contamination",
        required=False,
        default=False,
        action="store_true",
        help="Use Kraken2 to look for contamination in the input read (default : %(default)s)",
    )
    input_group.add_argument(
        "--illumina-fastq",
        required=False,
        default=None,
        nargs="+",
        metavar="<R1> [R2]",
        help="Files containing R1 & R2 Illumina reads",
    )
    input_group.add_argument(
        "--genome",
        required=False,
        default=None,
        metavar="<GENOME_FASTA>",
        help="A genome FASTA file to be used in place of assembly",
    )

    output_group = parser.add_argument_group("Output options")
    output_group.add_argument(
        "--output",
        required=False,
        default=None,
        metavar="<OUTPUT_DIR>",
        help="Output directory for the analysis",
    )
    output_group.add_argument(
        "--overwrite",
        required=False,
        default=False,
        action="store_true",
        help="Overwrite an existing output directory (default : %(default)s)",
    )
    output_group.add_argument(
        "--keep-intermediates",
        required=False,
        default=False,
        action="store_true",
        help="Keep all intermediate files (mostly retains large bam/mpileup files) (default : %(default)s)"
    )

    # Assembly options
    assembly_group = parser.add_argument_group("Assembly options")
    assembly_group.add_argument(
        "--assembler",
        required=False,
        default="flye",
        type=str,
        choices=["flye", "flye_sup", "raven"],
        help="Assembler to use. Choose 'flye_sup' if you are providing R10 fastq data base called with a super_accuracy model and wish to use flye (default : %(default)s)",
    )
    assembly_group.add_argument(
        "--genome-size",
        required=False,
        default="estimate",
        type=str,
        metavar="<GENOME_SIZE>",
        help="Genome size estimate for the assembly & downsampling, can use M,m,K,k shorthand (e.g. 5.5m). If unknown, you can let pima estimate the genome size using 'estimate' (default : %(default)s)",
    )
    assembly_group.add_argument(
        "--assembly-coverage",
        required=False,
        default=200,
        type=int,
        metavar="<X>",
        help="Downsample the provided reads to this coverage (default : %(default)sX)",
    )
    assembly_group.add_argument(
        "--no-medaka",
        required=False,
        default=False,
        action="store_true",
        help="Skip Medaka polising of the ONT assembly (faster) (default : %(default)s)",
    )
    assembly_group.add_argument(
        "--illumina-polisher",
        required=False,
        default="pilon",
        choices=["pilon", "polypolish", "skip"],
        help="Polish the genome assembly using short-reads, to skip use '--illumina-polisher skip' (default : %(default)s)",
    )
    assembly_group.add_argument(
        "--only-assemble",
        required=False,
        default=False,
        action="store_true",
        help="Only carry out assembly steps with the given data; no downstream analysis.",
    )
    assembly_group.add_argument(
        "--no-assembly",
        required=False,
        default=False,
        action="store_true",
        help="Don't attempt to assembly/polish a given genome/set of reads (default : %(default)s)",
    )

    # Database/download options
    download_group = parser.add_argument_group("Database downloading arguments")
    download_group.add_argument(
        "--download",
        required=False,
        default=False,
        action="store_true",
        help="Attempt to download Kraken/Plasmid databases if not found locally."
        + "Use witout other options.",
    )

    # Plasmid options
    plasmid_group = parser.add_argument_group("Plasmid and vector search options")
    plasmid_group.add_argument(
        "--plasmids",
        required=False,
        default=False,
        action="store_true",
        help="Do a plasmid search (default : %(default)s)",
    )
    plasmid_group.add_argument(
        "--plasmid-database",
        required=False,
        default=settings.plasmid_database_default_fasta,
        metavar="<PLASMID_FASTA>",
        help="Path to a FASTA file with reference plasmid sequences",
    )
    ##add mob_suite (recommended by Gulvik) - https://github.com/phac-nml/mob-suite

    # AMR gene options
    amr_group = parser.add_argument_group("AMR gene search options")
    amr_group.add_argument(
        "--amr-database",
        required=False,
        default=settings.amr_database_default,
        metavar="<AMR_FASTA>",
        help="Path to a FASTA file with AMR gene sequences (default : %(default)s)",
    )
    amr_group.add_argument(
        "--no-amr",
        required=False,
        default=False,
        action="store_true",
        help="Skip AMR search (default : %(default)s)",
    )

    # Inc group options
    inc_group = parser.add_argument_group("Incompatibility group search options")
    inc_group.add_argument(
        "--inc-database",
        required=False,
        default=settings.inc_database_default,
        metavar="<INC_FASTA>",
        help="Path to a FASTA file with incompatibility group sequences (default : %(default)s)",
    )
    inc_group.add_argument(
        "--no-inc",
        required=False,
        default=False,
        action="store_true",
        help="Skip incompatibility group search (default : %(default)s)",
    )

    # Pull in custom feature sets
    other_feature_group = parser.add_argument_group("Other feature search options")
    other_feature_group.add_argument(
        "--feature",
        required=False,
        default=None,
        metavar="<FEATURE_FASTA>",
        action="append",
        help="Path to a FASTA file with feature sequences",
    )

    # Drawing options
    drawing_group = parser.add_argument_group("Drawing options")
    drawing_group.add_argument(
        "--no-drawing",
        required=False,
        default=False,
        action="store_true",
        help="Skip drawing of contigs & Features (default : %(default)s)",
    )
    drawing_group.add_argument(
        "--self-circos",
        required=False,
        default=False,
        action="store_true",
        help="Use the assembled genome as the reference to draw circos images. NOT recommended for Illumina assemblies (default : %(default)s)",
    )

    # Options for comparing to a reference genome
    reference_group = parser.add_argument_group("Reference options")
    reference_group.add_argument(
        "--reference-dir",
        required=False,
        default=settings.reference_dir_default,
        metavar="<REFERNCE_DIR>",
        help="Directory containing refrence organisms (default : %(default)s)",
    )
    reference_group.add_argument(
        "--organism",
        required=False,
        default=None,
        metavar="<Genus_species>",
        help="Reference organism to compare against",
    )
    reference_group.add_argument(
        "--list-organisms",
        required=False,
        default=False,
        action="store_true",
        help="List the reference organisms available to this pipeline",
    )
    reference_group.add_argument(
        "--reference-genome",
        required=False,
        default=False,
        metavar="<GENOME_FASTA>",
        nargs="?",
        help="Reference genome to compare against (default : %(default)s)",
    )
    reference_group.add_argument(
        "--mutation-regions",
        required=False,
        default=False,
        metavar="<REGION_BED>",
        nargs="?",
        help="Regions in the reference genome to screen for mutations (default : %(default)s)",
    )

    # Other arguments
    other_group = parser.add_argument_group("Other options")
    other_group.add_argument(
        "--name",
        required=False,
        type=str,
        default="Genome",
        metavar="<NAME>",
        help="Name of this analysis for reporting.",
    )
    other_group.add_argument(
        "--threads",
        required=False,
        type=int,
        default=1,
        metavar="<NUM_THREADS>",
        help="Number of worker threads to use (default : %(default)s)",
    )
    other_group.add_argument(
        "--verbosity",
        required=False,
        type=int,
        default=1,
        metavar="<INT>",
        help="How much information to print as PiMA runs (default : %(default)s)",
    )
    other_group.add_argument(
        "--logfile",
        required=False,
        type=str,
        default="pipeline.log",
        metavar="<FILENAME>",
        help="Name of logfile written in output directory (default : %(default)s)",
    )
    other_group.add_argument(
        "--resume",
        required=False,
        action="store_true",
        help="Restart pipeline from last completed step (default : %(default)s)",
    )
    other_group.add_argument(
        "--bundle",
        required=False,
        type=str,
        default=None,
        metavar="<PATH>",
        help="Local Tectonic bundle (default : %(default)s)",
    )
    other_group.add_argument(
        "--fake-run",
        required=False,
        default=False,
        action="store_true",
        help="Don't actually run the pipeline, just pretend to (default : %(default)s)",
    )

    opts, unknown_args = parser.parse_known_args()

    if opts.organism:
        opts.organism = opts.organism.replace(" ", "_")

    if opts.help:
        print(Colors.HEADER)
        parser.print_help()
        print(Colors.ENDC)
        sys.exit(0)

    return opts, unknown_args
