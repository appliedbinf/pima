import Bio.SeqIO
import datetime

from Pima.pima_colors import Colors
import pandas as pd

class PimaData:
    def __init__(self, opts=None, unknown_args=None):
        # The actual steps to carry out in the analysis held as a list
        self.analysis = []
        
        # Verbosity levels and colors
        self.error_color = Colors.FAIL
        self.fail_verbosity = 1
        self.main_process_verbosity = 1
        self.warning_color = Colors.WARNING
        self.warning_verbosity = 1
        self.main_process_color = Colors.OKGREEN
        self.sub_process_verbosity = 2
        self.sub_process_color = Colors.OKBLUE
        self.command_verbosity = 3
        self.errors = []
        self.warnings = []

        # ONT FASTQ input
        self.ont_fastq = None
        self.ont_raw_fastq = self.ont_fastq
        self.ont_read_count = None
        self.ont_read_lengths = None
        self.will_have_ont_fastq = False
        self.ont_read_lengths = []

        # Read metadata
        self.read_metadata = pd.Series(dtype=object)

        # Demultiplexing
        self.multiplexed = None
        self.nextflow = None
        self.barcodes = None
        self.barcode_min_fraction = None
        self.barcode_summary = None

        # Contamination
        self.contam_check = False
        self.kraken_fracs = pd.Series(dtype=object)
        self.did_contamination_check = False

        # Genome FASTA input
        self.genome_fasta = None
        self.will_have_genome_fasta = False

        # Illumina FASTQ input
        self.illumina_fastq = None
        self.pilon_coverage_min = 25
        self.did_spades_illumina_fastq = False
        self.did_pilon_ont_assembly = False
        self.did_polypolish_ont_assembly = False

        # Output options
        self.output_dir = None
        self.overwrite = False
        self.resume = False
        self.keep_intermediates = False

        # Assembly options
        self.assembler = "flye"
        self.flye_sup = False
        self.genome_assembly_size = None
        self.genome_assembly_raw_size = None
        self.assembly_coverage = None
        self.no_medaka = False
        self.ont_n50 = None
        self.ont_n50_min = 2500
        self.ont_coverage_min = 30
        self.only_assemble = False
        self.no_assembly = False
        self.did_flye_ont_fastq = False
        self.did_raven_ont_fastq = False
        self.will_have_ont_assembly = False
        self.mean_coverage = dict()
        self.did_circos_plots = False
        
        # ONT polishing
        self.ont_model = None
        self.did_medaka_ont_assembly = False

        # Illumina polishing
        self.illumina_polisher = "pilon"

        # Feature options
        self.no_amr = False
        self.no_inc = False
        self.feature_fastas = None
        self.feature_hits = pd.Series(dtype=object)
        self.feature_plots = pd.Series(dtype=object)
        self.feature_dirs = []
        self.feature_names = []
        self.feature_colors = []
        self.did_blast_feature_sets = False

        # Download options
        self.download = False

        # Reference options
        self.reference_dir = None
        self.organism = None
        self.organism_dir = None
        self.list_organisms = False
        self.will_have_reference_fasta = False
        self.reference = None
        self.amr_mutations = pd.Series(dtype=object)
        self.mutation_regions = None
        self.amr_region_names = None
        self.virulence_genes_fp = None
        self.did_call_mutations = False
        self.amr_deletions = pd.DataFrame()
        self.did_call_large_indels = False
        self.reference_contig_order = None
        self.organism_amr_appendices = None
        # Files to remove when done
        self.files_to_clean = []

        # Plasmid  options
        self.plasmids = False
        self.did_call_plasmids = False

        # Notes for different sections of the analysis
        self.assembly_notes = pd.Series(dtype=object)
        self.alignment_notes = pd.Series(dtype=object)
        self.large_indel_notes = pd.Series(dtype=object)
        self.contig_alignment = pd.Series(dtype=object)
        self.versions = pd.Series(dtype=object)

        self.logging_handle = None
        self.fake_run = False

        self.bundle = None
        self.report = pd.Series(dtype=object)

        if opts is None or unknown_args is None:
            return

        # Date-time information
        self.start_time = datetime.datetime.now().strftime("%Y-%m-%d")

        # Logging information
        self.logging_file = None
        self.logging_handle = None

        # ONT FASTQ input
        self.ont_fastq = opts.ont_fastq
        self.ont_raw_fastq = self.ont_fastq

        # Demultiplexing
        self.multiplexed = opts.multiplexed
        self.nextflow = opts.nextflow
        self.barcodes = None
        self.barcode_min_fraction = opts.barcode_min_fraction

        # Contamination
        self.contam_check = opts.contamination
        self.did_contamination_check = False

        # Illumina FASTQ input
        self.illumina_fastq = opts.illumina_fastq

        # Genome FASTA input
        self.genome_fasta = opts.genome

        # Output options
        self.output_dir = opts.output
        self.overwrite = opts.overwrite
        self.resume = opts.resume
        self.keep_intermediates = opts.keep_intermediates

        # Assembly options
        self.assembler = opts.assembler
        self.genome_assembly_size = opts.genome_size
        self.assembly_coverage = opts.assembly_coverage
        self.only_assemble = opts.only_assemble
        self.no_assembly = opts.no_assembly

        # ONT polishing
        self.ont_model = opts.ont_model
        self.no_medaka = opts.no_medaka

        # Illumina polishing
        self.illumina_polisher = opts.illumina_polisher

        # Illumina metrics
        self.illumina_length_mean = None
        self.illumina_coverage_min = 30
        self.did_pilon_ont_assembly = False
        self.did_polypolish_ont_assembly = False

        # The assembly itself
        self.genome = None
        self.contig_info = None

        # Vs. reference options
        self.reference_identity_min = 98.0
        self.reference_alignment_min = 97.0
        self.query_alignment_min = 97.0

        #placeholders for comparison to given reference
        self.reference_identity = 0
        self.reference_aligned_bases = 0
        self.query_aligned_bases = 0
        self.reference_aligned_fraction = 0
        self.query_aligned_fraction = 0

        # Plasmid and feature options
        self.plasmids = opts.plasmids
        self.plasmid_database = opts.plasmid_database
        self.did_call_plasmids = False
        self.no_drawing = opts.no_drawing
        self.amr_database = opts.amr_database
        self.no_amr = opts.no_amr
        self.inc_database = opts.inc_database
        self.no_inc = opts.no_inc
        self.feature_fastas = opts.feature
        self.feature_hits = pd.Series(dtype="float64")
        self.feature_plots = pd.Series(dtype="float64")
        self.feature_dirs = []
        self.feature_names = []
        self.feature_colors = []
        self.download = opts.download

        # Reference options
        self.reference_dir = opts.reference_dir
        self.organism = opts.organism
        self.list_organisms = opts.list_organisms
        self.reference_fasta = opts.reference_genome
        self.mutation_region_bed = opts.mutation_regions
        self.self_circos = opts.self_circos
        self.threads = opts.threads

        # How much stuff to print
        self.verbosity = opts.verbosity

        # Files to remove when done
        self.files_to_clean = []

        # Don't actully run any commands
        self.fake_run = opts.fake_run

        # Reporting
        self.no_report = False
        self.bundle = opts.bundle
        self.analysis_name = opts.name
        self.mutation_title = "Mutations"
        self.report[self.mutation_title] = pd.Series(dtype="float64")
        self.large_indels = pd.Series(dtype="float64")
        self.plasmid_title = "Plasmid annotation"
        self.report[self.plasmid_title] = pd.Series(dtype="float64")
        self.amr_matrix_title = "AMR matrix"
        self.did_draw_amr_matrix = False
        self.report[self.amr_matrix_title] = pd.Series(dtype="float64")
        self.methods_title = "Methods summary"
        self.report[self.methods_title] = pd.Series(dtype="float64")
        self.basecalling_methods = "Basecalling & processing"
        self.report[self.methods_title][self.basecalling_methods] = pd.Series(
            dtype="float64"
        )
        self.assembly_methods = "Assembly & polishing"
        self.report[self.methods_title][self.assembly_methods] = pd.Series(
            dtype="float64"
        )
        self.mutation_methods = "Mutation screening "
        self.report[self.methods_title][self.mutation_methods] = pd.Series(
            dtype="float64"
        )
        self.plasmid_methods = "Plasmid annotation"
        self.report[self.methods_title][self.plasmid_methods] = pd.Series(
            dtype="float64"
        )
        self.meta_title = "PIMA meta-information"

        # See if we got any unknown args.  Not allowed.
        if len(unknown_args) != 0:
            self.errors = self.errors + [
                "Unknown argument: " + unknown for unknown in unknown_args
            ]

    def load_reference(self):
        self.reference = self.load_fasta(self.reference_fasta)
        self.will_have_reference_fasta = True

        self.reference_size = 0
        for i in self.reference:
            self.reference_size += len(i.seq)

    @staticmethod
    def load_fasta(fasta: str):
        sequence = pd.Series(dtype=object)
        for contig in Bio.SeqIO.parse(fasta, "fasta"):
            sequence[contig.id] = contig
        return sequence

    def load_genome(self):
        self.genome = self.load_fasta(self.genome_fasta)
        self.genome_size = 0
        for i in self.genome:
            self.genome_size += len(i.seq)
