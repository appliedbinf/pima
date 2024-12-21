from .download_references import (
    validate_download,
    validate_organism,
)
from .outdir import validate_output_dir
from .multiplexed import validate_multiplex_fastq, initialize_multiplex_analysis
from .fastq import (
    validate_ont_fastq, 
    info_given_ont_fastq,
    validate_illumina_fastq,
    info_illumina_fastq,
    validate_genome_estimate,
    estimate_genome_size,
)
from .check_contamination import validate_contamination_check, fastq_contamination

from .assembly import (
    validate_genome_fasta,
    validate_genome_assembly_size,
    validate_assembly_info,
    validate_assembler,
    flye_ont_fastq,
    raven_ont_fastq,
    spades_illumina_fastq,
    check_assembly_coverages,
)

from .ont_polishing import (
    validate_medaka,
    medaka_ont_assembly,
)

from .illumina_polishing import(
    validate_illumina_polish,
    pilon_assembly,
    polypolish_assembly,
)

from .evaluate_assembly import(
    validate_evaluate_assembly,
    check_for_small_contigs_and_fragmentation,
)

from .annotations import (
    validate_features,
    validate_blast,
    blast_feature_sets,
)

from .plasmids import (
    validate_plasmids,
    call_plasmids,
)

from .compare_to_ref import (
    validate_reference_fasta,
    validate_mutations,
    validate_quast,
    call_insertions,
    quast_genome,
    call_amr_mutations,
)

from .visualizations import (
    validate_draw_amr_matrix,
    validate_draw_features,
    validate_draw_circos,
    draw_features,
    draw_amr_matrix,
    draw_circos,
)

from .report import (
    validate_make_report,
    make_report,
)

__all__ = [
    # Validation
    "validate_download",
    "validate_organism",
    "validate_output_dir",
    "validate_multiplex_fastq",
    "validate_genome_estimate",
    "validate_ont_fastq",
    "validate_illumina_fastq",
    "validate_contamination_check",
    "validate_genome_fasta",
    "validate_genome_assembly_size",
    "validate_assembly_info", #check for coverages
    "validate_assembler",
    "validate_medaka",
    "validate_illumina_polish",
    "validate_evaluate_assembly", #check for contig size/# issues
    "validate_plasmids",
    "validate_features",
    "validate_blast",
    "validate_reference_fasta",
    "validate_quast",
    "validate_mutations",
    "validate_draw_amr_matrix",
    "validate_draw_features",
    "validate_draw_circos",
    "validate_make_report",
    # Analysis
    "initialize_multiplex_analysis",
    "estimate_genome_size",
    "info_given_ont_fastq",
    "info_illumina_fastq",
    "fastq_contamination",
    "flye_ont_fastq",
    "raven_ont_fastq",
    "spades_illumina_fastq",
    "medaka_ont_assembly",
    "pilon_assembly",
    "polypolish_assembly",
    "check_assembly_coverages",
    "check_for_small_contigs_and_fragmentation",
    "call_plasmids",
    "blast_feature_sets",
    "call_insertions",
    "quast_genome",
    "call_amr_mutations",
    "draw_features",
    "draw_amr_matrix",
    "draw_circos",
    "make_report",
]
