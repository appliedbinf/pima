
import pandas as pd

from Pima.pima_data import PimaData
from Pima.utils.utils import (
    print_and_run,
    print_and_log,
)

def validate_evaluate_assembly(pima_data: PimaData) :

    if pima_data.no_assembly:
        return
    
    if not pima_data.will_have_ont_assembly :
        return

    pima_data.analysis.append(['check_for_small_contigs_and_fragmentation', pima_data])

def check_for_small_contigs_and_fragmentation(pima_data: PimaData):
    """
    Will run on any assembly required, including denovo ONT, ONT polished, and/or Illumina polished assemblies
    Queue up for final step after all assembly manipulation steps are complete
    """
    print_and_log(
        pima_data,
        "Evaluating assembly",
        pima_data.main_process_verbosity,
        pima_data.main_process_color,
    )
    genome_sizes = pima_data.genome_fasta.replace(".fasta", ".sizes")
    command = " ".join(
        ["faidx -i chromsizes", pima_data.genome_fasta, ">", genome_sizes]
    )
    print_and_run(pima_data, command)

    assembly_info = pd.read_csv(genome_sizes, sep="\t", header=None)
    assembly_info.columns = ["contig", "length"]
    pima_data.contig_sizes = assembly_info

    # Take a look at the number of contigs, their sizes, and circularity.  Warn if things don't look good
    if assembly_info.shape[0] > 4:
        warning = f"Assembly produced {assembly_info.shape[0]} contigs, more than ususally expected; assembly may be fragmented."
        print_and_log(
            pima_data, warning, pima_data.warning_verbosity, pima_data.warning_color
        )
        pima_data.assembly_notes = pd.concat([pima_data.assembly_notes, pd.Series(warning, dtype='object')])
    small_contigs = assembly_info.loc[assembly_info["length"] <= 3000, :]

    if small_contigs.shape[0] > 0:
        warning = f"Assembly produced {small_contigs.shape[0]} small contigs ({', '.join(small_contigs['contig'])}); assembly may include spurious sequences."
        print_and_log(
            pima_data, warning, pima_data.warning_verbosity, pima_data.warning_color
        )
        pima_data.assembly_notes = pd.concat([pima_data.assembly_notes, pd.Series(warning, dtype='object')])