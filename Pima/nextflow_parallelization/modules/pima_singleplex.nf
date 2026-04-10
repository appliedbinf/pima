process PIMA_SINGLEPLEX {
    tag "$meta.sample_id"
    
    input:
    val meta

    output:
    path(meta.sample_id), emit: output_directory
    path("$meta.sample_id/report/report.pdf")

    script:
    def optional_arg_map = [
        ont_fastq : '--ont-fastq',
        illumina_fastq: '--illumina-fastq',
        genome    : '--genome',
        genome_size: '--genome-size',
        reference_organism: '--organism',
        reference_genome: '--reference-genome',
        reference_mutation_bed_file: '--mutation-regions',
    ]
    def optional_args = []

    optional_arg_map.each { csv_col, cli_flag -> 
        if (meta[csv_col]) {
            optional_args.add("${cli_flag} ${meta[csv_col]}")
        }
    }

    if (meta.self_circos == "yes") {
        optional_args.add("--self-circos")
    } 

    """
    python3 ${meta.pima_cmd} \\
        --threads ${task.cpus} \\
        --output ${meta.sample_id} \\
        ${optional_args.join(' ')}
    """
}