process PIMA_SINGLEPLEX {
    tag "$sample"
    
    input:
    // [sample, output, ont_fastq, pima_cmd]
    tuple val(sample), val(output), path(ont_fastq), val(pima_cmd)

    output:
    path(sample), emit: output_directory
    path("$sample/report/report.pdf")

    script:
    """
    python3 $pima_cmd --ont-fastq $ont_fastq --output $sample --threads ${task.cpus}
    """
}