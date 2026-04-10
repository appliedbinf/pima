#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { PIMA_SINGLEPLEX } from './modules/pima_singleplex.nf'
include { COPY_RESULTS } from './modules/copy_results.nf'

params.sample_sheet = null
params.output = null

workflow {
    outdir = file( params.output )

    samples = Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header:true, sep:",")
        .map { row -> 
            [ 
                sample_id                    : row.sample_name,
                output_directory             : row.output_directory,
                ont_fastq                    : row.ont_fastq,
                illumina_fastq               : row.illumina_r1,
                genome                       : row.genome,
                genome_size                  : row.genome_size,
                reference_organism           : row.reference_organism,
                reference_genome             : row.reference_genome,
                reference_mutation_bed_file  : row.reference_mutation_bed_file,
                self_circos                  : row.self_circos,
                pima_cmd                     : row.original_command
            ]
        }

    PIMA_SINGLEPLEX( samples )

    COPY_RESULTS( PIMA_SINGLEPLEX.out.output_directory, outdir )

}

workflow.onComplete {
    println "Project : $workflow.projectDir"
    println "Workdir : $workflow.workDir"
    println "homeDir : $workflow.homeDir"
    println "launchDir : $workflow.launchDir"
}