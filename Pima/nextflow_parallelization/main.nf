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
        .splitCsv(header:false, sep:",")

    PIMA_SINGLEPLEX( samples )

    COPY_RESULTS( PIMA_SINGLEPLEX.out.output_directory, outdir )

}

workflow.onComplete {
    println "Project : $workflow.projectDir"
    println "Workdir : $workflow.workDir"
    println "homeDir : $workflow.homeDir"
    println "launchDir : $workflow.launchDir"
}