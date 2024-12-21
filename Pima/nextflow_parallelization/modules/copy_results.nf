process COPY_RESULTS {
    input:
    path(sample_outdir)
    val(workflow_outdir)

    output:
    stdout

    script:
    """
    #move the pima results folder to the original output location
    cp -Lr $sample_outdir $workflow_outdir/
    """
}