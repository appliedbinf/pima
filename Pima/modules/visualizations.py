import re
import os

from importlib.metadata import version
import pandas as pd
from dna_features_viewer import GraphicFeature, GraphicRecord
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np 

from Pima.accessory_scripts.building_pycircos_figures import (
    BuildCircosPlots,
)

from Pima.pima_data import PimaData
from Pima.utils.settings import Settings
from Pima.utils.utils import (
    print_and_log,
    print_and_run,
    find_checkpoint,
    make_start_file,
    make_finish_file,
    add_warning,
)
from Pima.utils.mapping import (
    minimap_and_sort,
    mpileup_bam,
)
 
def validate_draw_amr_matrix(pima_data: PimaData):
    # skip conditions
    if pima_data.only_assemble:
        return
    if pima_data.no_drawing:
        return
    if not (pima_data.will_have_genome_fasta or pima_data.reference_fasta):
        return
    
    try:
        pima_data.versions['matplotlib'] = version('matplotlib')
    except:
        add_warning(pima_data, "Matplotlib version could not be determined")

    pima_data.analysis.append(["draw_amr_matrix", pima_data])


def validate_draw_features(pima_data: PimaData):
    # skip conditions
    if pima_data.only_assemble:
        return
    if pima_data.no_drawing:
        return
    if not pima_data.will_have_genome_fasta:
        return

    try:
        pima_data.versions['dna_features_viewer'] = version('dna_features_viewer')
    except:
        add_warning(pima_data, "dna_features_viewer version could not be determined")

    pima_data.analysis.append(["draw_features", pima_data])


def validate_draw_circos(pima_data: PimaData, settings: Settings):
    # skip conditions
    if pima_data.only_assemble:
        return
    if pima_data.no_drawing:
        return

    if pima_data.self_circos or pima_data.will_have_reference_fasta:
        pima_data.analysis.append(["draw_circos", pima_data, settings])


def draw_features(pima_data: PimaData):

    if pima_data.unique_hits.empty:
        return

    print_and_log(
        pima_data,
        'Drawing features', 
        pima_data.main_process_verbosity, 
        pima_data.main_process_color,
    )

    pima_data.drawing_dir = os.path.join(pima_data.output_dir, 'drawing')
    
    # Check if results exist that we can use
    if find_checkpoint(pima_data, pima_data.drawing_dir):
        feature_figs = [fig for fig in os.listdir(pima_data.drawing_dir) if fig.endswith('.svg') and not fig == "amr_matrix.svg"]
        for fig in feature_figs:
            pima_data.feature_plots[fig] = os.path.join(pima_data.drawing_dir, fig)
        return
    os.mkdir(pima_data.drawing_dir)
    make_start_file(pima_data, pima_data.drawing_dir)
    
    figure_width = 13

    # Draw one plot per contig for simplicity
    for contig in pima_data.genome:
        
        contig_plot_svg = os.path.join(pima_data.drawing_dir, f"{contig.id}.svg")
        #contig_plot_png = os.path.join(pima_data.drawing_dir, f"{contig.id}.png")
        
        contig_features = pima_data.unique_hits[pima_data.unique_hits['chrom'] == contig.id]
        if contig_features.empty:
            continue

        #contains amr & inc, nested in each is a list of graphic features to be drawn
        #amr: GF(tetL, start-end), GF(bla, start-end)...
        feature_sets_to_plot = pd.Series(dtype = object)

        print_and_log(
            pima_data,
            'Drawing features for ' + contig.id, 
            pima_data.sub_process_verbosity, 
            pima_data.sub_process_color,
        )
        
        for feature_type in contig_features['source'].unique(): # amr / inc
            feature_plots = []
            features = contig_features[contig_features['source'] == feature_type]
            
            if features.empty:
                continue

            for i in range(features.shape[0]):
                feature = features.iloc[i,]
                feature_plots.append(
                    GraphicFeature(
                        start = feature['start'], 
                        end = feature['end'], 
                        label=feature['gene'],
                        strand=f"{feature['strand']}1",
                        color = getattr(pima_data, f"{feature['source']}_color"),
                        box_color = getattr(pima_data, f"{feature['source']}_color"),
                    )
                )
            feature_sets_to_plot[feature_type] = feature_plots

        if len(feature_sets_to_plot) == 0:
            continue
            
        # Add blank feature sets for the header and ruler
        real_sets = feature_sets_to_plot.index.tolist() # = amr and/or inc
        empty_set = [GraphicFeature(start = 1, end = len(contig), color = '#FFFFFF')]

        # Figure out high each plot will be on its own for later scaling
        expected_plot_heights = []
        for i in range(len(feature_sets_to_plot)):
            record = GraphicRecord(sequence_length = len(contig), features = feature_sets_to_plot.iloc[i])

            with_ruler  = False
            if i == len(feature_sets_to_plot) - 1:
                with_ruler = True

            plot, _ = record.plot(figure_width = figure_width, with_ruler = with_ruler)
            expected_plot_heights += [plot.figure.get_size_inches()[1]]
            
        plot_height_sum = sum(expected_plot_heights)

        # Make a figure with separate plots for each feature class
        plots = plt.subplots(nrows = len(feature_sets_to_plot), ncols = 1, sharex = True,
                                figsize=(figure_width, plot_height_sum * .66666),
                                gridspec_kw={"height_ratios": expected_plot_heights})
        figure = plots[0]
        plots = plots[1]
        if len(feature_sets_to_plot) == 1: 
            plots = [plots]

        # Add each feature class's plot with the pre-determined height
        for i in range(len(feature_sets_to_plot)):
            record = GraphicRecord(sequence_length = len(contig), features = feature_sets_to_plot.iloc[i])

            with_ruler  = False
            if i == len(feature_sets_to_plot) - 1:
                with_ruler = True
                
            plot, _ = record.plot(ax = plots[i], with_ruler = with_ruler, figure_width = figure_width)
            ymin, ymax = plot.figure.axes[0].get_ylim()

            if i == 0: 
                plot.text(x = 0, y = ymax, s = contig.id) 
            
        figure.tight_layout()
        figure.savefig(contig_plot_svg)
        #figure.savefig(contig_plot_png, dpi=300)

        pima_data.feature_plots[contig.id] = contig_plot_svg

    make_finish_file(pima_data, pima_data.drawing_dir)


def draw_amr_matrix(pima_data: PimaData):

    have_gene_hits = not pima_data.unique_hits.empty #only run this is we have variant_data OR gene_data to draw
    if not (pima_data.did_call_mutations or have_gene_hits):
        return

    print_and_log(
        pima_data,
        "Drawing AMR matrix", 
        pima_data.main_process_verbosity, 
        pima_data.main_process_color,
    )

    pima_data.drawing_dir = os.path.join(pima_data.output_dir, 'drawing')
    if not os.path.isdir(pima_data.drawing_dir):
        os.mkdir(pima_data.drawing_dir)

    amr_to_draw = pd.DataFrame(columns =['gene', 'drug'])
    amr_df = pd.DataFrame(columns =['gene', 'drug'])
    # Roll up Resfinder / AMRFinder gene hits
    if not pima_data.unique_hits.empty:
        amr_df = pd.concat(
            [
                amr_df,
                pima_data.unique_hits[pima_data.unique_hits['source'] == 'amr'][['gene', 'class']].set_axis(['gene', 'drug'], axis=1)
            ]
        )

    # Roll up potentially resistance conferring mutations
    snp_df = pd.DataFrame(columns =['gene', 'drug'])
    indel_df = pd.DataFrame(columns =['gene', 'drug'])

    if pima_data.did_call_mutations and pima_data.amr_mutations.shape[0] > 0:
        mutations = pima_data.amr_mutations.copy(deep=True)
        ## Duplicated from the report.py PimaReport class, probably should consolidate into compare_to_ref, but wanted to keep all the data in the pima_data object for now
        #restrict to only known mutations if the input genome is too far from the reference genome
        if pima_data.reference_identity < pima_data.reference_identity_min and pima_data.reference_identity != 0:
            mutations = mutations.query('classification_type == "confirmed_location"')
            # hide homopolymer indels from ONT variant results UNLESS an R10 SUP model was used
        if not re.search(r"r1041.*sup.*", pima_data.ont_model):
            ##remove the +/- from the variant call & check if it is a homopolymer (same base), this does replace SNPs with NA, but we sort that out later
            mutations = mutations.copy(deep=False)
            mutations['var_mod'] = mutations['var'] \
                                        .astype(str) \
                                        .apply(lambda x: f'assembly_del{x}' if (x.isnumeric()) else x) \
                                        .apply(lambda x: ''.join(filter(str.isalpha, x))) \
                                        .apply(lambda x: x if (x != len(x)*x[0]) else 'na')

            mutations = mutations[~((mutations['var'].str.match(r"[+-]")) & (mutations['var_mod'] == 'na'))] #drop homopolymer indels
        column_names = ["gene", "Contig(ref)", "Pos(ref)", "drug"]
        #SNPs
        snp_df = mutations.query("var_type == 'snp'")
        snp_df = snp_df.assign(gene = snp_df['region_name'] + " " + snp_df['ref'] + " -> " + snp_df['var'])
        snp_df = snp_df[["gene", "GE", "loc", "amr_class"]]
        snp_df = snp_df.set_axis(column_names, axis=1)
        #indels [includes deletion-only indels]
        indel_df = mutations.query("var_type == 'indel' | var_type == 'large-indel'")
        indel_df = indel_df.assign(var = indel_df['var'].apply(lambda x: f"{x[0:5]}...{len(x)-5}" if len(x) > 5 else x))
        indel_df = indel_df.assign(gene = indel_df['region_name'] + " " + indel_df['ref'] + " -> " + indel_df['var'])
        indel_df = indel_df[["gene", "GE", "loc", "amr_class"]]
        indel_df = indel_df.set_axis(column_names, axis=1)

        #Build matrix
        amr_to_draw = pd.concat([
            amr_df[['gene', 'drug']],
            snp_df[['gene', 'drug']],
            indel_df[['gene', 'drug']],
        ])
        amr_to_draw.to_csv(
            os.path.join(
                pima_data.mutations_dir, 
                "amr_mutations_to_draw.csv",
            ), 
            index = False, 
        )
    elif amr_df.shape[0] > 0:
        #Build matrix
        amr_to_draw = amr_df[['gene', 'drug']]
        amr_to_draw.to_csv(
            os.path.join(
                pima_data.drawing_dir, 
                "amr_mutations_to_draw.csv",
            ), 
            index = False, 
        )
    # If there are no AMR hits, we can't draw the matrix
    if amr_to_draw.shape[0] < 1:
        return
    
    present_genes = amr_to_draw['gene'].unique()
    present_drugs = amr_to_draw['drug'].unique()
                        
    amr_matrix = pd.DataFrame(0, index = present_genes, columns = present_drugs)
    for hit_idx, hit in amr_to_draw.iterrows():
        amr_matrix.loc[hit['gene'], hit['drug']] = 1
    amr_matrix = amr_matrix.sort_values(list(present_drugs), ascending=False)  
    amr_matrix_fig = os.path.join(pima_data.drawing_dir, 'amr_matrix.svg')
    int_matrix = amr_matrix[amr_matrix.columns].astype(int)

    figure, axis = plt.subplots()
    c = mpl.colors.ListedColormap(['white', 'blue']) #setup colorscale
    n = mpl.colors.Normalize(vmin=0, vmax=1) #define color scale, correctly display if there are no 0s in int_matrix
    axis.pcolor(int_matrix, cmap = c, norm=n, linewidth = 0) #setup heatmap
    axis.invert_yaxis()
    axis.set_yticks(np.arange(0.5, len(amr_matrix.index)), minor = False)
    axis.set_yticklabels(int_matrix.index.values)

    axis.set_xticks(np.arange(0.5, len(amr_matrix.columns)), minor = False)
    axis.set_xticklabels(amr_matrix.columns.values, rotation = 90)
    axis.xaxis.tick_top()
    axis.xaxis.set_label_position('top')

    # scale for variable number of mutations
    num_mutations = axis.get_yticklabels()
    maxsize = max([i_mutation.get_window_extent().height for i_mutation in num_mutations])
    plot_height = maxsize/figure.dpi*len(amr_matrix.index)+2
    #margin = m/figure.get_size_inches()[1]

    ## If size is too big, drop the dpi
    if plot_height > 8:
        set_dpi = 8/plot_height*100
    else:
        set_dpi = 100

    figure.set_size_inches(figure.get_size_inches()[0], plot_height)
    figure.tight_layout()

    plt.savefig(amr_matrix_fig)
    ##figure.savefig(amr_matrix_png, dpi = set_dpi)
    pima_data.amr_matrix_fig = amr_matrix_fig
    pima_data.did_draw_amr_matrix = True 


def draw_circos(pima_data: PimaData, settings: Settings):

    #only plot virulence genes when the reference is Ba and it is present in the settings
    virulence_genes_fp = None

    if pima_data.self_circos:
        print_and_log(
            pima_data,
            f"Drawing Circos plot of denovo assembly",
            pima_data.main_process_verbosity, 
            pima_data.main_process_color,
        )
    #don't fail this if the scores are 0 since that just means we didn't measure them
    elif (pima_data.reference_identity < 95. and pima_data.reference_identity != 0) or (pima_data.reference_aligned_fraction < 70 and pima_data.reference_aligned_fraction != 0):
        message = (
            f"Circos plots supressed because genome is less than 95% nucleotide similarity to provided reference genome ({pima_data.reference_identity})."
            f" or genome shares less than 70% of the bases present in the reference genome ({pima_data.reference_aligned_fraction})"
        )
        add_warning(pima_data, message)
        return
    elif pima_data.organism == "Bacillus_anthracis" and settings.virulence_genes_fp is not None:
        print_and_log(
            pima_data,
            f"Drawing Circos plot using Bacillus anthracis, input read coverage mapped to B. anthracis, and the virulence genes",
            pima_data.main_process_verbosity, 
            pima_data.main_process_color,
        )
        virulence_genes_fp = settings.virulence_genes_fp
    elif not pima_data.genome_fasta:
        print_and_log(
            pima_data,
            f"Drawing Circos plot of reference and input read coverage mapped to reference",
            pima_data.main_process_verbosity, 
            pima_data.main_process_color,
        )
    else:
        print_and_log(
            pima_data,
            f"Drawing Circos plot of assembly vs. provided reference",
            pima_data.main_process_verbosity, 
            pima_data.main_process_color,
        )

    ont_reference_coverage_fp = None
    illumina_reference_coverage_fp = None

    # Make the directory for drawings
    pima_data.circos_dir = os.path.join(pima_data.output_dir, 'circos')
    if find_checkpoint(pima_data, pima_data.circos_dir):
        contig_dirs = [ dir.path for dir in os.scandir(pima_data.circos_dir) if dir.is_dir() ]
        if not pima_data.reference_contig_order:
            pima_data.reference_contig_order = [contig for contig in pd.read_csv(os.path.join(pima_data.circos_dir, 'genome.sizes'), sep="\t", header = None).iloc[:,0]]
        sorted_contig_dirs = [path for contig in pima_data.reference_contig_order for path in contig_dirs if contig in path] # or re.search(contig, path)
        for contig_dir in sorted_contig_dirs:
            contig = os.path.basename(contig_dir)
            circos_svg = os.path.join(contig_dir, 'circos.svg')
            pima_data.contig_alignment[contig] = circos_svg
        pima_data.did_circos_plots = True
        return
    
    os.makedirs(pima_data.circos_dir)
    make_start_file(pima_data, pima_data.circos_dir)
    
    # If no reference was given but user indicated circos plots should be generated
    if pima_data.self_circos:
        #specify our reference
        pima_data.reference_fasta = pima_data.genome_fasta
        #genome_sizes = re.sub('.fasta', '.sizes', pima_data.reference_fasta)
        reference_sizes_fp = os.path.join(pima_data.circos_dir, 'genome.sizes')
        command = " ".join(
            [
                'faidx -i chromsizes', 
                pima_data.reference_fasta, 
                '>', reference_sizes_fp
            ]
        )
        print_and_run(pima_data, command)

        contig_alignment_fp = None
        if pima_data.ont_fastq:
            pima_data.reference_mapping_ont_bam = os.path.join(pima_data.info_dir, 'ont_coverage.bam')
            ont_reference_coverage_fp = os.path.join(pima_data.circos_dir, 'ont_self_mapping.mpileup')
            mpileup_bam(
                pima_data,
                pima_data.reference_fasta,
                pima_data.reference_mapping_ont_bam, 
                ont_reference_coverage_fp, 
                pima_data.circos_dir)

        if pima_data.illumina_fastq: #We have the illumina mapped data already
            pima_data.reference_mapping_illumina_bam = os.path.join(pima_data.info_dir, 'illumina_coverage_unfilt.bam')
            illumina_reference_coverage_fp = os.path.join(pima_data.circos_dir, 'illumina_self_mapping_unfilt.mpileup')
            mpileup_bam(
                pima_data,
                pima_data.reference_fasta,
                pima_data.reference_mapping_illumina_bam, 
                illumina_reference_coverage_fp, 
                pima_data.circos_dir
            )

    elif not pima_data.genome_fasta: #denovo assembly turned off, just plot the reference & read mapping coverage
        genome_sizes = re.sub('.fasta', '.sizes', pima_data.reference_fasta)
        reference_sizes_fp = os.path.join(pima_data.circos_dir, os.path.basename(genome_sizes))
        command = " ".join(
            [
                'faidx -i chromsizes', 
                pima_data.reference_fasta, 
                '>', reference_sizes_fp
            ]
        )
        print_and_run(pima_data, command)
        contig_alignment_fp = None

        if pima_data.ont_fastq:
            pima_data.reference_mapping_ont_bam = os.path.join(pima_data.circos_dir, 'ont_coverage.bam')
            ont_reference_coverage_fp = os.path.join(pima_data.circos_dir, 'ont_self_mapping.mpileup')
            minimap_and_sort(
                    pima_data,
                    pima_data.reference_fasta,
                    pima_data.reference_mapping_ont_bam,
                    pima_data.ont_fastq,
                    ont = True,
                )
            mpileup_bam(
                pima_data,
                pima_data.reference_fasta,
                pima_data.reference_mapping_ont_bam, 
                ont_reference_coverage_fp, 
                pima_data.circos_dir)

        if pima_data.illumina_fastq and pima_data.mutation_region_bed: #We have the illumina mapped data already
            pima_data.reference_mapping_illumina_bam = os.path.join(pima_data.mutations_dir, 'reference_mapping_illumina_unfilt.bam')
            illumina_reference_coverage_fp = os.path.join(pima_data.circos_dir, 'illumina_ref_mapping_unfilt.mpileup')
            mpileup_bam(
                pima_data,
                pima_data.reference_fasta,
                pima_data.reference_mapping_illumina_bam, 
                illumina_reference_coverage_fp, 
                pima_data.circos_dir
            )

        elif pima_data.illumina_fastq and not pima_data.mutation_region_bed: #We don't have the Illumina data mapped
                pima_data.reference_mapping_illumina_bam = os.path.join(pima_data.circos_dir, 'reference_mapping_illumina_unfiltered.bam')
                minimap_and_sort(
                    pima_data,
                    pima_data.reference_fasta,
                    pima_data.reference_mapping_illumina_bam,
                    pima_data.illumina_fastq,
                    ont = False,
                )
            
                # Make an mpileup file out of the reads
                reference_mapping_mpileup = os.path.join(pima_data.circos_dir, 'reference_mapping_illumina.mpileup')
                mpileup_bam(
                    pima_data,
                    pima_data.reference_fasta,
                    pima_data.reference_mapping_illumina_bam, 
                    reference_mapping_mpileup, 
                    pima_data.circos_dir
                )
                illumina_reference_coverage_fp = os.path.join(pima_data.circos_dir, 'reference_mapping_illumina.mpileup')

    else: #standard path, use reference genome
        # Use the reference.sizes file to guide the cicros build
        reference_sizes_fp = os.path.join(pima_data.insertions_dir, 'reference.sizes')
        reference_1coords_fp = os.path.join(pima_data.insertions_dir, 'vs_reference.1coords')

        # if no mutation_regions.bed file provided a mutations_dir is not created and we do not have coverage data
        # however, if a reference genome was provided we can still draw the circos plots
        if not pima_data.mutation_region_bed:
            #have both files, but ONT has not been mapped to the reference yet
            if pima_data.ont_fastq and pima_data.illumina_fastq:
                #Illumina
                pima_data.reference_mapping_illumina_bam = os.path.join(pima_data.circos_dir, 'reference_mapping_illumina_unfiltered.bam')
                minimap_and_sort(
                    pima_data,
                    pima_data.reference_fasta,
                    pima_data.reference_mapping_illumina_bam,
                    pima_data.illumina_fastq,
                    ont = False,
                )
            
                # Make an mpileup file out of the reads
                reference_mapping_mpileup = os.path.join(pima_data.circos_dir, 'reference_mapping_illumina.mpileup')
                mpileup_bam(
                    pima_data,
                    pima_data.reference_fasta,
                    pima_data.reference_mapping_illumina_bam, 
                    reference_mapping_mpileup, 
                    pima_data.circos_dir
                )
                illumina_reference_coverage_fp = os.path.join(pima_data.circos_dir, 'reference_mapping_illumina.mpileup') 

                #ONT
                pima_data.reference_mapping_ont_bam = os.path.join(pima_data.circos_dir, 'reference_mapping_ont.bam')
                minimap_and_sort(
                    pima_data,
                    pima_data.reference_fasta, 
                    pima_data.reference_mapping_ont_bam,
                    pima_data.ont_fastq, 
                    ont = True,
                )

                # Make an mpileup file out of the reads
                reference_mapping_mpileup = os.path.join(pima_data.circos_dir, 'reference_mapping_ont.mpileup')
                mpileup_bam(
                    pima_data,
                    pima_data.reference_fasta, 
                    pima_data.reference_mapping_ont_bam, 
                    reference_mapping_mpileup, 
                    pima_data.circos_dir,
                )
                ont_reference_coverage_fp = os.path.join(pima_data.circos_dir, 'reference_mapping_ont.mpileup')
            
            #Illumina only                               
            elif pima_data.illumina_fastq:
                #We want to use the unfiltered bam file to prevent apparent coverage gaps due to repeats
                pima_data.reference_mapping_illumina_bam = os.path.join(pima_data.circos_dir, 'reference_mapping_illumina_unfiltered.bam')
                minimap_and_sort(
                    pima_data,
                    pima_data.reference_fasta,
                    pima_data.reference_mapping_illumina_bam,
                    pima_data.illumina_fastq,
                    ont = False,
                )
            
                # Make an mpileup file out of the reads
                reference_mapping_mpileup = os.path.join(pima_data.circos_dir, 'reference_mapping_illumina.mpileup')
                mpileup_bam(
                    pima_data,
                    pima_data.reference_fasta,
                    pima_data.reference_mapping_illumina_bam, 
                    reference_mapping_mpileup, 
                    pima_data.circos_dir
                )
                illumina_reference_coverage_fp = os.path.join(pima_data.circos_dir, 'reference_mapping_illumina.mpileup')

            #ONT only - but didn't already map to reference genome because no mutation file provided
            elif pima_data.ont_fastq:
                pima_data.reference_mapping_ont_bam = os.path.join(pima_data.circos_dir, 'reference_mapping_ont.bam')
                minimap_and_sort(
                    pima_data,
                    pima_data.reference_fasta, 
                    pima_data.reference_mapping_ont_bam,
                    pima_data.ont_fastq, 
                    ont = True,
                )

                # Make an mpileup file out of the reads
                reference_mapping_mpileup = os.path.join(pima_data.circos_dir, 'reference_mapping_ont.mpileup')
                mpileup_bam(
                    pima_data,
                    pima_data.reference_fasta, 
                    pima_data.reference_mapping_ont_bam, 
                    reference_mapping_mpileup, 
                    pima_data.circos_dir,
                )
                ont_reference_coverage_fp = os.path.join(pima_data.circos_dir, 'reference_mapping_ont.mpileup')
        
        #Have mutation file -> so ONT-only reads have been mapped to call variants
        else:
            ## ILLUMINA READS
            # Illumina reads used to call mutations if within the minimum perc identity
            ## however, if we did this then we didn't use the ONT reads to call mutations, so we'll need to do so
            ### Need to use the unfiltered bam file for the coverage plot so repeat regions don't look like missing data
            if pima_data.illumina_fastq and pima_data.reference_identity >= pima_data.reference_identity_min:
                illumina_unfiltered_bam = os.path.join(pima_data.mutations_dir, 'reference_mapping_illumina_unfilt.bam')
                illumina_reference_coverage_fp = os.path.join(pima_data.mutations_dir, 'reference_mapping_illumina_unfilt.mpileup')
                mpileup_bam(
                    pima_data,
                    pima_data.reference_fasta,
                    illumina_unfiltered_bam,
                    illumina_reference_coverage_fp,
                    pima_data.mutations_dir,
                )

            #we provided illumina data, but we didn't call mutations because the genomes were too distant
            elif (isinstance(pima_data.illumina_fastq, list) and pima_data.illumina_fastq) and (95.0 <= pima_data.reference_identity <= pima_data.reference_identity_min):
                illumina_unfiltered_bam = os.path.join(pima_data.circos_dir, 'reference_mapping_illumina_unfilt.bam')
                illumina_reference_coverage_fp = os.path.join(pima_data.circos_dir, 'reference_mapping_illumina.mpileup')
                minimap_and_sort(
                    pima_data,
                    pima_data.reference_fasta,
                    illumina_unfiltered_bam,
                    pima_data.illumina_fastq,
                    ont = False,
                )
                mpileup_bam(
                    pima_data,
                    reference_genome = pima_data.reference_fasta,
                    bam = illumina_unfiltered_bam, 
                    mpileup = illumina_reference_coverage_fp, 
                    output_dir = pima_data.circos_dir,
                )
                pima_data.files_to_clean.extend(illumina_unfiltered_bam)
                pima_data.files_to_clean.append(reference_mapping_mpileup)            
            
            ## ONT READS
            ## didn't call variants so now need to map coverage to reference genome, handled the illumina cov in above statement
            if pima_data.will_have_ont_fastq and pima_data.illumina_fastq:
                pima_data.reference_mapping_ont_bam = os.path.join(pima_data.circos_dir, 'reference_mapping_ont.bam')
                minimap_and_sort(
                    pima_data,
                    pima_data.reference_fasta, 
                    pima_data.reference_mapping_ont_bam,
                    pima_data.ont_fastq, 
                    ont = True,
                )

                # Make an mpileup file out of the reads
                reference_mapping_mpileup = os.path.join(pima_data.circos_dir, 'reference_mapping_ont.mpileup')
                mpileup_bam(
                    pima_data,
                    pima_data.reference_fasta, 
                    pima_data.reference_mapping_ont_bam, 
                    reference_mapping_mpileup, 
                    pima_data.circos_dir,
                )
                ont_reference_coverage_fp = os.path.join(pima_data.circos_dir, 'reference_mapping_ont.mpileup')       

            # if we have coverage data from ONT fastq reads, and we called amr mutations we can add coverage 
            elif pima_data.will_have_ont_fastq and pima_data.reference_identity >= pima_data.reference_identity_min:
                ont_reference_coverage_fp = os.path.join(pima_data.mutations_dir, 'reference_mapping_ont.mpileup')

            # genome is too divergent from reference to call mutations or we only used the illumina data, but still above 95%, so we can map
            elif pima_data.will_have_ont_fastq and (95.0 <= pima_data.reference_identity <= pima_data.reference_identity_min):
                pima_data.reference_mapping_ont_bam = os.path.join(pima_data.circos_dir, 'reference_mapping_ont.bam')
                minimap_and_sort(
                    pima_data,
                    pima_data.reference_fasta, 
                    pima_data.reference_mapping_ont_bam,
                    pima_data.ont_fastq, 
                    ont = True,
                )

                # Make an mpileup file out of the reads
                reference_mapping_mpileup = os.path.join(pima_data.circos_dir, 'reference_mapping_ont.mpileup')
                mpileup_bam(
                    pima_data,
                    pima_data.reference_fasta, 
                    pima_data.reference_mapping_ont_bam, 
                    reference_mapping_mpileup, 
                    pima_data.circos_dir,
                )
                ont_reference_coverage_fp = os.path.join(pima_data.circos_dir, 'reference_mapping_ont.mpileup')

    # Draw one circos plot for each of the contigs in the reference sequence
    with open(reference_sizes_fp, "r") as fin:
        for line in fin:
            contig, contig_size = line.rstrip().rsplit("\t")
            print_and_log(
                pima_data,
                f"Drawing Circos plot for {contig}",
                pima_data.sub_process_verbosity, 
                pima_data.sub_process_color,
            )

            contig_dir = os.path.join(pima_data.circos_dir, contig)
            os.makedirs(contig_dir)
            
            # Pull the aligned regions out of the dnadiff 1coords output
            if not (pima_data.self_circos or pima_data.genome_fasta is None):
                contig_alignment_fp = os.path.join(contig_dir, 'alignment.txt')
                command = " ".join(
                    [
                        'grep', contig, reference_1coords_fp,
                        '| awk \'{OFS = "\t";print $(NF - 1),$1,$2,$13}\'',
                        '| bedtools merge -d 25 -c 4 -o distinct -i -',
                        '1>', contig_alignment_fp, 
                        '2> /dev/null',
                    ]
                )
                print_and_run(pima_data, command)
        
            # Pull the coverage data from the mpileup file
            # Coverage files absent when no sequence data is provided (just a genome assembly)
            if (ont_reference_coverage_fp is None and illumina_reference_coverage_fp is None):
                contig_coverage_fp = None
                illumina_contig_coverage_fp = None

            # coverage files present with ONT fastq reads
            elif (ont_reference_coverage_fp is not None and illumina_reference_coverage_fp is None):
                contig_coverage_fp = os.path.join(contig_dir, 'coverage.mpileup')
                illumina_contig_coverage_fp = None
                command = " ".join(
                    [
                        'grep', contig, ont_reference_coverage_fp, 
                        '| awk \'{OFS = "\t"; print $1,$2,$4}\'',
                        '>', contig_coverage_fp,
                    ]
                )
                print_and_run(pima_data, command)

            # coverage files present for Illumina data only
            elif (ont_reference_coverage_fp is None and illumina_reference_coverage_fp is not None):
                contig_coverage_fp = None
                illumina_contig_coverage_fp = os.path.join(contig_dir, 'illumina_coverage.mpileup')
                command = " ".join(
                    [
                        'grep', contig, illumina_reference_coverage_fp, 
                        '| awk \'{OFS = "\t"; print $1,$2,$4}\'',
                        '>', illumina_contig_coverage_fp,
                    ]
                )
                print_and_run(pima_data, command)

            # coverage files present for ONT and Illumina data
            elif (ont_reference_coverage_fp is not None and illumina_reference_coverage_fp is not None):
                contig_coverage_fp = os.path.join(contig_dir, 'coverage.mpileup')
                command = " ".join(
                    [
                        'grep', contig, ont_reference_coverage_fp, 
                        '| awk \'{OFS = "\t"; print $1,$2,$4}\'',
                        '>', contig_coverage_fp, 
                    ]
                )
                print_and_run(pima_data, command)

                illumina_contig_coverage_fp = os.path.join(contig_dir, 'illumina_coverage.mpileup')
                command = " ".join(
                    [
                        'grep', contig, illumina_reference_coverage_fp, 
                        '| awk \'{OFS = "\t"; print $1,$2,$4}\'',
                        '>', illumina_contig_coverage_fp,
                    ]
                )
                print_and_run(pima_data, command)

            # Build circos plot with pycircos
            circos_elem = BuildCircosPlots(
                ge_name=contig, 
                ge_size=int(contig_size),
                aln_file=contig_alignment_fp,
                cov_file=contig_coverage_fp,
                illumina_cov_file=illumina_contig_coverage_fp,
                gene_file=virulence_genes_fp,
                outdir=contig_dir,
            )

            circos_fig = circos_elem.main()

            #draws a solid blue line for the reference backbone if you're using the denovo genome
            if pima_data.self_circos:
                circos_elem.ge_circle.barplot(
                    contig,
                    data = [1],
                    positions=[int(1)],
                    width=[int(contig_size)-1],
                    raxis_range=[600,700],
                    facecolor="#1f77b4",
                    linewidth=0,
                )
 
            #rasterize the coverage plots to shrink the svg
            # need to extract the matplotlib object that pycircos encapsulates
            # While doing so, lets try and squeeze as much whitespace as possible from the plots
            mod_fig = circos_fig.figure
            for ax in mod_fig.get_axes():
                ax.set_rasterization_zorder(2)
                ax.margins(0)
                ax.set_position([0,0,1,1])
                if hasattr(ax, 'collections') and ax.collections:
                    for collection in ax.collections:
                        collection.set_zorder(1)
                if hasattr(ax, 'patches') and ax.patches:
                    for patch in ax.patches:
                        patch.set_zorder(1)

            mod_fig.savefig(fname=os.path.join(contig_dir,'circos.svg'), format='svg', dpi=300, bbox_inches='tight', pad_inches=0)

            # Keep track of images for the report
            circos_svg = os.path.join(contig_dir, 'circos.svg')
            pima_data.contig_alignment[contig] = circos_svg
        pima_data.did_circos_plots = True
        make_finish_file(pima_data, pima_data.circos_dir)         