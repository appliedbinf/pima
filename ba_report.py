import pandas as pd
from si_prefix import si_format

from pylatex import Document, Section, \
    Command, \
    PageStyle, Head, Foot, MiniPage, VerticalSpace, \
    simple_page_number, NewPage, \
    LargeText, MediumText, LineBreak, NewLine, TextColor, \
    Subsection, Tabular, Math, TikZ, Axis, \
    Plot, Figure, Matrix, Alignat

from pylatex.position import FlushLeft

from pylatex.utils import italic, bold

import os

cdc_advisory = 'The analysis and report presented here should be treated as preliminary.  Please contact the CDC/BDRD with any results regarding Bacillus anthracis.'

class PimaReport :
    
    def __init__(self, analysis) :

        self.analysis = analysis
        self.report = analysis.report
        self.doc = None

        self.methods = pd.Series()
        self.summary_title = 'Summary'
        self.methods_title = 'Methods'
        self.basecalling_methods_title = 'Basecalling'
        self.contamination_methods_title = 'Contamination check'
        self.methods[self.contamination_methods_title] = pd.Series()
        self.assembly_methods_title = 'Assembly'
        self.methods[self.assembly_methods_title] = pd.Series()
        self.reference_methods_title = 'Reference comparison'
        self.methods[self.reference_methods_title] = pd.Series()
        self.mutation_methods_title = 'Mutation screening'
        self.methods[self.mutation_methods_title] = pd.Series()
        self.feature_methods_title = 'Feature annotation'
        self.methods[self.feature_methods_title] = pd.Series()
        self.plasmid_methods_title = 'Plasmid annotation'
        self.methods[self.plasmid_methods_title] = pd.Series()

        
        self.basecalling_title = 'Basecalling'
        self.assembly_notes_title = 'Assembly notes'
        self.alignment_title = 'Comparison with reference'
        self.alignment_notes_title = 'Alignment notes'
        self.contig_alignment_title = 'Alignment vs. reference contigs'
        self.large_indel_title = 'Large insertions & deletions'
        self.snp_indel_title = 'SNPs and small indels'
        self.feature_title = 'Features found in the assembly'
        self.feature_plot_title = 'Feature annotation plots'

        
    def start_doc(self) :

        geometry_options = {"margin": "0.75in"}
        self.doc = Document(geometry_options=geometry_options)
        self.doc.preamble.append(Command('usepackage{float}'))

        
    def add_header(self) :

        header_text = 'Analysis of ' + self.analysis.analysis_name
        
        header = PageStyle('header')
        with header.create(Head('L')) :
            header.append(header_text)
            header.append(LineBreak())

        with header.create(Foot('R')) :
            header.append(simple_page_number())

        self.doc.preamble.append(header)
        self.doc.change_document_style('header')


    def add_summary(self) :

        with self.doc.create(Section(self.summary_title, numbering = False)) :

            with self.doc.create(Subsection('CDC Advisory', numbering = False)) as subsection:
                self.doc.append(cdc_advisory)

            with self.doc.create(Subsection('Run information', numbering = False)) :
                with self.doc.create(Tabular('p{0.15\linewidth}p{0.65\linewidth}', width = 2)) as table :
                    table.add_row(('Date', self.analysis.start_time))
                    if self.analysis.ont_fast5 :
                        table.add_row(('ONT FAST5', self.analysis.ont_fast5))
                    if self.analysis.ont_raw_fastq :
                        table.add_row(('ONT FASTQ', self.analysis.ont_raw_fastq))
                    if self.analysis.illumina_fastq :
                        table.add_row(('Illumina FASTQ', ', '.join(self.analysis.illumina_fastq)))
                    if self.analysis.genome_fasta :
                        table.add_row(('Assembly', self.analysis.genome_fasta))
                    if self.analysis.reference_fasta :
                        table.add_row(('Reference', self.analysis.reference_fasta))
                        
                self.doc.append(VerticalSpace("10pt"))

            if not (self.analysis.ont_n50 is None ):
                with self.doc.create(Subsection('ONT library statistics', numbering = False)) :
                    with self.doc.create(Tabular('ll', width = 2)) as table :
                        table.add_row(('ONT N50', '{:,}'.format(self.analysis.ont_n50)))
                        table.add_row(('ONT reads', '{:,}'.format(self.analysis.ont_read_count)))
                        table.add_row(('ONT bases', '{:s}'.format(self.analysis.ont_bases)))
                    self.doc.append(VerticalSpace("10pt"))

            methods = []
            if self.analysis.did_guppy_ont_fast5 :
                methods += ['ONT reads were basecalled using guppy (v ' + self.analysis.versions['guppy'] + ').']
            if self.analysis.did_qcat_ont_fastq :
                methods += ['ONT reads were demultiplexed and trimmed using qcat (v ' + self.analysis.versions['qcat'] + ').']
            self.methods[self.basecalling_methods_title] = pd.Series(methods)

            if not (self.analysis.illumina_length_mean is None ):
                with self.doc.create(Subsection('Illumina library statistics', numbering = False)) :
                    with self.doc.create(Tabular('ll', width = 2)) as table :
                        table.add_row(('Illumina mean length', '{:.1f}'.format(self.analysis.illumina_length_mean)))
                        table.add_row(('Illumina reads', '{:,}'.format(self.analysis.illumina_read_count)))
                        table.add_row(('Illumina bases', '{:s}'.format(self.analysis.illumina_bases)))
                    self.doc.append(VerticalSpace("10pt"))
                    
            if not (self.analysis.genome is None) :
                with self.doc.create(Subsection('Assembly statistics', numbering = False)) :
                    with self.doc.create(Tabular('ll', width = 2)) as table :
                        table.add_row(('Contigs', len(self.analysis.genome)))
                        
                        genome_size = 0
                        for i in self.analysis.genome :
                            genome_size += len(i.seq)
                        genome_size = si_format(genome_size, precision = 1)
                        table.add_row(('Assembly size', genome_size))

                    self.doc.append(VerticalSpace("10pt"))

                if self.analysis.did_flye_ont_fastq :
                    method = 'ONT reads were assembled using Flye (v ' + self.analysis.versions['flye'] + ').'
                    self.methods[self.assembly_methods_title] = self.methods[self.assembly_methods_title].append(pd.Series(method))
                if self.analysis.did_medaka_ont_assembly :
                    method = 'The genome assembly was polished using ONT reads and Medaka (v ' +\
                        self.analysis.versions['medaka'] + ').'
                    self.methods[self.assembly_methods_title] = self.methods[self.assembly_methods_title].append(pd.Series(method))



            if len(self.analysis.assembly_notes) > 0 :
                with self.doc.create(Subsection(self.assembly_notes_title, numbering = False)) :
                        left = FlushLeft()
                        for note in self.analysis.assembly_notes :
                            left.append(note)
                            left.append(LineBreak())
                            self.doc.append(left)
                        self.doc.append(VerticalSpace("10pt"))
                
            if not (self.analysis.contig_info is None) :

                for method in ['ONT', 'Illumina'] :

                    if not method in self.analysis.contig_info.index :
                            continue

                    with self.doc.create(Subsection('Assembly coverage by ' + method, numbering = False)) :
                        
                        table_format = 'l' * self.analysis.contig_info[method].shape[1]

                        self.doc.append('')
                        
                        with self.doc.create(Tabular(table_format)) as table :
                            table.add_row(('Contig', 'Length (bp)', 'Coverage (X)'))
                            table.add_hline()
                            formatted = self.analysis.contig_info[method].copy()
                            formatted.iloc[:,1] = formatted.iloc[:,1].apply(lambda x: '{:,}'.format(x))
                            for i in range(self.analysis.contig_info[method].shape[0]) :
                                table.add_row(formatted.iloc[i, :].values.tolist())

                        self.doc.append(LineBreak())
                        self.doc.append(VerticalSpace("10pt"))


    def add_contamination(self) :

        if self.analysis.kraken_fracs is None :
            return

        self.doc.append(NewPage())

        with self.doc.create(Section('Contamination check', numbering = False)) :

            for read_type, kraken_fracs in self.analysis.kraken_fracs.iteritems() :

                left = FlushLeft()
                left.append(read_type + ' classifications')
                left.append(VerticalSpace('5pt'))
                self.doc.append(left)
                            
                with self.doc.create(Tabular(''.join(['l'] * kraken_fracs.shape[1]), width = kraken_fracs.shape[1])) as table :
                    table.add_row(('Percent of reads', 'Level', 'Label'))
                    table.add_hline()
                    for index, row in kraken_fracs.iterrows() :
                        table.add_row(row.tolist())
                        
                self.doc.append(LineBreak())
                self.doc.append(VerticalSpace('5pt'))

                if not self.contamination_methods_title in self.methods :
                    self.methods[self.contamination_methods_title] = ''

        method = 'Kraken2 (' + self.analysis.versions['kraken2'] + ') was used to assign the raw reads into taxa.'
        self.methods[self.contamination_methods_title] = self.methods[self.contamination_methods_title].append(pd.Series(method))

                
    def add_alignment(self) :

        if len(self.analysis.contig_alignment) > 0:
            alignments = self.analysis.contig_alignment
        else :
            return
            
        self.doc.append(NewPage())

        with self.doc.create(Section(self.alignment_title, numbering = False)) :

            with self.doc.create(Subsection(self.snp_indel_title, numbering = False)) :

                with self.doc.create(Tabular('ll', width = 2)) as table :
                    table.add_row(('SNPs', '{:,}'.format(self.analysis.quast_mismatches)))
                    table.add_row(('Small indels', '{:,}'.format(self.analysis.quast_indels)))
                self.doc.append(LineBreak())
                    
            if len(self.analysis.alignment_notes) > 0 :

                with self.doc.create(Subsection(self.alignment_notes_title, numbering = False)) :

                    left = FlushLeft()
                    for note in self.analysis.alignment_notes :
                        left.append(note)
                        left.append(LineBreak())
                    self.doc.append(left)
                        
            for contig in alignments.index.tolist() :
                                        
                contig_title = 'Alignment to ' + contig
                self.doc.append(Command('graphicspath{{../circos/' + contig + '/}}'))
                image_png = os.path.basename(alignments[contig])
                
                with self.doc.create(Subsection(contig_title, numbering = False)) :
                    with self.doc.create(Figure(position = 'H')) as figure :
                        figure.add_image(image_png, width = '5in')

                    
        method = 'The genome assembly was aligned against the reference sequencing using dnadiff (v ' \
            + self.analysis.versions['dnadiff'] + ').'
        self.methods[self.reference_methods_title] = self.methods[self.reference_methods_title].append(pd.Series(method))


    def add_features(self) :

        if len(self.analysis.feature_hits) == 0:
            return

        self.doc.append(NewPage())
        
        with self.doc.create(Section(self.feature_title, numbering = False)) :
            for feature_name in self.analysis.feature_hits.index.tolist() :

                features = self.analysis.feature_hits[feature_name].copy()
                if features.shape[0] == 0 :
                    continue
                    
                features.iloc[:, 1] = features.iloc[:, 1].apply(lambda x: '{:,}'.format(x))
                features.iloc[:, 2] = features.iloc[:, 2].apply(lambda x: '{:,}'.format(x))
                
                table_format = 'l' * (features.shape[1] - 1)

                with self.doc.create(Subsection(feature_name, numbering = False)) :
                    if (features.shape[0] == 0) :
                        self.doc.append('None')
                        continue

                    for contig in pd.unique(features.iloc[:, 0]) :

                        self.doc.append(contig)

                        contig_features = features.loc[(features.iloc[:, 0] == contig), :]
                        
                        with self.doc.create(Tabular(table_format)) as table :
                            table.add_row(('Start', 'Stop', 'Feature', 'Identity (%)', 'Strand'))
                            table.add_hline()
                            for i in range(contig_features.shape[0]) :
                                feature = contig_features.iloc[i, :].copy(deep = True)
                                feature[4] = '{:.3f}'.format(feature[4])
                                table.add_row(feature[1:].values.tolist())

                        self.doc.append(LineBreak())
                        self.doc.append(VerticalSpace("10pt"))

        method = 'The genome assembly was queried for features using blastn (v ' + self.analysis.versions['blastn'] + ').  ' + \
            'Feature hits were clustered using bedtools (v ' + self.analysis.versions['bedtools'] + ') ' + \
            'and the highest scoring hit for each cluster was reported.'
        self.methods[self.feature_methods_title] = self.methods[self.feature_methods_title].append(pd.Series(method))

                                        
    def add_feature_plots(self) :
        
        if len(self.analysis.feature_plots) == 0 :
            return

        self.doc.append(NewPage())

        self.doc.append(Command('graphicspath{{../drawing/}}'))
        
        with self.doc.create(Section('Feature plots', numbering = False)) :

            self.doc.append('Only contigs with features are shown')

            for contig in self.analysis.feature_plots.index.tolist() :

                image_png = os.path.basename(self.analysis.feature_plots[contig])
                
                with self.doc.create(Figure(position = 'h!')) as figure :
                    figure.add_image(image_png, width = '7in')

                                
    def add_mutations(self) :

        # Make sure we looked for mutations
        if len(self.report[self.analysis.mutation_title]) ==  0 :
            return
        
        mutations = self.report[self.analysis.mutation_title]
        
        self.doc.append(NewPage())

        table_format = 'p{0.1\linewidth}p{0.1\linewidth}p{0.12\linewidth}p{0.12\linewidth}p{0.12\linewidth}p{0.34\linewidth}'
        
        with self.doc.create(Section(self.analysis.mutation_title, numbering = False)) :

            for region_name in mutations.index.tolist() :

                region_mutations = mutations[region_name].copy()

                with self.doc.create(Subsection(region_name, numbering = False)) :
                    if (region_mutations.shape[0] == 0) :
                        self.doc.append('None')
                        continue

                    region_mutations.iloc[:, 1] = region_mutations.iloc[:, 1].apply(lambda x: '{:,}'.format(x))
                    
                    with self.doc.create(Tabular(table_format)) as table :
                        table.add_row(('Reference contig', 'Position', 'Reference', 'Alternate', 'Drug', 'Note'))
                        table.add_hline()
                        for i in range(region_mutations.shape[0]) :
                            table.add_row(region_mutations.iloc[i, [0,1,3,4,5,6]].values.tolist())

        method = self.analysis.mutations_read_type + ' reads were mapped to the reference sequence using minimap2 (v '\
            + self.analysis.versions['minimap2'] + ').'
        self.methods[self.mutation_methods_title] = self.methods[self.mutation_methods_title].append(pd.Series(method))

        method = ' '.join(['Mutations were identified using'
                           'samtools mpileup (v', self.analysis.versions['samtools'],  ')',
                           'and varscan (v', self.analysis.versions['varscan'], ').'])
        self.methods[self.mutation_methods_title] = self.methods[self.mutation_methods_title].append(pd.Series(method))

                            
    def add_amr_matrix(self) :

        # Make sure that we have an AMR matrix to plot
        amr_matrix = self.report[self.analysis.amr_matrix_title]

        if len(amr_matrix) == 0 :
            return
        
        if amr_matrix is None :
            return

        self.doc.append(NewPage())

        with self.doc.create(Section(self.analysis.amr_matrix_title, numbering = False)) :

            self.doc.append('AMR genes and mutations with their corresponding drugs.')

            self.doc.append(Command('graphicspath{{../}}'))
            
            with self.doc.create(Figure(position = 'h!')) as figure :
                figure.add_image(amr_matrix['png'], width = '7in')
                            

    def add_large_indels(self) :

        # Make sure we looked for mutations
        if len(self.analysis.large_indels) == 0 :
            return
        
        large_indels = self.analysis.large_indels

        self.doc.append(NewPage())

        with self.doc.create(Section(self.large_indel_title, numbering = False)) :
        
            for genome in ['Reference insertions', 'Query insertions'] :

                genome_indels = large_indels[genome].copy()
                
                with self.doc.create(Subsection(genome, numbering = False)) :
                    if (genome_indels.shape[0] == 0) :
                        self.doc.append('None')
                        continue

                    genome_indels.iloc[:, 1] = genome_indels.iloc[:, 1].apply(lambda x: '{:,}'.format(x))
                    genome_indels.iloc[:, 2] = genome_indels.iloc[:, 2].apply(lambda x: '{:,}'.format(x))
                    genome_indels.iloc[:, 3] = genome_indels.iloc[:, 3].apply(lambda x: '{:,}'.format(x))

                    table_format = 'l' * genome_indels.shape[1]
                    
                    with self.doc.create(Tabular(table_format)) as table :
                        table.add_row(('Reference contig', 'Start', 'Stop', 'Size (bp)'))
                        table.add_hline()
                        for i in range(genome_indels.shape[0]) :
                            table.add_row(genome_indels.iloc[i,:].values.tolist())

        method = 'Large insertions or deletions were found as the complement of aligned ' + \
            'regions using bedtools (v ' + self.analysis.versions['bedtools'] + ').'
        self.methods[self.reference_methods_title] = self.methods[self.reference_methods_title].append(pd.Series(method))

                            
    def add_plasmids(self) :

        if not self.analysis.did_call_plasmids :
            return
        
        # Make sure we looked for mutations
        plasmids = self.analysis.plasmids

        if plasmids is None :
            return

        plasmids = plasmids.copy()
        
        self.doc.append(NewPage())

        with self.doc.create(Section(self.analysis.plasmid_title, numbering = False)) :
        
            if (plasmids.shape[0] == 0) :
                self.doc.append('None')
                return

            plasmids.iloc[:, 3] = plasmids.iloc[:, 3].apply(lambda x: '{:,}'.format(x))
            plasmids.iloc[:, 4] = plasmids.iloc[:, 4].apply(lambda x: '{:,}'.format(x))
            plasmids.iloc[:, 5] = plasmids.iloc[:, 5].apply(lambda x: '{:,}'.format(x))
            
            table_format = 'p{0.1\linewidth}p{0.3\linewidth}p{0.15\linewidth}p{0.08\linewidth}p{0.08\linewidth}p{0.08\linewidth}'

            with self.doc.create(Tabular(table_format)) as table :
                table.add_row(('Genome contig', 'Plasmid hit', 'Plasmid acc.', 'Contig size', 'Aliged', 'Plasmid size'))
                table.add_hline()
                for i in range(plasmids.shape[0]) :
                    table.add_row(plasmids.iloc[i, 0:6].values.tolist())
                    
        method = ' '.join(['The plasmid reference database was queried against the genome assembly using minimap2 (v',
                           self.analysis.versions['minimap2'], ').'])
        self.methods[self.plasmid_methods_title] = self.methods[self.plasmid_methods_title].append(pd.Series(method))
        
        method = 'The resulting SAM was converted to a PSL using a custom version of sam2psl.'
        self.methods[self.plasmid_methods_title] = self.methods[self.plasmid_methods_title].append(pd.Series(method))

        method = 'Plasmid-to-genome hits were resolved using the pChunks algorithm.'
        self.methods[self.plasmid_methods_title]  = self.methods[self.plasmids_methods_title].append(pd.Series(method))


    def add_methods(self) :

        if len(self.methods) == 0 :
            return
        
        self.doc.append(NewPage())

        with self.doc.create(Section(self.methods_title, numbering = False)) :

            for methods_section in self.methods.index.tolist() :

                if len(self.methods[methods_section]) == 0 :
                    continue
                
                with self.doc.create(Subsection(methods_section, numbering = False)) :
                    self.doc.append(' '.join(self.methods[methods_section]))


    def make_tex(self) :

        self.doc.generate_tex(self.analysis.report_prefix)
                    
                    
    def make_report(self) :

        self.start_doc()
        self.add_header()
        self.add_summary()
        self.add_contamination()
        self.add_alignment()
        self.add_features()
        self.add_feature_plots()
        self.add_mutations()
        self.add_large_indels()
        self.add_plasmids()
        self.add_amr_matrix()
        #self.add_snps()
        self.add_methods()
        self.make_tex()
        

            


        


        
