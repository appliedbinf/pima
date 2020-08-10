

import pandas
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

class PimaReport :
    
    def __init__(self, analysis) :

        self.analysis = analysis
        self.report = analysis.report
        self.doc = None

        
    def start_doc(self) :

        geometry_options = {"margin": "0.75in"}
        self.doc = Document(geometry_options=geometry_options)
        self.doc.preamble.append(Command('usepackage{float}'))

        
    def add_header(self) :

        header_text = 'Analysis of ' + self.report['name']
        
        header = PageStyle('header')
        with header.create(Head('L')) :
            header.append(header_text)
            header.append(LineBreak())

        with header.create(Foot('R')) :
            header.append(simple_page_number())

        self.doc.preamble.append(header)
        self.doc.change_document_style('header')


    def add_summary(self) :

        with self.doc.create(Section(self.analysis.summary_title, numbering = False)) :

            with self.doc.create(Subsection('Run information', numbering = False)) :
                with self.doc.create(Tabular('lp{6cm}lp{20cm}', width = 2)) as table :
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

            if len(self.report[self.analysis.assembly_title]) > 0 :
                if len(self.report[self.analysis.assembly_title][self.analysis.assembly_notes_title]) > 0 :
                    with self.doc.create(Subsection(self.analysis.assembly_notes_title, numbering = False)) :
                        left = FlushLeft()
                        for note in self.report[self.analysis.assembly_title][self.analysis.assembly_notes_title] :
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

                        
    def add_alignment(self) :

        if len(self.report[self.analysis.alignment_title][self.analysis.contig_alignment_title]) > 0:
            alignments = self.report[self.analysis.alignment_title][self.analysis.contig_alignment_title]
        else :
            return
            
        self.doc.append(NewPage())

        with self.doc.create(Section(self.analysis.alignment_title, numbering = False)) :

            with self.doc.create(Subsection(self.analysis.snp_indel_title, numbering = False)) :

                with self.doc.create(Tabular('ll', width = 2)) as table :
                    table.add_row(('SNPs', '{:,}'.format(self.analysis.quast_mismatches)))
                    table.add_row(('Small indels', '{:,}'.format(self.analysis.quast_indels)))
                self.doc.append(LineBreak())
                    
            if len(self.report[self.analysis.alignment_title][self.analysis.alignment_notes_title]) > 0 :

                with self.doc.create(Subsection(self.analysis.alignment_notes_title, numbering = False)) :

                    left = FlushLeft()
                    for note in self.report[self.analysis.alignment_title][self.analysis.alignment_notes_title] :
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


    def add_features(self) :

        if len(self.report[self.analysis.feature_title]) == 0:
            return

        self.doc.append(NewPage())
        
        with self.doc.create(Section(self.analysis.feature_title, numbering = False)) :
            for feature_name in self.report[self.analysis.feature_title].index.tolist() :

                features = self.report[self.analysis.feature_title][feature_name].copy()
                if features.shape[0] == 0 :
                    continue
                features.iloc[:, 1] = features.iloc[:, 1].apply(lambda x: '{:,}'.format(x))
                features.iloc[:, 2] = features.iloc[:, 2].apply(lambda x: '{:,}'.format(x))
                
                table_format = 'l' * (features.shape[1] - 1)

                with self.doc.create(Subsection(feature_name, numbering = False)) :
                    if (features.shape[0] == 0) :
                        self.doc.append('None')
                        continue

                    for contig in pandas.unique(features.iloc[:, 0]) :

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
                                        
                                        
    def add_feature_plots(self) :
        
        if len(self.report[self.analysis.feature_plot_title]) == 0 :
            return

        self.doc.append(NewPage())

        self.doc.append(Command('graphicspath{{../drawing/}}'))
        
        with self.doc.create(Section('Feature plots', numbering = False)) :

            self.doc.append('Only contigs with features are shown')

            for contig in self.report[self.analysis.feature_plot_title].index.tolist() :

                image_png = os.path.basename(self.report[self.analysis.feature_plot_title][contig])
                
                with self.doc.create(Figure(position = 'h!')) as figure :
                    figure.add_image(image_png, width = '7in')

                                
    def add_mutations(self) :

        # Make sure we looked for mutations
        if len(self.report[self.analysis.mutation_title]) ==  0 :
            return
        
        mutations = self.report[self.analysis.mutation_title]
        
        self.doc.append(NewPage())

        table_format = 'l' * 4
        
        with self.doc.create(Section(self.analysis.mutation_title, numbering = False)) :

            for region_name in mutations.index.tolist() :

                region_mutations = mutations[region_name].copy()

                with self.doc.create(Subsection(region_name, numbering = False)) :
                    if (region_mutations.shape[0] == 0) :
                        self.doc.append('None')
                        continue

                    region_mutations.iloc[:, 1] = region_mutations.iloc[:, 1].apply(lambda x: '{:,}'.format(x))
                    
                    with self.doc.create(Tabular(table_format)) as table :
                        table.add_row(('Reference contig', 'Position', 'Reference', 'Alternate'))
                        table.add_hline()
                        for i in range(region_mutations.shape[0]) :
                            table.add_row(region_mutations.iloc[i, [0,1,3,4]].values.tolist())

                            
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
        if len(self.report[self.analysis.large_indel_title]) == 0 :
            return
        
        large_indels = self.report[self.analysis.large_indel_title]

        if large_indels is None :
            return
        
        self.doc.append(NewPage())

        with self.doc.create(Section(self.analysis.large_indel_title, numbering = False)) :
        
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

                            
    def add_plasmids(self) :

        if len(self.report[self.analysis.plasmid_title]) == 0 :
            return

        # Make sure we looked for mutations
        plasmids = self.report[self.analysis.plasmid_title].copy()

        if plasmids is None :
            return
        
        self.doc.append(NewPage())

        with self.doc.create(Section(self.analysis.plasmid_title, numbering = False)) :
        
            if (plasmids.shape[0] == 0) :
                self.doc.append('None')
                return

            plasmids.iloc[:, 3] = plasmids.iloc[:, 3].apply(lambda x: '{:,}'.format(x))
            plasmids.iloc[:, 4] = plasmids.iloc[:, 4].apply(lambda x: '{:,}'.format(x))
            plasmids.iloc[:, 5] = plasmids.iloc[:, 5].apply(lambda x: '{:,}'.format(x))
            
            table_format = 'p{0.15\linewidth}p{0.3\linewidth}p{0.1\linewidth}p{0.08\linewidth}p{0.08\linewidth}p{0.08\linewidth}'

            with self.doc.create(Tabular(table_format)) as table :
                table.add_row(('Genome contig', 'Plasmid hit', 'Plasmid acc.', 'Contig size', 'Aliged', 'Plasmid size'))
                table.add_hline()
                for i in range(plasmids.shape[0]) :
                    table.add_row(plasmids.iloc[i, 0:6].values.tolist())

                    
    def add_methods(self) :

        self.doc.append(NewPage())

        methods = self.report[self.analysis.methods_title]

        with self.doc.create(Section(self.analysis.methods_title, numbering = False)) :

            for methods_section in methods.index.tolist() :

                if len(methods[methods_section]) == 0 :
                    continue
                
                with self.doc.create(Subsection(methods_section, numbering = False)) :

                    self.doc.append('  '.join(methods[methods_section]))


    def make_tex(self) :

        self.doc.generate_tex(self.analysis.report_prefix)
                    
                    
    def make_report(self) :

        self.start_doc()
        self.add_header()
        self.add_summary()
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
        

            


        


        
