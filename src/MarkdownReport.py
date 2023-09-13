import pandas as pd
from mdutils.mdutils import MdUtils
from mdutils.tools import TableOfContents
from si_prefix import si_format

import os

cdc_advisory = 'The analysis and report presented here should be treated as preliminary.  Please contact the CDC/BDRD with any results regarding _Bacillus anthracis_.'

class PimaReport:

    def __init__(self, analysis):

        self.analysis = analysis
        self.doc = None

        self.summary_title = 'Summary'
        self.basecalling_title = 'Basecalling'
        self.assembly_notes_title = 'Assembly notes'
        self.alignment_title = 'Comparison with reference'
        self.alignment_notes_title = 'Alignment notes'
        self.contig_alignment_title = 'Alignment vs. reference contigs'
        self.large_indel_title = 'Large insertions & deletions'
        self.snp_indel_title = 'SNPs and small indels'
        self.feature_title = 'Features found in the assembly'
        self.feature_plot_title = 'Feature annotation plots'
        self.mutation_title = 'Mutations found in the sample'
        self.amr_matrix_title = 'AMR matrix'

        self.methods = pd.Series(dtype='float64')
        self.methods_title = 'Methods'
        self.basecalling_methods_title = 'Basecalling'
        self.contamination_methods_title = 'Contamination check'
        self.methods[self.contamination_methods_title] = pd.Series(dtype='float64')
        self.assembly_methods_title = 'Assembly'
        self.methods[self.assembly_methods_title] = pd.Series(dtype='float64')
        self.reference_methods_title = 'Reference comparison'
        self.methods[self.reference_methods_title] = pd.Series(dtype='float64')
        self.mutation_methods_title = 'Mutation screening'
        self.methods[self.mutation_methods_title] = pd.Series(dtype='float64')
        self.feature_methods_title = 'Feature annotation'
        self.methods[self.feature_methods_title] = pd.Series(dtype='float64')
        self.plasmid_methods_title = 'Plasmid annotation'
        self.methods[self.plasmid_methods_title] = pd.Series(dtype='float64')
        self.pilon_methods_title = 'Pilon assembly'
        self.methods[self.pilon_methods_title] = pd.Series(dtype='float64')

    def start_doc(self): #Converted
        header_text = 'Analysis of ' + self.analysis.analysis_name
        self.doc = MdUtils(file_name=self.analysis.report_md,title=header_text)

    def add_tableOfContents(self):
        self.doc.create_marker(text_marker="TableOfContents")
        self.doc.new_line()
        self.doc.new_line('<div style="page-break-after: always;"></div>')
        self.doc.new_line()

    def add_run_information(self): #Converted
        self.doc.new_line()
        self.doc.new_header(1,'Run Information')
        # Tables in md.utils are implemented as a wrapping function. Weird but ok.
        Table_list = [
            "Category",
            "Information",
            "Date",
            self.analysis.start_time,
            "ONT FAST5",
            self.wordwrap_markdown(self.analysis.ont_fast5),
            "ONT FASTQ",
            self.wordwrap_markdown(self.analysis.ont_raw_fastq),
            "Illumina FASTQ",
            self.wordwrap_markdown(self.analysis.illumina_fastq),
            "Assembly",
            self.wordwrap_markdown(self.analysis.genome_fasta),
            "Reference",
            self.wordwrap_markdown(self.analysis.reference_fasta),
        ]
        self.doc.new_table(columns=2,rows=7,text=Table_list,text_align='left')
        self.doc.new_line()
        self.add_tableOfContents()
        self.doc.new_line()

    def add_ont_library_information(self): #Converted

        if self.analysis.ont_n50 is None:
            return
        self.doc.new_line()
        self.doc.new_header(2, 'ONT library statistics')
        Table_List = [
            "Category",
            "Quantity",
            "ONT N50",
            '{:,}'.format(self.analysis.ont_n50),
            "ONT reads",
            '{:,}'.format(self.analysis.ont_read_count),
            "ONT bases",
            '{:s}'.format(self.analysis.ont_bases),
            "Illumina FASTQ",
            self.wordwrap_markdown(self.analysis.illumina_fastq),
            "Assembly",
            self.wordwrap_markdown(self.analysis.genome_fasta),
            "Reference",
            self.wordwrap_markdown(self.analysis.reference_fasta),
        ]
        self.doc.new_table(columns=2, rows=7, text=Table_List, text_align='left')
        self.doc.new_line()

    def add_illumina_library_information(self): #Converted
        if self.analysis.illumina_length_mean is None:
            return

        self.doc.new_line()
        self.doc.new_header(2, 'Illumina library statistics')
        Table_List = [
            "Illumina Info.",
            "Quantity",
            'Illumina mean length',
            '{:.1f}'.format(self.analysis.illumina_length_mean),
            'Illumina reads',
            '{:,}'.format(self.analysis.illumina_read_count),
            'Illumina bases',
            '{:s}'.format(self.analysis.illumina_bases)
        ]
        self.doc.new_table(columns=2, rows=4, text=Table_List, text_align='left')

    def add_assembly_information(self): #Converted
        if self.analysis.genome is None:
            return

        self.doc.new_line()
        self.doc.new_header(2, 'Assembly statistics')

        genome_size = si_format(sum([len(x) for x in self.analysis.genome]), # This sums up the lengths of teh genomes
                                precision=1)
        Table_List = [
            "Category",
            "Information",
            "Contigs",
            str(len(self.analysis.genome)),
            "Assembly size",
            genome_size
        ]
        self.doc.new_table(columns=2, rows=3, text=Table_List, text_align='left')

    def wordwrap_markdown(self,string):
        if string:
            if len(string) < 35:
                return(string)
            else:
                if '/' in string:
                    adjust = string.split('/')
                    out = ''
                    max = 35
                    for i in adjust:
                        out = out + '/' + i
                        if len(out) > max:
                            out += '<br>'
                            max += 35
                    return(out)
                else:
                    out = [string[i:i + 35] for i in range(0, len(string), 50)]
                    return('<br>'.join(out))
        else:
            return(string)

    def add_contig_info(self): #Converted

        if self.analysis.contig_info is None:
            return

        for method in ['ONT', 'Illumina']:

            if not method in self.analysis.contig_info.index:
                continue
            self.doc.new_line()
            self.doc.new_header(2, 'Assembly coverage by ' + method)
            Table_List = [
                "Contig",
                "Length (bp)",
                "Coverage (X)",
            ]
            formatted = self.analysis.contig_info[method].copy() #No idea what this does
            formatted.iloc[:, 1] = formatted.iloc[:, 1].apply(lambda x: '{:,}'.format(x))# Weird way to format a series but ok
            for i in range(self.analysis.contig_info[method].shape[0]):
                Table_List = Table_List + formatted.iloc[i, :].values.tolist()
            row_count = int(len(Table_List)/3)
            self.doc.new_table(columns=3, rows=row_count, text=Table_List, text_align='left')

    def add_assembly_notes(self): # Converted

        if len(self.analysis.assembly_notes) == 0:
            return

        self.doc.new_line()
        self.doc.new_line('<div style="page-break-after: always;"></div>')
        self.doc.new_line()
        self.doc.new_header(2, self.assembly_notes_title)

        for note in self.analysis.assembly_notes:
            self.doc.new_line(note)

    def add_contamination(self): #Converted

        if len(self.analysis.kraken_fracs) == 0:
            return
        self.doc.new_line()
        self.doc.new_header(2,'Contamination check')
        for read_type, kraken_fracs in self.analysis.kraken_fracs.iteritems():

            self.doc.new_line(read_type + ' classifications')
            self.doc.new_line()
            Table_List = [
                "Percent of Reads",
                "Reads",
                "Level",
                "Label"
            ]

            for index, row in kraken_fracs.iterrows():
                Table_List = Table_List + row.tolist()

            row_count = int(len(Table_List)/4)

            self.doc.new_table(columns=4, rows=row_count, text=Table_List, text_align='left')

            if not self.contamination_methods_title in self.methods:
                self.methods[self.contamination_methods_title] = ''

        method = 'Kraken2 (' + self.analysis.versions['kraken2'] + ') was used to assign the raw reads into taxa.'
        self.methods[self.contamination_methods_title] = pd.concat([self.methods[self.contamination_methods_title],
            pd.Series(method, dtype="object")])
        """ self.methods[self.contamination_methods_title] = self.methods[self.contamination_methods_title].append(
            pd.Series(method)) """

    def add_alignment(self): #Converted

        if len(self.analysis.contig_alignment) > 0:
            alignments = self.analysis.contig_alignment
        else:
            return
        self.doc.new_line()
        self.doc.new_header(level=2, title=self.alignment_title)
        self.doc.new_line()
        self.doc.new_header(level=3, title=self.snp_indel_title)

        Table_1 = [
            "Category",
            "Quantity",
            'SNPs',
            '{:,}'.format(self.analysis.quast_mismatches),
            'Small indels',
            '{:,}'.format(self.analysis.quast_indels)
        ]

        self.doc.new_table(columns=2, rows=3, text=Table_1, text_align='left')

        self.doc.new_line('<div style="page-break-after: always;"></div>')
        self.doc.new_line()

        if len(self.analysis.alignment_notes) > 0:
            self.doc.new_header(level=3, title=self.alignment_notes_title)
            for note in self.analysis.alignment_notes:
                self.doc.new_line(note)

        for contig in alignments.index.tolist():
            contig_title = 'Alignment to ' + contig
            image_png = alignments[contig]
            self.doc.new_line()
            self.doc.new_header(level=3,title=contig_title)
            self.doc.new_line()
            #self.doc.new_header(level=4,title="contig coverage:{0}".format(str(self.analysis.contig_info[0]["coverage"][0])))
            self.doc.write(
                self.doc.new_inline_image(
                    text='contig_title',
                    path=os.path.abspath(image_png)
                )
                ,wrap_width=0
            )
            self.doc.new_line('<div style="page-break-after: always;"></div>')
            self.doc.new_line()

        method = 'The genome assembly was aligned against the reference sequencing using dnadiff (v ' \
                 + self.analysis.versions['dnadiff'] + ').'
        
        self.methods[self.reference_methods_title] = pd.concat([self.methods[self.reference_methods_title], 
            pd.Series(method, dtype='object')])
        """self.methods[self.reference_methods_title] = self.methods[self.reference_methods_title].append(
            pd.Series(method))"""

# Only a few more to go

    def add_features(self): #Converted

        if len(self.analysis.feature_hits) == 0:
            return

        self.doc.new_line()
        self.doc.new_header(level=2,title=self.feature_title)

        for feature_name in self.analysis.feature_hits.index.tolist():

            features = self.analysis.feature_hits[feature_name].copy()
            if features.shape[0] == 0:
                continue

            features.iloc[:, 1] = features.iloc[:, 1].apply(lambda x: '{:,}'.format(x))
            features.iloc[:, 2] = features.iloc[:, 2].apply(lambda x: '{:,}'.format(x))

            self.doc.new_line()
            self.doc.new_header(level=3,title=feature_name)
            if (features.shape[0] == 0):
                continue

            for contig in pd.unique(features.iloc[:, 0]):

                self.doc.new_line(contig)

                contig_features = features.loc[(features.iloc[:, 0] == contig), :]

                Table_List = [
                    'Start', 'Stop', 'Feature', 'Identity (%)', 'Strand'
                ]

                for i in range(contig_features.shape[0]):
                    feature = contig_features.iloc[i, :].copy(deep=True)
                    feature[4] = '{:.3f}'.format(feature[4])
                    Table_List = Table_List + feature[1:].values.tolist()

                row_count = int(len(Table_List) / 5)
                self.doc.new_line()
                self.doc.new_table(columns=5, rows=row_count, text=Table_List, text_align='left')

        method = 'The genome assembly was queried for features using blastn (v ' + self.analysis.versions[
            'blastn'] + ').  ' + \
                 'Feature hits were clustered using bedtools (v ' + self.analysis.versions['bedtools'] + ') ' + \
                 'and the highest scoring hit for each cluster was reported.'
        self.methods[self.feature_methods_title] = pd.concat([self.methods[self.feature_methods_title], pd.Series(method, dtype='object')])
        """self.methods[self.feature_methods_title] = self.methods[self.feature_methods_title].append(pd.Series(method))"""

    def add_feature_plots(self): #Converted

        if len(self.analysis.feature_plots) == 0:
            return

        self.doc.new_line()
        self.doc.new_header(level=2,title='Feature Plots')
        self.doc.new_paragraph('Only contigs with features are shown')

        for contig in self.analysis.feature_plots.index.tolist():
            image_png = self.analysis.feature_plots[contig]
            self.doc.write(
                self.doc.new_inline_image(
                    text='Analysis',
                    path=os.path.abspath(image_png),
                )
                ,wrap_width=0
            )

    def add_mutations(self):

        # Make sure we looked for mutations
        if not getattr(self.analysis, 'did_call_amr_mutations', False):
            return

        mutations = self.analysis.amr_mutations

        self.doc.new_line()
        self.doc.new_header(level=2,title=self.mutation_title)
        for region_name in mutations.index.tolist():

            region_mutations = mutations[region_name].copy()
            self.doc.new_line()

            self.doc.new_header(level=3,title=region_name)
            if (region_mutations.shape[0] == 0):
                self.doc.append('None')
                continue
            region_mutations.iloc[:, 1] = region_mutations.iloc[:, 1].apply(lambda x: '{:,}'.format(x))

            Table_List = [
                'Reference contig', 'Position', 'Reference', 'Alternate', 'Drug', 'Note'
            ]
            for i in range(region_mutations.shape[0]):
                Table_List = Table_List + region_mutations.iloc[i, [0, 1, 3, 4, 5, 6]].values.tolist()
            row_count = int(len(Table_List)/6)
            self.doc.new_table(columns=6, rows=row_count, text=Table_List, text_align='left')

        method = self.analysis.mutations_read_type + ' reads were mapped to the reference sequence using minimap2 (v ' \
                 + self.analysis.versions['minimap2'] + ').'
        self.methods[self.mutation_methods_title] = pd.concat([self.methods[self.mutation_methods_title], pd.Series(method,dtype="object")])
        """ self.methods[self.mutation_methods_title] = self.methods[self.mutation_methods_title].append(pd.Series(method)) """

        method = ' '.join(['Mutations were identified using'
                           'samtools mpileup (v', self.analysis.versions['samtools'], ')',
                           'and varscan (v', self.analysis.versions['varscan'], ').'])
        self.methods[self.mutation_methods_title] = pd.concat([self.methods[self.mutation_methods_title], pd.Series(method,dtype="object")])
        """ self.methods[self.mutation_methods_title] = self.methods[self.mutation_methods_title].append(pd.Series(method)) """

    def add_amr_matrix(self): #Converted

        # Make sure that we have an AMR matrix to plot
        if not getattr(self.analysis, 'did_draw_amr_matrix', False):
            return

        amr_matrix_png = self.analysis.amr_matrix_png
        self.doc.new_line()
        self.doc.new_header(level=2,title=self.amr_matrix_title)
        self.doc.new_line('AMR genes and mutations with their corresponding drugs.')
        self.doc.write(
            self.doc.new_inline_image(
                text='AMR genes and mutations with their corresponding drugs',
                path=amr_matrix_png
            )
            ,wrap_width=0
        )

    def add_large_indels(self): #Converted

        # Make sure we looked for mutations
        if len(self.analysis.large_indels) == 0:
            return

        large_indels = self.analysis.large_indels

        self.doc.new_line()
        self.doc.new_header(level=2,title=self.large_indel_title)

        for genome in ['Reference insertions', 'Query insertions']:

            genome_indels = large_indels[genome].copy()
            self.doc.new_line()
            self.doc.new_header(level=3,title=genome)

            if (genome_indels.shape[0] == 0):
                continue

            genome_indels.iloc[:, 1] = genome_indels.iloc[:, 1].apply(lambda x: '{:,}'.format(x))
            genome_indels.iloc[:, 2] = genome_indels.iloc[:, 2].apply(lambda x: '{:,}'.format(x))
            genome_indels.iloc[:, 3] = genome_indels.iloc[:, 3].apply(lambda x: '{:,}'.format(x))

            Table_List = [
                'Reference contig', 'Start', 'Stop', 'Size (bp)'
            ]
            for i in range(genome_indels.shape[0]):
                Table_List = Table_List + genome_indels.iloc[i, :].values.tolist()

            row_count= int(len(Table_List)/4)
            self.doc.new_table(columns=4, rows=row_count, text=Table_List, text_align='left')

        method = 'Large insertions or deletions were found as the complement of aligned ' + \
                 'regions using bedtools (v ' + self.analysis.versions['bedtools'] + ').'
        self.methods[self.reference_methods_title] = pd.concat([self.methods[self.reference_methods_title], 
            pd.Series(method,dtype='object')])
        """self.methods[self.reference_methods_title] = self.methods[self.reference_methods_title].append(
            pd.Series(method))"""

        self.doc.new_line()
        self.doc.new_line('<div style="page-break-after: always;"></div>')
        self.doc.new_line()

    def add_plasmids(self): #Converted

        if not getattr(self.analysis, 'did_call_plasmids', False):
            return

        # Make sure we looked for mutations
        plasmids = self.analysis.plasmids

        if plasmids is None:
            return

        plasmids = plasmids.copy()

        self.doc.new_line()
        self.doc.new_header(level=2,title=self.analysis.plasmid_title)

        if (plasmids.shape[0] == 0):
            self.doc.new_line('None')
            return

        plasmids.iloc[:, 3] = plasmids.iloc[:, 3].apply(lambda x: '{:,}'.format(x))
        plasmids.iloc[:, 4] = plasmids.iloc[:, 4].apply(lambda x: '{:,}'.format(x))
        plasmids.iloc[:, 5] = plasmids.iloc[:, 5].apply(lambda x: '{:,}'.format(x))

        Table_List = [
            'Genome contig',
            'Plasmid hit',
            'Plasmid acc.',
            'Contig size',
            'Aliged',
            'Plasmid size'
        ]

        for i in range(plasmids.shape[0]):
            Table_List = Table_List + plasmids.iloc[i, 0:6].values.tolist()

        row_count = int(len(Table_List) / 6)

        self.doc.new_table(columns=6, rows=row_count, text=Table_List, text_align='left')

        method = ' '.join(['The plasmid reference database was queried against the genome assembly using minimap2 (v',
                           self.analysis.versions['minimap2'], ').'])
        self.methods[self.plasmid_methods_title] = pd.concat([self.methods[self.plasmid_methods_title],pd.Series(method,dtype="object")])
        """ self.methods[self.plasmid_methods_title] = self.methods[self.plasmid_methods_title].append(pd.Series(method)) """

        method = 'The resulting SAM was converted to a PSL using a custom version of sam2psl.'
        self.methods[self.plasmid_methods_title] = pd.concat([self.methods[self.plasmid_methods_title],pd.Series(method,dtype="object")])
        """ self.methods[self.plasmid_methods_title] = self.methods[self.plasmid_methods_title].append(pd.Series(method)) """

        method = 'Plasmid-to-genome hits were resolved using the pChunks algorithm.'
        self.methods[self.plasmid_methods_title] = pd.concat([self.methods[self.plasmid_methods_title],pd.Series(method,dtype="object")])
        """ self.methods[self.plasmid_methods_title] = self.methods[self.plasmid_methods_title].append(pd.Series(method)) """

    def add_methods(self): #Converted

        self.doc.new_line('<div style="page-break-after: always;"></div>')
        self.doc.new_line()

        if len(self.methods) == 0:
            return

        self.doc.new_line()
        self.doc.new_header(level=2, title=self.methods_title)


        for methods_section in self.methods.index.tolist():
            if len(self.methods[methods_section]) == 0:
                continue
            self.doc.new_line()
            self.doc.new_header(level=3,title=methods_section)
            self.doc.new_paragraph(' '.join(self.methods[methods_section]))

    def add_summary(self): #Converted
        # Add summary title
        # self.doc.new_header(level=1, title=self.summary_title)
        # First section of Summary
        self.doc.new_header(level=1, title='CDC Advisory')
        self.doc.new_paragraph(cdc_advisory)
        self.doc.new_line()
        self.add_run_information()
        self.add_ont_library_information()
        methods = []
        if self.analysis.did_guppy_ont_fast5:
            methods += ['ONT reads were basecalled using guppy (v ' + self.analysis.versions['guppy'] + ').']
        if self.analysis.did_qcat_ont_fastq:
            methods += [
                'ONT reads were demultiplexed and trimmed using qcat (v ' + self.analysis.versions['qcat'] + ').']
        self.methods[self.basecalling_methods_title] = pd.Series(methods, dtype='object')
        self.add_illumina_library_information()
        self.add_assembly_information()
        self.add_contig_info()
        self.add_assembly_notes()

        if self.analysis.did_flye_ont_fastq:
            method = 'ONT reads were assembled using Flye (v ' + self.analysis.versions['flye'] + ').'
            self.methods[self.assembly_methods_title] = pd.concat([self.methods[self.assembly_methods_title],
                pd.Series(method, dtype='object')])
            """ self.methods[self.assembly_methods_title] = self.methods[self.assembly_methods_title].append(
                pd.Series(method)) """
        if self.analysis.did_medaka_ont_assembly:
            method = 'The genome assembly was polished using ONT reads and Medaka (v ' + \
                     self.analysis.versions['medaka'] + ').'
            self.methods[self.assembly_methods_title] = pd.concat([self.methods[self.assembly_methods_title],
                pd.Series(method, dtype='object')])
            """ self.methods[self.assembly_methods_title] = self.methods[self.assembly_methods_title].append(
                pd.Series(method)) """
        if self.analysis.did_pilon_ont_assembly:
            if self.analysis.illumina_read_length_mean <= 50 :
                method = 'Illumina reads were mapped to the genome assembly using bwa aln.'
            else : # We have longer short reads
                method = 'Illumina reads were mapped to the genome assembly using minimap2 (v ' + \
                    str(self.analysis.versions['minimap2']) + ').'
            self.methods[self.pilon_methods_title] = pd.concat([self.methods[self.pilon_methods_title],
                pd.Series(method, dtype='object')])
            
            method = 'The Illumina mappings were then used to error-correct the assembly with Pilon (v ' + \
                str(self.analysis.versions['pilon']) + ').'
            self.methods[self.pilon_methods_title] = pd.concat([self.methods[self.pilon_methods_title],
                pd.Series(method, dtype='object')])

    def make_tex(self): #Converted
        self.doc.new_table_of_contents(table_title='Detailed Run Information', depth=2,marker="TableOfContents")
        text = self.doc.file_data_text
        text = text.replace("##--[","")
        text = text.replace("]--##","")
        self.doc.file_data_text = text
        self.doc.create_md_file()

    def make_report(self): # No need to Convert

        self.start_doc()
        self.add_summary()
        self.add_contamination()
        self.add_alignment()
        self.add_features()
        self.add_feature_plots()
        self.add_mutations()
        self.add_large_indels()
        self.add_plasmids()
        self.add_amr_matrix()
        # self.add_snps()
        self.add_methods()
        self.make_tex()