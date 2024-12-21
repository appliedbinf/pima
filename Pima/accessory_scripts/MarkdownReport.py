import re
import pandas as pd
from mdutils.mdutils import MdUtils
#from mdutils.tools import TableOfContents
from si_prefix import si_format

from Pima.utils.utils import format_kmg
import os

cdc_advisory = 'The analysis and report presented here should be treated as preliminary.  Please contact the CDC/BDRD with any results regarding _Bacillus anthracis_.'


class PimaReport:

    def __init__(self, analysis, settings):

        self.analysis = analysis
        self.doc = None
        self.pima_version = settings.pima_version

        self.summary_title = 'Summary'
        self.basecalling_title = 'Basecalling'
        self.assembly_notes_title = 'Assembly notes'
        self.alignment_title = 'Comparison with reference'
        self.reference_align_title = 'Reference sequences identified within query'
        self.query_align_title = 'Query sequences identified within reference'
        self.alignment_notes_title = 'Alignment notes'
        self.contig_alignment_title = 'Alignment vs. reference contigs'
        self.large_indel_title = 'Large insertions & deletions'
        self.large_indel_notes_title = 'Large insertions & deletions notes'
        self.snp_indel_title = 'SNPs and small indels'
        self.feature_title = 'Features found in the assembly'
        self.feature_plot_title = 'Feature annotation plots'
        self.mutation_title = 'AMR Conferring mutations found in the sample'
        self.amr_matrix_title = 'AMR matrix'
        self.appendix_title = "Appendices"

        self.methods = pd.Series(dtype='float64')
        self.methods_title = 'Methods'
        self.contamination_methods_title = 'Contamination check'
        self.methods[self.contamination_methods_title] = pd.Series(dtype='float64')
        self.assembly_methods_title = 'Assembly'
        self.methods[self.assembly_methods_title] = pd.Series(dtype='float64')
        self.illumina_polishing_methods_title = 'Polish assembly with Illumina'
        self.methods[self.illumina_polishing_methods_title] = pd.Series(dtype='float64')
        self.reference_methods_title = 'Reference comparison'
        self.methods[self.reference_methods_title] = pd.Series(dtype='float64')
        self.mutation_methods_title = 'Mutation screening'
        self.methods[self.mutation_methods_title] = pd.Series(dtype='float64')
        self.feature_methods_title = 'Feature annotation'
        self.methods[self.feature_methods_title] = pd.Series(dtype='float64')
        self.plasmid_methods_title = 'Plasmid annotation'
        self.methods[self.plasmid_methods_title] = pd.Series(dtype='float64')
        self.circos_methods_title = "Genome Visualization"
        self.methods[self.circos_methods_title] = pd.Series(dtype='float64')
        self.appendices = []

    def start_doc(self): #Converted
        header_text = 'Analysis of ' + self.analysis.analysis_name
        self.doc = MdUtils(file_name=self.analysis.report_md,title=header_text)
        self.doc.new_paragraph(f'PiMA Version: {self.pima_version}')

    def add_tableOfContents(self):
        self.doc.create_marker(text_marker="TableOfContents")
        self.doc.new_line()
        #self.doc.new_line('<div style="page-break-after: always;"></div>')
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
            "ONT FASTQ",
            self.wordwrap_markdown(self.analysis.ont_raw_fastq),
            "Illumina FASTQ",
            self.wordwrap_markdown(self.analysis.illumina_fastq),
            "Assembly",
            self.wordwrap_markdown(self.analysis.genome_fasta),
            "Reference",
            self.wordwrap_markdown(self.analysis.reference_fasta),
        ]
        self.doc.new_table(columns=2,rows=6,text=Table_list,text_align='left')
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
            f"{format_kmg(number = self.analysis.ont_raw_bases, decimals=1)}",
        ]
        self.doc.new_table(columns=2, rows=4, text=Table_List, text_align='left')
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

        genome_size = si_format(sum([len(x) for x in self.analysis.genome]), # This sums up the lengths of the genome
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

        ## Assembly related methods
        if self.analysis.did_flye_ont_fastq:
            method = 'ONT reads were assembled using Flye (v ' + self.analysis.versions['flye'] + ').'
            self.methods[self.assembly_methods_title] = pd.concat([self.methods[self.assembly_methods_title],
                pd.Series(method, dtype='object')])
        if self.analysis.did_raven_ont_fastq:
            method = 'ONT reads were assembled using Raven (v ' + self.analysis.versions['raven'] + ').'
            self.methods[self.assembly_methods_title] = pd.concat([self.methods[self.assembly_methods_title],
                pd.Series(method, dtype='object')])
        if self.analysis.did_medaka_ont_assembly:
            method = (f"The genome assembly was polished using ONT reads and Medaka (v '{self.analysis.versions['medaka']}')."
                      f"The basecalling model used was: {self.analysis.ont_model}")
            self.methods[self.assembly_methods_title] = pd.concat([self.methods[self.assembly_methods_title],
                pd.Series(method, dtype='object')])
        if self.analysis.did_spades_illumina_fastq:
            method = f"The genome was assembled using Illumina reads and SPAdes v( '{self.analysis.versions['spades']})."
            self.methods[self.assembly_methods_title] = pd.concat([self.methods[self.assembly_methods_title],
                pd.Series(method, dtype='object')])
        
        #Illumina polishing    
        if self.analysis.did_pilon_ont_assembly:
            if self.analysis.illumina_length_mean <= 50 :
                method = 'Illumina reads were mapped to the genome assembly using bwa aln.'
            else : # We have longer short reads
                method = 'Illumina reads were mapped to the genome assembly using minimap2 (v ' + \
                    str(self.analysis.versions['minimap2']) + ').'
            self.methods[self.illumina_polishing_methods_title] = pd.concat([self.methods[self.illumina_polishing_methods_title],
                pd.Series(method, dtype='object')])
            
            method = 'The Illumina mappings were then used to error-correct the assembly with Pilon (v ' + \
                str(self.analysis.versions['pilon']) + ').'
            self.methods[self.illumina_polishing_methods_title] = pd.concat([self.methods[self.illumina_polishing_methods_title],
                pd.Series(method, dtype='object')])
        if self.analysis.did_polypolish_ont_assembly:
            method = (
                f"Illumina reads were mapped to the genome using bwa mem (v{self.analysis.versions['bwa']})",
                f"and the assembly was corrected with polypolish (v{self.analysis.versions['polypolish']})",
            )
            self.methods[self.illumina_polishing_methods_title] = pd.concat([self.methods[self.illumina_polishing_methods_title],
                pd.Series(method, dtype='object')])

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
            formatted = self.analysis.contig_info[method].copy() #copy the contig information from the assembly (ONT or Illumina)
            formatted['size'] = formatted.apply(lambda x: "{:,}".format(x['size']), axis=1) #add commas to the size (e.g 50000 -> 50,000)
            for i in range(formatted.shape[0]):
                Table_List.extend(formatted.iloc[i, :].values.astype(str).tolist())
            row_count = int(len(Table_List)/3)
            self.doc.new_table(columns=3, rows=row_count, text=Table_List, text_align='left')

    def add_assembly_notes(self): # Converted

        if len(self.analysis.assembly_notes) == 0:
            return

        self.doc.new_line()
        #self.doc.new_line('<div style="page-break-after: always;"></div>')
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

        if self.analysis.genome_fasta is None:
            return
        
        if len(self.analysis.contig_alignment) == 0:
            return
        
        ##If a reference was specified
        if self.analysis.reference is not None:
        
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

            ### Add contig specific alignment stats for both query and reference
            if len(self.analysis.reference_alignment_stats) > 0:
                reference_alignments = self.analysis.reference_alignment_stats
            else:
                return
            
            self.doc.new_line()
            self.doc.new_header(level=3, title=self.reference_align_title)
            Table_List = [
                    "Reference Contig",
                    "Size (bp)",
                    "Bases Aligned to Query",
                    "Perc Align"
                ]
            for index, row in reference_alignments.iterrows():
                Table_List = Table_List + row.tolist()
            row_count = int(len(Table_List)/4)
            self.doc.new_table(columns=4, rows=row_count, text=Table_List, text_align='left')

            ### Add contig specific alignment stats for both query and reference
            if len(self.analysis.query_alignment_stats) > 0:
                query_alignments = self.analysis.query_alignment_stats
            else:
                return
            
            self.doc.new_line()
            self.doc.new_header(level=3, title=self.query_align_title)
            Table_List = [
                    "Query Contig",
                    "Size (bp)",
                    "Bases Aligned to Ref",
                    "Perc Align"
                ]
            #re-order query alignments based on contig size
            for index, row in query_alignments.sort_values(by="Size (bp)", key=lambda s: s.str[:].astype(float), ascending=False).iterrows():
                Table_List = Table_List + row.tolist()
            row_count = int(len(Table_List)/4)
            self.doc.new_table(columns=4, rows=row_count, text=Table_List, text_align='left')
            
            #self.doc.new_line('<div style="page-break-after: always;"></div>')
            self.doc.new_line()

            if len(self.analysis.alignment_notes) > 0:
                self.doc.new_header(level=3, title=self.alignment_notes_title)
                for note in self.analysis.alignment_notes:
                    self.doc.new_line(note)

            method = 'The genome assembly was aligned against the reference sequencing using dnadiff (v ' \
                    + self.analysis.versions['dnadiff'] + ').'
            
            self.methods[self.reference_methods_title] = pd.concat([self.methods[self.reference_methods_title], 
                pd.Series(method, dtype='object')])

    def add_circos(self):
        if self.analysis.did_circos_plots:
            if len(self.analysis.contig_alignment) > 0:
                alignments = self.analysis.contig_alignment
            else:
                return
            
            for contig in alignments.index.tolist():
                contig_title = 'Alignment to ' + contig
                image_png = alignments[contig]
                self.doc.new_line()
                self.doc.new_header(level=3,title=contig_title)
                self.doc.new_line()
                self.doc.write(
                    self.doc.new_inline_image(
                        text='contig_title',
                        path=os.path.abspath(image_png)
                    )
                    ,wrap_width=0
                )
                #self.doc.new_line('<div style="page-break-after: always;"></div>')
                self.doc.new_line()
    
            if self.analysis.self_circos:
                method = f"Sequence coverage plots were visualized using pycircos v0.3"
                self.methods[self.circos_methods_title] = pd.concat([self.methods[self.circos_methods_title], pd.Series(method,dtype="object")])

            else:
                method = f"Alignments of assembled genome and sequence coverage information to provided reference genome were visualized using pycircos v0.3"
                self.methods[self.circos_methods_title] = pd.concat([self.methods[self.circos_methods_title], pd.Series(method,dtype="object")])

    def add_features(self): #Converted

        #if we didn't produce an assembly, we don't search for genes and this series is empty
        if len(self.analysis.feature_hits) == 0:
            return
        #if we did produce an assembly we search for both amr and inc, so feature_hits has a length
        #check if both inc and amr dataframes are 0
        if all(size == 0 for size in (len(self.analysis.feature_hits['amr']), len(self.analysis.feature_hits['inc']))):
            return

        #Do we really want to include the inc hits within the features table? Would this be better split?
        self.doc.new_line()
        self.doc.new_header(level=2,title=self.feature_title)

        for feature_name in self.analysis.feature_hits.index.tolist():

            features = self.analysis.feature_hits[feature_name].copy()
            if features.shape[0] == 0:
                continue

            features[1] = features.apply(lambda x: '{:,}'.format(x[1]), axis=1) #start position of feature
            features[2] = features.apply(lambda x: '{:,}'.format(x[2]), axis=1) #end position

            self.doc.new_line()
            self.doc.new_header(level=3,title=feature_name)

            if (features.shape[0] == 0):
                continue

            for contig in pd.unique(features[0]): #contig name

                self.doc.new_line(contig)

                #subset features dataframe by the contig IDs
                contig_features = features.loc[(features[0] == contig)]
                Table_List = [
                    'Start', 'Stop', 'Feature', 'Identity (%)', 'Strand',
                ]

                for i in range(contig_features.shape[0]): #take each feature from the contig table
                    feature = contig_features.iloc[i, ].copy(deep=True) 
                    feature[3] = "_".join(feature[3].split("_")[:-1])
                    feature[4] = '{:.3f}'.format(feature[4]) #round the percID to 3 dec.
                    Table_List.extend(feature[1:].values.astype(str).tolist()) #pandas angry if different types

                row_count = int(len(Table_List) / 5)
                self.doc.new_line()
                self.doc.new_table(columns=5, rows=row_count, text=Table_List, text_align='left')

        method = 'The genome assembly was queried for features using blastn (v ' + self.analysis.versions[
            'blastn'] + ').  ' + \
                 'Feature hits were clustered using bedtools (v ' + self.analysis.versions['bedtools'] + ') ' + \
                 'and the highest scoring hit for each cluster was reported.'
        self.methods[self.feature_methods_title] = pd.concat([self.methods[self.feature_methods_title], pd.Series(method, dtype='object')])

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
        
        method = f"Detected features were visualized using dna_features_viewer v({self.analysis.versions['dna_features_viewer']})."
        self.methods[self.feature_methods_title] = pd.concat([self.methods[self.feature_methods_title], pd.Series(method, dtype='object')])

    def add_mutations(self):

        # Make sure we looked for mutations
        if not self.analysis.did_call_mutations:
            return

        mutations = self.analysis.amr_mutations
        #if the INDEL is too long, break it 
        mutations['var'] = mutations['var'].apply(lambda x: f"{x[0:8]}...{len(x)-8}" if len(x) > 8 else x)
        mutations['GE'] = mutations['GE'].apply(lambda x: f"{x[0:5]}" if len(x) > 5 else x)
        self.doc.new_line()
        self.doc.new_header(level=2,title=self.mutation_title)
        if mutations.size == 0:
            note = f"No mutations were confidently identified in any of the regions specified by the provided mutations bed file: {self.analysis.mutation_region_bed}"
            self.doc.new_paragraph(note)
            return
        else:
            # hide homopolymer indels from ONT variant results UNLESS an R10 SUP model was used
            if not re.search(r"r1041.*sup.*", self.analysis.ont_model):
                note = (
                    "The ONT basecaller was either: (1) not detected, (2) from an R9.4.1 flowcell or, (3) not in super accurarcy mode."
                    "Mutations in short homopolymer INDELs are therefore not shown due to the high chance of a false positive."
                    "If you suspect an important INDEL variant, the intermediate file: " 
                    f"`{os.path.join(self.analysis.mutations_dir, 'intersect_amr-regions_variants_deletions.tsv')}` contains the full output"
                )
                self.doc.new_paragraph(note)
                ##remove the +/- from the variant call & check if it is a homopolymer (same base), this does replace SNPs with NA, but we sort that out later
                mutations.loc[:, 'var_mod'] = mutations['var'] \
                                            .astype(str) \
                                            .apply(lambda x: f'assembly_del{x}' if (x.isnumeric()) else x) \
                                            .apply(lambda x: ''.join(filter(str.isalpha, x))) \
                                            .apply(lambda x: x if (x != len(x)*x[0]) else 'na')
                #mutations = mutations.query("var_type != 'indel' | var_mod != 'na'") #drop homopolymer indels in confirmed indels
                mutations = mutations[~((mutations['var'].str.match(r"[+-]")) & (mutations['var_mod'] == 'na'))] #drop homopolymer indels
            #restrict to only known mutations if the input genome is too far from the reference genome
            if self.analysis.reference_identity < self.analysis.reference_identity_min and self.analysis.reference_identity != 0:
                note1 = (
                    f"Due to < {self.analysis.reference_identity_min}% identity to the provided reference, "
                    "mutation detection within the regions specified is restricted to only those that have been experientally confirmed. "
                    "Futhermore, due to the distances, there is no guarantee that these mutations will confer the same resistances in this isolate "
                    "as has been observed for the reference organism. We recommend extreme caution when interpreting these results."
                )
                note2 = (
                    f"All detected mutations within the regions specified by the `{self.analysis.mutation_region_bed}` file can be found in the file "
                    f"`{os.path.join(self.analysis.mutations_dir, 'intersect_amr-regions_variants_deletions.tsv')}`. All detected mutations are reported " 
                    f"in the intermediate file `{os.path.join(self.analysis.mutations_dir,'varscan.vcf')}`."
                )
                self.doc.new_paragraph(note1)
                self.doc.new_paragraph(note2)
                mutations = mutations.query('classification_type == "confirmed_location"')

            drugs = set(mutations['amr_class'])
            mutations = mutations.assign(loc=mutations['loc'].astype(int).apply(lambda x: '{:,}'.format(x)))
            for amr_class in set(mutations['amr_class']):
                self.doc.new_header(level=3,title=amr_class.title())
                hits = mutations.query('amr_class == @amr_class')
                Table_List = [
                    'Region', 'Mutation Type', 'Position', 'Reference', 'Variant', 'Supporting Reads', 'Note',
                ]
                if hits.query('classification_type == "confirmed_location"').shape[0] > 0:
                    sub_hits = hits.query('classification_type == "confirmed_location"')
                    sub_hits.loc[:,'region_name'] = sub_hits['region_name'].astype(str).apply(lambda x: f"**{x}**")
                    Table_List.extend(sub_hits[['region_name', 'var_type', 'loc', 'ref', 'var', 'AF', 'note']].stack().values.tolist())

                if hits.query('classification_type == "potential_confer_amr"').shape[0] > 0:
                    sub_hits = hits.query('classification_type == "potential_confer_amr"')
                    Table_List.extend(sub_hits[['region_name', 'var_type', 'loc', 'ref', 'var', 'AF', 'note']].stack().values.tolist())
                    
                row_count = int(len(Table_List)/7)
                self.doc.new_table(columns=7, rows=row_count, text=Table_List, text_align='left')

            if any(mutations.classification_type == "confirmed_location") and any(mutations.classification_type == "potential_confer_amr"):
                message = (
                    "**Bold** gene names indicate mutations identified at locations that have been experimentally confirmed to confer resistance.\n"
                    "Unbolded gene names indicate mutations identified within regions that have been shown to confer resistance, but have not specifically been observed."
                )
            elif any(mutations.classification_type == "confirmed_location"):
                message = "**Bold** gene names indicate mutations identified at locations that have been experimentally confirmed to confer resistance"
            elif any(mutations.classification_type == "potential_confer_amr"):
                message = "Unbolded gene names indicate mutations identified within regions that have been shown to confer resistance, but have not specifically been observed"
            else:
                message = ""
                
            if len(message) > 0:
                self.doc.new_paragraph(message)

            if self.analysis.organism_amr_appendices is not None:
                self.appendices.extend(
                    [appendix for drug in drugs for appendix in self.analysis.organism_amr_appendices if drug in appendix]
                )
                note = (
                    "See the Appendix at the end of the report for more information about putative AMR conferring mutations. "
                    "This only concerns SNP/INDEL variants and not AMR conferring genes which will be reported in the 'Features' "
                    "section of the report"
                )
            self.doc.new_paragraph(note)

        method = self.analysis.mutations_read_type + ' reads were mapped to the reference sequence using minimap2 (v ' \
                 + self.analysis.versions['minimap2'] + ').'
        self.methods[self.mutation_methods_title] = pd.concat([self.methods[self.mutation_methods_title], pd.Series(method,dtype="object")])

        method = ' '.join(['Mutations were identified using '
                           'samtools mpileup (v', self.analysis.versions['samtools'], ')',
                           'and varscan (v', self.analysis.versions['varscan'], ').'])
        self.methods[self.mutation_methods_title] = pd.concat([self.methods[self.mutation_methods_title], pd.Series(method,dtype="object")])

    def add_appendix_title(self):
        if len(self.appendices) == 0:
            return

        else:
            self.doc.new_line('<div style="page-break-after: always;"></div>')
            self.doc.new_line()
            self.doc.new_header(2, self.appendix_title)

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
        method = (
            f"Detected SNPs or INDELs compared to the provided reference within regions specified by the mutations_regions.bed file were reported"
            f"and then visualized in a heatmap using matplotlib v({self.analysis.versions['matplotlib']})."
        )
        self.methods[self.mutation_methods_title] = pd.concat([self.methods[self.mutation_methods_title], pd.Series(method,dtype="object")])

    def add_large_indels(self): #Converted

        # Make sure we looked for mutations
        if len(self.analysis.large_indels) == 0:
            return
        
        large_indels = self.analysis.large_indels
        self.doc.new_line()
        self.doc.new_header(level=2,title=self.large_indel_title)

        if len(self.analysis.large_indel_notes) > 0:
            self.doc.new_header(level=3, title=self.large_indel_notes_title)
            for note in self.analysis.large_indel_notes:
                self.doc.new_line(note)

        for genome in ['Reference insertions', 'Query insertions']:

            genome_indels = large_indels[genome].copy()
            self.doc.new_line()
            self.doc.new_header(level=3,title=genome)
            if genome == 'Reference insertions':
                #self.doc.new_header(level=4, title="*(Deletions in query)*")
                self.doc.new_paragraph('*(Deletions in query)*')
                self.doc.new_line()
            if (genome_indels.shape[0] == 0):
                continue
            genome_indels['start'] = genome_indels['start'].astype('object').apply(lambda x: '{:,}'.format(x))
            genome_indels['stop'] = genome_indels['stop'].astype('object').apply(lambda x: '{:,}'.format(x))
            genome_indels['length'] = genome_indels['length'].astype('object').apply(lambda x: '{:,}'.format(x))

            Table_List = [
                'Reference contig', 'Start', 'Stop', 'Size (bp)'
            ]

            for i in range(genome_indels.shape[0]):
                Table_List.extend(genome_indels.iloc[i].values.tolist())

            row_count= int(len(Table_List)/4)
            self.doc.new_table(columns=4, rows=row_count, text=Table_List, text_align='left')

        method = 'Large insertions or deletions were found as the complement of aligned ' + \
                 'regions using bedtools (v ' + self.analysis.versions['bedtools'] + ').'
        self.methods[self.reference_methods_title] = pd.concat([self.methods[self.reference_methods_title], 
            pd.Series(method,dtype='object')])
        
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
        # First section of Summary
        self.doc.new_header(level=1, title='CDC Advisory')
        self.doc.new_paragraph(cdc_advisory)
        self.doc.new_line()
        self.add_run_information()
        self.add_ont_library_information()
        self.add_illumina_library_information()
        self.add_assembly_information()
        self.add_contig_info()
        self.add_assembly_notes()
            
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
        self.add_circos()
        self.add_features()
        self.add_feature_plots()
        self.add_mutations()
        self.add_amr_matrix()
        self.add_large_indels()
        self.add_plasmids()
        self.add_methods()
        self.add_appendix_title()
        self.make_tex()