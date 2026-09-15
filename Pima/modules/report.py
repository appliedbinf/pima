import os
import shutil
import re
import string
import pandas as pd
from mdutils.mdutils import MdUtils
from si_prefix import si_format
from pathlib import Path
from Pima.pima_data import PimaData
from Pima.utils.settings import Settings
from Pima.utils.utils import (
    format_kmg,
    print_and_log,
    print_and_run,
    validate_utility,
    find_checkpoint,
    validate_file_and_size_or_error,
    std_files,
    correct_nextflow_path,
)

class PimaReport:

    def __init__(self, pima_data: PimaData, settings: Settings):

        self.pima_data = pima_data
        self.doc = None
        self.pima_version = settings.pima_version

        #Report Sections ections
        self.run_summary = "Run Summary"
        self.run_metrics_title = "Run Metrics"
        self.feature_detection_title = "Feature Detection"
        self.reference_comparison_title = "Comparison with reference"
        self.methods_title = 'Methods'

        #Report layout
        self._run_metrics = [
            self.add_ont_library_information,
            self.add_illumina_library_information,
            self.add_assembly_information,
            self.add_contig_info,
            self.add_assembly_notes,
            self.add_contamination
        ]
        self._feature_detection = [
            self.add_amr_matrix,
            self.add_feature_plots,
            self.add_virulence_gene_hits,
            self.add_resfinder_amr,
            self.add_amrfinder_results,
            self.add_inc,            
            self.add_plasmids
        ]
        self._reference_comparison = [
            self.add_alignment,
            self.add_circos,
            self.add_mutations,
            self.add_large_indels
        ]

        #Track methods used
        self.methods = pd.Series(dtype='float64')
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
        self.feature_methods_title = 'AMR/INC Gene annotation'
        self.methods[self.feature_methods_title] = pd.Series(dtype='float64')
        self.plasmid_methods_title = 'Plasmid annotation'
        self.methods[self.plasmid_methods_title] = pd.Series(dtype='float64')
        self.circos_methods_title = "Genome Visualization"
        self.methods[self.circos_methods_title] = pd.Series(dtype='float64')
        self.appendices = []

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

    def start_doc(self): 
        header_text = f'Analysis of {self.pima_data.analysis_name}'
        self.doc = MdUtils(file_name=self.pima_data.report_md, title=header_text)
        self.doc.new_paragraph(f'PiMA Version: {self.pima_version}')

        cdc_advisory = (
            'The analysis and report presented here should be treated as preliminary.  '
            'Please contact the CDC/BDRD with any results regarding _Bacillus anthracis_.'
        )
        self.doc.new_header(level=1, title='CDC Advisory', add_table_of_contents='n')
        self.doc.new_paragraph(cdc_advisory)
        self.doc.new_line()

    def add_tableOfContents(self):
        self.doc.create_marker(text_marker="TableOfContents")
        self.doc.new_line()
        self.doc.new_line()

    def add_run_information(self): 
        self.doc.new_line()
        self.doc.new_header(level=1, title='Run Information', add_table_of_contents='n')

        if self.pima_data.organism:
            Table_list = [
                "Category",
                "Information",
                "Date",
                self.pima_data.start_time,
                "ONT FASTQ",
                self.wordwrap_markdown(self.pima_data.ont_raw_fastq or ''),
                "Illumina FASTQ",
                self.wordwrap_markdown(self.pima_data.illumina_fastq or ''),
                "Assembly",
                self.wordwrap_markdown(correct_nextflow_path(self.pima_data.genome_fasta) or ''),
                "Organism Version",
                self.wordwrap_markdown(self.pima_data.organism_version['version'] or ''),
                "Reference",
                self.wordwrap_markdown(self.pima_data.reference_fasta or ''),
                "Mutation Bed File",
                self.wordwrap_markdown(self.pima_data.mutation_region_bed or '')
            ]
            self.doc.new_table(columns=2,rows=8,text=Table_list,text_align='left')
            self.doc.new_line()
            self.add_tableOfContents()
            self.doc.new_line()

        else:
            Table_list = [
                "Category",
                "Information",
                "Date",
                self.pima_data.start_time,
                "ONT FASTQ",
                self.wordwrap_markdown(self.pima_data.ont_raw_fastq or ''),
                "Illumina FASTQ",
                self.wordwrap_markdown(self.pima_data.illumina_fastq or ''),
                "Assembly",
                self.wordwrap_markdown(self.pima_data.genome_fasta or ''),
                "Reference",
                self.wordwrap_markdown(self.pima_data.reference_fasta or ''),
                "Mutation Bed File",
                self.wordwrap_markdown(self.pima_data.mutation_region_bed or '')
            ]
            self.doc.new_table(columns=2,rows=7,text=Table_list,text_align='left')
            self.doc.new_line()
            self.add_tableOfContents()
            self.doc.new_line()            

# Section Headers & Info Orchestration
    def add_run_metrics(self):
        self.doc.new_line()
        self.doc.new_header(level=1, title=self.run_metrics_title)
        for metric_func in self._run_metrics:
            metric_func()

    def add_feature_detection(self):
        self.doc.new_line()
        self.doc.new_header(level=1, title=self.feature_detection_title)
        for metric_func in self._feature_detection:
            metric_func()

    def add_reference_comparison(self):
        self.doc.new_line()
        self.doc.new_header(level=1, title=self.reference_comparison_title)
        for metric_func in self._reference_comparison:
            metric_func()

# Subsections - info generators
    def add_ont_library_information(self): #done

        if self.pima_data.ont_n50 is None:
            return
        self.doc.new_line()
        self.doc.new_header(2, 'ONT library statistics')
        Table_List = [
            "Category",
            "Quantity",
            "ONT N50",
            '{:,}'.format(self.pima_data.ont_n50),
            "ONT reads",
            '{:,}'.format(self.pima_data.ont_read_count),
            "ONT bases",
            f"{format_kmg(number = self.pima_data.ont_raw_bases, decimals=1)}",
        ]
        self.doc.new_table(columns=2, rows=4, text=Table_List, text_align='left')
        self.doc.new_line()

    def add_illumina_library_information(self): #done
        if self.pima_data.illumina_length_mean is None:
            return

        self.doc.new_line()
        self.doc.new_header(2, 'Illumina library statistics')
        Table_List = [
            "Illumina Info.",
            "Quantity",
            'Illumina mean length',
            '{:.1f}'.format(self.pima_data.illumina_length_mean),
            'Illumina reads',
            '{:,}'.format(self.pima_data.illumina_read_count),
            'Illumina bases',
            '{:s}'.format(self.pima_data.illumina_bases)
        ]
        self.doc.new_table(columns=2, rows=4, text=Table_List, text_align='left')

    def add_assembly_information(self): #done
        if self.pima_data.genome is None:
            return

        self.doc.new_line()
        self.doc.new_header(2, 'Assembly statistics')

        genome_size = si_format(sum([len(x) for x in self.pima_data.genome]),
                                precision=1)
        Table_List = [
            "Category",
            "Information",
            "Contigs",
            str(len(self.pima_data.genome)),
            "Assembly size",
            genome_size
        ]
        self.doc.new_table(columns=2, rows=3, text=Table_List, text_align='left')

        ## Assembly related methods
        if self.pima_data.did_flye_ont_fastq:
            method = 'ONT reads were assembled using Flye (v ' + self.pima_data.versions['flye'] + ').'
            self.methods[self.assembly_methods_title] = pd.concat([self.methods[self.assembly_methods_title],
                pd.Series(method, dtype='object')])
        if self.pima_data.did_raven_ont_fastq:
            method = 'ONT reads were assembled using Raven (v ' + self.pima_data.versions['raven'] + ').'
            self.methods[self.assembly_methods_title] = pd.concat([self.methods[self.assembly_methods_title],
                pd.Series(method, dtype='object')])
        if self.pima_data.did_medaka_ont_assembly:
            method = (f"The genome assembly was polished using ONT reads and Medaka (v '{self.pima_data.versions['medaka']}')."
                      f"The basecalling model used was: {self.pima_data.ont_model}")
            self.methods[self.assembly_methods_title] = pd.concat([self.methods[self.assembly_methods_title],
                pd.Series(method, dtype='object')])
        if self.pima_data.did_spades_illumina_fastq:
            method = f"The genome was assembled using Illumina reads and SPAdes v( '{self.pima_data.versions['spades']})."
            self.methods[self.assembly_methods_title] = pd.concat([self.methods[self.assembly_methods_title],
                pd.Series(method, dtype='object')])
        
        #Illumina polishing    
        if self.pima_data.did_pilon_ont_assembly:
            if self.pima_data.illumina_length_mean <= 50 :
                method = 'Illumina reads were mapped to the genome assembly using bwa aln.'
            else : # We have longer short reads
                method = 'Illumina reads were mapped to the genome assembly using minimap2 (v ' + \
                    str(self.pima_data.versions['minimap2']) + ').'
            self.methods[self.illumina_polishing_methods_title] = pd.concat([self.methods[self.illumina_polishing_methods_title],
                pd.Series(method, dtype='object')])
            
            method = 'The Illumina mappings were then used to error-correct the assembly with Pilon (v ' + \
                str(self.pima_data.versions['pilon']) + ').'
            self.methods[self.illumina_polishing_methods_title] = pd.concat([self.methods[self.illumina_polishing_methods_title],
                pd.Series(method, dtype='object')])
        if self.pima_data.did_polypolish_ont_assembly:
            method = (
                f"Illumina reads were mapped to the genome using bwa mem (v{self.pima_data.versions['bwa']})",
                f"and the assembly was corrected with polypolish (v{self.pima_data.versions['polypolish']})",
            )
            self.methods[self.illumina_polishing_methods_title] = pd.concat([self.methods[self.illumina_polishing_methods_title],
                pd.Series(method, dtype='object')])

    def add_contig_info(self):

        if self.pima_data.contig_info is None:
            return

        for method in ['ONT', 'Illumina']:

            if not method in self.pima_data.contig_info.index:
                continue

            self.doc.new_line()
            self.doc.new_header(2, 'Assembly coverage by ' + method)
            Table_List = [
                "Contig",
                "Length (bp)",
                "Coverage (X)",
            ]
            formatted = self.pima_data.contig_info[method].copy() #copy the contig information from the assembly (ONT or Illumina)
            formatted['size'] = formatted.apply(lambda x: "{:,}".format(x['size']), axis=1) #add commas to the size (e.g 50000 -> 50,000)
            for i in range(formatted.shape[0]):
                Table_List.extend(formatted.iloc[i, :].values.astype(str).tolist())
            row_count = int(len(Table_List)/3)
            self.doc.new_table(columns=3, rows=row_count, text=Table_List, text_align='left')

    def add_assembly_notes(self):
        assembly_notes_title = 'Assembly notes'
        if len(self.pima_data.assembly_notes) == 0:
            return

        self.doc.new_line()
        self.doc.new_line()
        self.doc.new_header(2, assembly_notes_title)

        for note in self.pima_data.assembly_notes:
            self.doc.new_line(note)

    def add_contamination(self):

        if len(self.pima_data.kraken_fracs) == 0:
            return

        self.doc.new_line()
        self.doc.new_header(2,'Contamination check')
        for read_type, kraken_fracs in self.pima_data.kraken_fracs.items():

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

        method = 'Kraken2 (' + self.pima_data.versions['kraken2'] + ') was used to assign the raw reads into taxa.'
        self.methods[self.contamination_methods_title] = pd.concat([self.methods[self.contamination_methods_title],
            pd.Series(method, dtype="object")])

    def add_alignment(self):
        alignment_title = 'Comparison with reference'
        snp_indel_title = 'SNPs and small indels'
        reference_align_title = 'Reference sequences identified within query'
        query_align_title = 'Query sequences identified within reference'    
        alignment_notes_title = 'Alignment notes'

        if self.pima_data.genome_fasta is None:
            return
        
        if len(self.pima_data.contig_alignment) == 0:
            return
        
        ##If a reference was specified
        if self.pima_data.reference is not None:
        
            self.doc.new_line()
            self.doc.new_header(level=2, title=alignment_title)
            self.doc.new_line()
            self.doc.new_header(level=3, title=snp_indel_title, add_table_of_contents='n')

            Table_1 = [
                "Category",
                "Quantity",
                'SNPs',
                '{:,}'.format(self.pima_data.quast_mismatches),
                'Small indels',
                '{:,}'.format(self.pima_data.quast_indels)
            ]

            self.doc.new_table(columns=2, rows=3, text=Table_1, text_align='left')

            ### Add contig specific alignment stats for both query and reference
            if len(self.pima_data.reference_alignment_stats) > 0:
                reference_alignments = self.pima_data.reference_alignment_stats
            else:
                return
            
            self.doc.new_line()
            self.doc.new_header(level=3, title=reference_align_title, add_table_of_contents='n')
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
            if len(self.pima_data.query_alignment_stats) > 0:
                query_alignments = self.pima_data.query_alignment_stats
            else:
                return
            
            self.doc.new_line()
            self.doc.new_header(level=3, title=query_align_title, add_table_of_contents='n')
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
            self.doc.new_line()

            if len(self.pima_data.alignment_notes) > 0:
                self.doc.new_header(level=3, title=alignment_notes_title, add_table_of_contents='n')
                for note in self.pima_data.alignment_notes:
                    self.doc.new_line(note)

            method = 'The genome assembly was aligned against the reference sequencing using dnadiff (v ' \
                    + self.pima_data.versions['dnadiff'] + ').'
            
            self.methods[self.reference_methods_title] = pd.concat([self.methods[self.reference_methods_title], 
                pd.Series(method, dtype='object')])

    def add_circos(self):
        if self.pima_data.did_circos_plots:
            if len(self.pima_data.contig_alignment) > 0:
                alignments = self.pima_data.contig_alignment
            else:
                return
            
            for contig in alignments.index.tolist():
                contig_title = 'Alignment to ' + contig
                image = alignments[contig]
                self.doc.new_header(level=3, title=contig_title, add_table_of_contents='n')
                self.doc.write(
                    f"\n<img src='{image}' style='width: 3in; height: auto;' />\n"
                )
    
            if self.pima_data.self_circos:
                method = f"Sequence coverage plots were visualized using pycircos v0.3"
                self.methods[self.circos_methods_title] = pd.concat([self.methods[self.circos_methods_title], pd.Series(method,dtype="object")])

            else:
                method = f"Alignments of assembled genome and sequence coverage information to provided reference genome were visualized using pycircos v0.3"
                self.methods[self.circos_methods_title] = pd.concat([self.methods[self.circos_methods_title], pd.Series(method,dtype="object")])

    def add_resfinder_amr(self):
        resfinder_title = 'AMR genes from resfinder'
        if self.pima_data.resfinder_hits.empty:
            self.doc.new_line()
            self.doc.new_header(level=2,title=resfinder_title)
            self.doc.new_paragraph('No AMR genes detected above 90% identity')
            return
        
        self.doc.new_line()
        self.doc.new_header(level=2,title=resfinder_title)       

        features = self.pima_data.resfinder_hits.copy()
        features[1] = features.apply(lambda x: '{:,}'.format(x[1]), axis=1) #start position of feature
        features[2] = features.apply(lambda x: '{:,}'.format(x[2]), axis=1) #end position
        features[3] = features[3].str.split("_").str[:-1].str.join('_')
        features[4] = features.apply(lambda x: '{:.1f}'.format(x[4]*100), axis=1) #round the percID to 3 dec.
        
        row_count=len(features)+1
        Table_List = [
            'Contig', 'Start', 'Stop', 'Gene', 'Identity (%)', 'Abx Class'
        ]
        col_count=len(Table_List)
        Table_List.extend(features[[0,1,2,3,4,6]].values.flatten().tolist())
        self.doc.new_line()
        self.doc.new_table(columns=col_count, rows=row_count, text=Table_List, text_align='left')

        #Same method for both AMR/INC detection
        method = (
            f"The genome assembly was queried for features using blastn v{self.pima_data.versions['blastn']}. "
            f"Feature hits were clustered using bedtools v{self.pima_data.versions['bedtools']} "
                        f"and the highest scoring hit for each cluster was reported."
        )

        #if below exists skip
        if self.methods[self.feature_methods_title].empty:
            self.methods[self.feature_methods_title] = pd.concat([self.methods[self.feature_methods_title], pd.Series(method, dtype='object')])

    def add_inc(self):
        inc_title = 'Plasmid incompatibility genes'
        if self.pima_data.inc_hits.empty:
            self.doc.new_line()
            self.doc.new_header(level=2,title=inc_title)
            self.doc.new_paragraph('No INC genes detected above 95% identity')
            return
        
        self.doc.new_line()
        self.doc.new_header(level=2,title=inc_title)       

        features = self.pima_data.inc_hits.copy()
        features[1] = features.apply(lambda x: '{:,}'.format(x[1]), axis=1) #start position of feature
        features[2] = features.apply(lambda x: '{:,}'.format(x[2]), axis=1) #end position
        features[3] = features[3].str.split("_").str[:-1].str.join('_')
        features[4] = features.apply(lambda x: '{:.3f}'.format(x[4]), axis=1) #round the percID to 3 dec.


        row_count=len(features)+1
        Table_List = [
            'Contig', 'Start', 'Stop', 'Gene', 'Identity (%)',
        ]
        col_count=len(Table_List)
        Table_List.extend(features[[0,1,2,3,4]].values.flatten().tolist())
        self.doc.new_line()
        self.doc.new_table(columns=col_count, rows=row_count, text=Table_List, text_align='left')

        #Same method for both AMR/INC detection
        method = (
            f"The genome assembly was queried for features using blastn v{self.pima_data.versions['blastn']}. "
            f"Feature hits were clustered using bedtools v{self.pima_data.versions['bedtools']} "
                        f"and the highest scoring hit for each cluster was reported."
        )

        if self.methods[self.feature_methods_title].empty:
            self.methods[self.feature_methods_title] = pd.concat([self.methods[self.feature_methods_title], pd.Series(method, dtype='object')])

        # Report all notes for feature detection
        if not self.pima_data.annotation_notes.empty: 
            self.doc.new_header(level=2,title="Annotation warnings") 
            for note in self.pima_data.annotation_notes:
                self.doc.new_line(note)

    def add_amrfinder_results(self):
        amrfinder_title = 'AMR Finder results'
        if self.pima_data.amrfinder_results.empty:
            return
        
        self.doc.new_line()
        self.doc.new_header(level=2,title=amrfinder_title)
        amrfinder_results = self.pima_data.amrfinder_results.copy()
        amrfinder_results = amrfinder_results[["Contig id", "Start", "Stop", "% Identity to reference", "Element symbol", "Class"]]
        amrfinder_results['Start'] = amrfinder_results['Start'].apply(lambda x: f"{x:,}")
        amrfinder_results['Stop'] = amrfinder_results['Stop'].apply(lambda x: f"{x:,}")
        amrfinder_results['% Identity to reference'] = amrfinder_results['% Identity to reference'].apply(lambda x: f"{x:.3f}")
        amrfinder_results['Class'] = amrfinder_results['Class'].str.lower()

        row_count=len(amrfinder_results)+1
        Table_List = [
            'Contig', 'Start', 'Stop', 'Gene', 'Identity (%)', 'Abx Class',
        ]
        col_count=len(Table_List)
        Table_List.extend(amrfinder_results[['Contig id', 'Start', 'Stop', 'Element symbol', '% Identity to reference', 'Class']].values.flatten().tolist())
        self.doc.new_line()
        self.doc.new_table(columns=col_count, rows=row_count, text=Table_List, text_align='left')

        method = (
            f"The genome assembly was queried for features using amrfinder v.{self.pima_data.versions['amrfinder']}\n"
            f"AMRFinder was run with the database v.{self.pima_data.versions['amrfinder database']}"
        )
        self.methods[self.feature_methods_title] = pd.concat([self.methods[self.feature_methods_title], pd.Series(method, dtype='object')])

    def add_feature_plots(self):

        if len(self.pima_data.feature_plots) == 0:
            return

        self.doc.new_line()
        self.doc.new_header(level=2,title='Feature Plots')
        self.doc.new_paragraph('Only contigs with features are shown')

        for contig in self.pima_data.feature_plots.index.tolist():
            image = self.pima_data.feature_plots[contig]
            self.doc.write(
                    f"\n![Analysis]({os.path.abspath(image)})\n"
            )
            # self.doc.write(
            #     self.doc.new_inline_image(
            #         text='Analysis',
            #         path=os.path.abspath(image_png),
            #     )
            #     ,wrap_width=0
            # )
        
        method = f"Detected features were visualized using dna_features_viewer v({self.pima_data.versions['dna_features_viewer']})."
        self.methods[self.feature_methods_title] = pd.concat([self.methods[self.feature_methods_title], pd.Series(method, dtype='object')])

    def add_virulence_gene_hits(self):
        if not self.pima_data.will_have_genome_fasta:
            return
        if not self.pima_data.ba_virulence_genes:
            return

        organism_label = (self.pima_data.organism or "").replace("_", " ")
        self.doc.new_line()
        self.doc.new_header(level=2, title=f'{organism_label} virulence genes')
        if len(self.pima_data.ba_virulence_hits) == 0:
            self.doc.new_paragraph('No virulence genes detected above 90% identity and 90% coverage')
            return

        features = self.pima_data.ba_virulence_hits.copy()
        features[1] = features.apply(lambda x: '{:,}'.format(x[1]), axis=1) #start position of feature
        features[2] = features.apply(lambda x: '{:,}'.format(x[2]), axis=1) #end position

        for contig in pd.unique(features[0]): #contig name

            self.doc.new_line(f'contig ID: {contig}')
            #subset features dataframe by the contig IDs
            contig_features = features.loc[(features[0] == contig)]
            Table_List = [
                'Virulence Gene', 'Start', 'Stop', 'Identity (%)', 'Strand',
            ]

            for i in range(contig_features.shape[0]): #take each feature from the contig table
                feature = contig_features.iloc[i, ].copy(deep=True) 
                feature[4] = '{:.1f}'.format(feature[4]*100)
                Table_List.extend([str(val) for val in [feature[3],feature[1],feature[2],feature[4],feature[5]]])

            row_count = int(len(Table_List) / 5)
            self.doc.new_line()
            self.doc.new_table(columns=5, rows=row_count, text=Table_List, text_align='left')

    def add_mutations(self):
        mutation_title = 'AMR conferring variants detected relative to reference'

        if not self.pima_data.did_call_mutations:
            return

        mutations = self.pima_data.amr_mutations.copy()
        #if the INDEL is too long, break it 
        mutations['var'] = mutations['var'].apply(lambda x: f"{x[0:8]}...{len(x)-8}" if len(x) > 8 else x)
        mutations['GE'] = mutations['GE'].apply(lambda x: f"{x[0:5]}" if len(x) > 5 else x)
        self.doc.new_line()
        self.doc.new_header(level=2,title=mutation_title)
        if mutations.size == 0:
            note = f"No mutations were confidently identified in any of the regions specified by the provided mutations bed file: {self.pima_data.mutation_region_bed}"
            self.doc.new_paragraph(note)
            return
        else:
            # hide homopolymer indels from ONT variant results UNLESS an R10 SUP model was used
            # need to skip if (1) high qual ONT data, (2) illumina-only data, or (3) illumina polished data
            checks = [
                not re.search(r"r1041.*sup.*", self.pima_data.ont_model),
                self.pima_data.ont_fastq,
                not self.pima_data.did_illumina_polishing, 
            ]
            if all(checks):
                note = (
                    "The ONT basecaller was either: (1) not detected, (2) from an R9.4.1 flowcell or, (3) not in super accurarcy mode."
                    "Mutations in short homopolymer INDELs are therefore not shown due to the high chance of a false positive."
                    "If you suspect an important INDEL variant, the intermediate file: " 
                    f"`{os.path.join(self.pima_data.mutations_dir, 'intersect_amr-regions_variants_deletions.tsv')}` contains the full output"
                )
                self.doc.new_paragraph(note)
                ##remove the +/- from the variant call & check if it is a homopolymer (same base), this does replace SNPs with NA, but we sort that out later
                mutations['var_mod'] = mutations['var'] \
                                            .astype(str) \
                                            .apply(lambda x: f'assembly_del{x}' if (x.isnumeric()) else x) \
                                            .apply(lambda x: ''.join(filter(str.isalpha, x))) \
                                            .apply(lambda x: x if (x != len(x)*x[0]) else 'na')
                mutations = mutations[~((mutations['var'].str.match(r"[+-]")) & (mutations['var_mod'] == 'na'))] #drop homopolymer indels
            #restrict to only known mutations if the input genome is too far from the reference genome
            if self.pima_data.reference_identity < self.pima_data.reference_identity_min and self.pima_data.reference_identity != 0:
                note1 = (
                    f"Due to < {self.pima_data.reference_identity_min}% identity to the provided reference, "
                    "mutation detection within the regions specified is restricted to only those that have been experientally confirmed. "
                    "Futhermore, due to the distances, there is no guarantee that these mutations will confer the same resistances in this isolate "
                    "as has been observed for the reference organism. We recommend extreme caution when interpreting these results."
                )
                note2 = (
                    f"All detected mutations within the regions specified by the `{self.pima_data.mutation_region_bed}` file can be found in the file "
                    f"`{os.path.join(self.pima_data.mutations_dir, 'intersect_amr-regions_variants_deletions.tsv')}`. All detected mutations are reported " 
                    f"in the intermediate file `{os.path.join(self.pima_data.mutations_dir,'varscan.vcf')}`."
                )
                self.doc.new_paragraph(note1)
                self.doc.new_paragraph(note2)
                mutations = mutations[mutations["classification_type"] == "confirmed_location"].copy()

            drugs = set(mutations['amr_class'])
            mutations['loc'] = mutations['loc'].astype(int).apply(lambda x: '{:,}'.format(x))
            for amr_class in set(mutations['amr_class']):
                self.doc.new_header(level=3,title=amr_class.title(), add_table_of_contents='n')
                hits = mutations[mutations['amr_class'] == amr_class].copy()
                Table_List = [
                    'Region', 'Mutation Type', 'Position', 'Reference', 'Variant', 'Supporting Reads', 'Note',
                ]
                if not hits[hits['classification_type'] == "confirmed_location"].empty:
                    sub_hits = hits[hits['classification_type'] == "confirmed_location"].copy()
                    sub_hits['region_name'] = sub_hits['region_name'].astype(str).apply(lambda x: f"**{x}**")
                    Table_List.extend(sub_hits[['region_name', 'var_type', 'loc', 'ref', 'var', 'AF', 'note']].stack().values.tolist())

                if not hits[hits['classification_type'] == 'potential_confer_amr'].empty:
                    sub_hits = hits[hits['classification_type'] == 'potential_confer_amr'].copy()
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

            if self.pima_data.organism_amr_appendices is not None:
                self.appendices.extend(
                    [appendix for drug in drugs for appendix in self.pima_data.organism_amr_appendices if drug in appendix]
                )
                note = (
                    "See the Appendix at the end of the report for more information about putative AMR conferring mutations. "
                    "This only concerns SNP/INDEL variants and not AMR conferring genes which will be reported in the 'Features' "
                    "section of the report"
                )
            self.doc.new_paragraph(note)

        method = self.pima_data.mutations_read_type + ' reads were mapped to the reference sequence using minimap2 (v ' \
                 + self.pima_data.versions['minimap2'] + ').'
        self.methods[self.mutation_methods_title] = pd.concat([self.methods[self.mutation_methods_title], pd.Series(method,dtype="object")])

        method = ' '.join(['Mutations were identified using '
                           'samtools mpileup (v', self.pima_data.versions['samtools'], ')',
                           'and varscan (v', self.pima_data.versions['varscan'], ').'])
        self.methods[self.mutation_methods_title] = pd.concat([self.methods[self.mutation_methods_title], pd.Series(method,dtype="object")])

    def add_appendix_title(self):
        appendix_title = "Appendices"        
        if len(self.appendices) == 0:
            return

        else:
            self.doc.new_line('<div style="page-break-after: always;"></div>')
            self.doc.new_line()
            self.doc.new_header(2, appendix_title)

    def add_amr_matrix(self):
        amr_matrix_title = 'Summary of AMR Features Detected'
        if not getattr(self.pima_data, 'did_draw_amr_matrix', False):
            return

        #amr_matrix_png = self.pima_data.amr_matrix_png
        amr_matrix = self.pima_data.amr_matrix_fig
        self.doc.write(
                f"\n![{amr_matrix_title}]({amr_matrix})\n"
        )        
        # self.doc.new_line()
        # self.doc.new_header(level=2,title=amr_matrix_title)
        # self.doc.write(
        #     self.doc.new_inline_image(
        #         text='All AMR genes and variants detected',
        #         path=amr_matrix_png
        #     )
        #     ,wrap_width=0
        # )
        method = (
            f"Detected genes, SNPs or INDELs compared to the provided reference within regions specified by the mutations_regions.bed file were reported"
            f"and then visualized in a heatmap using matplotlib v({self.pima_data.versions['matplotlib']})."
        )
        self.methods[self.mutation_methods_title] = pd.concat([self.methods[self.mutation_methods_title], pd.Series(method,dtype="object")])


    def add_large_indels(self):
        large_indel_title = 'Large insertions & deletions'
        large_indel_notes_title = 'Large insertions & deletions notes'
        if len(self.pima_data.large_indels) == 0:
            return
        
        large_indels = self.pima_data.large_indels
        self.doc.new_line()
        self.doc.new_header(level=2,title=large_indel_title)

        if len(self.pima_data.large_indel_notes) > 0:
            self.doc.new_header(level=3, title=large_indel_notes_title, add_table_of_contents='n')
            for note in self.pima_data.large_indel_notes:
                self.doc.new_line(note)

        for genome in ['Reference insertions', 'Query insertions']:

            genome_indels = large_indels[genome].copy()
            self.doc.new_line()
            self.doc.new_header(level=3,title=genome, add_table_of_contents='n')
            if genome == 'Reference insertions':
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
                 'regions using bedtools (v ' + self.pima_data.versions['bedtools'] + ').'
        self.methods[self.reference_methods_title] = pd.concat([self.methods[self.reference_methods_title], 
            pd.Series(method,dtype='object')])
        
        self.doc.new_line()

    def add_plasmids(self):
        plasmid_title = "Plasmid annotation"
        plasmid_notes_title = 'Plasmid annotation notes'
        
        if not getattr(self.pima_data, 'did_call_plasmids', False):
            return

        plasmids = self.pima_data.plasmids

        # If contigs shared similarity with database
        if plasmids is not None:

            plasmids = plasmids.copy()
            self.doc.new_line()
            self.doc.new_header(level=2,title=plasmid_title)

            plasmids['query.size'] = plasmids['query.size'].apply(lambda x: '{:,}'.format(x))
            plasmids['aligned.bases'] = plasmids['aligned.bases'].apply(lambda x: '{:,}'.format(x))
            plasmids['plasmid.size'] = plasmids['plasmid.size'].apply(lambda x: '{:,}'.format(x))           
            
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

        #Plasmid Annotation Notes
        if len(self.pima_data.plasmid_notes) > 0:
            self.doc.new_header(level=3, title=plasmid_notes_title, add_table_of_contents='n')
            for note in self.pima_data.plasmid_notes:
                self.doc.new_line(note)

        #Shared methods
        method = ' '.join(['The plasmid reference database was queried against the genome assembly using minimap2 (v',
                           self.pima_data.versions['minimap2'], ').'])
        self.methods[self.plasmid_methods_title] = pd.concat([self.methods[self.plasmid_methods_title],pd.Series(method,dtype="object")])
        """ self.methods[self.plasmid_methods_title] = self.methods[self.plasmid_methods_title].append(pd.Series(method)) """

        method = 'The resulting SAM was converted to a PSL using a custom version of sam2psl.'
        self.methods[self.plasmid_methods_title] = pd.concat([self.methods[self.plasmid_methods_title],pd.Series(method,dtype="object")])
        """ self.methods[self.plasmid_methods_title] = self.methods[self.plasmid_methods_title].append(pd.Series(method)) """

        method = 'Plasmid-to-genome hits were resolved using the pChunks algorithm.'
        self.methods[self.plasmid_methods_title] = pd.concat([self.methods[self.plasmid_methods_title],pd.Series(method,dtype="object")])
        """ self.methods[self.plasmid_methods_title] = self.methods[self.plasmid_methods_title].append(pd.Series(method)) """

    def add_methods(self):

        self.doc.new_line('<div style="page-break-after: always;"></div>')
        self.doc.new_line()

        if len(self.methods) == 0:
            return

        self.doc.new_line()
        self.doc.new_header(level=1, title=self.methods_title)


        for methods_section in self.methods.index.tolist():
            if len(self.methods[methods_section]) == 0:
                continue
            self.doc.new_line()
            self.doc.new_header(level=2,title=methods_section)
            self.doc.new_paragraph(' '.join(self.methods[methods_section]))
            
    def make_tex(self):
        """Populate the ToC and create the md document"""
        self.doc.new_table_of_contents(table_title='Report Contents', depth=1, marker="TableOfContents")
        text = self.doc.file_data_text
        text = text.replace("##--[","")
        text = text.replace("]--##","")
        self.doc.file_data_text = text
        self.doc.create_md_file()

    def main(self):
        """Build the full report"""
        self.start_doc()
        self.add_run_information()
        self.add_run_metrics()
        self.add_feature_detection()
        self.add_reference_comparison()
        self.add_methods()
        self.add_appendix_title()
        self.make_tex()


def validate_make_report(pima_data: PimaData, settings: Settings):
    if pima_data.no_report:
        return

    print_and_log(
        pima_data,
        "Validating reporting utilities",
        pima_data.main_process_verbosity,
        pima_data.main_process_color,
    )

    validate_utility(
        pima_data, "pandoc", "pandoc is not on the PATH (required for reporting)."
    )
    
    pima_data.analysis.append(["make_report", pima_data, settings])

def make_report(pima_data: PimaData, settings: Settings):

    print_and_log(
        pima_data,
        "Making report", 
        pima_data.main_process_verbosity, 
        pima_data.main_process_color,
    )
    pima_data.report_dir = os.path.join(pima_data.output_dir, 'report')

    if find_checkpoint(pima_data, pima_data.report_dir):
        ##Always regenerate report in case downstream steps have changed the results
        shutil.rmtree(pima_data.report_dir) 
    os.mkdir(pima_data.report_dir)

    pima_data.report_prefix = os.path.join(pima_data.report_dir, 'report')
    pima_data.report_md = f'{pima_data.report_prefix}.md'

    pima_data.markdown_report = PimaReport(pima_data, settings)
    pima_data.markdown_report.main()

    ## Add appendices to the pima report
    if len(pima_data.markdown_report.appendices) > 0:
        #assign each new AMR class its own appendix ID (A-Z)
        for i, appendix in zip(string.ascii_uppercase, pima_data.markdown_report.appendices):
            png = re.sub(r"\.md", ".png", appendix)
            mod_appendix = os.path.join(pima_data.report_dir, os.path.basename(appendix))
            shutil.copyfile(png, os.path.join(pima_data.report_dir, os.path.basename(png)))
            with open(appendix, "rt") as fin:
                with open(mod_appendix, "wt") as fout:
                    for line in fin:
                        fout.write(line.replace("LETTER", i))

            # translate the appendix into the markdown report
            with open(pima_data.report_md, "a") as file:
                with open(mod_appendix, "r") as temp_file:
                    file.write(temp_file.read())
            #pima_data.files_to_clean.append(mod_appendix)
            #pima_data.files_to_clean.append(os.path.join(pima_data.report_dir, os.path.basename(png)))

    pima_data.report_pdf = f'{pima_data.report_prefix}.pdf'
    validate_file_and_size_or_error(pima_data, pima_data.report_md, 'Report MD', 'cannot be found', 'is empty')
    
    pandoc_stdout, pandoc_stderr = std_files(os.path.join(pima_data.report_dir, 'markdown2pdf'))
    command = ' '.join(
        [
            'pandoc -f gfm',
            pima_data.report_md,
            '-o', pima_data.report_pdf,
            '--pdf-engine=weasyprint',
            '--css ' + settings.pima_css,
            '1>', pandoc_stdout, 
            '2>', pandoc_stderr,
        ]
    )
    print_and_run(pima_data, command, change_exe_dir=pima_data.report_dir)
    validate_file_and_size_or_error(pima_data, pima_data.report_pdf, 'Report PDF', 'cannot be found', 'is empty')