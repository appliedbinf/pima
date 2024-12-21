import os
from pathlib import Path

class Settings():
    def __init__(self):
        self.data_dir = os.path.join(self.pima_path, 'data')
        self.amr_database_default = os.path.join(self.pima_path, 'data/amr.fasta')
        self.amr_gene_drug_tsv = os.path.join(self.pima_path, 'data/gene_drug.tsv')
        self.amr_default_color = '#FED976'
        self.inc_database_default = os.path.join(self.pima_path, 'data/inc.fasta')
        self.inc_default_color = '#0570B0'
        self.included_databases = [self.amr_database_default, self.inc_database_default]

        self.plasmid_database_default_fasta = os.path.join(self.pima_path, 'data/plasmids_and_vectors.fasta')
        self.kraken_database_default = os.path.join(self.pima_path, 'data/kraken2')
        self.reference_dir_default = os.path.join(self.pima_path, 'data/reference_sequences')
        self.pima_css = os.path.join(self.pima_path,'data/pima.css')
        self.virulence_genes_fp = os.path.join(self.data_dir, "reference_sequences/Bacillus_anthracis/ba_virulence_genes.bed")

        ## Docker specific paths
        self.DockerPathPlasmid = os.path.join('/home/DockerDir/Data/Temp_Data/plasmids_and_vectors.fasta')
        self.DockerPathKraken = os.path.join('/home/DockerDir/Data/Temp_Data/kraken2')
        
    @property
    def pima_path(self):
        # Is __name__ the most robust way to print the path of importing scripts, not this one?
        #return os.path.dirname(os.path.realpath(__name__))
        return Path(__file__).parent.parent
    @property
    def pima_version(self):
        with open(os.path.join(self.pima_path, "VERSION"), "r") as version_fp:
            VERSION = version_fp.read().strip()
        return VERSION