import os
from pathlib import Path

class Settings():
    def __init__(self):
        self.data_dir = os.path.join(self.pima_path, 'data')
        self.amr_database = os.path.join(self.pima_path, 'data/amr.fasta')
        self.amr_gene_drug_tsv = os.path.join(self.pima_path, 'data/gene_drug.tsv')
        self.inc_database = os.path.join(self.pima_path, 'data/inc.fasta')
        self.included_databases = [self.amr_database, self.inc_database]
        self.ba_virulence_genes = os.path.join(self.pima_path, 'data/ba_virulence_genes.fasta')
        self.plasmid_database_fasta = os.path.join(self.pima_path, 'data/plasmids_and_vectors.fasta')
        self.kraken_database = os.path.join(self.pima_path, 'data/kraken2')
        self.amrfinder_database = os.path.join(self.pima_path, "data/amrfinder_db")
        self.reference_dir = os.path.join(self.pima_path, 'data/reference_sequences')
        self.pima_css = os.path.join(self.pima_path,'data/pima.css')
        #draw on Ba plasmids
        self.virulence_genes_fp = os.path.join(self.data_dir, "reference_sequences/Bacillus_anthracis/ba_virulence_genes.bed")
       
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