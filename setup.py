import setuptools
import os
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open(os.path.join("Pima", "VERSION"), "r", encoding="utf-8") as version_fp:
    VERSION = version_fp.read().strip()

setuptools.setup(
    name="pima",
    version=VERSION,
    author="Applied Bioinformatics Laboratory",
    author_email="woverholt@asrtinc.com",
    description="Genomic characterization pipeline for Bacillus anthracis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/appliedbinf/pima",
    packages=setuptools.find_namespace_packages(),
    python_requires='>=3.8',
    package_data={
        "Pima": [
            "data/reference_sequences/**",
            "data/ames_single_copy_genes.fna",
            "data/amr.fasta",
            "data/ba_virulence_genes.fasta",
            "data/example_sample_sheet.tsv",
            "data/gene_drug.tsv",
            "data/inc.fasta",
            "data/pima.css",
            "VERSION",
            "nextflow_parallelization/**",
            "accessory_scripts/pChunks.R",
            "accessory_scripts/sam2psl.py",
        ],
    },
    scripts = ['Pima/pima.py', 'Pima/accessory_scripts/building_pycircos_figures.py'],
    zip_safe=False,
    include_package_data=True,
    entry_points={
        # Optional: specify any entry points for your package here
        'console_scripts': [
            'pima = Pima.pima:main',
        ],
    },
)