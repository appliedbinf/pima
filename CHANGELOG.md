# Changelog

All notable changes to this project from version 2.2.0 onwards are documented in this file.

> **Note:** Releases prior to 2.2.0 were tracked via internal `VERSION` bumps and are not documented here.

## Added
- Added check to ensure the provided genome is a fasta file.
- Implemented steps to dereplicate overlapping results from resfinder and amrfinder.
- Updated plasmidfinder and resfinder databases.
- Enabled multiplex mode without specifying genome size.
- Enhanced visualizations for AMR features.
- Added AMRFinder
- Added support for serial samplesheet mode.
- Added support for serial multiplex mode.

## Changed
- Updated the resfinder phenotype file.
- Modified pandas commands in the report to avoid warnings.
- Upgraded environments to use the latest medaka models.
- Revised the conda recipe.
- Updated the README.
- Restructured the pima_data class to remove redundant variables.
- Improved report creation to enhance PDF resolution.
- Ensured Python3 compatibility.

## Fixed
- Resolved the slicing bug in pChunks.
- Addressed the issue in multiplex mode where an empty illumina_fastq list triggered the wrong circos mode.
- Corrected report header formatting.
- Fixed the logic for plotting the AMR matrix.
- Resolved issues with plasmid reporting.
- Fixed the barcode_min_fraction bug.
- Cleaned up after kraken2.
- Corrected errors in passing illumina data to nextflow.
- Addressed issues in pChunks when no plasmids were detected and updated plasmid reporting.
- Ensured pima multiplex mode correctly handles cases where genome size is not specified.
- Resolved bugs in samplesheet and multiplexed mode.
- Addressed bugs related to plasmid reporting, report paths, and options in serial multiplex mode.
- Corrected issues in multiplex mode where it didn't skip functions if ONT data wasn't provided.
- Fixed issues in pChunks when no plasmids were detected and updated plasmid reporting.
- Addressed the nextflow config syntax to stop max-memory warnings.
- Resolved various other bugs related to plasmid reporting, report paths, and options in serial multiplex mode.

## Chores/Maintenance
- Incremented version numbers multiple times to reflect various changes and fixes.
- Updated build numbers to reflect incremental changes.
- Ensured compatibility with Python3.
- Cleaned up and removed unused CLI options.
- Updated conda environments and notes.
- Simplified CSS, updated database paths, and fixed issues with Nextflow and Illumina paths.
- Performed bug squashing and clarified data/settings classes.
- Updated the conda recipe.
- Removed all Docker-specific paths that are no longer used.
- Squashed various bugs and fixed plasmid reporting.
- Fixed report paths and options when run in serial multiplex mode.
- Removed the CLI option for a reference directory, allowing PiMA to manage it internally.
- Fixed bugs with serial samplesheet mode related to pima_data attributes not updating correctly.
- Fixed bugs in pima multiplex without specifying genome-size.
- Consolidated the 2 report scripts into the 'report.py' copy and removed accessory scripts.