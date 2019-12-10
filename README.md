# **pima**
Plasmid and antimicrobial resistance pipeline

### Example Installation
Pima is not yet packaged, but its primary dependencies are _mostly_ available. `guppy` must be manually installed with info [here](https://nanoporetech.com/community). If you don't yet have `conda` installed you'll need to do that [here](https://docs.anaconda.com/anaconda/install/). Once installed you can create an environment where `pima.py` will run.
```sh
conda create -n pima_env2 python=2.7.16
source activate pima_env2
# install python modules
conda install -c conda-forge joblib numpy pandas pathos -y
# install bioinfx packages
conda install -c bioconda bcftools bedtools blast flye mummer medaka miniasm minimap2 nanopolish parallel porechop qcat racon samtools spades wtdbg -y
git clone https://github.com/abconley/pima.git $HOME
python ~/pima/pima.py --help
```

##### Compatibility
Currently several unix utilities are used, so this is not compatible with Windows OS. CentOS and Ubuntu are both known working OSes.



##### Dependencies
| System Command | Package |
| -------------- | ------- |
| `bcftools` | BCFtools |
| `bedtools` | BEDTools |
| `blastn` | BLAST+ |
| `dnadiff` | MUMmer |
| `faidx` | SAMtools |
| `flye` | Flye |
| `parallel` | GNU parallel |
| `guppy_basecaller` | Guppy |
| `makeblastdb` | BLAST+ |
| `medaka_consensus` | Medaka |
| `miniasm` | miniasm |
| `minimap2` | minimap2 |
| `nanopolish` | Nanopolish |
| `nanopolish_makerange.py` | Nanopolish |
| `pblat` | pblat |
| `pChunks.R` | local R script |
| `porechop` | Porechop |
| `qcat` | Qcat |
| `racon` | racon |
| `read_fast5_basecaller.py` | Guppy |
| `samtools` | SAMtools |
| `spades.py` | SPAdes |
| `wtdbg2` | Wtdbg2 |
| `wtpoa-cns` | Wtdbg2 |
|  |  |
| `awk` | Unix |
| `cat` | Unix |
| `grep` | Unix |
| `head` | Unix |
| `ls` | Unix |
| `sed` | Unix |
| `sort` | Unix |

#### Literature Citations
- McLaughlin HP, Bugrysheva JV, Conley AB, Gulvik CA, Kolton CB, Marston C, Swaney E, Lonsway DR, Cherney B, Gargis AS, Kongphet-Tran T, Lascols C, Michel P, Villanueva J, Hoffmaster ER, Gee JE, Sue D. 2020. When minutes matter: rapid nanopore whole genome sequencing for anthrax emergency preparedness. Emerging Infectious Diseases
- [BCFtools and SAMtools](https://www.ncbi.nlm.nih.gov/pubmed/19505943)
- [BEDTools](https://www.ncbi.nlm.nih.gov/pubmed/25199790)
- [BLAST+](https://www.ncbi.nlm.nih.gov/pubmed/20003500)
- [GNU parallel](https://www.usenix.org/publications/login/february-2011-volume-36-number-1/gnu-parallel-command-line-power-tool)
- [Flye](https://www.ncbi.nlm.nih.gov/pubmed/30936562)
- [miniasm and minimap2](https://www.ncbi.nlm.nih.gov/pubmed/27153593)
- [MUMmer](https://www.ncbi.nlm.nih.gov/pubmed/14759262)
- [pblat](https://www.ncbi.nlm.nih.gov/pubmed/30646844)
- [racon](https://www.ncbi.nlm.nih.gov/pubmed/28100585)
- [SPAdes](https://www.ncbi.nlm.nih.gov/pubmed/22506599)
- [Wtdbg2](https://www.nature.com/articles/s41592-019-0669-3)
