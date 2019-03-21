# 2019 CasX model organism site catalog
Repository for scripts associated with manuscript.

Associated paper:

## Requirements
### General
* Linux + bash
* R
* R Studio
* Python

### Additional command-line
* make
* wget
* pigz (parallel implementation of gzip)

## Running analysis

### 1 - install R packages
* tidyverse - CRAN
* here - CRAN
* knitr - CRAN
* devtools - CRAN
* multidplyr - GitHub via devtools (devtools::install_github("hadley/multidplyr"))
* GenomicRanges - BioC
* GenomicFeatures - BioC

### 2 - install python requirements
Recommended that you do this in a python virtual environment

```bash
pip install motif-scraper
```

```bash
pip install pandas
```

```bash
pip install numpy
```

```bash
pip install scipy
```

### 3 - clone repository
```bash
git clone https://github.com/RobersonLab/2019CasXModelOrgCatalog.git
```

### 4 - change needed parameters
File: changeable_parameters/num_motif_cpus.txt
Number of threads used to run site identification with motif scraper and compress certain files with pigz.

File: changeable_parameters/num_r_cpus.txt
Number of CPUs used by multidplyr during annotation.

### 5 - Choose your adventure!
[Path A - Run only with Make.](#identify-and-annotate-with-makefile)
[Path B - Run initial identification with Make. Annotate on a PBS / Torque HPC cluster.](#identify-with-make-and-annotate-with-hpc)

#### Identify and annotate with makefile
```bash
nohup make -f casx_run_analysis.make 1>run.log 2>&1 &
```

#### Identify with Make and annotate with HPC
* In "all:" section of Makefile remove $(ROUTPUT) $(MOUSE_TARGETS) $(HAMMING).
* In PBS files (mplyr_run_r, pull_mouse_run, run_mouse_hamming), make sure PBS declarations at the top of the scripts are compatible with the equipment on your cluster.
* Also in run_mouse_hamming.pbs: a conda virtual environment is loaded. Edit / remove as needed to match your environment. This should be the same enviroment the Python dependencies were installed.

```bash
nohup make -f casx_run_analysis.make 1>run.log 2>&1 & 
```

After Makefile has run successfully:
```bash
ANNOTATEID=$(qsub mplyr_run_r.pbs)
```

```bash
MOUSEID=$(qsub -W depend=afterok:$ANNOTATEID pull_mouse_run.pbs)
```

```bash
qsub -W depend=afterok:$MOUSEID run_mouse_hamming.pbs
```

### 6 - knit markdowns with R Studio to generate figures
For the markdowns to work properly the site identification, site annotation, mouse site extraction, and hamming distance calculations must have completed successfully.

Knit the general analysis to generate the PAM figure and associated tables. Knit the mouse cor matrix markdown to generate the mouse strain Hamming figure.

## Data associated with the paper
[CasX site catalog on FigShare](https://figshare.com/projects/2019_CasX_genome_editing_site_annotations/61103)

