# metat
The pipeline for metatranscriptomics analysis

### Prerequirements

To reproduce the output, you need to use `Bioconda`.

Please follow the instruction [here](https://bioconda.github.io) to install `Bioconda`. 
And then you need to install `snakemake` and Python package `click` and `pandas`:

```shell
conda install snakemake=5.5.4
conda install Click=7.0
conda install pandas=0.25.0
```


After this has been done, download the pipeline onto your system:

```shell
git clone git@github.com:dawnmy/metat.git
```

### Modify the config file: `config/config.yaml`
All the paths must be either relative path to the parent directory of `config` folder or absolute path.

```yaml
dataset: mouse # name for the dataset
fq_dir: ../data/seq # dir of the raw FASTQ files
out_dir: ../outputs # dir to put the results
paired: true # is paried end reads?
suffix: # the suffixs of reads; the suffix is the comman suffix for all samples besides the sample name
  - _R1.fastq.gz # please keep the hyphen sign
  - _R2.fastq.gz
# host_ref: ../ref/mouse.fa
ref: ../ref/mouse_gut_gene_catalog.fa # The gene catalog for quantifying the expression
threads: 20
```


### Run the pipelines

#### Get the expression table for genes

```shell
snakemake -s metat.smk -j 20 --use-conda
```
`-s` to specify the pipeline file, and `-j` to set the number of threads to use and `--use-conda` to 
let the pipeline install required softwares with specified version. The `conda` ENVs will be created under 
the path of the program by default. The program may take ten minutes to create the ENV for the first time.
If you do not wish to create the conda ENV in the working directory, 
please use --conda-prefix parameter to specify the desired path to create the `conda` ENV. 


**If you use SGE for the job submission, you can use the following cmd:**

```shell
snakemake -s metat.smk --latency-wait 30 --use-conda -c "qsub -cwd -q <the job submission queue> \
 -pe multislot {threads} -i /dev/null -e <dir for std error logs> -o <dir for std output logs> \
 -v PATH" -j 2
```

### The output structure

```
outputs
└── mouse
    ├── data
    │   ├── bam
    │   └── qc_fq
    │       |── mrna
    │       └── rrna
    ├── reports
    │   ├── benchmarks
    │   ├── bwa
    │   ├── fastp
    │   └── samtools
    └── results
        └── count
```