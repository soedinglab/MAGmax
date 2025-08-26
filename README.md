# MAGmax
MAGmax is a dereplication tool designed to maximize the recovery of Metagenome-Assembled Genomes (MAGs) through bin Merging and reAssembly. It performs dereplication in three stages: (i) grouping bins based on average sequence identity, (ii) merging bins within each group, and (iii) reassembling the merged bins.
![Maxmax](https://github.com/user-attachments/assets/d5f01672-d894-4554-9329-ca5f324395f8)

## Inputs
MAGmax requires three input directories,
1. `<binsdir>`, directory containing bin files in FASTA format that need to be dereplicated. (e.g., output files from any metagenome binning tool)

2. `<readdir>`, directory containing read files in FASTQ format for each sample. 
   
3. `<mapid_dir>`, directory containing mapping files for each sample. Each file is a text file listing read IDs and the corresponding contig IDs they mapped to. These files are used to retrieve reads that map to each merged bin from the FASTQ files in `readdir` and to generate new bin-specific FASTQ files for reassembly.

## Outputs
An output directory named `mags_<x>comp_<y>purity` will be created, where `x` and `y` correspond to the user-specified completeness and purity thresholds used to select final bins. By default, MAGmax uses a percentage of 50 for completeness and 95 for purity.

The output directory contains dereplicated bins, and a text file listing the completeness and contamination scores for each bin as calculated by CheckM2.

## Example command line call

    magmax -b <binsdir> -r <readdir> -m <mapid_dir> -f fasta -t 24
    magmax -b <binsdir> -r <readdir> -m <mapid_dir> -f fasta -t 24 -q quality_report.tsv // if CheckM2 result is already available
    magmax -b <binsdir> -r <readdir> -m <mapid_dir> -f fasta -t 24 --split // if input bins are not already split by sample id 

## Installation
### Prerequisites

- **Rust**: Follow the instructions [here](https://www.rust-lang.org/tools/install) to install Rust.
- **Conda**: You can install Conda via [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/distribution).

### Dependencies

- **CheckM2**: Install [CheckM2](https://github.com/chklovski/CheckM2), download [checkm2 database](https://zenodo.org/api/records/14897628/files/checkm2_database.tar.gz/content) and set CHECKM2DB variable correctly. CheckM2 should already be installed and accessible in your PATH, regardless of the options used to install MAGmax.

Option 1: Use conda package

    conda install -c bioconda magmax
    or
    mamba install -c bioconda magmax # faster installation

Option 2: Use the pre-built executable.

    # For x86_64 Linux (glibc-based systems)
    wget https://github.com/soedinglab/MAGma/releases/download/v1.0.0/magmax-linux.tar.gz
    cd magmax-linux/bin
    chmod +x magmax
    ./magmax -h
    sudo cp magmax /usr/local/bin/ # to access globally

To use this option, in addition to CheckM2, [skani](https://github.com/bluenote-1577/skani), [SPAdes](https://github.com/ablab/spades), and [seqtk](https://github.com/lh3/seqtk), and [MEGAHIT](https://github.com/voutcn/megahit) (optional) must be installed already and available in your PATH. Alternatively, use environment.yml to create conda environment and activate it to run magmax.

    conda env create -f environment.yml
    conda activate magmax_env

Option 2: Build from source

    git clone https://github.com/soedinglab/MAGmax.git
    cd MAGmax
    conda env create -f environment.yml
    conda activate magmax_env
    cargo install --path .
    magmax -h

    
## Options
        -b, --bindir <BINDIR>
                Directory containing fasta files of bins
        -i, --ani <ANI>
                ANI for clustering bins (%) [default: 99]
        -c, --completeness <COMPLETENESS_CUTOFF>
                Minimum completeness of bins (%) [default: 50]
        -p, --purity <PURITY_CUTOFF>
                Mininum purity (1- contamination) of bins (%) [default: 95]
        -m, --mapdir <MAPDIR>
                Directory containing mapids files
        -r, --readdir <READDIR>
                Directory containing read files
        -f, --format <FORMAT>
                Bin file extension [default: fasta]
        -t, --threads <THREADS>
                Number of threads to use [default: 8]
            --split
                Split clusters into sample-wise bins before processing
        -q, --qual <QUAL>
                Quality file produced by CheckM2 (quality_report.tsv)
            --assembler <ASSEMBLER>
                assembler choice for reassembly step (spades|megahit) [default: spades, recommended]
        -h, --help
                Print help
        -V, --version
                Print version

## Test run using toy data
This example test run demonstrates dereplication of bins using the provided toy dataset. In the `test/bins` directory, example bins generated with MetaBAT2 are given. In the `test/reads` directory, paired-end read files for two samples are given and in the `test/mapids` directory, mapid files mapping reads to contigs for each sample are given. Precomputed CheckM2 quality scores for the input bins are given in the `test/quality_report.tsv`. Run the following command to execute the test:

    magmax -b test/bins -r test/reads -m test/mapids -t 24 -q test/quality_report.tsv
After running MAGmax, an output folder named `mags_50comp_95purity` will be created in the `test` directory. This folder contains the following files:

- `bins_checkm2_qualities.tsv` — Table summarizing the quality metrics of the dereplicated bins.  
- `sample_ERR3405607_metabat2_results.63.fasta` — Final bin obtained after dereplication of the input bins.

## Input specifications
1. Input contigs must have IDs prefixed with the sample ID, separated by a `C`. This is a common practice for both single- and multi-sample binning. Ensure mapping and binning are performed on contig files with these updated contig IDs.

2. Ensure that headers in the fastq files have read ID separated from sequencer details by a space or tab, not by a dot. This is important for `seqtk`, which is used by MAGmax, to fetch reads correctly.

    `Correct format: @SRR25448374.1 A00214R:157:HLMVMDSXY:1:1101:19868:1016:N:0.length=151#0/1`

    `Wrong format: @SRR25448374.1.A00214R:157:HLMVMDSXY:1:1101:19868:1016:N:0.length=151#0/1`

To fix, use the below bash command

    sed -i -E 's/^(@[^.]+\.[^.]+)\./\1 /' read.fastq

3. Mapid files can be created using `aligner2counts` (https://github.com/soedinglab/binning_benchmarking/tree/main/util#aligner2counts) with the `only-mapids` option. An example file format is given below,

    File name: `<sampleid>_mapids`
    ```
    read1_id    <sampleid>Ccontig1_id
    read2_id    <sampleid>Ccontig2_id
    read2_id    <sampleid>Ccontig4_id
    read3_id    <sampleid>Ccontig2_id
    read4_id    <sampleid>Ccontig3_id
    read4_id    <sampleid>Ccontig4_id
    ```

4. FASTQ and MAPID filenames must contain the sample ID (e.g., SRR25448374.fastq, SRR25448374_mapids). This is the default unless filenames are renamed manually.

## Notes
1. If input bins are not separated by sample IDs (e.g., when using MetaBAT2 or COMEBin on concatenated contigs), use the `--split` option to let MAGmax automatically separate bins by sample ID.

2. We recommend Spades for reassembly which produces bins with higher purity than bins assembled using Megahit.
