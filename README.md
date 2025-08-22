# MAGmax
MAGmax is a dereplication tool designed to maximize the recovery of Metagenome-Assembled Genomes (MAGs) through bin Merging and reAssembly. It performs dereplication in three stages: (i) grouping bins based on average sequence identity, (ii) merging bins within each group, and (iii) reassembling the merged bins.

## INPUTS
MAGmax requires three input directories,
1. `binsdir`, a directory containing bin files in FASTA format that need to be dereplicated. (e.g., output files from any metagenome binning tool)

2. `readdir`, a directory containing read files in FASTQ format for each sample. 
   
3. `mapid_dir`, a directory containing mapping files for each sample. Each file is a text file listing read IDs and the corresponding contig IDs they mapped to. These files are used to retrieve reads that map to each merged bin from the FASTQ files in `readdir` and to generate new bin-specific FASTQ files for reassembly.

## OUTPUTS
An output directory named `mags_<x>comp_<y>purity` will be created, where `x` and `y` correspond to the user-specified completeness and purity thresholds used to select final bins. By default, MAGmax uses a percentage of 50 for completeness and 95 for purity.   
The output directory contains dereplicated bins, and a text file listing the completeness and contamination scores for each bin as calculated by CheckM2.

### Example command line call

    magmax -b <binsdir> -r <readdir> -m <mapid_dir> -f fasta -t 24
    magmax -b <binsdir> -r <readdir> -m <mapid_dir> -f fasta -t 24 -q quality_report.tsv // if CheckM2 result is already available
    magmax -b <binsdir> -r <readdir> -m <mapid_dir> -f fasta -t 24 --split // if input bins are not already split by sample id 

## Install
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
    cd MAGma
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

### Test run using toy data
This example test run demonstrates dereplication of bins using the provided toy dataset. In the `test/bins` directory, example bins generated with MetaBAT2 are given. In the `test/reads` directory, paired-end read files for two samples are given and in the `test/mapids` directory, mapid files mapping reads to contigs for each sample are given. Precomputed CheckM2 quality scores for the input bins are given in the `test/quality_report.tsv`. Run the following command to execute the test:

    magmax -b test/bins -r test/reads -m test/mapids -t 24 -q test/quality_report.tsv


## Notes
1. Input contigs should have id prefixed with the sample ID, separated by 'C', as commonly practiced in the single-sample and multi-sample binning. Perform mapping and binning on contig files with these updated contig ids.
2. Mapid files can be generated using aligner2counts (https://github.com/soedinglab/binning_benchmarking/tree/main/util#aligner2counts) with `only-mapids` option.

    File name: `<sampleid>_mapids`
    ```
    read1_id    <sampleid>Ccontig1_id
    read2_id    <sampleid>Ccontig2_id
    read2_id    <sampleid>Ccontig4_id
    read3_id    <sampleid>Ccontig2_id
    read4_id    <sampleid>Ccontig3_id
    read4_id    <sampleid>Ccontig4_id
    ```

3. If input bins are not separated by sample IDs, such as when using MetaBAT2 or COMEBin on a concatenated set of contigs, use the `--split` option to automatically separate input bin by sample IDs.
4. Make sure that headers in the read fastq files have read_id separated by space/tab (not by `.`) from other sequencer details. This is important for `seqtk` to fetch reads correctly.

    `Correct format: @SRR25448374.1 A00214R:157:HLMVMDSXY:1:1101:19868:1016:N:0.length=151#0/1`

    `Wrong format: @SRR25448374.1.A00214R:157:HLMVMDSXY:1:1101:19868:1016:N:0.length=151#0/1`

When read ids are not seperated by space in the headers, run the below script and use the updated read file for mapping.
 
    sed -i -E 's/^(@[^.]+\.[^.]+)\./\1 /' read.fastq

MAGma works for paired-end (in separate files: SRR25448374_1.fastq and SRR25448374_2.fastq) and single-end read files.

5. Sample IDs must be in the file name of fastq and mapid files. (E.g., SRR25448374_1.fastq & SRR25448374_2.fastq or SRR25448374.fastq and SRR25448374_mapids)
6. We recommend Spades for reassembly which produces bins with higher purity than bins assembled using Megahit.
