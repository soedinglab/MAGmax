# MAGmax
MAGmax is a tool to maximize the yield of Metagenome-Assembled Genomes (MAGs) through bin Merging and reAssembly.

### Example run

    magmax -b binsdir -m mapid_dir -r readdir -f fasta -t 24
    magmax -b binsdir -m mapid_dir -r readdir -f fasta -t 24 -q quality_report.tsv // if CheckM2 result is already available
    magmax -b binsdir -m mapid_dir -r readdir -f fasta -t 24 --split // if input bins are not already split by sample id 

### Test run
    magmax -b test/bins -m test/mapids -r test/reads -t 24 -q test/quality_report.tsv

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


## Notes
1. Input contigs should have id prefixed with the sample ID, separated by 'C'. Perform mapping and binning on contig files with these updated contig ids.
2. Mapid files can be generated using aligner2counts (https://github.com/soedinglab/binning_benchmarking/tree/main/util#aligner2counts) with `only-mapids` option.

    File name: `<sampleid>_mapids`
    ```
    read1_id    sampleidCcontig1_id
    read2_id    sampleidCcontig2_id
    read2_id    sampleidCcontig4_id
    read3_id    sampleidCcontig2_id
    read4_id    sampleidCcontig3_id
    read4_id    sampleidCcontig4_id
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
