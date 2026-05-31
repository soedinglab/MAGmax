## Tutorial: generate custom species-level database

Genome dereplication is not always perfect due to inherent limitations of hierarchical clustering algorithms used in dereplication tools (dRep and Galah). Alternatively, taxonomic classification using GTDBtk followed by grouping genomes by taxonomy assignment is another option for dereplication, but it has limitations too: 1) ANI radius of under-represented species may be inaccurate, causing wrong taxonomy labeling; 2) novel species cannot be assigned. Combining dereplication and taxonomic classification can enhance the discovery of novel species with improved accuracy.

---

### Overview

The `magmax customdb` subcommand builds a species-level non-redundant genome database by combining two complementary strategies:

1. **GTDB-Tk-guided dereplication** — bins that are assigned to a known species (ANI ≥ species ANI radius and aligned fraction ≥ cutoff) are grouped by species name and one representative is chosen per species.
2. **ANI-based dereplication of remaining bins** — bins that are unclassified by GTDB-Tk or are not assigned to a species (s__) are clustered by ANI (default 95%), using the same graph-based approach as regular MAGmax dereplication.

The result is a combined set of representatives covering both known and novel species.


### Prerequisites

All tools required for regular MAGmax runs plus:
- **GTDB-Tk** — taxonomic classification of input bins (`gtdbtk classify_wf`).  
  Requires only the summary output file (`gtdbtk.bac120.summary.tsv` or equivalent).


### Typical workflow

#### Step 1 — Run GTDB-Tk on your bins

```bash
gtdbtk classify_wf \
    --genome_dir bins/ \
    --extension fasta \
    --out_dir gtdbtk_output/ \
    --cpus 32
```

The relevant summary file is one of:
- `gtdbtk_output/classify/gtdbtk.bac120.summary.tsv` (bacteria)
- `gtdbtk_output/classify/gtdbtk.ar53.summary.tsv`   (archaea)

For a combined database, concatenate both files (keeping only one header line):

```bash
head -1 gtdbtk.bac120.summary.tsv > gtdbtk.summary.tsv
tail -n +2 gtdbtk.bac120.summary.tsv >> gtdbtk.summary.tsv
tail -n +2 gtdbtk.ar53.summary.tsv  >> gtdbtk.summary.tsv
```

#### Step 2 — Run CheckM2 (optional but recommended)

Providing a pre-computed quality file avoids re-running CheckM2 inside `customdb`.

```bash
checkm2 predict --threads 32 --input bins/ --output-directory checkm2_output/ -x fasta
```

The output file is `checkm2_output/quality_report.tsv`.

#### Step 3 — Run `magmax customdb`

Minimal run (CheckM2 and skani are executed automatically):

```bash
magmax customdb \
    -g gtdbtk.summary.tsv \
    -b bins/ \
    -t 32
```

With pre-computed quality and ANI files to save time on re-runs:

```bash
magmax customdb \
    -g gtdbtk.summary.tsv \
    -b bins/ \
    -q checkm2_output/quality_report.tsv \
    --anifile ani_edges \
    -t 32
```

Including cultivated isolate genomes as priority representatives:

```bash
magmax customdb \
    -g gtdbtk.summary.tsv \
    -b bins/ \
    -q checkm2_output/quality_report.tsv \
    --isolate-genomes isolates.txt \
    -t 32
```

Using sensitive mode for representative selection among unclassified bins:

```bash
magmax customdb \
    -g gtdbtk.summary.tsv \
    -b bins/ \
    -q checkm2_output/quality_report.tsv \
    --sensitive \
    -t 32
```

### Input files

| Input | Flag | Required | Description |
|-------|------|----------|-------------|
| GTDB-Tk summary | `-g` | Yes | Tab-separated GTDB-Tk classification file (`gtdbtk.bac120.summary.tsv` or combined) |
| Bin directory | `-b` | Yes | Directory containing FASTA files of input bins |
| CheckM2 quality | `-q` | No | `quality_report.tsv` from CheckM2; computed automatically if omitted |
| Isolate genome list | `--isolate-genomes` | No | Plain text file, one genome name per line (with or without extension); these are prioritized as representatives |
| Pre-computed ANI | `--anifile` | No | Output of `skani triangle <bindir> -E -o <anifile>`; computed automatically if not given |

#### GTDB-Tk summary columns used

The parser reads the following columns (0-indexed):

| Column | Name in GTDB-Tk output | Used for |
|--------|------------------------|----------|
| 0 | `user_genome` | Bin identifier |
| 1 | `classification` | Full taxonomy string; species extracted from `s__` tag |
| 3 | `ani_radius` | Per-species ANI radius reported by GTDB-Tk |
| 5 | `closest_placement_ani` | ANI to closest reference genome |
| 6 | `closest_placement_af` | Aligned fraction to closest reference |

A bin is classified as **perfect** (confidently assigned to a known species) when:
- `closest_placement_ani` ≥ max(`ani_radius`, `--species-ani`) 
- `closest_placement_af` ≥ `--species-alignedfrac`
- Species field (`s__`) is non-empty

All other quality-passing bins are treated as **remaining** (unclassified or novel species).

#### Isolate genome list format

```
# Lines starting with '#' are ignored
isolate_genome_1           # with or without .fasta extension
path/to/isolate_genome_2   # path prefix is stripped; only basename is used
```


### Output files

Output is written to `specieslevel_customdb/` by default (use `-o` to override), created next to the bin directory.

| File | Description |
|------|-------------|
| `memberships.tsv` | All representatives and their cluster members (GTDB-Tk classified + unclassified). This is the final complete dereplication result. Tab-separated: representative, then a comma-separated member list |
| `bins_checkm2_qualities.tsv` | Completeness and contamination values of all final representatives. Columns: `#Bin`, `Completeness`, `Contamination` |
| `gtdbtk_species_representatives.tsv` | Representatives selected from GTDB-Tk-classified bins. Columns: `#gtdbtk_species_representative`, `species_name` |
| `unclassified_clusterrepresentatives_gtdbtkspecies_ani_connections.tsv` | ANI connections between representatives of novel-clusters and known species clusters that exceed the species ANI radius. Columns: `#unclassified_cluster_representative`, `gtdbtk_species_representative`, `ANI`, `species_ANI_radius`. It informs whether any unclassified representative might actually belong to a known species. |


### Options
    -b, --bindir <BINDIR>
            Directory containing fasta files of bins
    -g, --gtdbtk <GTDBTK>
            GTDB-Tk classification summary file
    -q, --qual <QUAL>
            Quality file produced by CheckM2 (quality_report.tsv)
        --isolate-genomes <ISOLATE_GENOMES>
            File listing isolate genomes in the input bins; these are prioritized as species representatives
        --sensitive
            Select representatives based on high connectivity. Bin merging and reassembly steps are disabled
        --species-ani <SPECIES_ANI>
            ANI for clustering bins (%), as per GTDB-Tk criteria [default: 95]
        --species-alignedfrac <SPECIES_ALIGNEDFRAC>
            Minimum aligned fraction (%) for species-level clustering, as per GTDB-Tk criteria [default: 50]
    -c, --completeness <COMPLETENESS_CUTOFF>
            Minimum completeness of bins (%) [default: 90]
    -p, --purity <PURITY_CUTOFF>
            Purity cutoff for custom database generation (%) [default: 5]
    -t, --threads <THREADS>
            Number of threads to use [default: 8]
        --split
            Split clusters into sample-wise bins before processing
    -f, --format <FORMAT>
            Bin file extension [default: fasta]
        --anifile <ANIFILE>
            ANI file produced by skani using command: skani triangle <bindir> -E -o <anifile>
    -o, --outdir <OUTPUT>
            Directory of output
    -h, --help
            Print help

### How representatives are selected

**GTDB-Tk-classified known species bins**

All bins assigned to the same GTDB-Tk species are grouped together. One representative is chosen per species:
1. If isolate genomes are present among the species members, the isolate with the lowest contamination is preferred.
2. Otherwise, the bin with the highest quality score (`completeness − 5 × contamination`) is selected.

**Unclassified or novel-species bins**

Remaining bins are clustered by pairwise ANI (default 95%, aligned fraction ≥ 50%) using single-linkage (connected components). Within each cluster:
- **Default mode**: selects the highest-quality bin (completeness ≥ 90% required; isolates are prioritized).
- **Sensitive mode** (`--sensitive`): selects the bin with the highest weighted ANI connectivity (`Σ max(0, ANI − threshold)` over neighbors), favoring bins that are more similar to a larger number of neighbors.


### Notes

1. **Quality thresholds for database creation are stricter than regular dereplication.** The defaults are completeness ≥ 90% and contamination ≤ 5%. Adjust with `-c` and `-p` if needed.

2. **Pre-computing ANI saves time on re-runs.** Generate the ANI file once with `skani triangle <bindir> -E -o ani_edges` and pass it via `--anifile` to skip re-computation. If not provided, ANI is computed among unclassified bins and representatives of known species clusters selected in the first step from GTDB-Tk classification. The result file is cached automatically as `<bindir>/subset_ani_edges`, which can be reused in future runs.

3. **The `--split` flag is needed when bins are not already separated by sample ID.** When running multi-sample binning on concatenated contigs (e.g., MetaBAT2 or COMEBin), use `--split` to let MAGmax separate bins by sample before processing. Isolate genomes will not be split by sample ID.

4. **Isolate genome names must match bin filenames.**

5. **The `unclassified_clusterrepresentatives_gtdbtkspecies_ani_connections.tsv` file is a diagnostic resource.** It lists novel-cluster representatives whose ANI to a known GTDB-Tk species representative meets or exceeds that species' ANI radius. This happens when unclassified cluster representatives have lower ANI to the GTDB reference species than the representatives selected from the user's input dataset.


### Building a unified species-level database: integrating MAGmax dereplication results with GTDB reference genomes

The `unifygtdb.sh` script combines magmax customdb output and GTDB reference genomes. This is useful when users wants to create a complete species-level genome reference database including all known species and unknown species covered in the input data.


    bash unifygtdb.sh gtdb_taxonomy.tsv gtdbtk_species_representatives.tsv memberships.tsv > unified.tsv

The script takes three inputs in this specific order:
1. **GTDB taxonomy file** — full GTDB genomes with their taxonomy assignments
   
    Example: `gtdb_taxonomy.tsv`
    RS_GCF_000016525.1	d__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanocatella;s__Methanocatella smithii
    GB_GCA_002686525.1	d__Archaea;p__Thermoplasmatota;c__Poseidoniia;o__Poseidoniales;f__Thalassarchaeaceae;g__MGIIb-O2;s__MGIIb-O2 sp002686525
    RS_GCF_000970205.1	d__Archaea;p__Halobacteriota;c__Methanosarcinia;o__Methanosarcinales;f__Methanosarcinaceae;g__Methanosarcina;s__Methanosarcina mazei

2. **GTDB-Tk species representatives selected from the user data** — the representatives from `magmax customdb` run that were classified by GTDB-Tk
   
    Example: `gtdbtk_species_representatives.tsv`
    #gtdbtk_species_representative	species_name
    38222_26_bin.5	d__Bacteria;p__Bacillota;c__Bacilli;o__Staphylococcales;f__Staphylococcaceae;g__Staphylococcus;s__Staphylococcus haemolyticus
    38222_27_bin.1	d__Archaea;p__Halobacteriota;c__Methanosarcinia;o__Methanosarcinales;f__Methanosarcinaceae;g__Methanosarcina;s__Methanosarcina mazei
    38222_27_bin.4	d__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__

3. **Membership file** — the complete representative list from your `magmax customdb` output
   
    Example:`memberships.tsv`
    #representative	member_genomes
    38222_26_bin.5	
    38222_27_bin.1	38510_111_bin.22,38510_42_bin.1,39907_150_bin.19,39923_2#63_bin.6,47680_139_bin.10,47681_116_bin.15
    38222_27_bin.4	38354_3_bin.7,38354_18_bin.10,38510_115_bin.7,38510_4_bin.3


### Output

A tab-separated file with two columns:
- **Column 1**: MAGmax customdb representative genomes + GTDB representative genomes whose species are not covered in the user input data
- **Column 2**: Matching GTDB reference genome for common species, or `unknown` if no match found

Example: `unified.tsv`
```
#user_representative	gtdb_representative
38222_26_bin.5	RS_GCF_006094395.1
38222_27_bin.1	RS_GCF_000970205.1
38222_27_bin.4	unknown
RS_GCF_000016525.1  RS_GCF_000016525.1
GB_GCA_002686525.1  GB_GCA_002686525.1
...
...
remaining GTDB reference genomes
```
