## Generate custom species-level database

Genome dereplication is not always perfect due to inherent limitations of hierarchical clustering algorithms used in dereplication tools (dRep and Galah). Alternatively, taxonomic classification using GTDBtk followed by grouping genomes by taxonomy assignment is another option for dereplication, but it has limitations too: 1) ANI radius of under-represented species may be inaccurate, causing wrong taxonomy labeling; 2) novel species cannot be assigned. Combining dereplication and taxonomic classification can enhance the discovery of novel species with improved accuracy.

---

### Overview

The `magmax customdb` subcommand builds a species-level non-redundant genome database by combining two complementary strategies:

1. **GTDB-Tk-guided dereplication** — bins that are confidently assigned to a known species (ANI ≥ species ANI radius and aligned fraction ≥ cutoff) are grouped by species name and one representative is chosen per species.
2. **ANI-based dereplication of remaining bins** — bins that are unclassified or do not confidently hit a known species are clustered by ANI (default 95%), using the same graph-based approach as regular MAGmax dereplication.

The result is a combined set of representatives covering both known and novel species, each supported by a high-quality genome.

---

### Prerequisites

All tools required for regular MAGmax runs plus:
- **GTDB-Tk** — taxonomic classification of input bins (`gtdbtk classify_wf`).  
  Requires only the summary output file (`gtdbtk.bac120.summary.tsv` or equivalent).
- **CheckM2** — genome quality assessment (run automatically if `-q` is not provided).
- **skani** — pairwise ANI computation (run automatically; or provide a pre-computed file with `--anifile`).

---

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

---

### Input files

| Input | Flag | Required | Description |
|-------|------|----------|-------------|
| GTDB-Tk summary | `-g` | Yes | Tab-separated GTDB-Tk classification file (`gtdbtk.bac120.summary.tsv` or combined) |
| Bin directory | `-b` | Yes | Directory containing FASTA files of input bins |
| CheckM2 quality | `-q` | No | `quality_report.tsv` from CheckM2; computed automatically if omitted |
| Isolate genome list | `--isolate-genomes` | No | Plain text file, one genome name per line (with or without extension); these are prioritized as representatives |
| Pre-computed ANI | `--anifile` | No | Output of `skani triangle <bindir> -E -o <anifile>`; computed automatically if omitted |

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

---

### Output files

Output is written to `specieslevel_customdb/` by default (use `-o` to override), created next to the bin directory.

| File | Description |
|------|-------------|
| `gtdbtk_species_representatives.tsv` | Representatives selected from GTDB-Tk-classified (perfect) bins. Columns: `#gtdbtk_species_representative`, `species_name` |
| `memberships.tsv` | All representatives and their cluster members (perfect + unclassified). Tab-separated: representative, then tab-separated member list |
| `bins_checkm2_qualities.tsv` | Completeness and contamination of all final representatives. Columns: `#Bin`, `Completeness`, `Contamination` |
| `unclassified_clusterrepresentatives_gtdbtkspecies_ani_connections.tsv` | ANI connections between novel-cluster representatives and GTDB-Tk species representatives that exceed the species ANI radius. Columns: `#unclassified_cluster_representative`, `gtdbtk_species_representative`, `ANI`, `species_ANI_radius`. Use this to check whether any unclassified representative might actually belong to a known species at higher ANI stringency. |

---

### Options reference

```
magmax customdb [OPTIONS] -g <GTDBTK>

Required:
  -g, --gtdbtk <GTDBTK>            GTDB-Tk classification summary file

Quality filtering:
  -c, --completeness <FLOAT>       Minimum completeness of bins (%) [default: 90]
  -p, --purity <FLOAT>             Maximum contamination of bins (%) [default: 5]
  -q, --qual <QUAL>                CheckM2 quality file (quality_report.tsv);
                                   run automatically if omitted

Species-level ANI criteria:
      --species-ani <FLOAT>        ANI threshold for species-level clustering (%) [default: 95]
      --species-alignedfrac <FLOAT>
                                   Minimum aligned fraction for species-level clustering (%) [default: 50]

Input:
  -b, --bindir <BINDIR>            Directory containing FASTA files of bins
  -f, --format <FORMAT>            Bin file extension [default: fasta]
      --split                      Split bins by sample ID before processing
      --isolate-genomes <FILE>     File listing isolate genomes; prioritized as representatives
      --anifile <ANIFILE>          Pre-computed skani ANI file
                                   (skani triangle <bindir> -E -o <anifile>)

Representative selection:
      --sensitive                  Select representatives based on high ANI connectivity
                                   instead of quality score (no reassembly)

Output:
  -o, --outdir <OUTPUT>            Output directory [default: specieslevel_customdb/]
  -t, --threads <THREADS>          Number of threads [default: 8]
```

---

### How representatives are selected

**GTDB-Tk-classified (perfect) bins**

All bins assigned to the same GTDB-Tk species are grouped together. One representative is chosen per species:
1. If isolate genomes are present among the species members, the isolate with the lowest contamination is preferred.
2. Otherwise, the bin with the highest quality score (`completeness − 5 × contamination`) is selected.

**Unclassified or novel-species (remaining) bins**

Remaining bins are clustered by pairwise ANI (default 95%, aligned fraction ≥ 50%) using single-linkage (connected components). Within each cluster:
- **Default mode**: selects the highest-quality bin (completeness ≥ 90% required; isolates are prioritized).
- **Sensitive mode** (`--sensitive`): selects the bin with the highest weighted ANI connectivity (`Σ max(0, ANI − threshold)` over neighbors), favoring bins that are more similar to a larger number of neighbors.

After initial selection, a deduplication step merges any two remaining representatives that are themselves redundant (ANI ≥ species cutoff), keeping the higher-quality one.

---

### Notes

1. **Quality thresholds for database creation are stricter than regular dereplication.** The defaults are completeness ≥ 90% and contamination ≤ 5%. Adjust with `-c` and `-p` if needed.

2. **Pre-computing ANI saves time on re-runs.** Generate the ANI file once with `skani triangle <bindir> -E -o ani_edges` and pass it via `--anifile` to skip re-computation. The file is also cached automatically as `<bindir>/subset_ani_edges` during the first run.

3. **The `--split` flag is needed when bins are not already separated by sample ID.** When running multi-sample binning on concatenated contigs (e.g., MetaBAT2 or COMEBin), use `--split` to let MAGmax separate bins by sample before processing.

4. **Isolate genome names must match bin filenames** (extension is stripped automatically; path prefixes are ignored). A warning is printed for any isolate genome that does not pass quality filtering.

5. **The `unclassified_clusterrepresentatives_gtdbtkspecies_ani_connections.tsv` file is a diagnostic resource.** It lists novel-cluster representatives whose ANI to a known GTDB-Tk species representative meets or exceeds that species' ANI radius, suggesting potential taxonomic placement at stricter thresholds.
