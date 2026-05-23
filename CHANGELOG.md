# Changelog

## [1.4.0] - 2026-05-22
### Added
  - customdb subcommand: incorporates GTDB taxonomy and isolate genome information to group bins and prioritize reference selection during dereplication at the species-level
  - support for gzip-compressed FASTQ input (.fastq.gz) with transparent on-the-fly decompression

### Changed
  - memory usage and hash lookup overhead in ANI data structures
  - mutex contention in parallel processing with lock-free concurrent maps
  - double for loop to parallelization of graph construction for large bin clusters
  - streaming I/O for FASTA files to reduce peak memory consumption

### Bug Fixes
  - explicit error message on CheckM2 execution failure
  - completeness threshold takes user-defined threshold (default 90%) while selecting high-quality representative bin
  - raise error when necessary info is missing in ANI input file

  Code Quality
  - Centralised quality score calculation to a single function, eliminating duplicated logic across modules

## [1.3.0] - 2026-02-06
### Added
- New **sensitive** mode (from `repsel_bymaxedges` branch) to select representative bins by maximum number of high-ANI connectivity
- Implemented using greedy maximum weighted dominating set algorithm
- Minor refactoring to improve code reability

## [1.2.0] - 2026-01-05
### Added
- New **parallel clique-based dereplication workflow** (from `parallel-cliques` branch) to improve scalability on large connected components
- `--anifile` command-line option to directly use existing skani output (avoids repeated skani runs if user wants to dereplicate genomes under different quality and ANI thresholds)
- `-a`/`--alignedfrac` command-line option to apply a fraction-aligned threshold when detecting redundant genomic pairs

### Changed
- Final dereplication step in **no-reassembly mode** now reuses the original skani output instead of rerunning skani
- Minor refactoring to improve code readability

### Performance
- Benchmarked MAGmax performance against *Galah*

## [1.1.0] - 2025-08-26
### Added
- New feature: `no-reassembly` command-line option.
- It enables dereplication of input bins without bin merging and reassembly