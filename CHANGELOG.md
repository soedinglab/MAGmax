# Changelog

## [1.3.0] - 2026-02-06
### Added
- New **sensitive** (from `repsel_bymaxedges` branch) to select representative bins by maximum number of high-ANI connectivity.
- Implemented using greedy maximum weighted dominating set algorithm
- Minor refactoring to improve code reability & updated README.md

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