# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [2.0.0-beta.7] - 2026-05-11

### Added
- Added version number to splash screen.

### Changed
- Changed all version displays to be the same format
- Updated the documentation

### Fixed
- Fixed multiple inclusions of the schema files
- Fixed some version messages being logged in the wrong format

## [2.0.0-beta.6] - 2026-04-30

### Fixed
- fixed pypi builds not including all relevant files

## [2.0.0-beta.5] - 2026-04-30

### Changed
- Switched to `pyproject.toml` for handling the packaging
- Enhanced JSON schema definitions for path validation in configuration files
- Unified log handling across all parts of the pipeline

## [2.0.0-beta.4] - 2026-03-24

### Added

- Expose feature prediction parameters to user via optional config parameters.
- Append genomic sequences to GFF file, following the Bakta/Prokka convention, allowing for greater downstream application compatibility.

### Changed

- Refactor code structure to add cli-package while keeping Snakemake workflow catalog compatibility.

## [2.0.0-beta.3] - 2026-03-17

### Added
- Added CRISPR prediction via Diced.
- Add RNA and CRISPR predictions to GFF and GBK output
- Add ARM64 support
### Changed
- Relaxed Snakemake minimum version requirement to v9.0.1

## [2.0.0-beta.2] - 2026-02-23

### Changed

- Updated the mapping file parsing to make it more resilient against unexpected formats and warn the user accordingly.
- Added option to join on accession or query_name, as some databases may use one or the other.

## [2.0.0-beta.1] - 2026-02-12

> ⚠️ **BREAKING CHANGES** ⚠️ — This is a complete rewrite of Prokanota with significant architectural improvements.

### Added

- Config-driven database management via `config/databases.yaml`.
- Dynamic rule generation — the pipeline automatically creates Snakemake rules for each enabled database.
- Standardized output format across all database types.
- Support for pyhmmer for HMM databases (faster than hmmsearch).
- Support for rpsblast for RPS-BLAST databases.
- Support for DIAMOND for BLAST-style searches.
- Support for mmseqs2 for ultra-fast sequence searches.
- Protein-only mode: annotate pre-computed protein sequences without gene prediction.
- Mixed sample support: combine genomic and protein-only samples in the same run.
- Flexible input handling for diverse annotation workflows.
- Comprehensive testing suite with automated CI/CD.
- Detailed wiki documentation with step-by-step database setup guides.

## [1.2.0] - 2025-12-05

### Added

- Automatic switching to pyrodigal meta mode if input sequence is too short (<20 kb). This action is logged in the log files.
- Config schema validation.
- Safer thread calculation.

### Fixed

- In rare cases, the parsing of the CDD database results could use wrong temp files when rule order was not correct, leading to wrong results. This has been fixed.

## [1.1.1] - 2025-12-04

### Changed

- Updated dependencies to newest versions.
- Added more documentation to tRNA-scan version parsing.

### Fixed

- When no tRNAs were found (e.g. when processing incomplete genomes), the pipeline would fail with an error. This was changed to return an empty `rna.tsv` file and write the info that no tRNAs were found to the log.
- Fixed minor typos.

## [1.1.0] - 2025-03-20

### Changed

- Renamed CDS prediction and all surrounding rules/scripts/files to "feature prediction", as this more accurately reflects what happens in this step. This also changes the output directory from `results/cds` to `results/features`.
- Added option for genomes to be output with contig headers replaced by internal ones for downstream applications.

### Fixed

- Fixed GFF3 spec compliance.

### Breaking Changes

- The output directory `results/cds` has been renamed to `results/features`.

## [1.0.0] - 2025-03-19

### Added

- Implemented tRNA prediction using tRNAscan-SE in `-G` mode.
- Added intensive logging to the main CDS and RNA prediction script.

### Fixed

- Fixed some write protection issues.

## [0.2.0] - 2025-03-18

### Added

- Implemented rRNA prediction using pybarrnap in `--accurate` mode, using CM models from Rfam (14.10).

### Fixed

- Fixed pipeline failure on some systems missing `curl`.

## [0.1.1] - 2025-03-12

### Changed

- Applied changes to comply with the [Snakemake Workflow Catalog rules](https://snakemake.github.io/snakemake-workflow-catalog/?rules=true).

## [0.1.0] - 2025-03-07

### Added

- Initial Snakemake Pipeline release.


---
[Unreleased]: https://github.com/richardstoeckl/prokanota/compare/v2.0.0-beta.7...HEAD
[2.0.0-beta.7]: https://github.com/richardstoeckl/prokanota/compare/v2.0.0-beta.6...v2.0.0-beta.7
[2.0.0-beta.6]: https://github.com/richardstoeckl/prokanota/compare/v2.0.0-beta.5...v2.0.0-beta.6
[2.0.0-beta.5]: https://github.com/richardstoeckl/prokanota/compare/v2.0.0-beta.4...v2.0.0-beta.5
[2.0.0-beta.4]: https://github.com/richardstoeckl/prokanota/compare/v2.0.0-beta.3...v2.0.0-beta.4
[2.0.0-beta.3]: https://github.com/richardstoeckl/prokanota/compare/v2.0.0-beta.2...v2.0.0-beta.3
[2.0.0-beta.2]: https://github.com/richardstoeckl/prokanota/compare/v2.0.0-beta.1...v2.0.0-beta.2
[2.0.0-beta.1]: https://github.com/richardstoeckl/prokanota/compare/v1.2.0...v2.0.0-beta.1
[1.2.0]: https://github.com/richardstoeckl/prokanota/compare/v1.1.1...v1.2.0
[1.1.1]: https://github.com/richardstoeckl/prokanota/compare/v1.1.0...v1.1.1
[1.1.0]: https://github.com/richardstoeckl/prokanota/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/richardstoeckl/prokanota/compare/v0.2.0...v1.0.0
[0.2.0]: https://github.com/richardstoeckl/prokanota/compare/v0.1.1...v0.2.0
[0.1.1]: https://github.com/richardstoeckl/prokanota/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/richardstoeckl/prokanota/releases/tag/v0.1.0
