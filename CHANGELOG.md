# Change Log
All notable changes to the misc tools will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/), 
and this project adheres to [Semantic Versioning](http://semver.org/).

## [0.2.3] - In progress
###
- Updated download-srr-files to use `misc.shell_utils`

## [0.2.2] - 2017-06-14
### Fixed
- References to old gtf functions in `misc.bio`

## [0.2.1] - 2017-06-14
### Fixed
- All references to `misc.bio` and `misc.bio_utils`. Please see
  [Issue #1](https://github.com/bmmalone/pybio-utils/issues/1) in the new repo
  for more details.

## [0.2.0] - 2017-05-31
This is a new version which moves the project from Bitbucket to GitHub.
Additionally, the other utilities (`misc.***`) have been completely
removed. They will be added to a new
[`pymisc-utils`](https://github.com/bmmalone/pymisc-utils) repository.

## [0.1.6] - 2017-05-10
### Fixed
- Missing import in counting alignments for bam files

## [0.1.5] - 2017-05-09
### Updated
- `get-read-length-distribution` script to handle bam, fasta and fastq files.
  See [the bio docs](docs/bio.md#get-read-length-distributions) for more
  details.

## [0.1.4] - 2017-05-09
### Removed
- bed, bam and gtf helpers from `bio.py`. These had already been deprecated for
  quite some time.

## [0.1.3] - 2017-03-30
### Added
- Script to remove duplicate entries from bed files. See
  [the bio docs](docs/bio.md#merge-bed12-files-and-remove-duplicate-entries)
  for more details.

## [0.1.2] and previous versions

The initial versions have not been documented in the change log.


