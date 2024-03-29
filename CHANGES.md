# CHANGES

## 1.12.0

- Add mnv flagging support.
  - Use MNVFLAGLIST to define flags to be run over MNVs in combination with mnvFlag in the original FLAGLIST section.

## 1.11.0

- flag.vcf.config.ini can now be passed as a .yaml or .yml file.
- Added a utility script to convert old config ini to new yaml format.

## 1.10.0

- add Circleci support
- Addition of optional WXS flags developed alongside Dermatlas project

## 1.9.3

- Fix bug in caluclating unmatched normal proportion when using VCF not BED.

## 1.9.2

- Correct bug in exome matched normal proportion settings

## 1.9.0

- Fix edgecase in overlapping read processing
  - Resulted in very small difference in passed variants at WGS level (+1)
  - See #36
- Significant performance improvements, see #37

## 1.8.9

- Speedups in callback methods
- Addition of IntervalTree use instead of tabix

## 1.8.7

- Reduce memory overhead when using BED as unmatched normal panel input

## 1.8.6

- Code modified for overlapping reads. Where reads overlap but carry the same base on each,
- they will be assigned alternately to each strand so as to ensure an even spread of strandedness.

### Behaviour change

**Where the proper pair filter flag is used, this code now checks that the paired-end orientation is also used.**
**This will mean that mate-pair orientation (F/F or R/R) will be rejected**
**If you wish to use mate-pair data, please use previous version**

- Where a proper pair filter is used, now check for the correct paired-end orientation of F/R.
- If this is not met the read is ignored.

## 1.8.5

- Added TARGETED to GRCh38 human config ini

## 1.8.4

- added 3 soft flags (CLPM, ASMD and ASRD) to all species builds and design_types (except RNAseq)

## 1.8.2

- Update Tabix `query` to `query_full` in order to cope with GRCh38 contig names (requires `Bio::DB::HTS` v2.10 to be installed by PCAP).

## 1.8.1

- `cgpFlagCaVEMan.pl` option `-sa` did not capture value

## 1.8.0

- Overlapping reads now handled
- Updates to licenese headers and contact information
- Support for csi and crai indexes
- README cleanup
- Ability to exclude contigs from analysis

## 1.7.3

- Change tabix->query to tabix->query_full

## 1.7.2

- Updating species config files to ensure they always point to gz files

## 1.7.1

- Convert warning to error where no flags were included in list
- Added species and build specific config files

## 1.7.0

- Added travis CI testing

## 1.6.9

- Correct sorting of flags when output to vcf file

## 1.6.7

- Add targeted to list of recognised types in cgpFlagCaVEMan.pl

## 1.6.6

- Add Amplicon to list of recognised types in cgpFlagCaVEMan.pl

## 1.6.5

- Correct bug where sequencing type names in script and config differ

## 1.6.3

- Removal of extra reference to Bio::DB::Bam::Alignment

## 1.6.1

- Replaced prerequestite test for Bio::DB::Sam with Bio::DB::HTS in setup.sh

## 1.6.0

- Now uses UCSC defined types as opposed to genome/exome etc. WXS,WGS,AMPLICON,RNASEQ defined so far
 -Uses Bio::DB::HTS rather than Bio::DB::Sam - cram can be used
