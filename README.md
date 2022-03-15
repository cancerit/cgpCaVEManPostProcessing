# cgpCaVEManPostProcessing

[![DOI][zenodo-badge]][zenodo-link]

cgpCaVEManPostProcessing is used to apply filtering on raw VCF calls generated using [CaVEMan][caveman].

For details of the underlying algorithm please see the [CaVEMan][caveman] site.

| Master                                        | Develop                                         |
| --------------------------------------------- | ----------------------------------------------- |
| [![Master Badge][circle-master]][circle-base] | [![Develop Badge][circle-develop]][circle-base] |

- [Docker, Singularity and Dockstore](#docker-singularity-and-dockstore)
- [Usage](#usage)
  - [Flagging CaVEMan Files](#flagging-caveman-files)
  - [Utility Scripts](#utility-scripts)
    - [cavemanPostProcessing_ini_to_yaml.pl](#cavemanpostprocessing_ini_to_yamlpl)
- [Dependencies/Install](#dependenciesinstall)
- [Creating a release](#creating-a-release)
  - [Preparation](#preparation)
  - [Cutting the release](#cutting-the-release)
- [LICENCE](#licence)

## Docker, Singularity and Dockstore

cgpCaVEManPostProcessing is available as a separate docker image on quay.io.

- [cgpcavemanpostprocessing][ds-cg-cpp]

And as part of pre-built full analysis images on quay.io.

- [dockstore-cgpwxs][ds-cgpwxs-git]
  - Contains tools specific to WXS analysis.
- [dockstore-cgpwgs][ds-cgpwgs-git]
  - Contains additional tools for WGS analysis.

These were primarily designed for use with dockstore.org but can be used as normal containers.

## Usage

More detailed instructions can be found on the [wiki]

As of version 1.10.0 cgpCaVEManPostProcessing has new WXS flags available that are not used by default.
These were developed in conjunction with the Dermatlas project and caution is advised when using them. 

- cavemanMatchNormalProportionFlag
- withinGapRangeFlag

Full flag definitions can be found [here](config/flag.to.vcf.convert.ini)

Flags can be tuned by modifying their parameters in the species ini file.
Human example is [here](config/Human/GRCh37d5/flag.vcf.config.ini).
The parameter names correspond to names in the flag descriptions.
Further details of flags and parameters are available in the [wiki].

### Flagging CaVEMan Files

Flag/Post Process CaVEMan files.

```bash
cgpFlagCaVEMan.pl [-h] -f vcfToFlag.vcf -o flaggedVCF.vcf -c
    configFile.yaml -s human -t pulldown -v vcfFlagNames.ini -n norm.bam -m
    tum.bam [-u unmatchedStore.tmp]

      General Options:

        --help                 (-h)       Brief documentation

        --version              (-version) Output the version number and exit

        --input                (-i)       The VCF input file to flag.

        --outFile              (-o)       The VCF output file to write.

        --species              (-s)       Species associated with this vcf file to use.

        --species-assembly     (-sa)      Species assembly for (output in VCF)

        --tumBam               (-m)       Tumour bam file

        --normBam              (-n)       Normal bam file

        --bedFileLoc           (-b)       Path to a folder containing the centromeric, snp, hi sequence depth,
                                          and simple repeat sorted (gzipped and tabixed) bed files (if required) i.e. the non annotation bed files.
                                          Names of files will be taken from the config file.

        --indelBed             (-g)       A bed file containing germline indels to filter on

        --unmatchedVCFLoc      (-umv)     Path to a directory containing the unmatched VCF normal files listed in the
                                          config file or unmatchedNormal.bed.gz (bed file is used in preference).

        --annoBedLoc           (-ab)      Path to bed files containing annotatable regions and coding regions.

        --reference            (-ref)     Reference index (fai) file corresponding to the mapping of the data being processed.
                                            (must have corresponding fasta file co-located)

        --index                (-idx)     Index of the job (to override LSB_JOBINDEX as used on LSF farms)

        --verbose

      OPTIONAL:

        --sampleToIgnoreInUnmatched    (-sp) Unmatched normal to ignore (to be used if the sample is one of those with a normal in the panel).

        --processid                    (-p)  Id anaylsis process to be added at a CGP specific header.

        --flagConfig                   (-c)  Config ini file to use for flag list and settings.

        --flagToVcfConfig             (-v)  Config::Inifiles style config file containing VCF flag code to flag name conversions see
                                             ../config/flag.to.vcf.convert.ini for example

        --studyType            (-t)       Study type, used to decide parameters in file (genome|genomic|WGS|pulldown|exome|WXS|followup|AMPLICON|targeted|RNA_seq).

      Examples:

        cgpFlagCaVEMan.pl [-h] -f vcfToFlag.vcf -o flaggedVCF.vcf -c configFile.ini -s human -t pulldown
```

### Utility Scripts

#### cavemanPostProcessing_ini_to_yaml.pl

Convert old .ini file to .yaml format.

```bash
cavemanPostProcessing_ini_to_yaml.pl [-h] -f flag_config.ini -o flag_config.yml 

  General Options:

    --help                 (-h)       Brief documentation

    --version              (-version) Output the version number and exit

    --input                (-i)       The VCF input file to flag.

    --outfile              (-o)       The VCF output file to write.

  Examples:

    cavemanPostProcessing_ini_to_yaml.pl [-h] -f flag_config.ini -o flag_config.yml
```

## Dependencies/Install

Please ensure the following packages are available. Alternatively use the Docker image.

- [cgpVcf][cgpvcf]
- [Bio::DB::HTS][bio-db-hts]
  - If you have an install of PCAP-core this is already available

Once complete please run:

```bash
./setup.sh /some/install/location
```

`setup.sh` will also install bedtools for you.

## Creating a release

### Preparation

- Commit/push all relevant changes.
- Pull a clean version of the repo and use this for the following steps.

### Cutting the release

1. Update `lib/Sanger/CGP/CavemanPostProcessor.pm` to the correct version.
2. Update `CHANGES.md`
3. Run `./prerelease.sh`
4. Check all tests and coverage reports are acceptable.
5. Commit the updated docs tree and updated module/version.
6. Push commits.
7. Use the GitHub tools to draft a release.

## LICENCE

```txt
Copyright (c) 2014-2018 Genome Research Ltd.

Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>

This file is part of cgpCaVEManPostProcessing.

cgpCaVEManPostProcessing is free software: you can redistribute it and/or modify it under
the terms of the GNU Affero General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option) any
later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

1. The usage of a range of years within a copyright statement contained within
this distribution should be interpreted as being equivalent to a list of years
including the first and last year specified and all consecutive years between
them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
2009, 2011-2012’ should be interpreted as being identical to a statement that
reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
2009, 2010, 2011, 2012’."
```

<!-- References -->
[caveman]: http://cancerit.github.io/CaVEMan
[cgpvcf]: https://github.com/cancerit/cgpVcf/releases
[bio-db-hts]: http://search.cpan.org/dist/Bio-DB-HTS
[ds-cgpwxs-git]: https://github.com/cancerit/dockstore-cgpwxs
[ds-cgpwgs-git]: https://github.com/cancerit/dockstore-cgpwgs
[ds-cg-cpp]: https://quay.io/repository/wtsicgp/cgpcavemanpostprocessing
[wiki]: https://github.com/cancerit/cgpCaVEManPostProcessing/wiki

<!-- Travis -->
[travis-base]: https://travis-ci.org/cancerit/cgpCaVEManPostProcessing
[travis-master]: https://travis-ci.org/cancerit/cgpCaVEManPostProcessing.svg?branch=master
[travis-develop]: https://travis-ci.org/cancerit/cgpCaVEManPostProcessing.svg?branch=dev

<!-- Circle-ci -->
[circle-base]: https://circleci.com/gh/cancerit/cgpCaVEManPostProcessing.svg?style=shield
[circle-master]: https://circleci.com/gh/cancerit/cgpCaVEManPostProcessing.svg?style=shield&branch=master;
[circle-develop]: https://circleci.com/gh/cancerit/cgpCaVEManPostProcessing.svg?style=shield&branch=dev;

<!-- Zenodo -->
[zenodo-badge]: https://zenodo.org/badge/doi/10.5281/zenodo.18404.svg
[zenodo-link]: http://dx.doi.org/10.5281/zenodo.18404
