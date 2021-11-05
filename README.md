# cgpCaVEManPostProcessing

[![DOI][zenodo-badge]][zenodo-link]

cgpCaVEManPostProcessing is used to apply filtering on raw VCF calls generated using [CaVEMan][caveman].

For details of the underlying algorithm please see the [CaVEMan][caveman] site.

| Master                                        | Develop                                         |
| --------------------------------------------- | ----------------------------------------------- |
| [![Master Badge][circle-master]][circle-base] | [![Develop Badge][circle-develop]][circle-base] |

- [Docker, Singularity and Dockstore](#docker-singularity-and-dockstore)
- [Usage](#usage)
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
