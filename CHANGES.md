## 1.7.1
*   Convert warning to error where no flags were included in list
*   Added species and build specific config files

## 1.7.0
*   Added travis CI testing

## 1.6.9
*   Correct sorting of flags when output to vcf file

## 1.6.7
*   Add targeted to list of recognised types in cgpFlagCaVEMan.pl

## 1.6.6
*   Add Amplicon to list of recognised types in cgpFlagCaVEMan.pl

## 1.6.5
*   Correct bug where sequencing type names in script and config differ

## 1.6.3
*   Removal of extra reference to Bio::DB::Bam::Alignment

## 1.6.1
*   Replaced prerequestite test for Bio::DB::Sam with Bio::DB::HTS in setup.sh

## 1.6.0
*   Now uses UCSC defined types as opposed to genome/exome etc. WXS,WGS,AMPLICON,RNASEQ defined so far
*   Uses Bio::DB::HTS rather than Bio::DB::Sam - cram can be used
