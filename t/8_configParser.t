# Copyright (c) 2014-2021
#
# Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
#
# This file is part of cgpCaVEManPostProcessing.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# 1. The usage of a range of years within a copyright statement contained within
# this distribution should be interpreted as being equivalent to a list of years
# including the first and last year specified and all consecutive years between
# them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
# 2009, 2011-2012’ should be interpreted as being identical to a statement that
# reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
# statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
# identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
# 2009, 2010, 2011, 2012’.
#

use strict;
use Sanger::CGP::CavemanPostProcessor::ConfigParser;
use Const::Fast qw(const);
use Data::Dumper;

use Test::More tests => 5;

use FindBin qw($Bin);
my $script_path = "$Bin/../bin/";
my $test_data_path = "$Bin/../testData/";
my $test_output_yaml = "$Bin/../testData/test_out.yml";
my $exp_out_yaml = "$Bin/../testData/expected_converted_config.yaml";

const my $YAML_CONFIG => $test_data_path.'test_config.yaml';
const my $INI_CONFIG => $test_data_path.'flag.vcf.config_test.ini';
const my $SPECIES => 'HUMAN';
const my $SEQ_TYPE => 'WGS';

const my %EXP_SECT_PARAMS => (
          'minNormMutAllelequal' => 15,
          'samePosMaxPercent' => 80,
          'minDepthQual' => 25,
          'vcfUnmatchedMinSamplePct' => 1,
          'maxNormIndelProportion' => 10,
          'pentamerMinPassAvgQual' => 20,
          'minPassPhaseQual' => 21,
          'minRdPosDepth' => 8,
          'readPosBeginningOfReadIgnoreProportion' => '0.08',
          'keepSW' => 1,
          'minAnalysedQual' => 11,
          'vcfUnmatchedMinMutAlleleCvg' => 3,
          'matchedNormalMaxMutProportion' => '0.2',
          'maxMatchedNormalAlleleProportion' => '0.05',
          'maxMatchedNormalAlleleHiCvgProportion' => '0.03',
          'readPosTwoThirdsOfReadExtendProportion' => '0.08',
          'maxTumIndelProportion' => 10,
          'maxPhasingMinorityStrandReadProportion' => '0.04',
          'matchedNormalAlleleHiCvgCutoff' => 2,
          'minSingleEndCoverage' => 10,
          'minPassAvgMapQual' => 21,
          );

const my @EXP_FLAGLIST => (
          'depthFlag',
          'readPositionFlag',
          'matchedNormalFlag',
          'pentamericMotifFlag',
          'avgMapQualFlag',
          'simpleRepeatFlag',
          'centromericRepeatFlag',
          'snpFlag',
          'phasingFlag',
          'hiSeqDepthFlag',
          'germlineIndelFlag',
          'unmatchedNormalVcfFlag',
          'singleEndFlag',
          'matchedNormalProportion',
          'alignmentScoreReadLengthAdjustedFlag',
          'clippingMedianFlag',
          'alnScoreMedianFlag'
          );

const my %EXP_BED_PARAMS => (
          'centromericRepeatBed' => 'centromeric_repeats.bed.gz',
          'highSeqDepthBed' => 'hi_seq_depth.bed.gz',
          'codingBed' => 'codingexon_regions.sub.bed.gz',
          'snpBed' => 'snps.bed.gz',
          'annotatableBed' => 'gene_regions.bed.gz',
          'germlineIndelBed' => undef,
          'simpleRepeatBed' => 'simple_repeats.bed.gz',
        );

const my %EXP_BED_PARAMS_INI => (
          'centromericRepeatBed' => 'centromeric_repeats.bed.gz',
          'highSeqDepthBed' => 'hi_seq_depth.bed.gz',
          'codingBed' => 'codingexon_regions.sub.bed.gz',
          'snpBed' => 'snps.bed.gz',
          'annotatableBed' => 'gene_regions.bed.gz',
          'germlineIndelBed' => '',
          'simpleRepeatBed' => 'simple_repeats.bed.gz',
        );


#----------------
#	Init tests
#----------------

subtest 'Initialise module' => sub {
  my $parser = new_ok('Sanger::CGP::CavemanPostProcessor::ConfigParser' => []);
  done_testing();
};

subtest '_get_yml_config_params' => sub {
  my ($sectParams, $flagList, $bedFileParams, $mnv_flaglist) = 
  Sanger::CGP::CavemanPostProcessor::ConfigParser::_get_yml_config_params($YAML_CONFIG, $SPECIES, $SEQ_TYPE);
  is_deeply($sectParams,\%EXP_SECT_PARAMS, '_get_yml_config_params'." sectParams");
  is_deeply($flagList,\@EXP_FLAGLIST, '_get_yml_config_params'." flagList");
  is_deeply($bedFileParams,\%EXP_BED_PARAMS, '_get_yml_config_params'." bedFileParams");
  done_testing();
};

subtest '_get_cfg_ini_config_params' => sub {
  my ($sectParams, $flagList, $bedFileParams, $mnv_flaglist) = 
  Sanger::CGP::CavemanPostProcessor::ConfigParser::_get_cfg_ini_config_params($INI_CONFIG, $SPECIES, $SEQ_TYPE, undef);
  is_deeply($sectParams,\%EXP_SECT_PARAMS, '_get_cfg_ini_config_params'." sectParams");
  is_deeply($flagList,\@EXP_FLAGLIST, '_get_cfg_ini_config_params'." flagList");
  is_deeply($bedFileParams,\%EXP_BED_PARAMS_INI, '_get_cfg_ini_config_params'." bedFileParams");
  done_testing();
};

subtest 'getConfigParams' => sub {
  my %opts = (
          'c' => $INI_CONFIG,
          's' => $SPECIES,
          't' => $SEQ_TYPE,
          );
  my ($sectParams, $flagList, $bedFileParams, $mnv_flaglist) = 
  Sanger::CGP::CavemanPostProcessor::ConfigParser::getConfigParams(\%opts);
  is_deeply($sectParams,\%EXP_SECT_PARAMS, 'getConfigParams_ini'." sectParams");
  is_deeply($flagList,\@EXP_FLAGLIST, 'getConfigParams_ini'." flagList");
  is_deeply($bedFileParams,\%EXP_BED_PARAMS_INI, 'getConfigParams_ini'." bedFileParams");

  $opts{'c'} = $YAML_CONFIG;

  ($sectParams, $flagList, $bedFileParams) = 
  Sanger::CGP::CavemanPostProcessor::ConfigParser::getConfigParams(\%opts);
  is_deeply($sectParams,\%EXP_SECT_PARAMS, 'getConfigParams_yaml'." sectParams");
  is_deeply($flagList,\@EXP_FLAGLIST, 'getConfigParams_yaml'." flagList");
  is_deeply($bedFileParams,\%EXP_BED_PARAMS, 'getConfigParams_yaml'." bedFileParams");
  done_testing();
};

subtest 'convert_ini_to_yaml' => sub {
  Sanger::CGP::CavemanPostProcessor::ConfigParser::convert_ini_to_yaml($INI_CONFIG, $test_output_yaml);
  compare_files($test_output_yaml, $exp_out_yaml);
  #Delete output file
  ok(unlink($test_output_yaml)==1);
  done_testing();
};

sub compare_files{
  my ($input_file, $exp_file) = @_;
  my @dat;
  my $exp;
  my $input;
  open(my $IN, '<', $input_file) || croak("Error reading test output $input_file:$!");
    @dat = <$IN>;
    $input = join("",@dat);
  close($IN);
  open(my $EXP,'<',$exp_file) || croak("Error reading expected output $exp_file:$!");
    @dat = <$EXP>;
    $exp = join("",@dat);
  close($EXP);
  if(!$exp eq $input){
    warn ("Compare files ".$input_file." - ".$exp_file." failed. \n$input \n!=\n$exp\n");
    ok(0==1);
  }else{
    ok(1==1, "Compare files ".$input_file." - ".$exp_file."\n");
  }
}
