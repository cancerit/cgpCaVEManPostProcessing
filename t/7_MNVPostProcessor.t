# Copyright (c) 2014-2022
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
use Sanger::CGP::CavemanPostProcessor::MNVPostProcessor;
use Data::Dumper;
use Bio::DB::HTS;
use Const::Fast qw(const);

use Test::More tests => 5;

use FindBin qw($Bin);
my $lib_path = "$Bin/../lib";
my $test_data_path = "$Bin/../testData/";

const my $T_BAM => $test_data_path.'test.bam';
const my $T_BAI => $test_data_path.'test.bam.bai';

const my $PENT_BAM => $test_data_path.'pent.bam';
const my $PENT_BAI => $test_data_path.'pent.bam.bai';

const my $GAP_N_BAM => $test_data_path.'gap_flag_test_normal.bam';
const my $GAP_N_BAI => $test_data_path.'gap_flag_test_normal.bam.bai';
const my $GAP_T_BAI => $test_data_path.'gap_flag_test_tumour.bam.bai';
const my $GAP_T_BAM => $test_data_path.'gap_flag_test_tumour.bam';

const my $CLIP_M_BAM => $test_data_path.'clip.m.bam';
const my $CLIP_M_BAI => $test_data_path.'clip.m.bam.bai';
const my $CLIP_N_BAM => $test_data_path.'clip.n.bam';
const my $CLIP_N_BAI => $test_data_path.'clip.n.bam.bai';

const my $MNV_VCF => $test_data_path.'mnv.vcf';
const my $MNV_T_BAM => $test_data_path.'test_MNV_tum.bam';
const my $MNV_N_BAM => $test_data_path.'test_MNV_norm.bam';

my $processor;

#----------------
#	Init tests
#----------------

subtest 'Initialise module (no params)' => sub {
  eval{Sanger::CGP::CavemanPostProcessor::MNVPostProcessor->new();};
  ok(defined($@),"Defined error from eval");
  ok($@ =~ m/tumBam and normBam are required for initialisation/, "Correct error thrown");
  done_testing();
};

subtest 'Initialise module (bam params 1)' => sub {
  my $processor = new_ok('Sanger::CGP::CavemanPostProcessor::MNVPostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
  isa_ok($processor->tumBam(), "Bio::DB::HTS", "Test tumour bam");
  isa_ok($processor->normBam(), "Bio::DB::HTS", "Test normal bam");
  done_testing();
};

subtest 'Initialise module (bam params)' => sub {
  my $processor = new_ok('Sanger::CGP::CavemanPostProcessor::MNVPostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
  isa_ok($processor->tumBam(), "Bio::DB::HTS", "Test tumour bam");
  isa_ok($processor->normBam(), "Bio::DB::HTS", "Test normal bam");
  ok($processor->minDepthQual == 25,"Min depth qual");
  ok($processor->depthCutoffProportion == (0.333333), "depthCutoffProportion");
  ok($processor->minNormalMutAlleleQual == 15,"Min normal mut allele qual");
  ok($processor->minAnalysedQual == 11,"Min analysed qualities");
  ok($processor->percentageSamePos == 80,"Same position max pct");
  ok($processor->maxTumIndelProportion == 10,"Max tumour indel proportion");
  ok($processor->pentamerMinPassAvgQual == 20,"Pentamer minimum pass average quality");
  ok($processor->maxNormIndelProportion == 10,"Max normal indel proportion");
  ok($processor->minPassAvgBaseQualPhasing == 21,"Min pass phase ");
  ok($processor->minPassAvgMapQual == 21,"Min pass avg map quality");
  ok($processor->keepSW == 0,"Keep SW off by default");
  ok($processor->maxPhasingMinorityStrandReadProportion == 0.04,"maxPhasingMinorityStrandReadProportion");
  ok($processor->maxMatchedNormalAlleleProportion == 0.05,"maxMatchedNormalAlleleProportion");
  ok($processor->readPosBeginningOfReadIgnoreProportion == 0.08,"readPosBeginningOfReadIgnoreProportion");
  ok($processor->readPosTwoThirdsOfReadExtendProportion == 0.08,"readPosTwoThirdsOfReadExtendProportion");
  ok($processor->minRdPosDepth == 8,"Min rd pos depth default");
  ok($processor->matchedNormalMaxMutProportion == 0.2,'minUnmatchNormMutAlleleProportion');
  ok($processor->minSingleEndCoverage == 10,'minSingleEndCoverage');
  ok($processor->minGapFlagDistEndOfReadPercent == 75,'minGapFlagDistEndOfReadPercent');
  ok($processor->maxGapFlagDistFromEndOfReadProp == 0.13, 'maxGapFlagDistFromEndOfReadProp');
  ok($processor->meanMapQualGapFlag == 10, 'meanMapQualGapFlag');
  ok($processor->minGapPresentInPercentReads == 30, 'minGapPresentInPercentReads');
  ok($processor->withinXBpOfDeletion == 10, 'withinXBpOfDeletion');
  ok($processor->maxCavemanMatchedNormalProportion == 0.2, 'maxCavemanMatchedNormalProportion');
  done_testing();
};

subtest 'Initialise module (ALL params)' => sub {
  my $processor = new_ok('Sanger::CGP::CavemanPostProcessor::MNVPostProcessor' => [
                                                                'tumBam' => $T_BAM,
                                                                'normBam' => $T_BAM,
                                                                'minDepthQual' => 2,
                                                                'depthCutoffProportion' => (1/2),
                                                                'minNormMutAllelequal' => 3,
                                                                'minAnalysedQual' => 5,
                                                                'samePosMaxPercent' => 8,
                                                                'keepSW' => 1,
                                                                'maxTumIndelProportion' => 9,
                                                                'maxNormIndelProportion' => 12 ,
                                                                'pentamerMinPassAvgQual'  => 11,
                                                                'minPassPhaseQual'=> 13,
                                                                'minPassAvgMapQual' => 14,
                                                                'maxPhasingMinorityStrandReadProportion' => 0.01,
                                                                'maxMatchedNormalAlleleProportion' => 0.08,
                                                                'readPosBeginningOfReadIgnoreProportion' => 0.03,
                                                                'readPosTwoThirdsOfReadExtendProportion' => 0.07,
                                                                'minRdPosDepth' => 10,
                                                                'matchedNormalMaxMutProportion' => 0.5,
                                                                'minSingleEndCoverage' => 5,
                                                                'minGapFlagDistEndOfReadPercent' => 65,
                                                                'maxGapFlagDistFromEndOfReadProp' => 0.05,
                                                                'minMeanMapQualGapFlag' => 7,
                                                                'minGapPresentInReads' => 25,
                                                                'withinXBpOfDeletion' => 3,
                                                                'maxCavemanMatchedNormalProportion' => 0.01]);

  ok($processor->maxGapFlagDistFromEndOfReadProp == 0.05, 'maxGapFlagDistFromEndOfReadProp');
  ok($processor->meanMapQualGapFlag == 7, 'minMeanMapQualGapFlag');
  ok($processor->minGapPresentInPercentReads == 25, 'minGapPresentInPercentReads');
  ok($processor->withinXBpOfDeletion == 3, 'withinXBpOfDeletion');
  ok($processor->maxCavemanMatchedNormalProportion == 0.01, 'maxCavemanMatchedNormalProportion');
  ok($processor->minGapFlagDistEndOfReadPercent == 65, 'minGapFlagDistEndOfReadPercent');
  ok($processor->minSingleEndCoverage == 5, 'minSingleEndCoverage');
  ok($processor->matchedNormalMaxMutProportion == 0.5, 'minUnmatchNormMutAlleleProportion');
  ok($processor->minDepthQual == 2, "Min depth qual");
  ok($processor->depthCutoffProportion == (1/2), "depthCutoffProportion");
  ok($processor->minNormalMutAlleleQual == 3, "Min normal mut allele qual");
  ok($processor->minAnalysedQual == 5, "Min analysed qualities");
  ok($processor->percentageSamePos == 8, "Same position max pct");
  ok($processor->maxTumIndelProportion == 9, "Max tumour indel proportion");
  ok($processor->pentamerMinPassAvgQual == 11, "Pentamer minimum pass average quality");
  ok($processor->maxNormIndelProportion == 12, "Max normal indel proportion");
  ok($processor->minPassAvgBaseQualPhasing == 13, "Min pass phase ");
  ok($processor->minPassAvgMapQual == 14, "Min pass avg map quality");
  ok($processor->keepSW == 1, "Keep smith waterman reads on");
  ok($processor->maxPhasingMinorityStrandReadProportion == 0.01, "maxPhasingMinorityStrandReadProportion reset");
  ok($processor->maxMatchedNormalAlleleProportion == 0.08, "maxMatchedNormalAlleleProportion reset");
  ok($processor->readPosBeginningOfReadIgnoreProportion == 0.03, "readPosBeginningOfReadIgnoreProportion reset");
  ok($processor->readPosTwoThirdsOfReadExtendProportion == 0.07, "readPosTwoThirdsOfReadExtendProportion reset");
  ok($processor->minRdPosDepth == 10, "Min rd pos depth changed");

  $processor = undef;
  $processor = new_ok('Sanger::CGP::CavemanPostProcessor::MNVPostProcessor' => [
                                      'tumBam' => $T_BAM,
                                      'normBam' => $T_BAM]);
  isa_ok($processor->tumBam(), "Bio::DB::HTS", "Test tumour bam");
  isa_ok($processor->normBam(), "Bio::DB::HTS", "Test normal bam");
  ok($processor->minDepthQual == 25,"Min depth qual");
  ok($processor->depthCutoffProportion == (0.333333), "depthCutoffProportion");
  ok($processor->minNormalMutAlleleQual == 15,"Min normal mut allele qual");
  ok($processor->minAnalysedQual == 11,"Min analysed qualities");
  ok($processor->percentageSamePos == 80,"Same position max pct");
  ok($processor->maxTumIndelProportion == 10,"Max tumour indel proportion");
  ok($processor->pentamerMinPassAvgQual == 20,"Pentamer minimum pass average quality");
  ok($processor->maxNormIndelProportion == 10,"Max normal indel proportion");
  ok($processor->minPassAvgBaseQualPhasing == 21,"Min pass phase ");
  ok($processor->minPassAvgMapQual == 21,"Min pass avg map quality");
  ok($processor->keepSW == 0,"Keep SW off by default");
  ok($processor->maxPhasingMinorityStrandReadProportion == 0.04,"maxPhasingMinorityStrandReadProportion");
  ok($processor->maxMatchedNormalAlleleProportion == 0.05,"maxMatchedNormalAlleleProportion");
  ok($processor->readPosBeginningOfReadIgnoreProportion == 0.08,"readPosBeginningOfReadIgnoreProportion");
  ok($processor->readPosTwoThirdsOfReadExtendProportion == 0.08,"readPosTwoThirdsOfReadExtendProportion");
  ok($processor->minRdPosDepth == 8,"Min rd pos depth default");
  ok($processor->matchedNormalMaxMutProportion == 0.2,'minUnmatchNormMutAlleleProportion');
  ok($processor->minSingleEndCoverage == 10,'minSingleEndCoverage');
  ok($processor->minGapFlagDistEndOfReadPercent == 75,'minGapFlagDistEndOfReadPercent');
  ok($processor->maxGapFlagDistFromEndOfReadProp == 0.13, 'maxGapFlagDistFromEndOfReadProp');
  ok($processor->meanMapQualGapFlag == 10, 'meanMapQualGapFlag');
  ok($processor->minGapPresentInPercentReads == 30, 'minGapPresentInPercentReads');
  ok($processor->withinXBpOfDeletion == 10, 'withinXBpOfDeletion');
  ok($processor->maxCavemanMatchedNormalProportion == 0.2, 'maxCavemanMatchedNormalProportion');
  done_testing();
};

subtest '_isCurrentPosCoveredFromAlignment' => sub {
  my @cig_array = ( ['M', 10], ['D', 2], ['M', 5], ['S', 5],);
  my $pos = 0;# 0 based leftmost position of read
  my $currentPos_start = 9;
  my $currentPos_end = 10;
  my $res = Sanger::CGP::CavemanPostProcessor::MNVPostProcessor::_isCurrentPosCoveredFromAlignment(
                                                      $pos, 
                                                      \@cig_array, 
                                                      $currentPos_start, 
                                                      $currentPos_end);
  ok($res==1, "Position is covered 1 before deletion $res != 1");
  $currentPos_start = 10;
  $currentPos_end = 11;
  $res = Sanger::CGP::CavemanPostProcessor::MNVPostProcessor::_isCurrentPosCoveredFromAlignment(
                                                      $pos, 
                                                      \@cig_array, 
                                                      $currentPos_start, 
                                                      $currentPos_end);
  ok($res==-1, "Position is not covered - overlaps deletion $res != 0");
  $currentPos_start = 11;
  $currentPos_end = 12;
  $res = Sanger::CGP::CavemanPostProcessor::MNVPostProcessor::_isCurrentPosCoveredFromAlignment(
                                                      $pos, 
                                                      \@cig_array, 
                                                      $currentPos_start, 
                                                      $currentPos_end);
  ok($res==-1, "Position is not covered in deletion $res != 0");
  $currentPos_start = 17;
  $currentPos_end = 18;
  $res = Sanger::CGP::CavemanPostProcessor::MNVPostProcessor::_isCurrentPosCoveredFromAlignment(
                                                      $pos, 
                                                      \@cig_array, 
                                                      $currentPos_start, 
                                                      $currentPos_end);
  ok($res==0, "Position is not covered in skipped region (soft clipped) $res != 0");
  done_testing();
};
