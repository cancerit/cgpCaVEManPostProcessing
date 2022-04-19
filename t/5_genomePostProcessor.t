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
use Sanger::CGP::CavemanPostProcessor::GenomePostProcessor;
use Data::Dumper;
use Const::Fast qw(const);

use Test::More tests => 8;

use FindBin qw($Bin);
my $lib_path = "$Bin/../lib";
my $test_data_path = "$Bin/../testData/";

const my $T_BAM => $test_data_path.'test.bam';
const my $T_BAI => $test_data_path.'test.bam.bai';

#----------------
#  Init tests for parent modules
#----------------
subtest 'Initialise module (no params)' => sub {
  eval{Sanger::CGP::CavemanPostProcessor::GenomePostProcessor->new();};
  ok(defined($@),"Defined error from eval");
  ok($@ =~ m/tumBam and normBam are required for initialisation/, "Correct error thrown");
  done_testing();
};

subtest 'Initialise module (bam params)' => sub {
  my $processor = new_ok('Sanger::CGP::CavemanPostProcessor::GenomePostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
  isa_ok($processor->tumBam(), "Bio::DB::HTS", "Test tumour bam");
  isa_ok($processor->normBam(), "Bio::DB::HTS", "Test normal bam");
  ok($processor->minAnalysedQual == 11,"Min analysed qualities");
  ok($processor->keepSW == 0,"Keep SW off by default");
  done_testing();
};

subtest 'Min analysed qual and keep SW getters/setters' => sub{
  my $processor = new_ok('Sanger::CGP::CavemanPostProcessor::GenomePostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
  ok($processor->minAnalysedQual == 11,"Min analysed qualities");
  ok($processor->keepSW == 0,"Keep SW off by default");
  $processor->minAnalysedQual(0);
  $processor->keepSW(1);
  ok($processor->minAnalysedQual == 0,"Min analysed quality changed");
  $processor->minAnalysedQual(11);
  ok($processor->minAnalysedQual == 11,"Min analysed qualities");
  ok($processor->keepSW == 1,"Keep SW on after switch");
  $processor->keepSW(0);
  ok($processor->keepSW == 0,"Keep SW off by default");
  done_testing();
};

subtest 'Initialise module (bam params)' => sub {
  my $processor = new_ok('Sanger::CGP::CavemanPostProcessor::GenomePostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
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
  done_testing();
};

subtest 'Test runProcess & related storage methods' => sub{
  my $processor = new_ok('Sanger::CGP::CavemanPostProcessor::GenomePostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
  my $chr = 1;
  my $pos = 10011533;
  my $ref = "G";
  my $mut = "T";
  $processor->runProcess($chr,$pos,$pos,$ref,$mut);
  ok($processor->_chromosome() eq $chr,"Chromosome correct");
  ok($processor->_currentPos() == $pos,"Current pos updated");
  ok($processor->_refBase() eq $ref,"Ref base changed");
  ok($processor->_mutBase() eq $mut,"Mut base changed");
  done_testing();
};

subtest 'Initialise module (ALL params)' => sub {
  my $processor = new_ok('Sanger::CGP::CavemanPostProcessor::GenomePostProcessor' => [
                                      'tumBam' => $T_BAM,
                                      'normBam' => $T_BAM,
                                      'minDepthQual' => 2,
                                      'depthCutoffProportion' => 0.5,
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
                                      'matchedNormalMaxMutProportion' => 0.3,
                                      'minSingleEndCoverage' => 12]);

  ok($processor->minDepthQual == 2,"Min depth qual");
  ok($processor->depthCutoffProportion == 0.5, "depthCutoffProportion");
  ok($processor->minNormalMutAlleleQual == 3,"Min normal mut allele qual");
  ok($processor->minAnalysedQual == 5,"Min analysed qualities");
  ok($processor->percentageSamePos == 8,"Same position max pct");
  ok($processor->maxTumIndelProportion == 9,"Max tumour indel proportion");
  ok($processor->pentamerMinPassAvgQual == 11,"Pentamer minimum pass average quality");
  ok($processor->maxNormIndelProportion == 12,"Max normal indel proportion");
  ok($processor->minPassAvgBaseQualPhasing == 13,"Min pass phase ");
  ok($processor->minPassAvgMapQual == 14,"Min pass avg map quality");
  ok($processor->keepSW == 1,"Keep smith waterman reads on");
  ok($processor->maxPhasingMinorityStrandReadProportion == 0.01,"maxPhasingMinorityStrandReadProportion reset");
  ok($processor->maxMatchedNormalAlleleProportion == 0.08,"maxMatchedNormalAlleleProportion reset");
  ok($processor->readPosBeginningOfReadIgnoreProportion == 0.03,"readPosBeginningOfReadIgnoreProportion reset");
  ok($processor->readPosTwoThirdsOfReadExtendProportion == 0.07,"readPosTwoThirdsOfReadExtendProportion reset");
  ok($processor->minRdPosDepth == 10,"Min rd pos depth changed");
  ok($processor->matchedNormalMaxMutProportion == 0.3,"matchedNormalMaxMutProportion changed");
  ok($processor->minSingleEndCoverage == 12, "minSingleEndCoverage changed");
  $processor = new_ok('Sanger::CGP::CavemanPostProcessor::GenomePostProcessor' => [
                                      'tumBam' => $T_BAM,
                                      'normBam' => $T_BAM]);
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
  ok($processor->matchedNormalMaxMutProportion == 0.2,"matchedNormalMaxMutProportion default");
  ok($processor->minSingleEndCoverage == 10,"minSingleEndCoverage default");
  done_testing();
};

#----------------
#  Getter setter tests
#----------------
subtest 'Test all getters/setters' => sub {
  my $processor = new_ok('Sanger::CGP::CavemanPostProcessor::GenomePostProcessor' => [
                                      'tumBam' => $T_BAM,
                                      'normBam' => $T_BAM]);

  ok($processor->pentamerMinPassAvgQual == 20,"Pentamer minimum pass average quality");
  $processor->pentamerMinPassAvgQual(2);
  ok($processor->pentamerMinPassAvgQual == 2,"Pentamer minimum pass average quality changed");
  ok($processor->percentageSamePos == 80,"Same position max pct");
  $processor->percentageSamePos(65);
  ok($processor->percentageSamePos == 65,"Same position max pct changed");
  ok($processor->maxTumIndelProportion == 10,"Max tumour indel proportion");
  $processor->maxTumIndelProportion(5);
  ok($processor->maxTumIndelProportion == 5,"Max tumour indel proportion changed");
  ok($processor->maxNormIndelProportion == 10,"Max normal indel proportion");
  $processor->maxNormIndelProportion(6);
  ok($processor->maxNormIndelProportion == 6,"Max normal indel proportion changed");
  ok($processor->minPassAvgMapQual == 21,"Min pass avg map quality");
  $processor->minPassAvgMapQual(23);
  ok($processor->minPassAvgMapQual == 23,"Min pass avg map quality changed");
  ok($processor->minPassAvgBaseQualPhasing == 21,"Min pass phase ");
  $processor->minPassAvgBaseQualPhasing(24);
  ok($processor->minPassAvgBaseQualPhasing == 24,"Min pass phase changed");
  ok($processor->minDepthQual == 25,"Min depth qual");
  $processor->minDepthQual(17);
  ok($processor->minDepthQual == 17,"Min depth qual changed");
    ok($processor->depthCutoffProportion == (0.333333),"depthCutoffProportion");
  $processor->depthCutoffProportion(1/2);
  ok($processor->depthCutoffProportion == (1/2),"depthCutoffProportion changed");
  ok($processor->minNormalMutAlleleQual == 15,"Min normal mut allele qual");
  $processor->minNormalMutAlleleQual(19);
  ok($processor->minNormalMutAlleleQual == 19,"Min normal mut allele qual change");
  ok($processor->maxPhasingMinorityStrandReadProportion == 0.04,"maxPhasingMinorityStrandReadProportion");
  $processor->maxPhasingMinorityStrandReadProportion(0.01);
  ok($processor->maxPhasingMinorityStrandReadProportion == 0.01,"maxPhasingMinorityStrandReadProportion");
  ok($processor->maxMatchedNormalAlleleProportion == 0.05,"maxMatchedNormalAlleleProportion");
  $processor->maxMatchedNormalAlleleProportion(0.02);
  ok($processor->maxMatchedNormalAlleleProportion == 0.02,"maxMatchedNormalAlleleProportion");
  ok($processor->readPosBeginningOfReadIgnoreProportion == 0.08,"readPosBeginningOfReadIgnoreProportion");
  $processor->readPosBeginningOfReadIgnoreProportion(0.03);
  ok($processor->readPosBeginningOfReadIgnoreProportion == 0.03,"readPosBeginningOfReadIgnoreProportion");
  ok($processor->readPosTwoThirdsOfReadExtendProportion == 0.08,"readPosTwoThirdsOfReadExtendProportion");
  $processor->readPosTwoThirdsOfReadExtendProportion(0.04);
  ok($processor->readPosTwoThirdsOfReadExtendProportion == 0.04,"readPosTwoThirdsOfReadExtendProportion");
  ok($processor->minRdPosDepth == 8,"Min rd pos depth default");
  $processor->minRdPosDepth(6);
  ok($processor->minRdPosDepth == 6,"Min rd pos depth changed");
  ok($processor->matchedNormalMaxMutProportion == 0.2,"matchedNormalMaxMutProportion default");
  $processor->matchedNormalMaxMutProportion(0.4);
  ok($processor->matchedNormalMaxMutProportion == 0.4,"matchedNormalMaxMutProportion changed");
  ok($processor->minSingleEndCoverage == 10,"minSingleEndCoverage default");
  $processor->minSingleEndCoverage(12);
  ok($processor->minSingleEndCoverage == 12,"minSingleEndCoverage changed");
  done_testing();
};

#----------------
#  Genome specific post processing tests
#----------------
subtest 'Test genome specific post processing subroutines' => sub {
  my $processor = new_ok('Sanger::CGP::CavemanPostProcessor::GenomePostProcessor' => [
                                      'tumBam' => $T_BAM,
                                      'normBam' => $T_BAM]);

  $processor->_muts->{'nqs'} = [12,11,20,21];
  $processor->_muts->{'totalNCoverage'} = 50;
  ok($processor->getNormMutsAllelesResult == 1,"Pass matched normal");
  $processor->clearResults();
  $processor->_muts->{'totalNCoverage'} = 70;
  $processor->_muts->{'nqs'} = [12,11,20,21];
  ok($processor->getNormMutsAllelesResult == 1,"Pass matched normal");
  $processor->clearResults();
  $processor->_muts->{'nqs'} = [12,11,20,21,22];
  $processor->_muts->{'totalNCoverage'} = 70;
  ok($processor->getNormMutsAllelesResult == 0,"Fail matched normal");
  $processor->clearResults();
  $processor->_muts->{'nqs'} = [];
  $processor->_muts->{'tqs'} = [12,15,20,25,40];
  $processor->_muts->{'normcvg'} = 20;
  $processor->_muts->{'tumcvg'} = 20;
  ok($processor->getMatchedNormalProportionResult == 1,"Easy pass matched normal proportion");
  $processor->clearResults();
  $processor->_muts->{'nqs'} = [40,45];
  $processor->_muts->{'tqs'} = [12,15,20,25,40];
  $processor->_muts->{'normcvg'} = 50;
  $processor->_muts->{'tumcvg'} = 20;
  ok($processor->getMatchedNormalProportionResult == 1,"Close pass matched normal proportion");
  $processor->clearResults();
  $processor->clearResults();
  $processor->_muts->{'nqs'} = [12,15,20,25,40];
  $processor->_muts->{'tqs'} = [12,15,20,25,40];
  $processor->_muts->{'normcvg'} = 20;
  $processor->_muts->{'tumcvg'} = 20;
  ok($processor->getMatchedNormalProportionResult == 0,"Fail matched normal proportion");
  $processor->clearResults();
  $processor->_muts->{'pcvg'} = 9;
  $processor->_muts->{'ncvg'} = 12;
  $processor->_muts->{'tstr'} = [1,-1];
  ok($processor->getSingleEndResult() == 1,"Single end pass low coverage +ve strand");
  $processor->clearResults();
  $processor->_muts->{'pcvg'} = 12;
  $processor->_muts->{'ncvg'} = 9;
  $processor->_muts->{'tstr'} = [1,-1];
  ok($processor->getSingleEndResult() == 1,"Single end pass low coverage -ve strand");
  $processor->clearResults();
  $processor->_muts->{'pcvg'} = 12;
  $processor->_muts->{'ncvg'} = 12;
  $processor->_muts->{'tstr'} = [1,-1];
  ok($processor->getSingleEndResult() == 1,"Single end pass mut alleles present");
  $processor->clearResults();
  $processor->_muts->{'pcvg'} = 12;
  $processor->_muts->{'ncvg'} = 12;
  $processor->_muts->{'tstr'} = [+1];
  ok($processor->getSingleEndResult() == 0,"Single end fail only +ve");
  $processor->clearResults();
  $processor->_muts->{'pcvg'} = 12;
  $processor->_muts->{'ncvg'} = 12;
  $processor->_muts->{'tstr'} = [-1];
  ok($processor->getSingleEndResult() == 0,"Single end fail only -ve");
  $processor->clearResults();
  done_testing();
};
