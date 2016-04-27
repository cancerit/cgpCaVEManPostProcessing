##########LICENCE##########
# Copyright (c) 2014 Genome Research Ltd.
#
# Author: Cancer Genome Project cgpit@sanger.ac.uk
#
# This file is part of cgpCaVEManPostProcessing.
#
# cgpCaVEManPostProcessing is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
##########LICENCE##########

use strict;
use Sanger::CGP::CavemanPostProcessor::PostProcessor;
use Data::Dumper;
use Bio::DB::HTS;
use Const::Fast qw(const);

use Test::More tests => 20;

use FindBin qw($Bin);
my $lib_path = "$Bin/../lib";
my $test_data_path = "$Bin/../testData/";

const my $T_BAM => $test_data_path.'test.bam';
const my $T_BAI => $test_data_path.'test.bam.bai';

const my $PENT_BAM => $test_data_path.'pent.bam';
const my $PENT_BAI => $test_data_path.'pent.bam.bai';

const my $CLIP_M_BAM => $test_data_path.'clip.m.bam';
const my $CLIP_M_BAI => $test_data_path.'clip.m.bam.bai';
const my $CLIP_N_BAM => $test_data_path.'clip.n.bam';
const my $CLIP_N_BAI => $test_data_path.'clip.n.bam.bai';

my $processor;


#----------------
#	Init tests
#----------------

subtest 'Initialise module (no params)' => sub {
  eval{Sanger::CGP::CavemanPostProcessor::PostProcessor->new();};
  ok(defined($@),"Defined error from eval");
  ok($@ =~ m/tumBam and normBam are required for initialisation/, "Correct error thrown");
  done_testing();
};

subtest 'Initialise module (bam params)' => sub {
	my $processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
	isa_ok($processor->tumBam(), "Bio::DB::HTS", "Test tumour bam");
	isa_ok($processor->normBam(), "Bio::DB::HTS", "Test normal bam");
	ok($processor->minAnalysedQual == 11,"Min analysed qualities");
	ok($processor->keepSW == 0,"Keep SW off by default");
  done_testing();
};

subtest 'Min analysed qual and keep SW getters/setters' => sub{
	my $processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
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
	my $processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
	isa_ok($processor->tumBam(), "Bio::DB::HTS", "Test tumour bam");
	isa_ok($processor->normBam(), "Bio::DB::HTS", "Test normal bam");
	ok($processor->minDepthQual == 25,"Min depth qual");
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
  done_testing();
};

subtest 'Test runProcess & related storage methods' => sub{
	my $processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
	my $chr = 1;
	my $pos = 10011533;
	my $ref = "G";
	my $mut = "T";
	$processor->runProcess($chr,$pos,$pos,$ref,$mut);
	ok($processor->_chromosome eq $chr,"Chromosome correct");
	ok($processor->_currentPos == $pos,"Current pos updated");
	ok($processor->_refBase eq $ref,"Ref base changed");
	ok($processor->_mutBase eq $mut,"Mut base changed");
	done_testing();
};

subtest 'Initialise module (ALL params)' => sub {
	my $processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [
																			'tumBam' => $T_BAM,
																			'normBam' => $T_BAM,
																			'minDepthQual' => 2,
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
																			'minRdPosDepth' => 10]);

	ok($processor->minDepthQual == 2,"Min depth qual");
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
	$processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [
																			'tumBam' => $T_BAM,
																			'normBam' => $T_BAM]);
	ok($processor->minDepthQual == 25,"Min depth qual");
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
	done_testing();
};

#----------------
#	Getter setter tests
#----------------
subtest 'Test all getters/setters' => sub {
	my $processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [
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
	done_testing();
};

#--------------------------------------
#	Post processing utility method tests
#--------------------------------------
subtest '_calcualteMeanBaseQualAfterMotif' => sub {
	my $processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [
																			'tumBam' => $T_BAM,
																			'normBam' => $T_BAM]);
	my $startPos = 2;
	my @quals = (2,3,3,3,3,3,3,3,3,3);
	ok($processor->_calcualteMeanBaseQualAfterMotif($startPos,\@quals) == 3,"Test avg base qual 3");
	$startPos = 1;
	ok($processor->_calcualteMeanBaseQualAfterMotif($startPos,\@quals) == 2.9,"Test avg base qual 2.9");
	done_testing();
};

#---------------------------------------------------------------------------------
#	Post processing filter tests - fake data. Real data is used in the DEV pipeline.
#---------------------------------------------------------------------------------


subtest 'getTumIndelReadDepthResult' => sub{
	my $processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
	$processor->_muts->{'indelTCount'} = 9;
	$processor->_muts->{'totalTCoverage'} = 100;
	#Will pass
	ok($processor->getTumIndelReadDepthResult == 1,"Pass tum indel filter 9%");
	#borderline pass
	$processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
	$processor->_muts->{'indelTCount'} = 10;
	$processor->_muts->{'totalTCoverage'} = 100;
	ok($processor->getTumIndelReadDepthResult == 1,"Pass (borderline) tum indel filter 10%");
	#Will fail
	$processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
	$processor->_muts->{'indelTCount'} = 11;
	$processor->_muts->{'totalTCoverage'} = 100;
	ok($processor->getTumIndelReadDepthResult == 0,"Fail tum indel filter 11%");

	#Check clear results and change of cutoff works.
	$processor->clearResults;
	$processor->maxTumIndelProportion(11);
	ok($processor->getTumIndelReadDepthResult == 1,"Pass tum indel filter 11% changed proportion");
	done_testing();
};

subtest 'getNormIndelReadDepthResult' => sub{
	my $processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
	$processor->_muts->{'indelNCount'} = 9;
	$processor->_muts->{'totalNCoverage'} = 100;
	#Will pass
	ok($processor->getNormIndelReadDepthResult == 1,"Pass norm indel filter 9%");
	#borderline pass
	$processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
	$processor->_muts->{'indelNCount'} = 10;
	$processor->_muts->{'totalNCoverage'} = 100;
	ok($processor->getNormIndelReadDepthResult == 1,"Pass (borderline) norm indel filter 10%");
	#Will fail
	$processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
	$processor->_muts->{'indelNCount'} = 11;
	$processor->_muts->{'totalNCoverage'} = 100;
	ok($processor->getNormIndelReadDepthResult == 0,"Fail norm indel filter 11%");

	#Check clear results and change of cutoff works.
	$processor->clearResults;
	$processor->maxNormIndelProportion(11);
	ok($processor->getTumIndelReadDepthResult == 1,"Pass norm indel filter 11% changed proportion");
	done_testing();
};

subtest 'getPhasingResult' => sub{
	my $processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
	$processor->_muts->{'allTumBaseQuals'} = [9,10,17,25,12,13];
	$processor->_muts->{'allTumStrands'} = [1,1,1,1,1,1];
	#FWD fail
	ok($processor->getPhasingResult == 0, "Fail fwd phasing");

	$processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
	$processor->_muts->{'allTumBaseQuals'} = [9,10,17,25,12,13];
	$processor->_muts->{'allTumStrands'} = [-1,-1,-1,-1,-1,-1];
	#REV fail
	ok($processor->getPhasingResult == 0, "Fail rev phasing");

	$processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
	$processor->_muts->{'allTumBaseQuals'} = [21,21,21,21,21,21];
	$processor->_muts->{'allTumStrands'} = [1,1,1,1,1,1];
	#PASS FWD
	ok($processor->getPhasingResult == 1, "Pass fwd phasing");

	$processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
	$processor->_muts->{'allTumBaseQuals'} = [21,21,21,21,21,21];
	$processor->_muts->{'allTumStrands'} = [-1,-1,-1,-1,-1,-1];
	#PASS REV
	ok($processor->getPhasingResult == 1, "Pass rev phasing");

	#Check clear results and change of cutoff works.
	$processor->clearResults;

	$processor->minPassAvgBaseQualPhasing(22);
	#PASS REV
	ok($processor->getPhasingResult == 0, "Fail rev changed av BQ phasing");

	#FWD and REV with opposite strand below and above new cutoffs.

	#change strand %age and check (and clear results).

	done_testing();
};

subtest 'getDepthResult' => sub{
	#Pass
	my $processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
	$processor->_muts->{'tqs'} = [9,8,25,25,25,25,25,25,25];
	ok($processor->getDepthResult == 1,"Pass depth check");
	#Fail
	$processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
	$processor->_muts->{'tqs'} = [9,8,10,21,21,21,21,21,21];
	ok($processor->getDepthResult == 0,"Fail depth check");
	#Change minDepthQual
	$processor->clearResults;
	$processor->minDepthQual(8);
	ok($processor->getDepthResult == 1,"Pass depth check, changed min depth quality.");
	done_testing();
};

subtest 'getDifferingReadPositionResult' => sub{
	#Fail same read pos
	my $processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
	$processor->_muts->{'trp'} = [51,75,75,75,75,75,75,75,75,75];
	ok($processor->getDifferingReadPositionResult == 0,"Fail same read pos");
	#Pass same read pos
	$processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
	$processor->_muts->{'trp'} = [51,45,55,22,75,75,75,75,75,75];
	ok($processor->getDifferingReadPositionResult == 1,"Pass same read pos");

	#Pass same read pos change % same pos
	$processor->clearResults;
	$processor->percentageSamePos(50);
	ok($processor->getDifferingReadPositionResult == 0,"Fail same read pos, changed %age");
	done_testing();
};

subtest 'getAvgMapQualResult' => sub{
	#Fail (21)
	my $processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
	$processor->_muts->{'tmq'} = [20,19,17,21];
	ok($processor->getAvgMapQualResult == 0,"Fail avg map qual");
	#Pass (21)
	$processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
	$processor->_muts->{'tmq'} = [21,22,23,24];
	ok($processor->getAvgMapQualResult == 1,"Pass avg map qual");
	#Fail - changed avg map qual pass
	$processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM,'minPassAvgMapQual' => 30]);
	$processor->_muts->{'tmq'} = [21,22,23,24,30];
	ok($processor->getAvgMapQualResult == 0,"Fail avg map qual - 30");
	$processor->minPassAvgMapQual(15);
	$processor->clearResults();
	$processor->_muts->{'tmq'} = [10,12,15,16];
	ok($processor->getAvgMapQualResult == 0,"Fail avg map qual - 15");
	#Pass - changed avg map qual pass
	$processor->minPassAvgMapQual(13);
	$processor->clearResults();
	ok($processor->getAvgMapQualResult == 1,"Pass avg map qual changed & clear results");
	done_testing();
};


subtest 'getReadPositionResult' => sub{
	#Fail, pos 1 allowed
	my $processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
	$processor->_muts->{'trp'} = [75,75,75,75,75,75,75,75];
	$processor->_muts->{'trl'} = [75,75,75,75,75,75,75,75];
	ok($processor->getReadPositionResult == 0,"Fail rd pos, 1st 8% not allowed, extra 8% - end of read, depth 8");
	#Pass as depth > 8
	$processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
	$processor->_muts->{'trp'} = [75,75,75,75,75,75,75,75,75];
	$processor->_muts->{'trl'} = [75,75,75,75,75,75,75,75,75];
	ok($processor->getReadPositionResult == 1,"Pass rd pos, depth > 8");
	#Fail after pos 8% only
	$processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
	$processor->_muts->{'trp'} = [8,75,75,75,75,75,75,75];
	$processor->_muts->{'trl'} = [75,75,75,75,75,75,75,75];
	ok($processor->getReadPositionResult == 1,"Pass rd pos, 1st 8% not allowed, extra 8%, depth 8");
	$processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
	$processor->_muts->{'trp'} = [1,75,75,75,75,75,75,75];
	$processor->_muts->{'trl'} = [75,75,75,75,75,75,75,75];
	ok($processor->getReadPositionResult == 0,"Fail rd pos, 1st 8% not allowed, extra 8% - beginning of read, depth 8");
	#Change %age each side.
	##Fail 1st 8% excluded.
	$processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
	$processor->_muts->{'trp'} = [5,75,75,75,75,75,75,75];
	$processor->_muts->{'trl'} = [75,75,75,75,75,75,75,75];
	ok($processor->getReadPositionResult == 0,"Fail rd pos, 1st 8% not allowed, extra 8%, depth 8");
	#pass 1st 6% excluded.
	$processor->readPosBeginningOfReadIgnoreProportion(0.06);
	$processor->clearResults();
	$processor->_muts->{'trp'} = [5,75,75,75,75,75,75,75];
	$processor->_muts->{'trl'} = [75,75,75,75,75,75,75,75];
	ok($processor->getReadPositionResult == 1,"Pass rd pos, 1st 6% not allowed, extra 8%, depth 8");
	#Pass extra 8% included.
	$processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
	$processor->_muts->{'trp'} = [56,75,75,75,75,75,75,75];
	$processor->_muts->{'trl'} = [75,75,75,75,75,75,75,75];
	ok($processor->getReadPositionResult == 1,"Pass rd pos, 1st 8% not allowed, extra 8% (6 bases), depth 8");
	#Fail extra 6% included.
	$processor->clearResults();
	$processor->readPosTwoThirdsOfReadExtendProportion(0.06);
	$processor->_muts->{'trp'} = [56,75,75,75,75,75,75,75];
	$processor->_muts->{'trl'} = [75,75,75,75,75,75,75,75];
	ok($processor->getReadPositionResult == 0,"Fail rd pos, 1st 8% not allowed, extra 6% (4.5 bases), depth 8");
	done_testing();
};

subtest 'getNormMutsAllelesResult' => sub{
	#Fail as more than proportion
	my $processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
	$processor->_muts->{'nqs'} = [1,3,8,10,21,21];
	$processor->_muts->{'totalNCoverage'} = 30;
	ok($processor->getNormMutsAllelesResult == 0,"Fail matched normal mut allele check");
	#Pass as less than proportion
	$processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM]);
	$processor->_muts->{'nqs'} = [1,3,8,10,21,12];
	$processor->_muts->{'totalNCoverage'} = 41;
	ok($processor->getNormMutsAllelesResult == 1,"Pass matched normal mut allele check");
	#Fail minNormalMutAlleleQual
	$processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM, 'minNormMutAllelequal' => 25]);
	$processor->_muts->{'nqs'} = [1,3,8,10,25,25];
	$processor->_muts->{'totalNCoverage'} = 30;
	ok($processor->getNormMutsAllelesResult == 0,"Fail matched normal mut allele check, minNormalMutAlleleQual = 25");
	#Pass minNormalMutAlleleQual
	$processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $T_BAM, normBam => $T_BAM, 'minNormMutAlleleQual' => 25]);
	$processor->_muts->{'nqs'} = [1,3,8,10,25,25];
	$processor->_muts->{'totalNCoverage'} = 30;
	ok($processor->getNormMutsAllelesResult == 0,"Fail matched normal mut allele check, minNormalMutAlleleQual = 25");
	$processor->minNormalMutAlleleQual(26);
	$processor->clearResults();
	$processor->_muts->{'nqs'} = [1,3,8,10,25,25];
	$processor->_muts->{'totalNCoverage'} = 30;
	ok($processor->getNormMutsAllelesResult == 1,"Fail matched normal mut allele check, minNormalMutAlleleQual = 25, clear results");
	done_testing();
};

subtest 'getPentamerResult' => sub{
	my $processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => ['tumBam' => $PENT_BAM, 'normBam' => $PENT_BAM]);


	#--Fails
	#6	138186703	A > G
	$processor->runProcess('6',138186703,138186703,"A","G");
	ok($processor->getPentamerResult == 0,"Fail pentamer check, 6:138186703");
	$processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => ['tumBam' => $PENT_BAM, 'normBam' => $PENT_BAM]);
	$processor->pentamerMinPassAvgQual(16);
	$processor->runProcess('6',138186703,138186703,"A","G");
	ok($processor->getPentamerResult == 1,"Pass pentamer check, 6:138186703, change mean base qual");

	#--Passes
	#X	48092088	C > T
	$processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => ['tumBam' => $PENT_BAM, 'normBam' => $PENT_BAM]);
	$processor->runProcess('X',48092088,48092088,"C","T");
	ok($processor->getPentamerResult == 1,"Pass pentamer check, X:48092088");
	done_testing();
};

subtest 'Initialise module (bam clip params)' => sub {
	my $processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $CLIP_M_BAM, normBam => $CLIP_N_BAM]);
	isa_ok($processor->tumBam(), "Bio::DB::HTS", "Test tumour bam");
	isa_ok($processor->normBam(), "Bio::DB::HTS", "Test normal bam");
	ok($processor->minAnalysedQual == 11,"Min analysed qualities");
	ok($processor->keepSW == 0,"Keep SW off by default");
  done_testing();
};

subtest 'Clipped Read tests' => sub {
	my $processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $CLIP_M_BAM, normBam => $CLIP_N_BAM]);
	my $chr = 1;
	my $pos = 10437;
	my $ref = "T";
	my $mut = "C";
	$processor->runProcess($chr,$pos,$pos,$ref,$mut);
	ok($processor->_chromosome eq $chr,"Chromosome correct");
	ok($processor->_currentPos == $pos,"Current pos updated");
	ok($processor->_refBase eq $ref,"Ref base changed");
	ok($processor->_mutBase eq $mut,"Mut base changed");

	#Manually set counts
	my $exp_sclp = [1,2,3,4,5,6,7,8,9,10];
	$processor->_muts->{'sclp'} = $exp_sclp;
	is_deeply($processor->_muts->{'sclp'}, [1,2,3,4,5,6,7,8,9,10], "softclipcounts");

	#Check manually set count results
	#getClipMedianResult
	my $exp_res = sprintf('%.2f',5.5);
	is($processor->getClipMedianResult, $exp_res,"getClipMedianResult");

	#Reset and use real data
	$processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $CLIP_M_BAM, normBam => $CLIP_N_BAM]);
	$processor->runProcess($chr,$pos,$pos,$ref,$mut);
	ok($processor->_chromosome eq $chr,"Chromosome correct");
	ok($processor->_currentPos == $pos,"Current pos updated");
	ok($processor->_refBase eq $ref,"Ref base changed");
	ok($processor->_mutBase eq $mut,"Mut base changed");

	#Check counts have been filled correctly.
	$exp_sclp = [0,9,78,78,0,83,49,104,38,35,0,20,0,0,0,66,0,0,0];
	is_deeply($processor->_muts->{'sclp'}, [0,9,78,78,0,83,49,104,38,35,0,20,0,0,0,66,0,0,0], "softclipcounts");

	#Check real data count results
	$exp_res = sprintf('%.2f',9);
	is($processor->getClipMedianResult, $exp_res,"getClipMedianResult");
  done_testing();
};

subtest 'Alignment score tests' => sub {
	my $processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $CLIP_M_BAM, normBam => $CLIP_N_BAM]);
	my $chr = 1;
	my $pos = 10437;
	my $ref = "T";
	my $mut = "C";
	$processor->runProcess($chr,$pos,$pos,$ref,$mut);
	ok($processor->_chromosome eq $chr,"Chromosome correct");
	ok($processor->_currentPos == $pos,"Current pos updated");
	ok($processor->_refBase eq $ref,"Ref base changed");
	ok($processor->_mutBase eq $mut,"Mut base changed");

	my $exp_prim = [100,80,50];
	my $rdlen = [150,150,150];

	$processor->_muts->{'alnp'} = $exp_prim;
	$processor->_muts->{'trl'} = $rdlen;
	is_deeply($processor->_muts->{'alnp'}, [100,80,50], "primary alignment scores");
	is_deeply($processor->_muts->{'trl'}, [150,150,150], "read lengths");

	#getAlnScoreMedianReadAdjusted
	my $exp_res = sprintf('%.2f',0.53333333333333);
	is($processor->getAlignmentScoreMedianReadAdjusted, $exp_res,"getAlignmentScoreMedianReadAdjusted");
	#getAlignmentScoreMedian
	$exp_res = sprintf('%.2f',80);
	is($processor->getAlignmentScoreMedian, $exp_res,"getAlignmentScoreMedian");

	$processor = new_ok('Sanger::CGP::CavemanPostProcessor::PostProcessor' => [tumBam => $CLIP_M_BAM, normBam => $CLIP_N_BAM]);
	$chr = 1;
	$pos = 10437;
	$ref = "T";
	$mut = "C";
	$processor->runProcess($chr,$pos,$pos,$ref,$mut);
	ok($processor->_chromosome eq $chr,"Chromosome correct");
	ok($processor->_currentPos == $pos,"Current pos updated");
	ok($processor->_refBase eq $ref,"Ref base changed");
	ok($processor->_mutBase eq $mut,"Mut base changed");
	is_deeply($processor->_muts->{'alnp'}, [66,102,51,56,110,63,61,44,82,80,123,88,123,76,87,60,118,93,139], "primary alignment scores");
	is_deeply($processor->_muts->{'trl'}, [151,151,151,151,151,151,151,151,151,151,151,151,151,151,151,151,151,151,151], "tumor read lengths");

	#getAlnScoreMedianReadAdjusted
	$exp_res = sprintf('%.2f',0.543046358);
	is($processor->getAlignmentScoreMedianReadAdjusted, $exp_res,"getAlignmentScoreMedianReadAdjusted");
	#getAlignmentScoreMedian
	$exp_res = sprintf('%.2f',82);
	is($processor->getAlignmentScoreMedian, $exp_res,"getAlignmentScoreMedian");

  done_testing();
};

