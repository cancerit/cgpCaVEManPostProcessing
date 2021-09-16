##########LICENCE##########
# Copyright (c) 2014-2018 Genome Research Ltd.
#
# Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
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
use Sanger::CGP::CavemanPostProcessing::Flagger;
use Sanger::CGP::CavemanPostProcessing::Config;
use Data::Dumper;
use Const::Fast qw(const);

use Test::More tests => 20;

use FindBin qw($Bin);
const my $lib_path => "$Bin/../lib";
const my $test_data_path => "$Bin/../testData/";

const my $CONFIG => $test_data_path."flag.vcf.config.new.ini";
const my $FLAG_CFG => $test_data_path."flag.to.vcf.convert.ini";
const my $SPECIES => 'HUMAN';
const my $TYPE => 'WGS';
const my $GERMLINE_INDEL => $test_data_path.'germline_indel.bed.gz';
const my $BED_LOC => $test_data_path;

const my $T_BAM => $test_data_path.'test.bam';
const my $T_BAI => $test_data_path.'test.bam.bai';

const my $PENT_BAM => $test_data_path.'pent.bam';
const my $PENT_BAI => $test_data_path.'pent.bam.bai';

const my $CLIP_M_BAM => $test_data_path.'clip.m.bam';
const my $CLIP_M_BAI => $test_data_path.'clip.m.bam.bai';
const my $CLIP_N_BAM => $test_data_path.'clip.n.bam';
const my $CLIP_N_BAI => $test_data_path.'clip.n.bam.bai';

#----------------
#	Init tests
#----------------

subtest 'Initialise module (new)' => sub {
  my $processor = new_ok('Sanger::CGP::CavemanPostProcessing::Flagger');
  done_testing();
};

subtest 'Initialise module (init)' => sub {
  my $processor = new_ok('Sanger::CGP::CavemanPostProcessing::Flagger');
  eval{$processor->init(undef,undef,undef);};
  ok(defined($@),"Defined error from eval");
  ok($@ =~ m/Tumour HTS file wasn't defined or an existing file./, "Correct tum file error thrown");
  eval{$processor->init(undef,undef,undef);};
  eval{$processor->init($T_BAM,undef,undef);};
  ok($@ =~ m/Normal HTS file wasn't defined or an existing file./, "Correct norm file error thrown");
  my $cfg = Sanger::CGP::CavemanPostProcessing::Config->new();
  my $opts;
  $opts->{'c'} = $CONFIG;
  $opts->{'v'} = $FLAG_CFG;
  $opts->{'s'} = $SPECIES;
  $opts->{'t'} = $TYPE;
  $opts->{'g'} = $GERMLINE_INDEL;
  $opts->{'b'} = $BED_LOC;
  $cfg->init_config($opts);
  $processor->init($T_BAM,$T_BAM,$cfg);
  done_testing();
};

subtest 'set_position' => sub {
  my $processor = new_ok('Sanger::CGP::CavemanPostProcessing::Flagger');
  my $cfg = Sanger::CGP::CavemanPostProcessing::Config->new();
  my $opts;
  $opts->{'c'} = $CONFIG;
  $opts->{'v'} = $FLAG_CFG;
  $opts->{'s'} = $SPECIES;
  $opts->{'t'} = $TYPE;
  $opts->{'g'} = $GERMLINE_INDEL;
  $opts->{'b'} = $BED_LOC;
  $cfg->init_config($opts);
  $processor->init($T_BAM,$T_BAM,$cfg);

  my $chr = 1;
	my $pos = 10011533;
	my $ref = "G";
	my $mut = "T";
  eval{$processor->set_position($chr,$pos,$ref,$mut);};
  ok($@ =~ m//,"No error from eval.");
  done_testing();
};

subtest 'matchedNomralProportion' => sub {
  my $processor = new_ok('Sanger::CGP::CavemanPostProcessing::Flagger');
  my $cfg = Sanger::CGP::CavemanPostProcessing::Config->new();
  my $opts;
  $opts->{'c'} = $CONFIG;
  $opts->{'v'} = $FLAG_CFG;
  $opts->{'s'} = $SPECIES;
  $opts->{'t'} = $TYPE;
  $opts->{'g'} = $GERMLINE_INDEL;
  $opts->{'b'} = $BED_LOC;
  $cfg->init_config($opts);
  $processor->init($T_BAM,$T_BAM,$cfg);

  my $chr = 1;
	my $pos = 10011533;
	my $ref = "G";
	my $mut = "T";
  eval{$processor->set_position($chr,$pos,$ref,$mut);};
  ok($@ =~ m//,"No error from eval.");

  $processor->_muts->{'nqs'} = [];
	$processor->_muts->{'tqs'} = [12,15,20,25,40];
	$processor->_muts->{'normcvg'} = 20;
	$processor->_muts->{'tumcvg'} = 20;
	ok($processor->matchedNormalProportion == 0,"Easy pass matched normal proportion");
	$processor->_muts->{'nqs'} = [40,45];
	$processor->_muts->{'tqs'} = [12,15,20,25,40];
	$processor->_muts->{'normcvg'} = 50;
	$processor->_muts->{'tumcvg'} = 20;
	ok($processor->matchedNormalProportion == 0,"Close pass matched normal proportion");
	$processor->_muts->{'nqs'} = [12,15,20,25,40];
	$processor->_muts->{'tqs'} = [12,15,20,25,40];
	$processor->_muts->{'normcvg'} = 20;
	$processor->_muts->{'tumcvg'} = 20;
	ok($processor->matchedNormalProportion == 1,"Fail matched normal proportion");
	done_testing();
};

subtest 'avgMapQualFlag' => sub {
  my $processor = new_ok('Sanger::CGP::CavemanPostProcessing::Flagger');
  my $cfg = Sanger::CGP::CavemanPostProcessing::Config->new();
  my $opts;
  $opts->{'c'} = $CONFIG;
  $opts->{'v'} = $FLAG_CFG;
  $opts->{'s'} = $SPECIES;
  $opts->{'t'} = $TYPE;
  $opts->{'g'} = $GERMLINE_INDEL;
  $opts->{'b'} = $BED_LOC;
  $cfg->init_config($opts);
  $processor->init($T_BAM,$T_BAM,$cfg);

  my $chr = 1;
	my $pos = 10011533;
	my $ref = "G";
	my $mut = "T";
  eval{$processor->set_position($chr,$pos,$ref,$mut);};
  ok($@ =~ m//,"No error from eval.");

  $processor->_muts->{'tmq'} = [20,19,17,21];
	ok($processor->avgMapQualFlag == 1,"Fail avg map qual");
	#Pass (21)
	$processor->_muts->{'tmq'} = [21,22,23,24];
	ok($processor->avgMapQualFlag == 0,"Pass avg map qual");
	#Fail - changed avg map qual pass
	$processor->_prms->{'minPassAvgMapQual'} = 30;
	$processor->_muts->{'tmq'} = [21,22,23,24,30];
	ok($processor->avgMapQualFlag == 1,"Fail avg map qual - 30");


	$processor->_prms->{'minPassAvgMapQual'} = 15;
	$processor->_muts->{'tmq'} = [10,12,15,16];
	ok($processor->avgMapQualFlag == 1,"Fail avg map qual - 15");

	#Pass - changed avg map qual pass
	$processor->_prms->{'minPassAvgMapQual'} = 13;
	ok($processor->avgMapQualFlag == 0,"Pass avg map qual changed & clear results");

  done_testing();
};

subtest 'depthFlag' => sub{
  my $processor = new_ok('Sanger::CGP::CavemanPostProcessing::Flagger');
  my $cfg = Sanger::CGP::CavemanPostProcessing::Config->new();
  my $opts;
  $opts->{'c'} = $CONFIG;
  $opts->{'v'} = $FLAG_CFG;
  $opts->{'s'} = $SPECIES;
  $opts->{'t'} = $TYPE;
  $opts->{'g'} = $GERMLINE_INDEL;
  $opts->{'b'} = $BED_LOC;
  $cfg->init_config($opts);
  $processor->init($T_BAM,$T_BAM,$cfg);

  my $chr = 1;
	my $pos = 10011533;
	my $ref = "G";
	my $mut = "T";
  eval{$processor->set_position($chr,$pos,$ref,$mut);};
  ok($@ =~ m//,"No error from eval.");

  $processor->_muts->{'tqs'} = [9,8,25,25,25,25,25,25,25];
	ok($processor->depthFlag == 0,"Pass depth check");
	#Fail
	$processor->_muts->{'tqs'} = [9,8,10,21,21,21,21,21,21];
	ok($processor->depthFlag == 1,"Fail depth check");
	#Change depthCutoffProportion

    $processor->_muts->{'tqs'} = [9,8,10,20,20,21,21,21,21];
	$processor->_prms->{'depthCutoffProportion'} = 0.5;
	ok($processor->depthFlag == 1,"Fail depth check, changed min depthCutoffProportion.");
    $processor->_prms->{'minDepthQual'} = 21;
    $processor->_muts->{'tqs'} = [9,8,10,20,21,21,21,21,21];
    ok($processor->depthFlag == 0,"Fail depth check, changed min depthCutoffProportion.");

  my $chr = 1;
  my $pos = 10011533;
  my $ref = "G";
  my $mut = "T";
  eval{$processor->set_position($chr,$pos,$ref,$mut);};
  ok($@ =~ m//,"No error from eval.");

  $processor->_prms->{'minDepthQual'} = 25;
  $processor->_prms->{'depthCutoffProportion'} = 0.333333;
  $processor->_muts->{'tqs'} = [9,8,25,25,25,25,25,25,25];
  ok($processor->depthFlag == 0,"Pass depth check");
  #Fail
  $processor->_muts->{'tqs'} = [9,8,10,21,21,21,21,21,21];
  ok($processor->depthFlag == 1,"Fail depth check");
  #Change depthCutoffProportion

  $processor->_prms->{'minDepthQual'} = 21;
  ok($processor->depthFlag == 0,"Pass depth check, changed min depth quality.");

  $processor->_muts->{'tqs'} = [9,8,25,25,25,25,25,25,25];
  ok($processor->depthFlag == 0,"Pass depth check");
  


  done_testing();
};

subtest 'readPositionFlag' => sub {
  my $processor = new_ok('Sanger::CGP::CavemanPostProcessing::Flagger');
  my $cfg = Sanger::CGP::CavemanPostProcessing::Config->new();
  my $opts;
  $opts->{'c'} = $CONFIG;
  $opts->{'v'} = $FLAG_CFG;
  $opts->{'s'} = $SPECIES;
  $opts->{'t'} = $TYPE;
  $opts->{'g'} = $GERMLINE_INDEL;
  $opts->{'b'} = $BED_LOC;
  $cfg->init_config($opts);
  $processor->init($T_BAM,$T_BAM,$cfg);

  my $chr = 1;
	my $pos = 10011533;
	my $ref = "G";
	my $mut = "T";
  eval{$processor->set_position($chr,$pos,$ref,$mut);};
  ok($@ =~ m//,"No error from eval.");

  $processor->_muts->{'trp'} = [75,75,75,75,75,75,75,75];
	$processor->_muts->{'trl'} = [75,75,75,75,75,75,75,75];
	ok($processor->readPositionFlag == 1,"Fail rd pos, 1st 8% not allowed, extra 8% - end of read, depth 8");
	#Pass as depth > 8
	$processor->_muts->{'trp'} = [75,75,75,75,75,75,75,75,75];
	$processor->_muts->{'trl'} = [75,75,75,75,75,75,75,75,75];
	ok($processor->readPositionFlag == 0,"Pass rd pos, depth > 8");
	#Fail after pos 8% only
	$processor->_muts->{'trp'} = [8,75,75,75,75,75,75,75];
	$processor->_muts->{'trl'} = [75,75,75,75,75,75,75,75];
	ok($processor->readPositionFlag == 0,"Pass rd pos, 1st 8% not allowed, extra 8%, depth 8");
	$processor->_muts->{'trp'} = [1,75,75,75,75,75,75,75];
	$processor->_muts->{'trl'} = [75,75,75,75,75,75,75,75];
	ok($processor->readPositionFlag == 1,"Fail rd pos, 1st 8% not allowed, extra 8% - beginning of read, depth 8");
	#Change %age each side.
	##Fail 1st 8% excluded.
	$processor->_muts->{'trp'} = [5,75,75,75,75,75,75,75];
	$processor->_muts->{'trl'} = [75,75,75,75,75,75,75,75];
	ok($processor->readPositionFlag == 1,"Fail rd pos, 1st 8% not allowed, extra 8%, depth 8");
	#pass 1st 6% excluded.
	$processor->_prms->{'readPosBeginningOfReadIgnoreProportion'} = 0.06;
	$processor->_muts->{'trp'} = [5,75,75,75,75,75,75,75];
	$processor->_muts->{'trl'} = [75,75,75,75,75,75,75,75];
	ok($processor->readPositionFlag == 0,"Pass rd pos, 1st 6% not allowed, extra 8%, depth 8");
	#Pass extra 8% included.
	$processor->_muts->{'trp'} = [56,75,75,75,75,75,75,75];
	$processor->_muts->{'trl'} = [75,75,75,75,75,75,75,75];
	ok($processor->readPositionFlag == 0,"Pass rd pos, 1st 8% not allowed, extra 8% (6 bases), depth 8");
	#Fail extra 6% included.
	$processor->_prms->{'readPosTwoThirdsOfReadExtendProportion'} = 0.06;
	$processor->_muts->{'trp'} = [56,75,75,75,75,75,75,75];
	$processor->_muts->{'trl'} = [75,75,75,75,75,75,75,75];
	ok($processor->readPositionFlag == 1,"Fail rd pos, 1st 8% not allowed, extra 6% (4.5 bases), depth 8");

  done_testing();
};

subtest 'matchedNormalFlag' => sub {
  my $processor = new_ok('Sanger::CGP::CavemanPostProcessing::Flagger');
  my $cfg = Sanger::CGP::CavemanPostProcessing::Config->new();
  my $opts;
  $opts->{'c'} = $CONFIG;
  $opts->{'v'} = $FLAG_CFG;
  $opts->{'s'} = $SPECIES;
  $opts->{'t'} = $TYPE;
  $opts->{'g'} = $GERMLINE_INDEL;
  $opts->{'b'} = $BED_LOC;
  $cfg->init_config($opts);
  $processor->init($T_BAM,$T_BAM,$cfg);

  my $chr = 1;
	my $pos = 10011533;
	my $ref = "G";
	my $mut = "T";
  eval{$processor->set_position($chr,$pos,$ref,$mut);};
  ok($@ =~ m//,"No error from eval.");

  $processor->_muts->{'nqs'} = [1,3,8,10,21,21];
	$processor->_muts->{'totalNCoverage'} = 30;
	ok($processor->matchedNormalFlag == 1,"Fail matched normal mut allele check");
	#Pass as less than proportion
	$processor->_muts->{'nqs'} = [1,3,8,10,21,12];
	$processor->_muts->{'totalNCoverage'} = 41;
	ok($processor->matchedNormalFlag == 0,"Pass matched normal mut allele check");
	#Fail minNormalMutAlleleQual
	$processor->_prms->{'minNormMutAllelequal'} = 25;
	$processor->_muts->{'nqs'} = [1,3,8,10,25,25];
	$processor->_muts->{'totalNCoverage'} = 30;
	ok($processor->matchedNormalFlag == 1,"Fail matched normal mut allele check, minNormalMutAlleleQual = 25");
	#Pass minNormalMutAlleleQual
	$processor->_prms->{'minNormMutAllelequal'} = 25;
	$processor->_muts->{'nqs'} = [1,3,8,10,25,25];
	$processor->_muts->{'totalNCoverage'} = 30;
	ok($processor->matchedNormalFlag == 1,"Fail matched normal mut allele check, minNormalMutAlleleQual = 25");
	$processor->_prms->{'minNormalMutAlleleQual'} = 26;
	$processor->_muts->{'nqs'} = [1,3,8,10,25,25];
	$processor->_muts->{'totalNCoverage'} = 30;
	ok($processor->matchedNormalFlag == 1,"Fail matched normal mut allele check, minNormalMutAlleleQual = 25, clear results");

  done_testing();
};

subtest 'pentamericMotifFlag' => sub {
  my $processor = new_ok('Sanger::CGP::CavemanPostProcessing::Flagger');
  my $cfg = Sanger::CGP::CavemanPostProcessing::Config->new();
  my $opts;
  $opts->{'c'} = $CONFIG;
  $opts->{'v'} = $FLAG_CFG;
  $opts->{'s'} = $SPECIES;
  $opts->{'t'} = $TYPE;
  $opts->{'g'} = $GERMLINE_INDEL;
  $opts->{'b'} = $BED_LOC;
  $cfg->init_config($opts);
  $processor->init($PENT_BAM,$PENT_BAM,$cfg);

  #--Fails
  my $chr = '6';
	my $pos = 138186703;
	my $ref = "A";
	my $mut = "G";
	$processor->_prms->{'keepSW'} = 0;
  $processor->set_position($chr,$pos,$ref,$mut);
  ok($processor->pentamericMotifFlag == 1,"Fail pentamer check, $chr:$pos");

  $processor->_prms->{'pentamerMinPassAvgQual'} = 16;
  ok($processor->pentamericMotifFlag == 0,"Pass pentamer check, $chr:$pos, change mean base qual");

	#--Passes
	$chr = 'X';
	$pos = 48092088;
	$ref = "C";
	$mut = "T";
	eval{$processor->set_position($chr,$pos,$ref,$mut);};
  ok($@ =~ m//,"No error from eval.");
	ok($processor->pentamericMotifFlag == 0,"Pass pentamer check, $chr:$pos");

  done_testing();
};

subtest 'phasingFlag' => sub {
  my $processor = new_ok('Sanger::CGP::CavemanPostProcessing::Flagger');
  my $cfg = Sanger::CGP::CavemanPostProcessing::Config->new();
  my $opts;
  $opts->{'c'} = $CONFIG;
  $opts->{'v'} = $FLAG_CFG;
  $opts->{'s'} = $SPECIES;
  $opts->{'t'} = $TYPE;
  $opts->{'g'} = $GERMLINE_INDEL;
  $opts->{'b'} = $BED_LOC;
  $cfg->init_config($opts);
  $processor->init($T_BAM,$T_BAM,$cfg);

  $processor->_muts->{'allTumBaseQuals'} = [9,10,17,25,12,13];
	$processor->_muts->{'allTumStrands'} = [1,1,1,1,1,1];
	#FWD fail
	ok($processor->phasingFlag == 1, "Fail fwd phasing");

	$processor->_muts->{'allTumBaseQuals'} = [9,10,17,25,12,13];
	$processor->_muts->{'allTumStrands'} = [-1,-1,-1,-1,-1,-1];
	#REV fail
	ok($processor->phasingFlag == 1, "Fail rev phasing");

	$processor->_muts->{'allTumBaseQuals'} = [21,21,21,21,21,21];
	$processor->_muts->{'allTumStrands'} = [1,1,1,1,1,1];
	#PASS FWD
	ok($processor->phasingFlag == 0, "Pass fwd phasing");

	$processor->_muts->{'allTumBaseQuals'} = [21,21,21,21,21,21];
	$processor->_muts->{'allTumStrands'} = [-1,-1,-1,-1,-1,-1];
	#PASS REV
	ok($processor->phasingFlag == 0, "Pass rev phasing");

	#Check clear results and change of cutoff works.

	$processor->_prms->{'minPassPhaseQual'} = 22;
	#PASS REV
	ok($processor->phasingFlag == 1, "Fail rev changed av BQ phasing");
	#FWD and REV with opposite strand below and above new cutoffs.
	#change strand %age and check (and clear results).

  done_testing();
};

subtest 'getBEDUnmatchedFlag' => sub {
  my $processor = new_ok('Sanger::CGP::CavemanPostProcessing::Flagger');
  my $cfg = Sanger::CGP::CavemanPostProcessing::Config->new();
  my $opts;
  $opts->{'c'} = $CONFIG;
  $opts->{'v'} = $FLAG_CFG;
  $opts->{'s'} = $SPECIES;
  $opts->{'t'} = $TYPE;
  $opts->{'g'} = $GERMLINE_INDEL;
  $opts->{'b'} = $BED_LOC;
  $cfg->init_config($opts);
  $processor->init($T_BAM,$T_BAM,$cfg);

  #Pass, existing position, not allele
  my $chr = '21';
	my $pos = 45251767;
	my $ref = "A";
	my $mut = "C";
	eval{$processor->set_position($chr,$pos,$ref,$mut);};
  ok($@ =~ m//,"No error from eval.");
	my $line = "21\t45251766\t45251767\t5\tT";
  ok($processor->getBEDUnmatchedFlag($line) == 0, "Pass as C allele not present");

  #Fail
  $mut = "T";
  eval{$processor->set_position($chr,$pos,$ref,$mut);};
  ok($@ =~ m//,"No error from eval.");
  ok($processor->getBEDUnmatchedFlag($line) == 1, "Fail as allele present");

  done_testing();
};

subtest 'getVCFUnmatchedFlag' => sub {
  my $processor = new_ok('Sanger::CGP::CavemanPostProcessing::Flagger');
  my $cfg = Sanger::CGP::CavemanPostProcessing::Config->new();
  my $opts;
  $opts->{'c'} = $CONFIG;
  $opts->{'v'} = $FLAG_CFG;
  $opts->{'s'} = $SPECIES;
  $opts->{'t'} = $TYPE;
  $opts->{'g'} = $GERMLINE_INDEL;
  $opts->{'b'} = $BED_LOC;
  $cfg->init_config($opts);
  $processor->init($T_BAM,$T_BAM,$cfg);

  #Pass, existing position, not allele
  my $chr = '21';
	my $pos = 45251767;
	my $ref = "A";
	my $mut = "C";
	eval{$processor->set_position($chr,$pos,$ref,$mut);};
  ok($@ =~ m//,"No error from eval.");
	my $line = "21\t45251767\t.\tA\t.\t.\t.\tAC=0;AN=4;DP=51\tGT:AZ:CZ:GZ:TZ\t0|0:41:0:0:0\t0|0:41:0:0:0\t0|0:41:0:0:0\t0|0:41:0:0:0\t0|0:41:0:0:0\t0|0:21:0:0:20\t0|0:21:0:0:20\t0|0:21:0:0:20\t0|0:21:0:0:20\t0|0:21:0:0:20";
  ok($processor->getVCFUnmatchedFlag($line) == 0, "Pass as C allele not present");

  #Fail
  $mut = "T";
  eval{$processor->set_position($chr,$pos,$ref,$mut);};
  ok($@ =~ m//,"No error from eval.");
  ok($processor->getVCFUnmatchedFlag($line) == 1, "Fail as allele present");

  done_testing();
};

subtest 'singleEndFlag' => sub {
  my $processor = new_ok('Sanger::CGP::CavemanPostProcessing::Flagger');
  my $cfg = Sanger::CGP::CavemanPostProcessing::Config->new();
  my $opts;
  $opts->{'c'} = $CONFIG;
  $opts->{'v'} = $FLAG_CFG;
  $opts->{'s'} = $SPECIES;
  $opts->{'t'} = $TYPE;
  $opts->{'g'} = $GERMLINE_INDEL;
  $opts->{'b'} = $BED_LOC;
  $cfg->init_config($opts);
  $processor->init($T_BAM,$T_BAM,$cfg);
  my $chr = 1;
	my $pos = 10011533;
	my $ref = "G";
	my $mut = "T";
  eval{$processor->set_position($chr,$pos,$ref,$mut);};
  ok($@ =~ m//,"No error from eval.");


  $processor->_muts->{'pcvg'} = 9;
	$processor->_muts->{'ncvg'} = 12;
	$processor->_muts->{'tstr'} = [1,-1];
	ok($processor->singleEndFlag() == 0,"Single end pass low coverage +ve strand");
	$processor->_muts->{'pcvg'} = 12;
	$processor->_muts->{'ncvg'} = 9;
	$processor->_muts->{'tstr'} = [1,-1];
	ok($processor->singleEndFlag() == 0,"Single end pass low coverage -ve strand");
	$processor->_muts->{'pcvg'} = 12;
	$processor->_muts->{'ncvg'} = 12;
	$processor->_muts->{'tstr'} = [1,-1];
	ok($processor->singleEndFlag() == 0,"Single end pass mut alleles present");
	$processor->_muts->{'pcvg'} = 12;
	$processor->_muts->{'ncvg'} = 12;
	$processor->_muts->{'tstr'} = [+1];
	ok($processor->singleEndFlag() == 1,"Single end fail only +ve");
	$processor->_muts->{'pcvg'} = 12;
	$processor->_muts->{'ncvg'} = 12;
	$processor->_muts->{'tstr'} = [-1];
	ok($processor->singleEndFlag() == 1,"Single end fail only -ve");

  done_testing();
};

subtest 'alignmentScoreReadLengthAdjustedFlag' => sub {
  my $processor = new_ok('Sanger::CGP::CavemanPostProcessing::Flagger');
  my $cfg = Sanger::CGP::CavemanPostProcessing::Config->new();
  my $opts;
  $opts->{'c'} = $CONFIG;
  $opts->{'v'} = $FLAG_CFG;
  $opts->{'s'} = $SPECIES;
  $opts->{'t'} = $TYPE;
  $opts->{'g'} = $GERMLINE_INDEL;
  $opts->{'b'} = $BED_LOC;
  $cfg->init_config($opts);
  $processor->init($T_BAM,$T_BAM,$cfg);
  my $chr = 1;
	my $pos = 10011533;
	my $ref = "G";
	my $mut = "T";
  eval{$processor->set_position($chr,$pos,$ref,$mut);};
  ok($@ =~ m//,"No error from eval.");


  done_testing();
};


subtest 'Alignment score tests' => sub {
  my $processor = new_ok('Sanger::CGP::CavemanPostProcessing::Flagger');
  my $cfg = Sanger::CGP::CavemanPostProcessing::Config->new();
  my $opts;
  $opts->{'c'} = $CONFIG;
  $opts->{'v'} = $FLAG_CFG;
  $opts->{'s'} = $SPECIES;
  $opts->{'t'} = $TYPE;
  $opts->{'g'} = $GERMLINE_INDEL;
  $opts->{'b'} = $BED_LOC;
  $cfg->init_config($opts);
  $processor->init($CLIP_M_BAM,$CLIP_N_BAM,$cfg);
  my $chr = 1;
	my $pos = 10437;
	my $ref = "T";
	my $mut = "C";
	eval{$processor->set_position($chr,$pos,$ref,$mut);};
  ok($@ =~ m//,"No error from eval.");


	my $exp_prim = [100,80,50];
	my $rdlen = [150,150,150];

	$processor->_muts->{'alnp'} = $exp_prim;
	$processor->_muts->{'trl'} = $rdlen;
	is_deeply($processor->_muts->{'alnp'}, [100,80,50], "primary alignment scores");
	is_deeply($processor->_muts->{'trl'}, [150,150,150], "read lengths");

	#getAlnScoreMedianReadAdjusted
	my $exp_res = sprintf('%.2f',0.53333333333333);
	is($processor->alignmentScoreReadLengthAdjustedFlag, $exp_res,"alignmentScoreReadLengthAdjustedFlag");
	#getAlignmentScoreMedian
	$exp_res = sprintf('%.2f',80);
	is($processor->alnScoreMedianFlag, $exp_res,"alnScoreMedianFlag");

  eval{$processor->set_position($chr,$pos,$ref,$mut);};

 # warn Dumper ($processor->_muts);

  ok($@ =~ m//,"No error from eval.");
	is_deeply($processor->_muts->{'alnp'}, [66,102,51,56,110,63,61,44,82,80,123,88,123,76,87,60,118,93,139], "primary alignment scores");
	is_deeply($processor->_muts->{'trl'}, [151,151,151,151,151,151,151,151,151,151,151,151,151,151,151,151,151,151,151], "tumor read lengths");

	#getAlnScoreMedianReadAdjusted
	$exp_res = sprintf('%.2f',0.543046358);
	is($processor->alignmentScoreReadLengthAdjustedFlag, $exp_res,"alignmentScoreReadLengthAdjustedFlag");
	#getAlignmentScoreMedian
	$exp_res = sprintf('%.2f',82);
	is($processor->alnScoreMedianFlag, $exp_res,"alnScoreMedianFlag");

  done_testing();
};

subtest 'Clipped Read tests' => sub {
	my $processor = new_ok('Sanger::CGP::CavemanPostProcessing::Flagger');
  my $cfg = Sanger::CGP::CavemanPostProcessing::Config->new();
  my $opts;
  $opts->{'c'} = $CONFIG;
  $opts->{'v'} = $FLAG_CFG;
  $opts->{'s'} = $SPECIES;
  $opts->{'t'} = $TYPE;
  $opts->{'g'} = $GERMLINE_INDEL;
  $opts->{'b'} = $BED_LOC;
  $cfg->init_config($opts);
  $processor->init($CLIP_M_BAM,$CLIP_N_BAM,$cfg);
  my $chr = 1;
	my $pos = 10437;
	my $ref = "T";
	my $mut = "C";
	eval{$processor->set_position($chr,$pos,$ref,$mut);};
  ok($@ =~ m//,"No error from eval.");

	#Manually set counts
	my $exp_sclp = [1,2,3,4,5,6,7,8,9,10];
	$processor->_muts->{'sclp'} = $exp_sclp;
	is_deeply($processor->_muts->{'sclp'}, [1,2,3,4,5,6,7,8,9,10], "softclipcounts");

	#Check manually set count results
	#getClipMedianResult
	my $exp_res = sprintf('%.2f',5.5);
	is($processor->clippingMedianFlag, $exp_res,"clippingMedianFlag - Fake");

	#Reset and use real data
	eval{$processor->set_position($chr,$pos,$ref,$mut);};
  ok($@ =~ m//,"No error from eval.");

	#Check counts have been filled correctly.
	$exp_sclp = [0,9,78,78,0,83,49,104,38,35,0,20,0,0,0,66,0,0,0];
	is_deeply($processor->_muts->{'sclp'}, [0,9,78,78,0,83,49,104,38,35,0,20,0,0,0,66,0,0,0], "softclipcounts");

	#Check real data count results
	$exp_res = sprintf('%.2f',9);
	is($processor->clippingMedianFlag, $exp_res,"clippingMedianFlag - Real");
  done_testing();
};


subtest 'tumIndelDepthFlag' => sub {
  my $processor = new_ok('Sanger::CGP::CavemanPostProcessing::Flagger');
  my $cfg = Sanger::CGP::CavemanPostProcessing::Config->new();
  my $opts;
  $opts->{'c'} = $CONFIG;
  $opts->{'v'} = $FLAG_CFG;
  $opts->{'s'} = $SPECIES;
  $opts->{'t'} = $TYPE;
  $opts->{'g'} = $GERMLINE_INDEL;
  $opts->{'b'} = $BED_LOC;
  $cfg->init_config($opts);
  $processor->init($T_BAM,$T_BAM,$cfg);
  my $chr = 1;
	my $pos = 10011533;
	my $ref = "G";
	my $mut = "T";
	eval{$processor->set_position($chr,$pos,$ref,$mut);};
  ok($@ =~ m//,"No error from eval.");

  $processor->_muts->{'indelTCount'} = 9;
	$processor->_muts->{'totalTCoverage'} = 100;
	#Will pass
	ok($processor->tumIndelDepthFlag == 0,"Pass tum indel filter 9%");

	#borderline pass
	$processor->_muts->{'indelTCount'} = 10;
	$processor->_muts->{'totalTCoverage'} = 100;
	ok($processor->tumIndelDepthFlag == 0,"Pass (borderline) tum indel filter 10%");
	#Will fail
	$processor->_muts->{'indelTCount'} = 11;
	$processor->_muts->{'totalTCoverage'} = 100;
	ok($processor->tumIndelDepthFlag == 1,"Fail tum indel filter 11%");

	#Check clear results and change of cutoff works.
	$processor->_prms->{'maxTumIndelProportion'} = 11;
	ok($processor->tumIndelDepthFlag == 0,"Pass tum indel filter 11% changed proportion");
  done_testing();
};

subtest 'sameReadPosFlag' => sub {
  my $processor = new_ok('Sanger::CGP::CavemanPostProcessing::Flagger');
  my $cfg = Sanger::CGP::CavemanPostProcessing::Config->new();
  my $opts;
  $opts->{'c'} = $CONFIG;
  $opts->{'v'} = $FLAG_CFG;
  $opts->{'s'} = $SPECIES;
  $opts->{'t'} = $TYPE;
  $opts->{'g'} = $GERMLINE_INDEL;
  $opts->{'b'} = $BED_LOC;
  $cfg->init_config($opts);
  $processor->init($T_BAM,$T_BAM,$cfg);
  my $chr = 1;
	my $pos = 10011533;
	my $ref = "G";
	my $mut = "T";
	eval{$processor->set_position($chr,$pos,$ref,$mut);};
  ok($@ =~ m//,"No error from eval.");

  $processor->_muts->{'trp'} = [51,75,75,75,75,75,75,75,75,75];
	ok($processor->sameReadPosFlag == 1,"Fail same read pos");
	#Pass same read pos
	$processor->_muts->{'trp'} = [51,45,55,22,75,75,75,75,75,75];
	ok($processor->sameReadPosFlag == 0,"Pass same read pos");

	#Pass same read pos change % same pos
	$processor->_prms->{'samePosMaxPercent'} = 50;
	ok($processor->sameReadPosFlag == 1,"Fail same read pos, changed %age");
  done_testing();
};

subtest 'cavemanMatchNormalProportionFlag' => sub {
  my $processor = new_ok('Sanger::CGP::CavemanPostProcessing::Flagger');
  my $cfg = Sanger::CGP::CavemanPostProcessing::Config->new();
  my $opts;
  $opts->{'c'} = $CONFIG;
  $opts->{'v'} = $FLAG_CFG;
  $opts->{'s'} = $SPECIES;
  $opts->{'t'} = $TYPE;
  $opts->{'g'} = $GERMLINE_INDEL;
  $opts->{'b'} = $BED_LOC;
  $cfg->init_config($opts);
  $processor->init($T_BAM,$T_BAM,$cfg);
  my $chr = 1;
  my $pos = 10011533;
  my $ref = "A";
  my $mut = "T";
  eval{$processor->set_position($chr,$pos,$ref,$mut);};
  ok($@ =~ m//,"No error from eval.");
  my $normal_col = '0/0:90:0:0:10:0.1'; #0.1
  my $normal_col_fail = '0/0:60:0:0:40:0.4'; #0.4
  my $tumcol = '1/0:50:0:0:50:0.5'; #0.5
  my $oldformat = 'GT:AA:CA:GA:TA:PM';
  ok($processor->_getCavemanMatchNormalProportionFlag($normal_col,$tumcol,$oldformat)==0,"Pass caveman matched normal check old format");
  ok($processor->_getCavemanMatchNormalProportionFlag($normal_col_fail,$tumcol,$oldformat)==1,"Fail caveman matched normal check old format");
  my $newformat = 'GT:FAZ:FCZ:FGZ:FTZ:RAZ:RCZ:RGZ:RTZ:PM';
  $normal_col = '0/0:45:0:0:5:45:0:0:5:0.1'; #0.1
  $normal_col_fail = '0/0:30:0:0:20:30:0:0:20:0.4'; #0.4
  $tumcol = '1/0:25:0:0:25:25:0:0:25:0.5'; #0.5
  ok($processor->_getCavemanMatchNormalProportionFlag($normal_col,$tumcol,$newformat)==0,"Pass caveman matched normal check new format");
  ok($processor->_getCavemanMatchNormalProportionFlag($normal_col_fail,$tumcol,$newformat)==1,"Fail caveman matched normal check new format");
};

subtest 'Read Gap Tests' => sub {
  #Fail as more than proportion
  my $processor = new_ok('Sanger::CGP::CavemanPostProcessing::Flagger');
  my $cfg = Sanger::CGP::CavemanPostProcessing::Config->new();
  my $opts;
  $opts->{'c'} = $CONFIG;
  $opts->{'v'} = $FLAG_CFG;
  $opts->{'s'} = $SPECIES;
  $opts->{'t'} = $TYPE;
  $opts->{'g'} = $GERMLINE_INDEL;
  $opts->{'b'} = $BED_LOC;
  $cfg->init_config($opts);
  $processor->init($T_BAM,$T_BAM,$cfg);
  my $chr = 1;
  my $pos = 10011533;
  my $ref = "A";
  my $mut = "T";
  eval{$processor->set_position($chr,$pos,$ref,$mut);};
  ok($@ =~ m//,"No error from eval.");

  ok($processor-> _prms->{'minGapPresentInReads'}==30,"Correct minGapPresentInPercentReads");
  ok($processor-> _prms->{'minMeanMapQualGapFlag'}==10,"Correct minMeanMapQualGapFlag");
  ok($processor-> _prms->{'withinXBpOfDeletion'}==10,"Correct withinXBpOfDeletion");

  ok($processor->withinGapRangeFlag==0,"Initially passes read gap flag");

  #Change to fail 
  my $allTumMapQuals = [60,29,60,60];
  my $allTumBases = ['T','G','G','T'];
  my $allMinGapDistances = [-1,2,2,-1];
  $processor->_muts->{'allTumMapQuals'} = $allTumMapQuals;
  $processor->_muts->{'allTumBases'} = $allTumBases;
  $processor->_muts->{'allMinGapDistances'} = $allMinGapDistances;
  ok($processor->withinGapRangeFlag==1,"Fail read gap.");

  #Change to pass on map qualities
  eval{$processor->set_position($chr,$pos,$ref,$mut);};
  ok($@ =~ m//,"No error from eval.");
  $allTumMapQuals = [5,5,5,6];
  $processor->_muts->{'allTumMapQuals'} = $allTumMapQuals;
  ok($processor->withinGapRangeFlag==0,"Pass read gap on map qualities.");

  #Change to fail
  $allTumMapQuals = [60,29,60,60];
  $allTumBases = ['T','G','G','T'];
  $allMinGapDistances = [-1,2,2,-1];
  $processor->_muts->{'allTumMapQuals'} = $allTumMapQuals;
  $processor->_muts->{'allTumBases'} = $allTumBases;
  $processor->_muts->{'allMinGapDistances'} = $allMinGapDistances;
  ok($processor->withinGapRangeFlag==1,"Fail read gap 2.");

  #Change to pass on distance from deletion
  $allTumMapQuals = [60,29,60,60];
  $allTumBases = ['T','G','G','T'];
  $allMinGapDistances = [-1,11,11,-1];
  $processor->_muts->{'allTumMapQuals'} = $allTumMapQuals;
  $processor->_muts->{'allTumBases'} = $allTumBases;
  $processor->_muts->{'allMinGapDistances'} = $allMinGapDistances;
  ok($processor->withinGapRangeFlag==0,"Pass on deletion distance.");

  #Change to fail
  $allTumMapQuals = [60,29,60,60];
  $allTumBases = ['T','G','G','T'];
  $allMinGapDistances = [-1,2,2,-1];
  $processor->_muts->{'allTumMapQuals'} = $allTumMapQuals;
  $processor->_muts->{'allTumBases'} = $allTumBases;
  $processor->_muts->{'allMinGapDistances'} = $allMinGapDistances;
  ok($processor->withinGapRangeFlag==1,"Fail read gap 3.");

  #Change to pass on percentage reads
  $allTumMapQuals = [60,29,60,60,60];
  $allTumBases = ['T','G','G','G','T'];
  $allMinGapDistances = [-1,11,11,2,-1];
  $processor->_muts->{'allTumMapQuals'} = $allTumMapQuals;
  $processor->_muts->{'allTumBases'} = $allTumBases;
  $processor->_muts->{'allMinGapDistances'} = $allMinGapDistances;
  ok($processor->withinGapRangeFlag==0,"Pass on percentage reads.");
  
  done_testing();
};
