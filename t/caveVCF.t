##########LICENCE##########
# Copyright (c) 2014-2016 Genome Research Ltd.
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
use Sanger::CGP::CavemanPostProcessing::CaveVCF;
use Sanger::CGP::CavemanPostProcessing::Config;

use Data::Dumper;
use Const::Fast qw(const);
use Carp;
use Capture::Tiny ':all';

use Test::More tests => 6;

use FindBin qw($Bin);
const my $lib_path => "$Bin/../lib";
const my $test_data_path => "$Bin/../testData/";
const my $GOOD_INPUT_VCF => $test_data_path."input_FlagCaVEManVCF.vcf";
const my $OLD_VER_INPUT => $test_data_path."test.oldversion.vcf";
const my $EXP_NEW_OUTPUT => $test_data_path."exp_new_output.vcf";
const my $OUT_VCF => $test_data_path."test_tmp_out.vcf";
const my $CONFIG => $test_data_path."flag.vcf.config.new.ini";
const my $FLAG_CFG => $test_data_path."flag.to.vcf.convert.ini";
const my $SPECIES => 'HUMAN';
const my $TYPE => 'WGS';
const my $GERMLINE_INDEL => $test_data_path.'germline_indel.bed';
const my $BED_LOC => $test_data_path;

subtest 'CaveVCF new' => sub {
  my $cfg = new_ok('Sanger::CGP::CavemanPostProcessing::CaveVCF');
  done_testing();
};

subtest 'CaveVCF init_vcfs' => sub {
  my $cfg = new_ok('Sanger::CGP::CavemanPostProcessing::CaveVCF');
  eval{$cfg->init_vcfs($OLD_VER_INPUT,$OUT_VCF);};
  ok($@ =~ m/Attempting to use a new version of CaVEManPostProcessing with an incompatible CaVEMan version/,"Correct error for old version.");
  $cfg->init_vcfs($GOOD_INPUT_VCF,$OUT_VCF);
  ok(-e $OUT_VCF,"Output file created.");
  ok(unlink($OUT_VCF),"Remove output VCF");
  done_testing();
};

subtest 'CaveVCF output' => sub {
  my $caveVCF = new_ok('Sanger::CGP::CavemanPostProcessing::CaveVCF');
  eval{$caveVCF->init_vcfs($OLD_VER_INPUT,$OUT_VCF);};
  ok($@ =~ m/Attempting to use a new version of CaVEManPostProcessing with an incompatible CaVEMan version/,"Correct error for old version.");
  $caveVCF->init_vcfs($GOOD_INPUT_VCF,$OUT_VCF);
  ok(-e $OUT_VCF,"Output file created.");
  my $cfg = Sanger::CGP::CavemanPostProcessing::Config->new();
  my $opts;
  $opts->{'c'} = $CONFIG;
  $opts->{'v'} = $FLAG_CFG;
  $opts->{'s'} = $SPECIES;
  $opts->{'t'} = $TYPE;
  $opts->{'g'} = $GERMLINE_INDEL;
  $opts->{'b'} = $BED_LOC;
  $opts->{'p'} = 12345;
  $opts->{'f'} = $GOOD_INPUT_VCF;
  $cfg->init_config($opts);
  $caveVCF->output_header($opts,$cfg);

  #Check header is as expected.
  my $exp_header = read_header($EXP_NEW_OUTPUT);
  my $got_header = read_header($OUT_VCF);
  is_deeply(_tab_data_to_object($got_header), _tab_data_to_object($exp_header),'compare');
  $caveVCF->close_vcfs;
  ok(unlink($OUT_VCF),"Remove output VCF");
  done_testing();
};

subtest 'get/set methods' => sub {
  my $caveVCF = new_ok('Sanger::CGP::CavemanPostProcessing::CaveVCF');
  my $test = "TEST";
  $caveVCF->caveman_version($test);
  ok($caveVCF->caveman_version eq $test,"Check caveman_version");
  $test = "TEST_out";
  $caveVCF->output_vcf($test);
  ok($caveVCF->output_vcf eq $test,"Check output_vcf");
  $test = "TEST_in";
  $caveVCF->input_vcf($test);
  ok($caveVCF->input_vcf eq $test,"Check input_vcf");
  done_testing();
};

subtest 'CaVEMan version testing' => sub {
  my $caveVCF = new_ok('Sanger::CGP::CavemanPostProcessing::CaveVCF');
  my $test_bad_ver = "0_3_0";
  ok($caveVCF->isOldVersionOfCaVEMan($test_bad_ver)==1,"Check old version of CaVEMan");
  my $test_good_ver = "1.8.0";
  ok($caveVCF->isOldVersionOfCaVEMan($test_good_ver)==0,"Check good version of CaVEMan");
  done_testing();
};

subtest 'output_vcf_lines' => sub {
  my $caveVCF = new_ok('Sanger::CGP::CavemanPostProcessing::CaveVCF');
  $caveVCF->init_vcfs($GOOD_INPUT_VCF,$OUT_VCF);
  ok(-e $OUT_VCF,"Output file created.");
  my $cfg = Sanger::CGP::CavemanPostProcessing::Config->new();
  my $opts;
  $opts->{'c'} = $CONFIG;
  $opts->{'v'} = $FLAG_CFG;
  $opts->{'s'} = $SPECIES;
  $opts->{'t'} = $TYPE;
  $opts->{'g'} = $GERMLINE_INDEL;
  $opts->{'b'} = $BED_LOC;
  $opts->{'p'} = 12345;
  $opts->{'f'} = $GOOD_INPUT_VCF;
  $cfg->init_config($opts);
  $caveVCF->output_header($opts,$cfg);
  my $test_line = "15\t22095504\tID_36061\tG\tA\t.\tDP\tDP=22;MP=8.4e-01;GP=1.5e-01;TG=GG/AG;TP=8.1e-01;SG=AG/GG;SP=7.8e-02\tGT:AA:CA:GA:TA:PM\t0|0:0:0:9:0:0.0e+00\t0|1:4:0:9:0:3.1e-01\n";
  my @lines = ($test_line);
  $caveVCF->output_vcf_lines(\@lines);
  $caveVCF->close_vcfs;
  my $exp_header = read_header($EXP_NEW_OUTPUT);
  my $got_header = read_header($OUT_VCF);
  is_deeply(_tab_data_to_object($got_header), _tab_data_to_object($exp_header),'compare headers');
  my $exp_body = read_non_header($EXP_NEW_OUTPUT);
  my $got_body = read_non_header($OUT_VCF);
  is_deeply(_tab_data_to_object($got_body), _tab_data_to_object($exp_body),'compare non-header');

  ok(unlink($OUT_VCF),"Remove output VCF");
  done_testing();
};

sub read_non_header{
  my ($file) = @_;
  my $cmd = "grep -ve '^#' $file";
  my ($stdout, $stderr, $exit) = capture{
    system($cmd);
  };
  croak("Error reading non-header from $file: $exit. STDERR: $stderr") unless($exit==0);
  my $body = $stdout;
  return $body;

}

sub read_header{
  my ($file) = @_;
  my $cmd = "grep -e '^#' $file";
  my ($stdout, $stderr, $exit) = capture{
    system($cmd);
  };
  croak("Error reading header from $file: $exit. STDERR: $stderr") unless($exit==0);
  my $head = $stdout;
  $head =~ s/InputVCFParam=<(.+)>>/InputVCFParam=<>>/g;
	$head =~ s/InputVCFVer=<(.+)>>/InputVCFVer=<>>/g;
	$head =~ s/cgpAnalysisProc_\d+\./cgpAnalysisProc_./g;
	$head =~ s/vcfProcessLog_\d+\./vcfProcessLog_./g;
  return $head;
}

sub _tab_data_to_object {
  my @lines = split /\n/, shift;
  my @final;
  for(@lines) {
    if($_ =~ m/^#/) {
      push @final, $_;
      next;
    }
    my @tmp = split /\t/, $_;
    $tmp[6] = [sort (split /;/,$tmp[6])];
    $tmp[7] = [sort (split /;/,$tmp[7])];
    push @final, \@tmp;
  }
  return \@final;
}