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
use warnings FATAL => 'all';
use autodie;
use Data::Dumper;
use Const::Fast qw(const);
use Capture::Tiny qw(capture);
use Carp;

use Test::More tests => 9;

use FindBin qw($Bin);
const my $lib_path => "$Bin/../lib";
const my $script_path => "$Bin/../bin/";
const my $test_data_path => "$Bin/../testData/";

const my $test_out_file => $test_data_path.'actual_output_flagVCF.vcf';
const my $test_input_file => $test_data_path.'input_FlagCaVEManVCF.vcf.1';
const my $expect_output => $test_data_path.'expected_output_FlagCaVEManVCF_new.vcf';
const my $test_mut_bam => $test_data_path.'testflagmt.bam';
const my $test_norm_bam => $test_data_path.'testflagwt.bam';
const my $test_germ_indel_bed => $test_data_path.'germline_indel.bed.gz';
const my $test_snp_bed => $test_data_path.'';
const my $test_simple_repeat_bed => $test_data_path.'simple_repeats.bed.gz';
const my $test_centro_repeat_bed => $test_data_path.'centromeric_repeats.bed.gz';
const my $test_exp_new_config => $test_data_path.'flag.vcf.config.new.ini';
const my $test_config => $test_data_path.'flag.vcf.config.tmp.ini';
const my $test_vcf_convert_config => $test_data_path.'flag.to.vcf.convert.ini';
const my $script => $script_path.'flagCaVEMan.pl';
const my $oldVersionNormBam => $test_data_path.'oldversion.wt.bam';
const my $oldVersionMutBam => $test_data_path.'oldversion.wt.bam';
const my $oldVersionVcf => $test_data_path.'test.oldversion.vcf';
const my $test_out_oldVersionVCF => $test_data_path.'test.oldversion.flagged.vcf';
const my $expected_out_oldVersionVCF => $test_data_path.'expected_output_oldVersion_flagged.vcf';
const my $unmatched_vcf => $test_data_path.'unmatchedVCF/';
const my $ref => $test_data_path.'genome.fa.fai';

my $perl_exec = "$^X -I ".(join ':', @INC);

const my $id_analysis_proc => 123456;

main(@ARGV);

sub main{
  my (@ARGV) = @_;
  setup_config();
  run_flag();
  checkSupportedButNoFlags();
  checkUnsupportedOption();
  run_old_flag();
  remove_created_files();
  done_testing();
}

sub run_flag{
	my $cmd = "$perl_exec $script -i $test_input_file -o $test_out_file -c $test_config".
			" -s human -t genome -m $test_mut_bam -n $test_norm_bam".
			" -g $test_germ_indel_bed -v $test_vcf_convert_config".
			" -b $test_snp_bed -ab $test_snp_bed -p $id_analysis_proc".
			" -u $unmatched_vcf -ref $ref";
	my ($out, $err, $exit) = capture{ system($cmd) };
	is($exit, 0,'Flagging ran correctly');
  #Check headers
  my $exp_header = read_header($expect_output);
  my $got_header = read_header($test_out_file);
  is_deeply(_tab_data_to_object($got_header), _tab_data_to_object($exp_header),'compare headers');
  #Check non headers
  my $exp_body = read_non_header($expect_output);
  my $got_body = read_non_header($test_out_file);
  is_deeply(_tab_data_to_object($got_body), _tab_data_to_object($exp_body),'compare non-header');
}

sub run_old_flag{
	my $cmd = "$perl_exec $script -i $oldVersionVcf -o $test_out_oldVersionVCF -c $test_config".
			" -s human -t genome -m $oldVersionMutBam -n $oldVersionNormBam".
			" -g $test_germ_indel_bed -v $test_vcf_convert_config ".
			" -b $test_snp_bed -ab $test_snp_bed -p $id_analysis_proc ".
			" -u $unmatched_vcf -ref $ref";
	my ($out, $err, $exit) = capture{ system($cmd) };
	isnt($exit, 0, 'Old version flagging ran correctly');
}

sub checkSupportedButNoFlags{
	my $cmd = "$perl_exec $script -i $test_input_file -o $test_out_file -c $test_config".
			" -s human -t pulldown -m $test_mut_bam -n $test_norm_bam".
			" -v $test_vcf_convert_config".
			" -b $test_snp_bed -ab $test_snp_bed -p $id_analysis_proc ".
			"-u $unmatched_vcf -ref $ref";
	my ($out, $err, $exit) = capture{ system($cmd) };
	like($err, qr/No flag list found in .+ for section. No flagging will be done./, 'Correctly runs with warnings as no flags available (check message).');
	is($exit, 0, 'Correctly runs with warnings as no flags available (check exitcode).');
}

sub checkUnsupportedOption{
	my $cmd = "$perl_exec $script -i $test_input_file -o $test_out_file -c $test_config".
			" -s human -t pulldown -m $test_mut_bam -n $test_norm_bam".
			" -v $test_vcf_convert_config".
			" -b $test_snp_bed -ab $test_snp_bed -p $id_analysis_proc -x".
			" -u $unmatched_vcf -ref $ref";
	my ($out, $err, $exit) = capture{ system($cmd) };
	isnt($exit, 0, 'Correctly fails for not supported CL option (check exitcode).');
}


sub setup_config{
	my $FH;
	open($FH, '<', $test_exp_new_config) || croak("Error trying to open config file $test_exp_new_config: $!");
		my @data = <$FH>;
		my $dat = join("",@data);
	close($FH);

	my $OUT;
	open($OUT, '>', $test_config) || croak("Error trying to open config file for write $test_config: $!");
		print $OUT $dat;
	close($OUT);
}

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

sub remove_created_files{
	#vcf
	ok(unlink($test_out_file));
	#config created
  ok(unlink($test_config));
}

