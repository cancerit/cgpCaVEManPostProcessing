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
use Test::More tests => 13;
use strict;
use warnings FATAL => 'all';
use Carp;
use Data::Dumper;
use Capture::Tiny qw(capture);

use FindBin qw($Bin);
my $script_path = "$Bin/../bin/";
my $test_data_path = "$Bin/../testData/";

my $index = 1;
my $test_out_file = $test_data_path.'actual_output_flagVCF.vcf';
my $test_input_file = $test_data_path.'input_FlagCaVEManVCF.vcf';
my $expect_output = $test_data_path.'expected_output_FlagCaVEManVCF.vcf';
my $test_mut_bam = $test_data_path.'testflagmt.bam';
my $test_norm_bam = $test_data_path.'testflagwt.bam';
my $test_germ_indel_bed = $test_data_path.'germline_indel.bed';
my $test_snp_bed = $test_data_path.'';
my $test_simple_repeat_bed = $test_data_path.'simple_repeats.bed';
my $test_centro_repeat_bed = $test_data_path.'centromeric_repeats.bed';
my $test_exp_config = $test_data_path.'flag.vcf.config.ini';
my $test_config = $test_data_path.'flag.vcf.config.tmp.ini';
my $test_vcf_convert_config = $test_data_path.'flag.to.vcf.convert.ini';
my $script = $script_path.'cgpFlagCaVEMan.pl';
my $oldVersionNormBam = $test_data_path.'oldversion.wt.bam';
my $oldVersionMutBam = $test_data_path.'oldversion.wt.bam';
my $oldVersionVcf = $test_data_path.'test.oldversion.vcf';
my $test_out_oldVersionVCF = $test_data_path.'test.oldversion.flagged.vcf';
my $expected_out_oldVersionVCF = $test_data_path.'expected_output_oldVersion_flagged.vcf';
my $unmatched_vcf = $test_data_path.'unmatchedVCF/';
my $ref = $test_data_path.'genome.fa.fai';

my $id_analysis_proc = 123456;

main(@ARGV);

sub main{
	setup_config();
	run_flag();
	compare($test_out_file, $expect_output,$index);
	checkSupportedButNoFlags();
	checkUnsupportedOption();
	run_old_flag();
	#check_old_version_correct
	compare($test_out_oldVersionVCF, $expected_out_oldVersionVCF);
	remove_created_files();
}

sub setup_config{
	my $FH;
	open($FH, '<', $test_exp_config) || croak("Error trying to open config file $test_exp_config: $!");
		my @data = <$FH>;
		my $dat = join("",@data);
	close($FH);

	my $OUT;
	open($OUT, '>', $test_config) || croak("Error trying to open config file for write $test_config: $!");
		print $OUT $dat;
	close($OUT);
}

sub run_flag{
	my $cmd = "perl $script -i $test_input_file -o $test_out_file -c $test_config".
			" -s human -t genome -m $test_mut_bam -n $test_norm_bam".
			" -g $test_germ_indel_bed -v $test_vcf_convert_config ".
			"--index $index -b $test_snp_bed -ab $test_snp_bed -p $id_analysis_proc ".
			"-u $unmatched_vcf -ref $ref";
	my ($out, $err, $exit) = capture{ system($cmd) };
	is($exit, 0,'Flagging ran correctly');
}

sub run_old_flag{
	my $cmd = "perl $script -i $oldVersionVcf -o $test_out_oldVersionVCF -c $test_config".
			" -s human -t genome -m $oldVersionMutBam -n $oldVersionNormBam".
			" -g $test_germ_indel_bed -v $test_vcf_convert_config ".
			" -b $test_snp_bed -ab $test_snp_bed -p $id_analysis_proc ".
			" -u $unmatched_vcf -ref $ref";
	my ($out, $err, $exit) = capture{ system($cmd) };
	is($exit, 0, 'Old version flagging ran correctly');
}

sub checkSupportedButNoFlags{
	my $cmd = "perl $script -i $test_input_file -o $test_out_file -c $test_config".
			" -s human -t pulldown -m $test_mut_bam -n $test_norm_bam".
			" -v $test_vcf_convert_config".
			" --index $index -b $test_snp_bed -ab $test_snp_bed -p $id_analysis_proc ".
			"-u $unmatched_vcf -ref $ref";
	my ($out, $err, $exit) = capture{ system($cmd) };
	like($err, qr/No flagList found in .+ for section HUMAN_PULLDOWN FLAGLIST. No flagging will be done./, 'Correctly runs with warnings as no flags available (check message).');
	is($exit, 0, 'Correctly runs with warnings as no flags available (check exitcode).');
}

sub checkUnsupportedOption{
	my $cmd = "perl $script -i $test_input_file -o $test_out_file -c $test_config".
			" -s human -t pulldown -m $test_mut_bam -n $test_norm_bam".
			" -v $test_vcf_convert_config".
			" --index $index -b $test_snp_bed -ab $test_snp_bed -p $id_analysis_proc -x".
			" -u $unmatched_vcf -ref $ref";
	my ($out, $err, $exit) = capture{ system($cmd) };
	like($err, qr/Unknown parameter: -x at .+ line \d+./, 'Correctly fails for not supported CL option (check message).');
	isnt($exit, 0, 'Correctly fails for not supported CL option (check exitcode).');
}

sub compare{
	my ($test_out_file, $expect_output,$idx) = @_;
	#Read new file
	my $NW;
	my $new;
	if(defined($idx)){
		$test_out_file .= ".$index"
	}
	open($NW,'<',$test_out_file) || croak("Error reading new output $test_out_file : $!");
		my @dat = <$NW>;
		$new = join("",@dat);
		close($NW);
	#Read expected file
	my $XP;
	my $exp;
	open($XP,'<',$expect_output) || croak("Error reading expected output $expect_output:$!");
		@dat = <$XP>;
		$exp = join("",@dat);
	close($XP);
	#Make modifications for vcfProcessLog entries etc

	my ($current_processLog_date) = $new =~ /##vcfProcessLog_([0-9]+)/;# this is a hack and will replace all occurences with the first date it comes across...
	$exp =~ s/##vcfProcessLog_[0-9]+/##vcfProcessLog_$current_processLog_date/g;
	$new =~ s/##vcfProcessLog_[0-9]+/##vcfProcessLog_$current_processLog_date/g;

	my ($current_analysis_proc_date) = $new =~ /##cgpAnalysisProc_([0-9]+)/;# this is a hack and will replace all occurences with the first date it comes across...
	$exp =~ s/##cgpAnalysisProc_[0-9]+/##cgpAnalysisProc_$current_processLog_date/g;
	$new =~ s/##cgpAnalysisProc_[0-9]+/##cgpAnalysisProc_$current_processLog_date/g;
	######################



	######################
	## Becuase the params go in in random order we have to take them out and test them seperately.
	#my ($vcf_test_log_loc,$output_test_string_proces_log_params) = $exp =~ /##vcfProcessLog_([0-9]+\.[0-9]+).+InputVCFParam=<(.+)>>/g;
	while($exp =~ /##vcfProcessLog_([0-9]+\.[0-9]+)=<InputVCF=<([^>]+)>.+InputVCFParam=<(.+)>>/g){
		my $vcf_test_log_loc = $1;
		my $input = $2;
		my $output_test_string_proces_log_params = $3;
		my ($vcf_log_loc, $output_proces_log_params) = $new =~ /##vcfProcessLog_($vcf_test_log_loc)=<InputVCF=<$input>.+InputVCFParam=<(.+)>>/g;
		my $output_test_string_proces_log_params_hash = {};
		my $output_proces_log_params_hash = {};

		foreach my $pair(split(",",$output_test_string_proces_log_params)){
			my ($key,$value) = split("=",$pair);
			$output_test_string_proces_log_params_hash->{$key} = $value;
		}

		foreach my $pair(split(",",$output_proces_log_params)){
			my ($key,$value) = split("=",$pair);
			$output_proces_log_params_hash->{$key} = $value;
		}
		is_deeply($output_proces_log_params_hash,$output_test_string_proces_log_params_hash,"convert:: process log params");
	}

	$exp =~ s/InputVCFParam=<(.+)>>/InputVCFParam=<>>/g;
	$new =~ s/InputVCFParam=<(.+)>>/InputVCFParam=<>>/g;
	$exp =~ s/InputVCFVer=<(.+)>>/InputVCFVer=<>>/g;
	$new =~ s/InputVCFVer=<(.+)>>/InputVCFVer=<>>/g;
	#check they equate.
	is_deeply(_tab_data_to_object($new), _tab_data_to_object($exp),'compare');
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
	unlink($test_out_file.".$index");
	#config created
	unlink($test_config);
	unlink($test_out_oldVersionVCF);
}
