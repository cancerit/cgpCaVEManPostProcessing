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
package Sanger::CGP::CavemanPostProcessing::CaveVCF;

use strict;
use warnings FATAL=>'all';
use autodie;
use Const::Fast;
use Carp qw(croak cluck);
use Vcf;

use Sanger::CGP::CavemanPostProcessing;

const my $GOOD_CAVE_VER => '0_3_1';

sub new {
	my ($proto) = @_;
	my $class = ref($proto) || $proto;
 	my $self = {};
  bless($self, $class);
  return $self;
}

sub init_vcfs{
  my ($self,$inputVCF,$outputVCF) = @_;
  return $self->_setup($inputVCF,$outputVCF);
}

sub _setup{
  my ($self,$inputVCF,$outputVCF) = @_;
  my $vcf = Vcf->new(file=>$inputVCF, version=>'4.1');
  $self->input_vcf($vcf);
  $self->input_vcf->parse_header();
  $self->caveman_version($self->getCavemanVersionFromVCF($self->input_vcf));
  croak("Attempting to use a new version of CaVEManPostProcessing with an incompatible CaVEMan version '".$self->caveman_version."'") if($self->isOldVersionOfCaVEMan($self->caveman_version));
  $self->_open_output($outputVCF);
  return;
}

sub _open_output{
  my ($self,$file) = @_;
  my $FH;
  open($FH, '>', $file);
  $self->output_vcf($FH);
  return;
}

sub output_vcf_lines{
  my ($self,$lines) = @_;
  foreach my $line(@$lines){
		my @all_fields = split /\t/, $line;
		my @flags = split /;/,$all_fields[6];
		$all_fields[6] = join ";" , sort {$a cmp $b} @flags;
		$line = join "\t", @all_fields;
    print {$self->output_vcf} $line or croak("Error trying to write VCF line '".$line."' to output file: $!");
  }
  return;
}

sub output_header{
  my ($self,$opts,$cfg) = @_;
  #If we have analysis_process
  if(exists($opts->{'p'}) && defined($opts->{'p'})){
		$self->input_vcf->add_header_line({key=>'cgpAnalysisProc',value=>$opts->{'p'}}, append=>1);
	}
  foreach my $flag(@{$cfg->flags}){
    if($flag->is_info){
      $self->input_vcf->add_header_line({key=>'INFO', ID=>$flag->id,Type=>$flag->type,Number=>$flag->value,Description=>$flag->description});
    }else{
      $self->input_vcf->add_header_line({key=>'FILTER', ID=>$flag->id,Description=>$flag->description});
    }
  }
  $self->add_trimmed_opts($opts);
  print {$self->output_vcf} $self->input_vcf->format_header() or croak ("Error printing header to output VCF file: $!");
}

sub add_trimmed_opts{
  my ($self,$opts) = @_;
  my $outputOpts;
  foreach my $key(keys %$opts){
    if(defined($opts->{$key})){
      if(-f $opts->{$key}){
        $outputOpts->{$key} = _trim_file_path($opts->{$key});
      }else{
        $outputOpts->{$key} = $opts->{$key};
      }
    }else{
			$outputOpts->{$key} = '.';
		}
  }
  $self->input_vcf->add_header_line({
                    key=>'vcfProcessLog',
                    'InputVCF'=>''.$outputOpts->{'f'}.'',
                    'InputVCFSource'=>'flagCaVEMan.pl',
                    'InputVCFVer'=>''.$VERSION.'',
                    'InputVCFParam'=>$outputOpts
                    } ,append => 1);
	return;
}

sub _trim_file_path{
  my ($string) = @_;
	my @bits = (split("/", $string));
	return pop @bits;
}

sub caveman_version{
  my ($self,$ver) = @_;
  if(defined($ver)){
    $self->{'cave_v'} = $ver;
  }
  return $self->{'cave_v'};
}

sub output_vcf{
  my ($self,$fh) = @_;
  if(defined($fh)){
    $self->{'out_vcf'} = $fh;
  }
  return $self->{'out_vcf'};
}

sub input_vcf{
  my ($self,$vcf) = @_;
  if(defined($vcf)){
    $self->{'in_vcf'} = $vcf;
  }
  return $self->{'in_vcf'};
}

sub close_vcfs{
  my ($self,) = @_;
  $self->input_vcf->close();
  close($self->output_vcf);
  return;
}

sub isOldVersionOfCaVEMan{
	my ($self,$ver) = @_;
	return 1 if(!defined($ver) || $ver eq "");
	my ($good_main,$good_feat,$good_bug) = split(/[_\.]/,$GOOD_CAVE_VER);
	$ver =~ s/^[^\|]+\|//g;
	$ver =~ s/\|[^\|]+$//g;
	my ($main,$feat,$bug) = split(/[\._]/,$ver);
	return 1 if( $main < $good_main || ($main == $good_main && $feat < $good_feat) || ($main == $good_main && $feat == $good_feat && $bug < $good_bug));
	return 0;
}

sub getCavemanVersionFromVCF{
	my ($self,$vcf) = @_;
	my $caveVer = undef;
	my $verLine = $vcf->get_header_line(key=>'cavemanVersion');
	$caveVer = $verLine->[0]->[0]->{'value'};
	#This might be an original CaVEMan VCF or an old run before the headers were fixed
	if(!defined($caveVer) || $caveVer eq ""){
		my $plLines = $vcf->get_header_line(key=>'vcfProcessLog');
		foreach my $plln(@$plLines){
			if($plln->[0]->{'InputVCFSource'} eq 'CaVEMan'){
				$caveVer = $plln->[0]->{'InputVCFVer'};
			}
		}
	}
	#Lastly, perhaps the VCF version is old so we need to cycle through each header ourselves to find the right one :-/.
	if(!defined($caveVer) || $caveVer eq ""){
		my $head = $vcf->format_header();
		my @plLines = split(/\n/,$head);
		foreach my $plln(@plLines){
			next unless($plln =~ m/^\s*#+vcfProcessLog.+InputVCFSource=<CaVEMan>.+/);
			$plln =~ s/^\s*#+//g;
			my @ln = split(/=/,$plln,2);
			my $ky = $ln[0];
			my @rest = split(/>,/,$ln[1]);
			#Look through the body to find the version
			foreach my $rt(@rest){
				next unless $rt =~ m/InputVCFVer/;
				$rt =~ s/<//g;
				my $tmp;
				(undef,$caveVer) = split(/=/,$rt);
				last;
			}
		}
	}
	if(!defined($caveVer) || $caveVer eq ""){
		$caveVer = ".";
	}	return $caveVer;
}

1;
