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

package Sanger::CGP::CavemanPostProcessor::GenomePostProcessor;

use strict;
use Carp;
use Const::Fast qw(const);

use Sanger::CGP::CavemanPostProcessor;
use Sanger::CGP::CavemanPostProcessor::Constants;

use parent qw(Sanger::CGP::CavemanPostProcessor::PostProcessor);

my $const = 'Sanger::CGP::CavemanPostProcessor::Constants';

sub new{
  my ($proto) = shift;
  my $class = ref($proto) || $proto;
  my %inputs = @_;
  my $self = $class->SUPER::new(%inputs);
  bless $self, $class;
  return $self;
}

#---------------
#	Init methods
#---------------
sub _init{
	my ($self,$inputs) = @_;
  my ($self,$inputs) = @_;
  if(!defined($inputs->{'tumBam'}) || !defined($inputs->{'normBam'})){
    croak("tumBam and normBam are required for initialisation.\n");
  }
  $self->tumBam($inputs->{'tumBam'}, $inputs->{'ref'});
  $self->normBam($inputs->{'normBam'}, $inputs->{'ref'});
  $self->keepSW($inputs->{'keepSW'}) if exists $inputs->{'keepSW'};
  $self->minAnalysedQual($inputs->{'minAnalysedQual'}) if exists $inputs->{'minAnalysedQual'};
  $self->matchedNormalAlleleHiCvgCutoff($inputs->{'matchedNormalAlleleHiCvgCutoff'});
  $self->maxMatchedNormalAlleleHiCvgProportion($inputs->{'maxMatchedNormalAlleleHiCvgProportion'});
  $self->SUPER::_init($inputs);
  return $self;
}

sub clearResults{
	my ($self) = @_;
	$self->SUPER::clearResults();
	return 1;
}

#-----------------------------
#	Post processing filter methods
#-----------------------------
#Overrides the postProcessor.pm subroutine as this is slightly scarier!
sub _checkNormMuts{
	my ($self) = @_;
	my $qualCnt = 0;
	my $call = "";
	foreach my $q(@{$self->_muts->{'nqs'}}){
		if($q >= $self->minNormalMutAlleleQual()){
			$qualCnt++;
		}
		my $proportion = $qualCnt / $self->_muts->{'totalNCoverage'};
		if($qualCnt > $self->matchedNormalAlleleHiCvgCutoff()){
			if($proportion > $self->maxMatchedNormalAlleleHiCvgProportion()){
				return 0;
			}
		}else{
			if($proportion > $self->maxMatchedNormalAlleleProportion()){
				return 0;
			}
		}

	}
	return 1;
}

#-----------------
#	Getters/setters
#-----------------
sub matchedNormalAlleleHiCvgCutoff{
	my ($self,$p) = @_;
	if(defined($p)){
		 $self->{'mnahcc'} = $p;
	}else{
		if(!defined($self->{'mnahcc'})){
			$self->{'mnahcc'} = $const->default_flag_values('MATCHED_NORMAL_ALLELE_HICVG_CUTOFF');
		}
	}
	return $self->{'mnahcc'};
}

sub maxMatchedNormalAlleleHiCvgProportion{
	my ($self,$p) = @_;
	if(defined($p)){
		 $self->{'mmnahcvp'} = $p;
	}else{
		if(!defined($self->{'mmnahcvp'})){
			$self->{'mmnahcvp'} = $const->default_flag_values('MAX_MATCHED_NORMAL_ALLELE_HICVG_PROPORTION');
		}
	}
	return $self->{'mmnahcvp'};
}

#----------
#	DESTROY!
#----------

sub DESTROY{
	my $self = shift;
	#warn "GenomePostProcessor::DESTROY\n";
	$self->SUPER::DESTROY;
}

return 1;

=head1 NAME

CavemanPostProcessor - Perl module for post processing CaVEMan data.

=head1 SYNOPSIS

  use CavemanPostProcessor;
  my $processor = CavemanPostProcessor->new(tumBam => 'tumBamPath', normBam => 'normBamPath'); #Required

  #Optional...
                        'minDepthQual' => 25,
                        'depthCutoffProportion' => 0.333333
                        'minNormMutAllelequal' => 20,
                        'maxNormalMutAlleleCount' => 1,
                        'minAnalysedQual' => 10,
                        'samePosMaxPercent' => 80,
                        'keepSW' => 1,
                        'maxTumIndelProportion' => 10,
                        'maxNormIndelProportion' => 10 ,
                        'pentamerMinPassAvgQual'  => 20,
                        'minPassPhaseQual'=> 21,
                        'minPassAvgMapQual' =>	,
                        'maxMatchedNormalAlleleProportion' => 0.05

=head1 DESCRIPTION

This module checks a mutant position in a sample against matched and unmatched samples, and returns a pass or fail (boolean) decision for various checks.
these include: depth check, normal sample mutant allele, read position and  unmatched normal mutant allele.
For details of each check see the relevant method documentation.

=head2 Methods

=over 4

=item * CavemanPostProcessor->new($tumBamFile,$normBamFile,$minDepthQual,$minNormMutAlleleQual,$maxNormalAlleleCount,
											$minAnalysedQual,$maxUnmatchedNormQualCutoff,$maxUnmatchedNormAlleleCount,$csvNormalList,$keepSW)

Creates and returns a new CavemanPostProcessor object. $tumBamFile,$normBamFile must be passed. The remainder can be null and will revert to their default values.
$tumBamFile - path to the tumour bam file.
$normBamFile - path to the normal bam file.

=item * $object->runProcess{$chr,$start,$stop,$refBase,$mutBase);

Returns 1. When passed a chromosome (in the same format as was used in mapping the bam file.),
genomic start, genomic stop, reference allele and mutant allele this method sets up the initial pile up required
for post processing.

=back

=head1 AUTHOR

David Jones (drj@sanger.ac.uk)

=head1 COPYRIGHT

=head1 SEE ALSO

perl(1).

=cut
