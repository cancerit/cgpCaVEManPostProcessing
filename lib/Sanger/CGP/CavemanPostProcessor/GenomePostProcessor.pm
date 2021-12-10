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

package Sanger::CGP::CavemanPostProcessor::GenomePostProcessor;

use strict;
use Carp;
use Const::Fast qw(const);

use Sanger::CGP::CavemanPostProcessor;
our $VERSION = Sanger::CGP::CavemanPostProcessor->VERSION;

use base qw(Sanger::CGP::CavemanPostProcessor::PostProcessor);

const my $MAX_MATCHED_NORM_MUT_ALLELE_PROP => 0.05;
const my $MAX_PHASING_MINORITY_STRAND_PROP => 0.04;
const my $RD_POS_BEGINNING_OF_RD_PROP => 0.08;
const my $RD_POS_END_OF_TWOTHIRDS_EXTEND_PROP => 0.08;
const my $MIN_PASS_AVG_QUAL_PENTAMER => 20;
const my $SAME_RD_POS_PERCENT => 80;
const my $MAX_TUM_INDEL_PROP => 10;
const my $MAX_NORM_INDEL_PROP => 10;
const my $MIN_AVG_MAP_QUAL => 21;
const my $MIN_AVG_PHASING_BASE_QUAL => 21;
const my $MIN_DEPTH_QUAL => 25;
const my $MIN_NORM_MUT_ALLELE_BASE_QUAL => 15;
const my $MIN_RD_POS_DEPTH => 8;
const my $MATCHED_NORMAL_ALLELE_HICVG_CUTOFF => 2;
const my $MAX_MATCHED_NORMAL_ALLELE_HICVG_PROPORTION => 0.03;

#---------------
#	Init methods
#---------------

sub _init{
	my ($self,$inputs) = @_;
	$self->matchedNormalAlleleHiCvgCutoff($inputs->{'matchedNormalAlleleHiCvgCutoff'});
	$self->maxMatchedNormalAlleleHiCvgProportion($inputs->{'maxMatchedNormalAlleleHiCvgProportion'});
	$self->minSingleEndCoverage($inputs->{'minSingleEndCoverage'});
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
			$self->{'mnahcc'} = $MATCHED_NORMAL_ALLELE_HICVG_CUTOFF;
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
			$self->{'mmnahcvp'} = $MAX_MATCHED_NORMAL_ALLELE_HICVG_PROPORTION;
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

	foreach (chromosome){
		foreach(mutant position){
			$processor->runProcess($chr,$start,$stop,$refBase,$mutBase);
			if($processor->getDepthResult == 1 &&
						$processor->getReadPositionResult == 1 &&
						$processor->getNormMutsAllelesResult == 1 &&
						$postProcessor->getUnmatchedNormalResult == 1){
				$pass = 1;
			}
		}
	}



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
