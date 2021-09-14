##########LICENCE##########
# Copyright (c) 2014-2019 Genome Research Ltd.
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

package Sanger::CGP::CavemanPostProcessing;

use strict;
use Bio::DB::HTS;
use Bio::DB::HTS::Constants;
use Bio::DB::HTS::Alignment;
use POSIX qw(strftime);
use Carp;
use Const::Fast qw(const);
use Attribute::Abstract;
use Data::Dumper;
use List::Util qw(min max);
use base 'Exporter';

our $VERSION = '1.9.3';
our @EXPORT = qw($VERSION);

const my $MATCH_CIG => 'M';
const my $SKIP_CIG => 'N';
const my $INS_CIG => 'I';
const my $DEL_CIG => 'D';
const my $SOFT_CLIP_CIG => 'S';
const my $HARD_CLIP_CIG => 'H';

const my $MIN_SINGLE_END_CVG => 10;
const my $MATCHED_NORMAL_MAX_MUT_PROP => 0.2;
const my $CAVEMAN_MATCHED_NORMAL_MAX_MUT_PROP => 0.2;

my $muts;
my $norms;
my $muts_rds;
my $norms_rds;
my $currentPos;
my $refBase;
my $mutBase;
my $keepSW = 0;
my $minAnalysedQual = 11;

sub new {
	my ($proto) = shift;
	my %inputs = @_;
	my $class = ref($proto) || $proto;
 	my $self = {};
  bless($self, $class);
	$self->_init_base(\%inputs);
	$self->_init(\%inputs);
	return $self;
}


sub _init_base{
	my ($self,$inputs) = @_;
	if(!defined($inputs->{'tumBam'}) || !defined($inputs->{'normBam'})){
		croak("tumBam and normBam are required for initialisation.\n");
	}
	$self->tumBam($inputs->{'tumBam'});
	$self->normBam($inputs->{'normBam'});
	$self->keepSW($inputs->{'keepSW'});
	$self->minAnalysedQual($inputs->{'minAnalysedQual'});
	return $self;
}


=item _init
	_init method required by inheriting classes.
=cut
sub _init : Abstract;

=item clearResults
	Clears last positions worth of data from the stored results.
=cut
sub clearResults : Abstract;

=item keepSW
  	Whether to include Smith-Waterman aligned reads in post processing
=cut
sub keepSW{
	my ($self,$keep) = @_;
	if(defined($keep) && ($keep == 1 || $keep == 0)){
		$keepSW = $keep;
	}
	return $keepSW;
}

sub getSingleEndResult{
	my ($self) = @_;
	if(!defined($self->{'single'})){
		$self->{'single'} = $self->_calculateSingleEndResult();
	}
	return $self->{'single'};
}

sub _calculateSingleEndResult{
	my ($self) = @_;
	return 1 if($self->_muts->{'pcvg'} < $self->minSingleEndCoverage() || $self->_muts->{'ncvg'} < $self->minSingleEndCoverage());
	my $hasPos = 0;
	my $hasNeg = 0;
	foreach my $str(@{$self->_muts->{'tstr'}}){
		if($str == -1){
			$hasNeg++;
		}elsif($str == 1){
			$hasPos++;
		}
		return 1 if($hasNeg > 0 && $hasPos > 0);
	}
	if($hasNeg == 0 || $hasPos == 0){
		return 0;
	}
	return 1;
}

sub minSingleEndCoverage{
	my ($self,$p) = @_;
	if(defined($p)){
		 $self->{'sec'} = $p;
	}else{
		if(!defined($self->{'sec'})){
			$self->{'sec'} = $MIN_SINGLE_END_CVG;
		}
	}
	return $self->{'sec'};
}

sub getMatchedNormalProportionResult{
	my ($self) = @_;
	if(!defined($self->{'umpropres'})){
		$self->{'umpropres'} = $self->_calculateMatchedNormalProportion();
	}
	return $self->{'umpropres'};
}

sub _calculateMatchedNormalProportion{
	my ($self) = @_;
	#Calculate tumour proportion of mut allele
	my $tumProp = 0;
	if(scalar(@{$self->_muts->{'tqs'}}) > 0){
		$tumProp = scalar(@{$self->_muts->{'tqs'}})/$self->_muts->{'tumcvg'};
	}
	#Calculate normal proportion of mut allele
	my $normProp = 0;
	if(exists($self->_muts->{'nqs'}) && scalar(@{$self->_muts->{'nqs'}}) > 0){
		$normProp = scalar(@{$self->_muts->{'nqs'}})/$self->_muts->{'normcvg'};
	}
	#Fail if the difference is less than the given proportion/percentage
	return 0 if($normProp > 0 && ($tumProp - $normProp) < $self->matchedNormalMaxMutProportion());
	return 1;
}

sub maxCavemanMatchedNormalProportion{
    my ($self, $val) = @_;
    if(defined($val)){
		 $self->{'cmnmmp'} = $val;
	}else{
		if(!defined($self->{'cmnmmp'})){
			$self->{'cmnmmp'} = $CAVEMAN_MATCHED_NORMAL_MAX_MUT_PROP;
		}
	}
	return $self->{'cmnmmp'};
}

sub matchedNormalMaxMutProportion{
	my ($self,$p) = @_;
	if(defined($p)){
		 $self->{'mnmmp'} = $p;
	}else{
		if(!defined($self->{'mnmmp'})){
			$self->{'mnmmp'} = $MATCHED_NORMAL_MAX_MUT_PROP;
		}
	}
	return $self->{'mnmmp'};
}

=item minAnalysedQual
	Holds the minimum base quality of reads to be used in analysis.
=cut
sub minAnalysedQual{
	my ($self,$q) = @_;
	if(defined($q)){
		$minAnalysedQual = $q;
	}else{
		if(!defined($minAnalysedQual)){
			$minAnalysedQual = 11;
		}
	}
	return $minAnalysedQual;
}

sub _currentPos{
	my ($self,$pos) = @_;
	if(defined($pos)){
		$currentPos = $pos;
	}
	return $currentPos;
}

sub _chromosome{
	my ($self,$c) = @_;
	if(defined($c)){
		$self->{'chro'} = $c;
	}
	return $self->{'chro'};
}

sub _refBase{
	my ($self,$b) = @_;
	if(defined($b)){
		$refBase = $b;
	}
	return $refBase;
}

sub _mutBase{
	my ($self,$b) = @_;
	if(defined($b)){
		$mutBase = $b;
	}
	return $mutBase;
}

sub tumBam{
	my ($self,$bam) = @_;
	if(defined($bam)){
		$self->{'tb'} = Bio::DB::HTS->new(-bam=>$bam);
	}
	return $self->{'tb'};
}

sub normBam{
	my ($self,$bam) = @_;
	if(defined($bam)){
		$self->{'nb'} = Bio::DB::HTS->new(-bam=>$bam);
	}
	return $self->{'nb'};
}

sub _muts{
	if(!defined($muts)){
		$muts = {};
	}
	return $muts;
}

sub _norms{
	my ($self,$new) = @_;
	if(defined($new)){
		$norms = $new;
	}
	if(!defined($norms)){
		$norms = {};
	}
	return $norms;
}

sub _get_soft_clip_count_from_cigar{
	my ($cig_arr) = @_;
	my $count = 0;
	foreach my $cigentry(@$cig_arr){
		if($cigentry->[0] eq $SOFT_CLIP_CIG){
			$count += $cigentry->[1];
		}
	}
	return $count;
}

sub _getDistanceFromGapInRead{
  my ($cigar_array,$rdPosIndexOfInterest) = @_;
  my $min_gap_dist = -1;
  my $currentRp = 0;
  foreach my $cigSect(@{$cigar_array}){
    if($cigSect->[0] eq $MATCH_CIG || $cigSect->[0] eq $SKIP_CIG ||
          $cigSect->[0] eq $INS_CIG || $cigSect->[0] eq $SOFT_CLIP_CIG){
      $currentRp+=$cigSect->[1];
    }elsif($cigSect->[0] eq $DEL_CIG){
      my $dp_start = $currentRp+1;
      my $dp_end = $currentRp+$cigSect->[1];
      my $tmp_dist = max(abs($rdPosIndexOfInterest-$dp_start),abs($dp_end-$rdPosIndexOfInterest));
      if($tmp_dist < $min_gap_dist || $min_gap_dist == -1){
        $min_gap_dist = $tmp_dist;
      }
    }
  }
  return $min_gap_dist;
}

sub _getReadPositionFromAlignment{
	my ($algn) = @_;
	my $rdPosIndexOfInterest = 0;
    	my $currentRefPos = $algn->start -1;
		foreach my $cigSect(@{$algn->cigar_array}){
			if($cigSect->[0] eq $MATCH_CIG){
				if($currentRefPos <= $currentPos && ($currentRefPos+$cigSect->[1]) >= $currentPos){
					for(my $i=0;$i<$cigSect->[1];$i++){
						$rdPosIndexOfInterest++;
						$currentRefPos++;
						if($currentPos == $currentRefPos){
							return ($rdPosIndexOfInterest,$currentRefPos);
						}
					}
				}else{
					$rdPosIndexOfInterest += $cigSect->[1];
					$currentRefPos += $cigSect->[1];
				}
			}elsif($cigSect->[0] eq $DEL_CIG || $cigSect->[0] eq $SKIP_CIG){
				$currentRefPos += $cigSect->[1];
			}elsif($cigSect->[0] eq $INS_CIG || $cigSect->[0] eq $SOFT_CLIP_CIG){
				$rdPosIndexOfInterest += $cigSect->[1];
			}
		}
}

sub _isCurrentPosCoveredFromAlignment{
	my ($aln) = @_;
	my $pos = $aln->start - 1;
	foreach my $cigSect(@{$aln->cigar_array}){

		if($cigSect->[0] eq $MATCH_CIG){
			if($pos <= $currentPos && ($pos+$cigSect->[1]) >= $currentPos){
				return 1;
			}
			$pos+= $cigSect->[1];
		}elsif($cigSect->[0] eq $DEL_CIG || $cigSect->[0] eq $SKIP_CIG){
			if($pos <= $currentPos && ($pos+$cigSect->[1]) > $currentPos){
				return 0;
			}
			$pos+= $cigSect->[1];
		}
	}
	return 0;
}

sub DESTROY{
	my $self = shift;
	$currentPos = undef;
	$refBase = undef;
	$mutBase = undef;
	$keepSW = 0;
	$minAnalysedQual = 11;
	$muts = undef;
	$norms = undef;
	#warn "Base::DESTROY\n";
}

return 1;
