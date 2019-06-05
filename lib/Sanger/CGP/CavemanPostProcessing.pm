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
use base 'Exporter';

our $VERSION = '1.8.8';
our @EXPORT = qw($VERSION);

const my $MATCH_CIG => 'M';
const my $SKIP_CIG => 'N';
const my $INS_CIG => 'I';
const my $DEL_CIG => 'D';
const my $SOFT_CLIP_CIG => 'S';
const my $HARD_CLIP_CIG => 'H';

const my $MIN_SINGLE_END_CVG => 10;
const my $MATCHED_NORMAL_MAX_MUT_PROP => 0.2;

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

sub runProcess{
	my ($self,$chr,$start,$stop,$refBase,$mutBase) = @_;
	$muts = undef;
	$norms = undef;
	$muts_rds = {};
	$norms_rds = {};
	$self->clearResults();
	$self->_chromosome($chr);
	$self->_currentPos($start);
	$self->_refBase($refBase);
	$self->_mutBase($mutBase);
	$self->{'tb'}->fetch($chr.':'.$start.'-'.$stop,\&_callbackTumFetch);
	$self->{'nb'}->fetch($chr.':'.$start.'-'.$stop,\&_callbackMatchedNormFetch);
	return 1;
}

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

sub _callbackTumFetch{
	my ($algn) = @_;
	my $flagValue = $algn->flag;
	#Check read and mate are mapped. If not return.
	return if((int($flagValue) & 8) != 0);
	return if((int($flagValue) & 4) != 0);
	#Check for duplicate status
	return if((int($flagValue) & 256) != 0);
	return if((int($flagValue) & 512) != 0);
	return if((int($flagValue) & 1024) != 0);
	return if((int($flagValue) & 2048) != 0); #Exclude supplementary alignments
	#Quick check that were covering the base with this read (skips/indels are ignored)
	if(_isCurrentPosCoveredFromAlignment($algn) == 1){
		#Get the correct read position.
		my ($rdPosIndexOfInterest,$currentRefPos) = _getReadPositionFromAlignment($algn);
		#print $rdPosIndexOfInterest,"\n";
		#print $algn->qseq,"\n";
		my @splt = split(//,$algn->qseq);
  		#Calculate other stuff
		my $totalPCovg = 0;
		my $totalNCovg = 0;
		my $indelRdCount = 0;
		my $nom = $algn->qname;
		my $start = $algn->start;
		#Read strand
		my $str = 1;
		if($algn->reversed){
			$str = -1;
		}
		return unless ($algn->proper_pair == 1);
    # Ensure that we keep
    return if((int($flagValue) & 16) != 0 && (int($flagValue) & 32) != 0);
    return if((int($flagValue) & 16) == 0 && (int($flagValue) & 32) == 0);

		$muts->{'totalTCoverage'} += 1;
		if($str == 1){
			$muts->{'totalTCoveragePos'} += 1;
		}else{
			$muts->{'totalTCoverageNeg'} += 1;
		}
		my $xt = $algn->aux_get('XT');
		if($algn->cigar_str =~ m/[ID]/){
			$muts->{'indelTCount'} += 1;
		}

		#Read base
		my $qbase = $splt[$rdPosIndexOfInterest-1];

		#Base quality
		my $qscore = $algn->qscore->[$rdPosIndexOfInterest-1];

		push(@{$muts->{'completeMutStrands'}},$str);

		push(@{$muts->{'allTumBases'}},$qbase);

		push(@{$muts->{'allTumBaseQuals'}},$qscore);

		push(@{$muts->{'allTumStrands'}},$str);

		#return if(uc($qbase) ne uc($mutBase));

		return if ($keepSW == 0 && defined($xt) && $xt eq 'M');

		return if($qscore < $minAnalysedQual);

		my $rdPos = $rdPosIndexOfInterest;
		my $ln = $algn->l_qseq;

		$muts->{'tumcvg'} += 1;

		#return if(uc($qbase) ne uc($mutBase));

		if($str == 1){
			$totalPCovg++;
		}else{
			$totalNCovg++;
			$rdPos = ($ln - $rdPos) + 1;
		}

		$muts->{'pcvg'} += $totalPCovg;

		$muts->{'ncvg'} += $totalNCovg;

		my $rdName = $algn->qname;

		return if(uc($qbase) ne uc($mutBase));

		my $softclipcount = _get_soft_clip_count_from_cigar($algn->cigar_array);
		my $primaryalnscore = $algn->get_tag_values('AS');

		#Tum quals
		push(@{$muts->{'tqs'}},$qscore);

		#Tum Rd Pos
		push(@{$muts->{'trp'}},$rdPos);

		#Tum rd length
		push(@{$muts->{'trl'}},$ln);

		#Tum XT tags
		push(@{$muts->{'txt'}},$xt);

		#Tum rd start
		push(@{$muts->{'trdst'}},$start);

		#Strands
		push(@{$muts->{'tstr'}},$str);

		#RdNames
		push(@{$muts->{'trn'}},$rdName);

		#Mapping quals
		push(@{$muts->{'tmq'}},$algn->qual);

		#AlnScoresPrm
		push(@{$muts->{'alnp'}},$primaryalnscore);

		#Softclipping
		push(@{$muts->{'sclp'}},$softclipcount);

		#print Dumper($a);
	}
	return 1;
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

sub _callbackMatchedNormFetch{
		my ($algn) = @_;
	my $flagValue = $algn->flag;
	#Check read and mate are mapped.
	return if((int($flagValue) & 8) != 0);
	return if((int($flagValue) & 4) != 0);
	#Check for duplicate status
	return if((int($flagValue) & 256) != 0);
	return if((int($flagValue) & 512) != 0);
	return if((int($flagValue) & 1024) != 0);
	return if((int($flagValue) & 2048) != 0); #Exclude supplementary alignments
	#Quick check that were covering the base with this read (skips/indels are ignored)
	if(_isCurrentPosCoveredFromAlignment($algn) == 1){
		#Get the correct read position.
		my ($rdPosIndexOfInterest,$currentRefPos) = _getReadPositionFromAlignment($algn,$currentPos);
		#print $rdPosIndexOfInterest,"\n";
		#print $algn->qseq,"\n";
		my @splt = split(//,$algn->qseq);
  		#Calculate other stuff
		my $totalPCovg = 0;
		my $totalNCovg = 0;
		my $indelRdCount = 0;
		my $nom = $algn->qname;
		return unless ($algn->proper_pair == 1);
    # Ensure that we keep
    return if((int($flagValue) & 16) != 0 && (int($flagValue) & 32) != 0);
    return if((int($flagValue) & 16) == 0 && (int($flagValue) & 32) == 0);

		if(!defined($muts->{'totalNCoverage'})){
			$muts->{'totalNCoverage'} = 0;
		}
		$muts->{'totalNCoverage'} += 1;
		my $xt = $algn->aux_get('XT');
		#Read base
		my $qbase = $splt[$rdPosIndexOfInterest-1];

		#Read strand
		my $str = 1;
		if($algn->reversed){
			$str = -1;
		}
		#Base quality
		my $qscore = $algn->qscore->[$rdPosIndexOfInterest-1];

		if(!defined($muts->{'allNormBases'})){
			$muts->{'allNormBases'} = [];
		}
		push(@{$muts->{'allNormBases'}},$qbase);

		if(!defined($muts->{'allNormBaseQuals'})){
			$muts->{'allNormBaseQuals'} = [];
		}
		push(@{$muts->{'allNormBaseQuals'}},$qscore);

		if(!defined($muts->{'allNormStrands'})){
			$muts->{'allNormStrands'} = [];
		}
		push(@{$muts->{'allNormStrands'}},$str);

		return if ($keepSW == 0 && defined $xt && $xt eq 'M');

		return if($qscore < $minAnalysedQual);

		if(!defined($muts->{'normcvg'})){
			$muts->{'normcvg'} = 0;
		}
		$muts->{'normcvg'} += 1;


		#return if(uc($qbase) ne uc($mutBase));

		my $rdPos = $rdPosIndexOfInterest;

		my $ln = length($algn->qseq);

		if($str == +1){
			$totalPCovg++;
		}else{
			$totalNCovg++;
			$rdPos = ($ln - $rdPos) + 1;
		}
		my $rdName = $algn->qname;
		return if(uc($qbase) ne uc($mutBase));
		#Tum quals
		if(!defined($muts->{'nqs'})){
			my @empty = ();
			$muts->{'nqs'} = \@empty;
		}
		push(@{$muts->{'nqs'}},$qscore);

		#Tum Rd Pos
		if(!defined($muts->{'nrp'})){
			my @empty = ();
			$muts->{'nrp'} = \@empty;
		}
		push(@{$muts->{'nrp'}},$rdPos);

		#Tum rd length
		if(!defined($muts->{'nrl'})){
			my @empty = ();
			$muts->{'nrl'} = \@empty;
		}
		push(@{$muts->{'nrl'}},$ln);

		if(!defined($muts->{'npcvg'})){
			$muts->{'npcvg'} = 0;
		}
		$muts->{'npcvg'} += $totalPCovg;

		if(!defined($muts->{'nncvg'} )){
			$muts->{'nncvg'} = 0;
		}
		$muts->{'nncvg'} += $totalNCovg;


	}
	return 1;
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
