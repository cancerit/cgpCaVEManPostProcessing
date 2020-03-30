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
package Sanger::CGP::CavemanPostProcessor::PostProcessor;

use strict;
use Bio::DB::HTS;
use Bio::DB::HTS::Constants;
use Bio::DB::HTS::Alignment;
use POSIX qw(strftime ceil);
use List::Util qw (sum);
use Carp;
use Const::Fast qw(const);

use Data::Dumper;

use Sanger::CGP::CavemanPostProcessor;
our $VERSION = Sanger::CGP::CavemanPostProcessor->VERSION;

use base qw(Sanger::CGP::CavemanPostProcessor);

#Defaults for this post processing module
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

sub _init{
	my ($self,$inputs) = @_;

	if(!defined($inputs->{'tumBam'}) || !defined($inputs->{'normBam'})){
		croak("tumBam and normBam are required for initialisation.\n");
	}

	$self->minDepthQual($inputs->{'minDepthQual'});
	$self->minNormalMutAlleleQual($inputs->{'minNormMutAllelequal'});
	$self->percentageSamePos($inputs->{'samePosMaxPercent'});
	$self->maxTumIndelProportion($inputs->{'maxTumIndelProportion'});
	$self->maxNormIndelProportion($inputs->{'maxNormIndelProportion'});
	$self->minPassAvgMapQual($inputs->{'minPassAvgMapQual'});
	$self->pentamerMinPassAvgQual($inputs->{'pentamerMinPassAvgQual'});
	$self->minPassAvgBaseQualPhasing($inputs->{'minPassPhaseQual'});
	$self->maxPhasingMinorityStrandReadProportion($inputs->{'maxPhasingMinorityStrandReadProportion'});
	$self->maxMatchedNormalAlleleProportion($inputs->{'maxMatchedNormalAlleleProportion'});
	$self->readPosBeginningOfReadIgnoreProportion($inputs->{'readPosBeginningOfReadIgnoreProportion'});
	$self->readPosTwoThirdsOfReadExtendProportion($inputs->{'readPosTwoThirdsOfReadExtendProportion'});
	$self->minRdPosDepth($inputs->{'minRdPosDepth'});

	return $self;
}

sub clearResults{
	my ($self) = @_;
	$self->{'depth'} = undef;
	$self->{'pos'} = undef;
	$self->{'norm'} = undef;
	$self->{'other'} = undef;
	$self->{'dpos'} = undef;
	$self->{'indelTum'} = undef;
	$self->{'indelNorm'} = undef;
	$self->{'motif'} = undef;
	$self->{'mapQ'} = undef;
	$self->{'otherOLD'} = undef;
	$self->{'phase'} = undef;
	$self->{'XTCheck'} = undef;
	$self->{'single'} = undef;
	$self->{'umpropres'} = undef;
	$self->{'clipmed'} = undef;
	$self->{'alnmedrd'} = undef;
	$self->{'algnmed'} = undef;
	return 1;
}

#-----------------------------
#	Getters/setters
#-----------------------------

sub maxMatchedNormalAlleleProportion{
	my ($self,$p) = @_;
	if(defined($p)){
		 $self->{'maxMatchedNormalAlleleProportion'} = $p;
	}else{
		if(!defined($self->{'maxMatchedNormalAlleleProportion'})){
			$self->{'maxMatchedNormalAlleleProportion'} = $MAX_MATCHED_NORM_MUT_ALLELE_PROP;
		}
	}
	return $self->{'maxMatchedNormalAlleleProportion'};
}

sub maxPhasingMinorityStrandReadProportion{
	my ($self,$p) = @_;
	if(defined($p)){
		$self->{'maxPhasingMinorityStrandReadProportion'} = $p;
	}else{
		if(!defined($self->{'maxPhasingMinorityStrandReadProportion'})){
			$self->{'maxPhasingMinorityStrandReadProportion'} = $MAX_PHASING_MINORITY_STRAND_PROP;
		}
	}
	return $self->{'maxPhasingMinorityStrandReadProportion'};
}

sub readPosBeginningOfReadIgnoreProportion{
	my ($self,$p) = @_;
	if(defined($p)){
		 $self->{'readPosBeginningOfReadIgnoreProportion'} = $p;
	}else{
		if(!defined($self->{'readPosBeginningOfReadIgnoreProportion'})){
			$self->{'readPosBeginningOfReadIgnoreProportion'} = $RD_POS_BEGINNING_OF_RD_PROP;
		}
	}
	return $self->{'readPosBeginningOfReadIgnoreProportion'};
}

sub readPosTwoThirdsOfReadExtendProportion{
	my ($self,$p) = @_;
	if(defined($p)){
		 $self->{'readPosTwoThirdsOfReadExtendProportion'} = $p;
	}else{
		if(!defined($self->{'readPosTwoThirdsOfReadExtendProportion'})){
			$self->{'readPosTwoThirdsOfReadExtendProportion'} = $RD_POS_END_OF_TWOTHIRDS_EXTEND_PROP;
		}
	}
	return $self->{'readPosTwoThirdsOfReadExtendProportion'};
}


sub pentamerMinPassAvgQual{
	my ($self,$p) = @_;
	if(defined($p)){
		 $self->{'pentamerMinPassAvgQual'} = $p;
	}else{
		if(!defined($self->{'pentamerMinPassAvgQual'})){
			$self->{'pentamerMinPassAvgQual'} = $MIN_PASS_AVG_QUAL_PENTAMER;
		}
	}
	return $self->{'pentamerMinPassAvgQual'};
}

sub percentageSamePos{
	my ($self,$p) = @_;
	if(defined($p)){
		$self->{'pctSamePos'} = $p;
	}else{
		if(!defined($self->{'pctSamePos'})){
			$self->{'pctSamePos'} = $SAME_RD_POS_PERCENT;
		}
	}
	return $self->{'pctSamePos'};
}

sub maxTumIndelProportion{
	my ($self,$p) = @_;
	if(defined($p)){
		$self->{'maxTumIndelProp'} = $p;
	}else{
		if(!defined($self->{'maxTumIndelProp'})){
			$self->{'maxTumIndelProp'} = $MAX_TUM_INDEL_PROP;
		}
	}
	return $self->{'maxTumIndelProp'};
}

sub maxNormIndelProportion{
	my ($self,$p) = @_;
	if(defined($p)){
		$self->{'maxNormIndelProp'} = $p;
	}else{
		if(!defined($self->{'maxNormIndelProp'})){
			$self->{'maxNormIndelProp'} = $MAX_NORM_INDEL_PROP;
		}
	}
	return $self->{'maxNormIndelProp'};
}

sub minPassAvgMapQual{
	my ($self,$p) = @_;
	if(defined($p)){
		$self->{'minAvgMq'} = $p;
	}else{
		if(!defined($self->{'minAvgMq'})){
			$self->{'minAvgMq'} = $MIN_AVG_MAP_QUAL;
		}
	}
	return $self->{'minAvgMq'};
}

sub minPassAvgBaseQualPhasing{
	my ($self,$bq) = @_;
	if(defined($bq)){
		$self->{'phaseQual'} = $bq;
	}elsif(!defined($self->{'phaseQual'})){
		$self->{'phaseQual'} = $MIN_AVG_PHASING_BASE_QUAL;
	}
	return $self->{'phaseQual'};
}

sub minDepthQual{
	my ($self,$q) = @_;
	if(defined($q)){
		$self->{'d'} = $q;
	}else{
		if(!defined($self->{'d'})){
			$self->{'d'} = $MIN_DEPTH_QUAL;
		}
	}
	return $self->{'d'};
}

sub minNormalMutAlleleQual{
	my ($self,$q) = @_;
	if(defined($q)){
		$self->{'q'} = $q;
	}else{
		if(!defined($self->{'q'})){
			$self->{'q'} = $MIN_NORM_MUT_ALLELE_BASE_QUAL;
		}
	}
	return $self->{'q'};
}

sub minRdPosDepth{
	my ($self,$q) = @_;
	if(defined($q)){
		$self->{'minRdPsDpth'} = $q;
	}else{
		if(!defined($self->{'minRdPsDpth'})){
			$self->{'minRdPsDpth'} = $MIN_RD_POS_DEPTH;
		}
	}
	return $self->{'minRdPsDpth'};
}

#-----------------------------
#	Post processing tests
#-----------------------------

sub getTumIndelReadDepthResult{
	my ($self) = @_;
	if(!defined($self->{'indelTum'})){
		$self->{'indelTum'} = $self->_mutIndelCheck();
	}
	return $self->{'indelTum'};
}

sub getPentamerResult{
	my ($self) = @_;
	if(!defined($self->{'motif'})){
		$self->{'motif'} = $self->_evaluatePentamerCheck();
	}
	return $self->{'motif'};
}

sub getPhasingResult{
	my ($self) = @_;
	if(!defined($self->{'phase'})){
		$self->{'phase'} = $self->_runPhasingCheck();
	}
	return $self->{'phase'};
}

sub _runPhasingCheck{
	my ($self) = @_;
	#Count mut bases on each strand.
	my ($f_q_sum, $r_q_sum, $f_count, $r_count) = (0,0,0,0);
	for(my $i=0;$i<scalar(@{$self->_muts->{'allTumStrands'}});$i++){
		if($self->_muts->{'allTumStrands'}->[$i] == 1){
			$f_count++;
	    	$f_q_sum += $self->_muts->{'allTumBaseQuals'}->[$i];
	  	}else {
	   	$r_count++;
	    	$r_q_sum += $self->_muts->{'allTumBaseQuals'}->[$i];
	  	}
	}
	#Calcualte the average qualities of mutant bases.
	my ($f_av_qual,$r_av_qual) = ('.','.');
	if($f_count > 0) {
	  $f_av_qual = $f_q_sum / $f_count;
	}
	if($r_count > 0) {
	  $r_av_qual = $r_q_sum / $r_count;
	}

	#Get proportions for strand
	#Use the counts and avg quals to evaluate phasing
	#Only on fwd strand
	if($f_count > 0 && ($r_count == 0 || ($r_count / ($self->_getTotalReadsOnStrandCount(-1))) <= $self->maxPhasingMinorityStrandReadProportion())){
		if($f_av_qual < $self->minPassAvgBaseQualPhasing){
			return 0;
		}
	}elsif($r_count > 0 && ($f_count == 0 || ($f_count / $self->_getTotalReadsOnStrandCount(1)) <= $self->maxPhasingMinorityStrandReadProportion())){#Only on rev strand
		if($r_av_qual < $self->minPassAvgBaseQualPhasing){
			return 0;
		}
	}
	return 1;
}

sub _getTotalReadsOnStrandCount{
	my ($self,$strand) = @_;
	my $strandCount = 0;
	foreach my $str(@{$self->_muts->{'completeMutStrands'}}){
		if($strand == $str){
			$strandCount++;
		}
	}
	return $strandCount;
}

sub getNormIndelReadDepthResult{
	my ($self) = @_;
	if(!defined($self->{'indelNorm'})){
		$self->{'indelNorm'} = $self->_normIndelCheck();
	}
	return $self->{'indelNorm'};
}

sub _normIndelCheck{
	my ($self) = @_;
	if($self->_muts->{'totalNCoverage'} == 0){
		return 1;
	}
	my $prop = ($self->_muts->{'indelNCount'} / $self->_muts->{'totalNCoverage'}) * 100;
	if($prop <= $self->maxNormIndelProportion){
		return 1;
	}
	return 0;

}

sub _mutIndelCheck{
	my ($self) = @_;
	if(!exists $self->_muts->{'totalTCoverage'} || $self->_muts->{'totalTCoverage'} == 0){
		return 1;
	}
	my $prop = ($self->_muts->{'indelTCount'} / $self->_muts->{'totalTCoverage'}) * 100;
	if($prop <= $self->maxTumIndelProportion){
		return 1;
	}
	return 0;

}

sub getDepthResult{
	my ($self) = @_;
	if(!defined($self->{'depth'})){
		$self->{'depth'} = $self->_checkDepth();
	}
	return $self->{'depth'};
}

sub getDifferingReadPositionResult{
	my ($self) = @_;
	if(!defined($self->{'dpos'})){
		$self->{'dpos'} = $self->_checkDiffReadPos();
	}
	return $self->{'dpos'};
}

sub getAvgMapQualResult{
	my ($self) = @_;
	if(!defined($self->{'mapQ'})){
		$self->{'mapQ'} = $self->_checkAvgMapQual();
	}
	return $self->{'mapQ'};
}

sub _checkAvgMapQual{
	my ($self) = @_;
	my $total = 0;
	my $noOfMQs = 0;
	foreach my $mq(@{$self->_muts->{'tmq'}}){
		$total += $mq;
		$noOfMQs++;
	}
	if($noOfMQs == 0){
		#warn ("checkAvgMapQual: Should not have encountered 0 mutant allele reads in the tumour for ",$self->_chromosome,":",$self->_currentPos,", returning pass.\n");
		return 0;
	}
	#If the mean is less than the min mq for a pass we fail.
	return 0 if(($total / $noOfMQs) < $self->minPassAvgMapQual);
	return 1;
}

sub _evaluatePentamerCheck{
	my ($self) = @_;
	my $minus = 0;
	my $plus = 0;
	#Check strands first.
	foreach my $str(@{$self->_muts->{'tstr'}}){
		if($str == 1){
			$plus++;
		}else{
			$minus++;
		}
		return 1 if($plus > 1 && $minus > 1);
	}

	#If all or (all - 1) mut allele reads are on one strand continue.
	my $sz = scalar(@{$self->_muts->{'tstr'}});
	if(!(($minus >= ($sz -1) && $plus <= 1) || ($plus >= ($sz -1) && $minus <= 1)) ){
		return 1;
	}

	my $avgBQOverall = 0;
	my $cntAvg = 0;
	my @trp = @{$self->_muts->{'trp'}};
	my @trl = @{$self->_muts->{'trl'}};
	my @tstr = @{$self->_muts->{'tstr'}};
	my @trn = @{$self->_muts->{'trn'}};
	#We've checked 3rds of read and the strand,
	#if we got this far we have to check the entire mutant read.
	#Read names are handily stored in $self->_muts->{'trn'} so we can look at each in turn.
	my $avgBQ = 0;
	my $analysed = 0;

	my %want_rds = map { $_ => 1 } @{$self->_muts->{'trn'}};
	my %aligns;
	#Fetch the reads we want
	my $t_bam = $self->{'tb'};
	$t_bam->hts_index->fetch(
		$t_bam->hts_file,
		$t_bam->header->parse_region($self->_chromosome.":".$self->_currentPos."-".$self->_currentPos),
		sub{
			my $a = shift;
			if(exists $want_rds{$a->qname}) {
				$aligns{$a->qname} = $a
			}
		}
	);

	#Iterate through each read name
	for(my $i=0;$i<$sz;$i++){
		my $pos = $self->_muts->{'trp'}->[$i];
		my $rdName = $self->_muts->{'trn'}->[$i];
		#If position is not in last 3rd we can return now.
		if($pos < ($self->_muts->{'trl'}->[$i] / 2)){ # trl is read length
			return 1;
		}

		my $rd = $aligns{$rdName};

		my $isReversed = $rd->reversed;
		my $readPosOfMut = ($self->_currentPos - ($rd->pos + 1)) + 1;
		my $seq = $rd->qseq;
		my $quals = $rd->qscore();
		my @matches = ();
		if($isReversed == 1){
			$seq =~ tr/acgtnACGTN/tgcanTGCAN/;
			$seq = reverse($seq);
			#Reverse the quality array
			@$quals = reverse(@$quals);
			$readPosOfMut = (length($seq) - $readPosOfMut) + 1;
		}
		#Check for motif match in second half of the read (rev comped as required).
		#No motif, so skip
		next if($seq !~ m/GGC[AT]G/);
		while($seq =~ m/GGC[AT]G/g){
			push(@matches,length($`).",".length($&)."");
		}

		my $lastPos = 0;

		my $halfLength = (length($seq)/2);
		foreach my $match(@matches){
			my @split = split(/,/,$match);
			next if($split[0] < $halfLength);
			last if($split[0]+ $split[1] >= $readPosOfMut);
			#This will eventually give us the last occurrence of the motif before the mutation
			if($split[0]+ $split[1] < $readPosOfMut){
				#Last pos of the match is the start pos + the length - 1
				$lastPos = ($split[0] + $split[1]);
			}
		}

		next if($lastPos <= 0);
		#Finally we check the average base quality after the motif.
		my $avg = $self->_calcualteMeanBaseQualAfterMotif($lastPos+1,$quals);
		$avgBQOverall += $avg;
		$cntAvg++;
	}#End of iteration through each read name
	return 1 if($cntAvg < 1);
	if(($avgBQOverall / $cntAvg) < $self->pentamerMinPassAvgQual()){
		return 0;
	}
	return 1;

}

sub _calcualteMeanBaseQualAfterMotif{
	my ($self,$startPos,$quals) = @_;
	my $qSum = 0;
	my $qCount = 0;
	for(my $i=$startPos-1;$i<scalar(@$quals);$i++){
		$qSum += $quals->[$i];
		$qCount++;
	}
	if($qCount == 0){
		return 0;
	}
	return ($qSum / $qCount);
}

sub _checkDiffReadPos{
	my ($self) = @_;
	my $tmp;
	my $noOfReads = scalar(@{$self->_muts->{'trp'}});
	my $permittedNo = $noOfReads  * ( $self->percentageSamePos / 100 );
	for(my $i=0;$i<scalar(@{$self->_muts->{'trp'}});$i++){
		my $rdPos = $self->_muts->{'trp'}->[$i];
		if(!defined($tmp->{$rdPos})){
			$tmp->{$rdPos} = 0;
		}
		$tmp->{$rdPos}++;
		if($tmp->{$rdPos} > $permittedNo){
			return 0;
		}
	}
	return 1;
}

sub _checkDepth{
	my ($self) = @_;
	my $depth = scalar(@{$self->_muts->{'tqs'}});
	my $overCutoff = 0;
	foreach my $q(@{$self->_muts->{'tqs'}}){
		if($q >= $self->minDepthQual){
			$overCutoff++;
		}
		if($overCutoff >= ($depth / 3)){
			return 1;
		}
	}
	return 0;
}

sub getReadPositionResult{
	my ($self) = @_;
	if(!defined($self->{'pos'})){
		$self->{'pos'} = $self->_checkReadPos();
	}
	return $self->{'pos'};
}

sub _checkReadPos{
	my ($self) = @_;
	my $sec3rd = 0;
	my $first3rd = 0;
	if(scalar(@{$self->_muts->{'trp'}}) > $self->{'minRdPsDpth'}){
		return 1;
	}
	for(my $i=0;$i<scalar(@{$self->_muts->{'trp'}});$i++){
		my $rdLn = $self->_muts->{'trl'}->[$i];
		my $thirds = $rdLn / 3;
		my $rdPos = $self->_muts->{'trp'}->[$i];
		my $prop = ($rdLn * $self->{'readPosBeginningOfReadIgnoreProportion'});
		my $halfprop = ($rdLn * $self->{'readPosTwoThirdsOfReadExtendProportion'});
		#In first or second third, but not first n%, extending by m%
		if($rdPos <= (($thirds * 2)+$halfprop) && $rdPos > $prop){
			return 1;
		}
	}
	return 0;

}

sub getNormMutsAllelesResult{
	my ($self) = @_;
	if(!defined($self->{'norm'})){
		$self->{'norm'} = $self->_checkNormMuts();
	}
	return $self->{'norm'};
}

sub _checkNormMuts{
	my ($self) = @_;
	my $qualCnt = 0;
	my $call = "";
	foreach my $q(@{$self->_muts->{'nqs'}}){
		if($q >= $self->minNormalMutAlleleQual()){
			$qualCnt++;
		}
		my $proportion = $qualCnt / $self->_muts->{'totalNCoverage'};
		if($proportion > $self->maxMatchedNormalAlleleProportion()){
			return 0;
		}

	}
	return 1;
}

sub median{
	return (sum( ( sort { $a <=> $b } @_ )[ int( $#_/2 ), ceil( $#_/2 ) ] )/2);
}

sub _checkMedianClipping{
	my ($self) = @_;
	return sprintf('%.2f',median(@{$self->_muts->{'sclp'}}));
}

sub _calcPrimAlignmentScoreReadAdjustedMedian{
  my ($self) = @_;
  my $adj;
  for(my $i=0; $i<scalar(@{$self->_muts->{'alnp'}}); $i++){
    $adj->[$i] = $self->_muts->{'alnp'}->[$i] / $self->_muts->{'trl'}->[$i];
  }
  return sprintf('%.2f',median(@{$adj}));
}

sub getClipMedianResult{
	my ($self) = @_;
	if(!defined($self->{'clipmed'})){
		$self->{'clipmed'} = $self->_checkMedianClipping();
	}
	return $self->{'clipmed'};
}

sub getAlignmentScoreMedianReadAdjusted{
  my ($self) = @_;
	if(!defined($self->{'alnmedrd'})){
    $self->{'alnmedrd'} = $self->_calcPrimAlignmentScoreReadAdjustedMedian();
  }
  return $self->{'alnmedrd'};
}

sub getAlignmentScoreMedian{
  my ($self) = @_;
  if(!defined($self->{'algnmed'})){
    $self->{'algnmed'} = $self->_calcPrimAlignmentScoreMedian();
  }
  return $self->{'algnmed'};
}

sub _calcPrimAlignmentScoreMedian{
  my ($self) = @_;
	return sprintf('%.2f',median(@{$self->_muts->{'alnp'}}));
}

#----------
#	DESTROY
#----------

sub DESTROY{
	my $self = shift;
	#warn "PostProcessor::DESTROY\n";
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
  											'minNormMutAllelequal' => 20,
  											'maxNormalMutAlleleCount' => 1,
  											'minAnalysedQual' => 10,
  											'minUnmatchedNormalQual' => 20,
  											'maxUnmatchedNormAlleleCount' => 2,
  											'unmatchedNormalSampleList' => 'sample1.bam,sample2.bam,sample3.bam',
  											'samePosMaxPercent' => 80,
  											'keepSW' => 1,
  											'maxTumIndelProportion' => 10,
  											'maxNormIndelProportion' => 10 ,
  											'pentamerMinPassAvgQual'  => 20,
  											'minPassPhaseQual'=> 21,
  											'minPassAvgMapQual' =>

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

=item * $object->clearResults()

Returns 1. Removes any results from the previous analysis. This is automatically run every time runProcess is called, so is only there for convenience.


=item * $object->tumBam($bam)

Returns the file location of the tumour sample bam file.
If $bam is passed the location is set to this value.
Returns undef if unset.

=item * $object->normBam($bam)

Returns the file location of the normal sample bam file.
If $bam is passed the location is set to this value.
Returns undef if unset.

=item * $object->unmatchedNormalSampleList($list)

Returns an array reference of Bio::DB::HTS objects, one pointed at each of the passed files.
If $list is passed (as a comma separated string list) this method will create a Bio::DB::HTS object
for each file in the UNMATCHED normal list. If unset, returns an empty array pointer.

=item * $object->minDepthQual(quality)

Returns the minimum quality permitted before a mutant allele is counted as high quality in the depth check.
Returns the set value of quality. Returns the default 25 if previously unset.

=item * $object->getDepthResult()

Returns 1 (pass)
IF
	at least 1/3 of mutant alleles are of base quality >= minDepthQual
OTHERWISE return 0

=item * $object->getReadPositionResult()

Returns 1 (pass)
IF
	coverage >= 10 and at least one mutant allele in the tumour is in the middle third of a read
OR
	coverage < 10 and at least one mutant allele in the tumour is in the first or middle third of a read
OTHERWISE
	Returns 0 (fail)

=item * $object->minAnalysedQual(qual)

Sets and returns the minimum analysed quality base. If unset or qual is not passed returns the default 11.

=item * $object->maxNormalMutAlleleCount(count)

Sets and returns the maximum number of high quality (see minNormalMutAlleleQual) mutant alleles permitted in the matched normal.
If unset returns the default 1

=item * $object->minNormalMutAlleleQual(qual)

Sets and returns the minimum base quality a matched normal mutant allele must have to be counted in the normal allele check.
If unset returns the default 20

=item * $object->getNormMutsAllelesResult()

Returns 1 (pass)
IF
	there are no more than $maxNormalMutAlleleCount mutant alleles of quality >= $minNormalMutAlleleQual in the matched normal.
OTHERWISE
	Returns 0 (fail)

=item * $object->getPhasingResult()

Returns 1 (pass)
IF
	The mutant is found on both strands,
OR
	The mutant is found on one strand, but the average mutant base quality > $object->minPassAvgBaseQualPhasing()
OTHERWISE
	Returns 0 (fail)


=item * $object->minPassAvgBaseQualPhasing()
Sets and returns the minimum average (mean) mutant base quality required for a single stranded mutant to pass the phasing flag.
Default = 21


=item * minPassAvgMapQual
Sets and returns the minimum average mapping qual to pass the map qual flag.
Default = 21

=item * minRdPosStart
Sets and returns the minimum read pos start for position on read.


=back

=head1 AUTHOR

David Jones (drj@sanger.ac.uk)

=head1 COPYRIGHT

=head1 SEE ALSO

perl(1).
Bio::DB::SAM

=cut
