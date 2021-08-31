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
package Sanger::CGP::CavemanPostProcessing::Flagger;

use strict;
use warnings FATAL=>'all';
use autodie;
use Carp;
use Const::Fast;
use POSIX qw(strftime ceil);
use List::Util qw (sum zip max);


use Bio::DB::HTS;
use Bio::DB::HTS::Constants;

use Data::Dumper;

use Sanger::CGP::CavemanPostProcessing;

const my $BAD_BITS => hex RFLAGS->{UNMAPPED}|RFLAGS->{M_UNMAPPED}|RFLAGS->{NOT_PRIMARY}|RFLAGS->{QC_FAILED}|RFLAGS->{DUPLICATE}|RFLAGS->{SUPPLEMENTARY};
const my $GOOD_BITS => hex RFLAGS->{MAP_PAIR};

const my $VCF_COLUMN_NORMAL => 'NORMAL';
const my $VCF_COLUMN_TUMOUR => 'TUMOUR';
const my $VCF_COLUMN_FORMAT => 'FORMAT';
const my $OLD_ALLELE_VCF_FORMAT => 'GT:AA:CA:GA:TA:PM';
const my $NEW_ALLELE_VCF_FORMAT => 'GT:FAZ:FCZ:FGZ:FTZ:RAZ:RCZ:RGZ:RTZ:PM';
const my %OLD_ALLELE_VCF_FORMAT_INDEX_HASH => ('A' => 1, 'C' => 2, 'G' => 3, 'T' => 4, );
const my %NEW_ALLELE_VCF_FORMAT_INDEX_HASH => ('A'=>[1,5], 'C' =>[2,6], 'G'=>[3,7], 'T'=>[4,8], );

const my $MATCH_CIG => 'M';
const my $SKIP_CIG => 'N';
const my $INS_CIG => 'I';
const my $DEL_CIG => 'D';
const my $SOFT_CLIP_CIG => 'S';
const my $HARD_CLIP_CIG => 'H';

my $ref;
my $mut;
my $pos;
my $chr;
my $prms = {};
my $muts = {};
my $is_stranded_format = 1;

sub new {
	my ($proto) = @_;
	my $class = ref($proto) || $proto;
 	my $self = {};
  bless($self, $class);
  return $self;
}

sub init{
  my ($self,$tum,$norm,$cfg) = @_;
  if(!defined($tum) || ! -e $tum){
    croak("Tumour HTS file wasn't defined or an existing file.");
  }
  if(!defined($norm) || ! -e $norm){
    croak("Normal HTS file wasn't defined or an existing file.");
  }
  $self->{'thts'} = Bio::DB::HTS->new(-bam=>$tum);
  $self->{'nhts'} = Bio::DB::HTS->new(-bam=>$norm);
  $self->{'cfg_params'} = _setupCfgParams($cfg);
  $prms = $self->{'cfg_params'};
  return;
}

sub _muts{
  return $muts;
}

sub _prms{
  return $prms;
}

sub _setupCfgParams{
  my ($cfg) = @_;
  my $params;
  foreach my $param_name(@{$cfg->params}){
    $params->{$param_name} = $cfg->param_val($param_name);
  }
  return $params;
}

sub set_position{
  my ($self,$c,$pt,$wt,$mt) = @_;
  $self->{'chr'} = $c;
  $chr = $c;
  $self->{'pos'} = $pt;
  $pos = $pt;
  $self->{'ref'} = $wt;
  $ref = $wt;
  $self->{'mut'} = $mt;
  $mut = $mt;
  $self->_populate_tum_data();
  $self->_populate_norm_data();
}

sub _populate_tum_data{
  my ($self) = @_;
  $muts = undef;
  $self->{'thts'}->fetch($self->{'chr'}.':'.$self->{'pos'}.'-'.$self->{'pos'},\&_tumFetch);
}

sub _isCurrentPosCoveredFromAlignment{
	my ($aln) = @_;
	my $cpos = $aln->start - 1;
	foreach my $cigSect(@{$aln->cigar_array}){

		if($cigSect->[0] eq CIGAR_SYMBOLS->[BAM_CMATCH]){
			if($cpos <= $pos && ($pos+$cigSect->[1]) >= $pos){
				return 1;
			}
			$cpos+= $cigSect->[1];
		}elsif($cigSect->[0] eq CIGAR_SYMBOLS->[BAM_CDEL] || $cigSect->[0] eq CIGAR_SYMBOLS->[BAM_CREF_SKIP]){
			if($cpos <= $pos && ($cpos+$cigSect->[1]) > $pos){
				return 0;
			}
			$cpos+= $cigSect->[1];
		}
	}
	return 0;
}

sub _getReadPositionFromAlignment{
	my ($algn) = @_;
	my $rdPosIndexOfInterest = 0;
  my $currentRefPos = $algn->start -1;
  foreach my $cigSect(@{$algn->cigar_array}){
    if($cigSect->[0] eq CIGAR_SYMBOLS->[BAM_CMATCH]){
      if($currentRefPos <= $pos && ($currentRefPos+$cigSect->[1]) >= $pos){
        for(my $i=0;$i<$cigSect->[1];$i++){
          $rdPosIndexOfInterest++;
          $currentRefPos++;
          if($pos == $currentRefPos){
            return ($rdPosIndexOfInterest,$currentRefPos);
          }
        }
      }else{
        $rdPosIndexOfInterest += $cigSect->[1];
        $currentRefPos += $cigSect->[1];
      }
    }elsif($cigSect->[0] eq CIGAR_SYMBOLS->[BAM_CDEL] || $cigSect->[0] eq CIGAR_SYMBOLS->[BAM_CREF_SKIP]){
      $currentRefPos += $cigSect->[1];
    }elsif($cigSect->[0] eq CIGAR_SYMBOLS->[BAM_CINS] || $cigSect->[0] eq CIGAR_SYMBOLS->[BAM_CSOFT_CLIP]){
      $rdPosIndexOfInterest += $cigSect->[1];
    }
  }
}

sub _get_soft_clip_count_from_cigar{
	my ($cig_arr) = @_;
	my $count = 0;
	foreach my $cigentry(@$cig_arr){
		if($cigentry->[0] eq CIGAR_SYMBOLS->[BAM_CSOFT_CLIP]){
			$count += $cigentry->[1];
		}
	}
	return $count;
}

sub _tumFetch{
  my ($algn) = @_;
  my $flagValue = $algn->flag;
  return if((int($flagValue) & $BAD_BITS) != 0);
  return if((int($flagValue) & $GOOD_BITS) != $GOOD_BITS);
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
		my $start = $algn->start;
		#Read strand
		my $str = 1;
		if($algn->reversed){
			$str = -1;
		}
		$muts->{'totalTCoverage'} += 1;
		if($str == 1){
			$muts->{'totalTCoveragePos'} += 1;
		}else{
			$muts->{'totalTCoverageNeg'} += 1;
		}
		if($algn->cigar_str =~ m/[ID]/){
			$muts->{'indelTCount'} += 1;
		}

		#Read base
		my $qbase = $splt[$rdPosIndexOfInterest-1];

		#Base quality
		my $qscore = $algn->qscore->[$rdPosIndexOfInterest-1];

    my $gapDist = _getDistanceFromGapInRead($algn->cigar_array,$rdPosIndexOfInterest);

		push(@{$muts->{'completeMutStrands'}},$str);

		push(@{$muts->{'allTumBases'}},$qbase);

		push(@{$muts->{'allTumBaseQuals'}},$qscore);

		push(@{$muts->{'allTumStrands'}},$str);

    push(@{$muts->{'allTumMapQuals'}},$algn->qual);

    push(@{$muts->{'allMinGapDistances'}},$gapDist);

		#return if(uc($qbase) ne uc($mutBase));
    my $xt = $algn->aux_get('XT');
		return if (get_param('keepSW') == 0 && defined($xt) && $xt eq 'M');

		return if($qscore < get_param('minAnalysedQual'));

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

		return if(uc($qbase) ne uc($mut));

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
		push(@{$muts->{'alnp'}},$primaryalnscore) if(defined($primaryalnscore));

		#Softclipping
		push(@{$muts->{'sclp'}},$softclipcount);
	}
	return 1;
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

sub _matchedNormFetch{
  my ($algn) = @_;
  my $flagValue = $algn->flag;
  return if((int($flagValue) & $BAD_BITS) != 0);
  return if((int($flagValue) & $GOOD_BITS) != $GOOD_BITS);
  if(_isCurrentPosCoveredFromAlignment($algn) == 1){
		#Get the correct read position.
		my ($rdPosIndexOfInterest,$currentRefPos) = _getReadPositionFromAlignment($algn,$pos);
		#print $rdPosIndexOfInterest,"\n";
		#print $algn->qseq,"\n";
		my @splt = split(//,$algn->qseq);
  		#Calculate other stuff
		my $totalPCovg = 0;
		my $totalNCovg = 0;
		my $indelRdCount = 0;
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

		return if (get_param('keepSW') == 0 && defined $xt && $xt eq 'M');

		return if($qscore < get_param('minAnalysedQual'));

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
		return if(uc($qbase) ne uc($mut));
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

sub _populate_norm_data{
  my ($self) = @_;
  $self->{'nhts'}->fetch($self->{'chr'}.':'.$self->{'pos'}.'-'.$self->{'pos'},\&_matchedNormFetch);
}

sub get_param{
  my ($pname) = @_;
  croak("Parameter '$pname' wasn't found. Is it in the config?") unless (exists $prms->{$pname} && defined $prms->{$pname});
  return $prms->{$pname};
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

sub _getTotalReadsOnStrandCount{
	my ($self,$strand) = @_;
	my $strandCount = 0;
	foreach my $str(@{$muts->{'completeMutStrands'}}){
		if($strand == $str){
			$strandCount++;
		}
	}
	return $strandCount;
}

sub median{
	return (sum( ( sort { $a <=> $b } @_ )[ int( $#_/2 ), ceil( $#_/2 ) ] )/2);
}

###### Below are the default flags provided by this package

sub depthFlag{
	my ($self) = @_;
	my $depth = scalar(@{$muts->{'tqs'}});
	my $overCutoff = 0;
	foreach my $q(@{$muts->{'tqs'}}){
		if($q >= get_param('minDepthQual')){
			$overCutoff++;
		}
		if($overCutoff >= ($depth * get_param('depthCutoffProportion'))){
			return 0; #Pass
		}
	}
	return 1; #Fail
}

sub cavemanMatchNormalProportionFlag{
    my ($self, $vcf, $x) = @_;
    my $normal_col = $vcf->get_column($x,$VCF_COLUMN_NORMAL);
    my $tumour_col = $vcf->get_column($x,$VCF_COLUMN_TUMOUR);
    my $format = $vcf->get_column($x,$VCF_COLUMN_FORMAT);
    return $self->_getCavemanMatchNormalProportionFlag($normal_col,$tumour_col,$format);
}

sub _getCavemanMatchNormalProportionFlag{
    my ($self,$normal_col,$tumour_col,$format) = @_;
    my @splitnorm = split(/:/,$normal_col);
    my @splittum = split(/:/,$tumour_col);
    my @splitformat = split(/:/,$format);
    $is_stranded_format = 0 if($format =~ m/$OLD_ALLELE_VCF_FORMAT/);
    my $total_norm_cvg = 0;
    my $mut_allele_cvg = 0;
    my $total_tumm_cvg = 0;
    my $mut_allele_tum_cvg = 0;
    my %decode_hash = %OLD_ALLELE_VCF_FORMAT_INDEX_HASH;
    if($is_stranded_format==1){
      %decode_hash = %NEW_ALLELE_VCF_FORMAT_INDEX_HASH;
      $total_norm_cvg = sum(@splitnorm[1..8]);
      $total_tumm_cvg = sum(@splittum[1..8]);
    }else{
      $total_norm_cvg = sum(@splitnorm[1..4]);
      $total_tumm_cvg = sum(@splittum[1..4]);
    }
    my $mutbase = $self->{'mut'};
    $mut_allele_cvg = sum(@splitnorm[$decode_hash{$mutbase}]);
    $mut_allele_tum_cvg = sum(@splittum[$decode_hash{$mutbase}]);
    my $norm_prop = $mut_allele_cvg/$total_norm_cvg;
    my $tum_prop = $mut_allele_tum_cvg/$total_tumm_cvg;
    return 1 if($norm_prop > 0 && (($tum_prop - $norm_prop) < get_param('maxCavemanMatchedNormalProportion'))); #Fail
    return 0; #Pass
}

sub readPositionFlag{
  my ($self) = @_;
  my $sec3rd = 0;
	my $first3rd = 0;
	if(scalar(@{$muts->{'trp'}}) > get_param('minRdPosDepth')){
		return 0; #Pass
	}
	for(my $i=0;$i<scalar(@{$muts->{'trp'}});$i++){
		my $rdLn = $muts->{'trl'}->[$i];
		my $thirds = $rdLn / 3;
		my $rdPos = $muts->{'trp'}->[$i];
		my $prop = ($rdLn * get_param('readPosBeginningOfReadIgnoreProportion'));
		my $halfprop = ($rdLn * get_param('readPosTwoThirdsOfReadExtendProportion'));
		#In first or second third, but not first n%, extending by m%
		if($rdPos <= (($thirds * 2)+$halfprop) && $rdPos > $prop){
			return 0; #Pass
		}
	}
	return 1; #Fail
}

sub matchedNormalFlag{
  my ($self) = @_;
  my $qualCnt = 0;
	my $call = "";
	foreach my $q(@{$muts->{'nqs'}}){
		if($q >= (get_param('minNormMutAllelequal'))){
			$qualCnt++;
			my $proportion = $qualCnt / $muts->{'totalNCoverage'};
      if($qualCnt > (get_param('matchedNormalAlleleHiCvgCutoff'))){
        if($proportion > (get_param('maxMatchedNormalAlleleHiCvgProportion'))){
          return 1;
        }
      }else{
        if($proportion > (get_param('maxMatchedNormalAlleleProportion'))){
          return 1;
        }
      }
		}
	}
	return 0;
}

sub pentamericMotifFlag{
  my ($self) = @_;
  my $minus = 0;
	my $plus = 0;
	my $lastThird = 0;
	#Check strands first.
	foreach my $str(@{$muts->{'tstr'}}){
		if($str == 1){
			$plus++;
		}else{
			$minus++;
		}
		return 0 if($plus > 1 && $minus > 1);
	}

	#If all or (all - 1) mut allele reads are on one strand continue.
	my $sz = scalar(@{$muts->{'tstr'}});
	if(!(($minus >= ($sz -1) && $plus <= 1) || ($plus >= ($sz -1) && $minus <= 1)) ){
		#Passes on strand rule
		return 0;
	}

	my $avgBQOverall = 0;
	my $cntAvg = 0;
	my @tqs = @{$muts->{'tqs'}};
	my @trp = @{$muts->{'trp'}};
	my @trl = @{$muts->{'trl'}};
	my @tstr = @{$muts->{'tstr'}};
	my @trn = @{$muts->{'trn'}};
	#We've checked 3rds of read and the strand,
	#if we got this far we have to check the entire mutant read.
	#Read names are handily stored in $muts->{'trn'} so we can look at each in turn.
	my $avgBQ = 0;
	my $analysed = 0;



  #Iterate through each read name
	for(my $i=0;$i<$sz;$i++){
		my $thispos = $muts->{'trp'}->[$i];
		my $length = $muts->{'trl'}->[$i];
		my $rdName = $muts->{'trn'}->[$i];
		my $start = $muts->{'trdst'}->[$i];
		#If position is not in last 3rd we can return now.
		#print "$rdName\n";
		if($thispos < ($length / 2)){
		  #Passed on last third of read $rdName $muts->{'tstr'}->[$i]\n";
			return 0;
		}
		my $rd;
		#Fetch the read we want
		$self->{'thts'}->fetch($chr.":".$pos."-".$pos,
						sub{
							my $a = shift;
							if($a->qname eq $rdName){# && $start == $a->start){
								$rd = $a;
								return;
							}
						});
		my $isReversed = $rd->reversed;
		my $readPosOfMut = ($pos - ($rd->pos + 1)) + 1;
		my $seq = $rd->qseq;
		my $quals = $rd->qscore();
		my $posMatch = 0;
		my $lastPos = 0;
		my @matches = ();
		if($isReversed == 1){
			$seq =~ tr/acgtnACGTN/tgcanTGCAN/;
			$seq = reverse($seq);
			#Reverse the quality array
			@$quals = reverse(@$quals);
			$readPosOfMut = (length($seq) - $readPosOfMut) + 1;
		}
		#croak("Error RDPosFirst $pos != RDPosNew $readPosOfMut : \t $rdName\n") if($pos != $readPosOfMut);
		#print "mut rd pos $readPosOfMut\n";
		#Check for motif match in second half of the read (rev comped as required).
		#No motif, so skip
		next if($seq !~ m/GGC[AT]G/);
			#print "No motif matches\n";
		while($seq =~ m/GGC[AT]G/g){
			#$posMatch = length($`);
			push(@matches,length($`).",".length($&)."");
		}
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
  return 0 if($cntAvg < 1); #Pass
	if(($avgBQOverall / $cntAvg) < get_param('pentamerMinPassAvgQual')){
		return 1; #Fail
	}
	return 0; #Pass

}

sub avgMapQualFlag{
  my ($self) = @_;
  my $total = 0;
	my $noOfMQs = 0;
	foreach my $mq(@{$muts->{'tmq'}}){
		$total += $mq;
		$noOfMQs++;
	}
	if($noOfMQs == 0){
		#warn ("checkAvgMapQual: Should not have encountered 0 mutant allele reads in the tumour for ",$self->_chromosome,":",$self->_currentPos,", returning pass.\n");
		return 1; #Fail
	}
	#If the mean is less than the min mq for a pass we fail.
	return 1 if(($total / $noOfMQs) < get_param('minPassAvgMapQual')); #Fail
	return 0; #Pass
}

sub phasingFlag{
  my ($self) = @_;
	#Count mut bases on each strand.
	my ($f_q_sum, $r_q_sum, $f_count, $r_count) = (0,0,0,0);
	for(my $i=0;$i<scalar(@{$muts->{'allTumStrands'}});$i++){
		if($muts->{'allTumStrands'}->[$i] == 1){
			$f_count++;
	    	$f_q_sum += $muts->{'allTumBaseQuals'}->[$i];
	  	}else {
	   	$r_count++;
	    	$r_q_sum += $muts->{'allTumBaseQuals'}->[$i];
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
	if($f_count > 0 && ($r_count == 0 || ($r_count / ($self->_getTotalReadsOnStrandCount(-1))) <= get_param('maxPhasingMinorityStrandReadProportion'))){
		if($f_av_qual < get_param('minPassPhaseQual')){
			return 1;
		}
	}elsif($r_count > 0 && ($f_count == 0 || ($f_count / $self->_getTotalReadsOnStrandCount(1)) <= get_param('maxPhasingMinorityStrandReadProportion'))){#Only on rev strand
		if($r_av_qual < get_param('minPassPhaseQual')){
			return 1;
		}
	}
	return 0;
}

sub getBEDUnmatchedFlag{
  my ($self,$line) = @_;
  my ($chr,$start,$stop,$score,$mutlist) = split(/\t/,$line);
  return 1 if($mutlist =~ m/$mut/i);
  return 0;
}

sub getVCFUnmatchedFlag{
  my ($self,$line) = @_;
  my ($ch,$po,$ident,$refAll,$altAll,$quality,$filts,$info,$format,@samples) = split(/\t/,$line);
  #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	PD4106b	PD4107b	PD4108b	PD4109b	PD4110b	PD4111b	PD4112b	PD4113b	PD4114b	PD4115b
  my ($geno,@formats) = split(':',$format);
  my $sampleHitCount = 0;
  my $totalSampleCnt = 0;
  foreach my $sampData(@samples){
    $totalSampleCnt++;
    next if($sampData eq q{0:.:.:.:.} || $sampData eq q{-} || $sampData eq q{.});
    #GT:GF:CF:TF:AF	0|0:41:0:0:0
    my ($gentype,@data) = split(':',$sampData);
    my $totalCvg = 0;
    my $mutAlleleCvg = 0;
    my $proportion = 0;
    for (my $i=0;$i<scalar(@data);$i++){
      next unless ($formats[$i] =~ m/^[ACGT]{1}[A-Z]$/);
      $totalCvg += $data[$i];
      if(substr($formats[$i],0,1) eq $mut){
        $mutAlleleCvg = $data[$i];
      }
    }
    if($mutAlleleCvg >= get_param('vcfUnmatchedMinMutAlleleCvg')){
      $sampleHitCount++;
    }
    if((($sampleHitCount/scalar($totalSampleCnt))*100) >= get_param('vcfUnmatchedMinSamplePct')){
      return 1;
    }
  }
  return 0 if($totalSampleCnt == 0);
  return 0;
}

sub singleEndFlag{
  my ($self) = @_;
  return 0 if($muts->{'pcvg'} < get_param('minSingleEndCoverage') || $muts->{'ncvg'} < get_param('minSingleEndCoverage'));
	my $hasPos = 0;
	my $hasNeg = 0;
	foreach my $str(@{$muts->{'tstr'}}){
		if($str == -1){
			$hasNeg++;
		}elsif($str == 1){
			$hasPos++;
		}
		return 0 if($hasNeg > 0 && $hasPos > 0);
	}
	if($hasNeg == 0 || $hasPos == 0){
		return 1;
	}
	return 0;
}

sub matchedNormalProportion{
  my ($self) = @_;
	#Calculate tumour proportion of mut allele
	my $tumProp = 0;
	if(scalar(@{$muts->{'tqs'}}) > 0){
		$tumProp = scalar(@{$muts->{'tqs'}})/$muts->{'tumcvg'};
	}
	#Calculate normal proportion of mut allele
	my $normProp = 0;
	if(exists($muts->{'nqs'}) && scalar(@{$muts->{'nqs'}}) > 0){
		$normProp = scalar(@{$muts->{'nqs'}})/$muts->{'normcvg'};
	}
	#Fail if the difference is less than the given proportion/percentage
	return 1 if($normProp > 0 && ($tumProp - $normProp) < get_param('matchedNormalMaxMutProportion'));
	return 0;
}

sub alignmentScoreReadLengthAdjustedFlag{
  my ($self) = @_;
  my $adj;
  return sprintf('%.2f',0) if(!exists($muts->{'alnp'}) || !defined($muts->{'alnp'}) || scalar(@{$muts->{'alnp'}}) == 0);
  for(my $i=0; $i<scalar(@{$muts->{'alnp'}}); $i++){
    $adj->[$i] = $muts->{'alnp'}->[$i] / $muts->{'trl'}->[$i];
  }
  return sprintf('%.2f',median(@{$adj}));
}

sub clippingMedianFlag{
  my ($self) = @_;
  return sprintf('%.2f',0) if(!defined($muts->{'sclp'}) || scalar(@{$muts->{'sclp'}}) == 0);
	return sprintf('%.2f',median(@{$muts->{'sclp'}}));
}

sub alnScoreMedianFlag{
  my ($self) = @_;
  return sprintf('%.2f',0) if(!defined($muts->{'alnp'}) || scalar(@{$muts->{'alnp'}}) == 0);
  return sprintf('%.2f',median(@{$self->_muts->{'alnp'}}));;
}

sub tumIndelDepthFlag{
	my ($self) = @_;
	if(!exists $self->_muts->{'totalTCoverage'} || $self->_muts->{'totalTCoverage'} == 0){
		return 0;
	}
	my $prop = ($self->_muts->{'indelTCount'} / $self->_muts->{'totalTCoverage'}) * 100;
	if($prop <= get_param('maxTumIndelProportion')){
		return 0;
	}
	return 1;

}

sub sameReadPosFlag{
  my ($self) = @_;
	my $tmp;
	my $noOfReads = scalar(@{$self->_muts->{'trp'}});
	my $permittedNo = $noOfReads  * ( get_param('samePosMaxPercent') / 100 );
	for(my $i=0;$i<scalar(@{$self->_muts->{'trp'}});$i++){
		my $rdPos = $self->_muts->{'trp'}->[$i];
		if(!defined($tmp->{$rdPos})){
			$tmp->{$rdPos} = 0;
		}
		$tmp->{$rdPos}++;
		if($tmp->{$rdPos} > $permittedNo){
			return 1;
		}
	}
	return 0;
}

sub withinGapRangeFlag{
  my ($self) = @_;
  my $meanMapQ = sum($self->_muts->{'allTumMapQuals'})/scalar(@{$self->_muts->{'allTumMapQuals'}});
  return 0 if($meanMapQ < get_param('minMeanMapQualGapFlag')); #Pass as likely mismapping
  my $total_reads = scalar(@{$self->_muts->{'allTumMapQuals'}});
  my @non_tum_base_dist = [];
  my $norm_base_dist_count = 0;
  foreach (zip($self->_muts->{'allTumBases'},$self->_muts->{'allMinGapDistances'})){
    my ($base, $distance) = @{$_};
    if($base eq $self->{'mut'}){
      return 0 if($distance != -1);
    }else{
      if($distance != -1 && $distance <= get_param('withinXBpOfDeletion')){
        push(@non_tum_base_dist, $distance);
        $norm_base_dist_count++;
      }
    }
  }
  return 0 if($norm_base_dist_count==0); #Pass if zero reference reads with gap
  my $percentage_reads_present = ($norm_base_dist_count/$total_reads) * 100;
  return 1 if($percentage_reads_present >= get_param('minGapPresentInReads'));
  return 0;
}

1;
