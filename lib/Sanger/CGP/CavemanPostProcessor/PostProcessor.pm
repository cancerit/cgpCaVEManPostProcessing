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

package Sanger::CGP::CavemanPostProcessor::PostProcessor;

use strict;
use Bio::DB::HTS;
use Bio::DB::HTS::Constants;
use Bio::DB::HTS::Alignment;
use POSIX qw(strftime ceil);
use List::Util qw (sum zip min max);
use Carp;
use Const::Fast qw(const);

use Data::Dumper;

use Sanger::CGP::CavemanPostProcessor;
use Sanger::CGP::CavemanPostProcessor::Constants;

use parent qw(Sanger::CGP::CavemanPostProcessor);

my $const = 'Sanger::CGP::CavemanPostProcessor::Constants';
my $muts;
my $norms;
my $muts_rds;
my $norms_rds;
my $refBase;
my $mutBase;
my $chromosome;
my $tum_readnames;
my $norm_readnames;
my $tum_readnames_arr;
my $norm_readnames_arr;
my $tum_readnames_hash;
my $norm_readnames_hash;

#Defaults for this post processing module

my $is_stranded_format = 1;

sub new{
  my ($proto) = shift;
  my $class = ref($proto) || $proto;
  my %inputs = @_;
  my $self = $class->SUPER::new(%inputs);
  bless $self, $class;
  return $self;
}

sub _init{
  my ($self,$inputs) = @_;
  
  return $self;
}

sub _init{
  my ($self,$inputs) = @_;
  if(!defined($inputs->{'tumBam'}) || !defined($inputs->{'normBam'})){
    croak("tumBam and normBam are required for initialisation.\n");
  }
  $self->tumBam($inputs->{'tumBam'}, $inputs->{'ref'});
  $self->normBam($inputs->{'normBam'}, $inputs->{'ref'});
  $self->keepSW($inputs->{'keepSW'}) if exists $inputs->{'keepSW'};
  $self->minAnalysedQual($inputs->{'minAnalysedQual'}) if exists $inputs->{'minAnalysedQual'};
  $self->minDepthQual($inputs->{'minDepthQual'});
  $self->depthCutoffProportion($inputs->{'depthCutoffProportion'});
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
  $self->matchedNormalMaxMutProportion($inputs->{'matchedNormalMaxMutProportion'});
  $self->maxCavemanMatchedNormalProportion($inputs->{'maxCavemanMatchedNormalProportion'});
  $self->withinXBpOfDeletion($inputs->{'withinXBpOfDeletion'});
  $self->minGapPresentInPercentReads($inputs->{'minGapPresentInReads'});
  $self->meanMapQualGapFlag($inputs->{'minMeanMapQualGapFlag'});
  $self->maxGapFlagDistFromEndOfReadProp($inputs->{'maxGapFlagDistFromEndOfReadProp'});
  $self->minGapFlagDistEndOfReadPercent($inputs->{'minGapFlagDistEndOfReadPercent'});
  $self->minSingleEndCoverage($inputs->{'minSingleEndCoverage'});

  return $self;
}

sub runProcess{
  my ($self,$chr,$start,$stop,$refBase,$mBase) = @_;
  $muts = undef;
  $norms = undef;
  $muts_rds = {};
  $norms_rds = {};
  $self->clearResults();
  $self->_chromosome($chr);
  $self->_currentPos($start);
  $self->_refBase($refBase);
  $self->_mutBase($mBase);
  $tum_readnames = undef;
  $tum_readnames_arr = undef;
  $tum_readnames_hash = undef;
  my $hts = $self->{'tb'};
  $hts->hts_index->fetch(
    $hts->hts_file,
    $hts->header->parse_region($chr.':'.$start.'-'.$stop),
    sub{$self->_callbackTumFetch(@_);},
    $hts
  );
  $self->process_hashed_reads(sub{$self->populate_muts(@_);}, $tum_readnames, $tum_readnames_arr);
  $norm_readnames = undef;
  $norm_readnames_arr = undef;
  $norm_readnames_hash = undef;
  $hts = $self->{'nb'};
  $hts->hts_index->fetch(
    $hts->hts_file,
    $hts->header->parse_region($chr.':'.$start.'-'.$stop),
    sub{$self->_callbackMatchedNormFetch(@_);},
    $hts
  );
  $self->process_hashed_reads(sub{$self->populate_norms(@_);}, $norm_readnames, $norm_readnames_arr);
  return 1;
}

=item keepSW
    Whether to include Smith-Waterman aligned reads in post processing
=cut
sub keepSW{
  my ($self,$keep) = @_;
  if(defined($keep) && ($keep == 1 || $keep == 0)){
    $self->{'keepSW'} = $keep;
  }
  if(! defined $self->{'keepSW'}){
    $self->{'keepSW'} = $const->default_flag_values('KEEPSW');
  }
  return $self->{'keepSW'};
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
      $self->{'sec'} = $const->default_flag_values('MIN_SINGLE_END_CVG');
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
      $self->{'mnmmp'} = $const->default_flag_values('MATCHED_NORMAL_MAX_MUT_PROP');
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
    $self->{'minaq'} = $q;
  }else{
    if(!defined($self->{'minaq'})){
      $self->{'minaq'} = $const->default_flag_values('MIN_ANALYSED_QUAL');
    }
  }
  return $self->{'minaq'};
}

sub _currentPos{
  my ($self,$pos) = @_;
  if(defined($pos)){
    $self->{'currPos'} = $pos;
  }
  return $self->{'currPos'};
}

sub _chromosome{
  my ($self,$c) = @_;
  if(defined($c)){
    $chromosome = $c;
  }
  return $chromosome;
}

sub _refBase{
  my ($self,$b) = @_;
  if(defined($b)){
    $refBase = $b;
  }
  return $refBase;
}

sub refBase{
  return $refBase;
}

sub mutBase{
  return $mutBase;
}

sub _mutBase{
  my ($self,$b) = @_;
  if(defined($b)){
    $mutBase = $b;
  }
  return $mutBase;
}

sub _tum_readnames_hash{
  my ($self,$h) = @_;
  if(defined($h)){
    $tum_readnames_hash = $h;
  }
  return $tum_readnames_hash;
}

sub _tum_readnames_arr{
  my ($self,$a) = @_;
  if(defined($a)){
    $tum_readnames_arr = $a;
  }
  return $tum_readnames_arr;
}

sub _tum_readnames{
  my ($self,$r) = @_;
  if(defined($r)){
    $tum_readnames = $r
  }
  return $tum_readnames;
}

sub _norm_readnames_hash{
  my ($self,$h) = @_;
  if(defined($h)){
    $norm_readnames_hash = $h;
  }
  return $norm_readnames_hash;
}

sub _norm_readnames_arr{
  my ($self,$a) = @_;
  if(defined($a)){
    $norm_readnames_arr = $a;
  }
  return $norm_readnames_arr;
}

sub _norm_readnames{
  my ($self,$r) = @_;
  if(defined($r)){
    $norm_readnames = $r
  }
  return $norm_readnames;
}

sub tumBam{
  my ($self,$bam,$fasta) = @_;
  if(defined($bam)){
    $self->{'tb'} = Bio::DB::HTS->new(-bam=>$bam, -fasta=>$fasta);
  }
  return $self->{'tb'};
}

sub normBam{
  my ($self,$bam,$fasta) = @_;
  if(defined($bam)){
    $self->{'nb'} = Bio::DB::HTS->new(-bam=>$bam, -fasta=>$fasta);
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

sub process_hashed_reads{
  my ($self, $code, $hashed_reads, $readname_arr) = @_;

  my %loc_counts = (
    1  => {A => 0, C => 0, G => 0, T => 0},
    -1 => {A => 0, C => 0, G => 0, T => 0},
  );

  foreach my $readnom(@$readname_arr){
    my $read = $hashed_reads->{$readnom};
    my ($read_to_use, $read_to_use_2);
    if(exists $read->{1} && exists $read->{-1}) {
      # if change is different have to include both versions
      if( $read->{1}->{qbase} ne $read->{-1}->{qbase}) {
        $read_to_use = $read->{1};
        $read_to_use_2 = $read->{-1};
        $loc_counts{1}{$read->{1}->{qbase}}++;
        $loc_counts{-1}{$read->{-1}->{qbase}}++;
      }
      # score is same check loc counts for that allele to determine where to put
      elsif($read->{1}->{qscore} == $read->{-1}->{qscore}) {
        if($loc_counts{1}{$read->{1}->{qbase}} <= $loc_counts{-1}{$read->{-1}->{qbase}}) {
          $read_to_use = $read->{1};
          $loc_counts{1}{$read->{1}->{qbase}}++;
        }
        else {
          $read_to_use = $read->{-1};
          $loc_counts{-1}{$read->{-1}->{qbase}}++;
        }
      }
      elsif($read->{1}->{qscore} > $read->{-1}->{qscore}) {
        $read_to_use = $read->{1};
        $loc_counts{1}{$read->{1}->{qbase}}++;
      }
      else{
        $read_to_use = $read->{-1};
        $loc_counts{-1}{$read->{-1}->{qbase}}++;
      }
    }
    else {
      if(exists $read->{1}) {
        $read_to_use = $read->{1};
        $loc_counts{1}{$read->{1}->{qbase}}++;
      }
      else {
        $read_to_use = $read->{-1};
        $loc_counts{-1}{$read->{-1}->{qbase}}++;
      }
    }
    &$code($read_to_use);
    if(defined $read_to_use_2){
      &$code($read_to_use_2);
    }
    delete $hashed_reads->{$readnom};
  }
}

sub _callbackTumFetch{
  my ($self, $a, $hts) = @_;
  my $flagValue = $a->flag;
  #Check read and mate are mapped. If not return.
  return if((int($flagValue) & 2) != 2); # Proper pair
  return if((int($flagValue) & 3852) != 0);
  # Ensure that we keep
  return if((int($flagValue) & 16) != 0 && (int($flagValue) & 32) != 0);
  return if((int($flagValue) & 16) == 0 && (int($flagValue) & 32) == 0);

  my $algn = Bio::DB::HTS::AlignWrapper->new($a,$hts);
  my $pos = $a->pos;
  my $cigar_array = $algn->cigar_array; # expensive and reused so save to variable
  #Quick check that were covering the base with this read (skips/indels are ignored)
  my $is_covered = _isCurrentPosCoveredFromAlignment($pos, $cigar_array, $self->_currentPos()); #1 is covered, -1 is covered but within indel
  if($is_covered != 0){
    my $this_read;

    my $rdname = $a->qname;
    #Read strand, faster than using $a->strand
    my $str = 1;
    if($algn->reversed){
      $str = -1;
    }
    my $cig_str = $algn->cigar_str;
    #Read base
    $this_read->{str} = $str;


    $this_read->{gapDist} = 0;
    $this_read->{ln} = $a->l_qseq;
    $this_read->{distFromEnd}=-1;
    if($is_covered == 1){
      #Get the correct read position.
      my ($rdPosIndexOfInterest,$currentRefPos) = $self->_getReadPositionFromAlignment($pos, $cigar_array);
      $this_read->{qbase} = substr $a->qseq, $rdPosIndexOfInterest-1, length($self->_mutBase());
      $this_read->{qscore} = unpack('C*', substr($a->_qscore, $rdPosIndexOfInterest-1, length($self->_mutBase())));
      $this_read->{rdPos} = $rdPosIndexOfInterest;
      $this_read->{gapDist} = _getDistanceFromGapInRead($algn->cigar_array,$rdPosIndexOfInterest);
      $this_read->{distFromEnd} = min(($rdPosIndexOfInterest/$this_read->{ln}),(($this_read->{ln}-$rdPosIndexOfInterest)/$this_read->{ln}));
    }
    $this_read->{matchesindel} = ($cig_str =~ m/[ID]/);
    $this_read->{xt} = $a->aux_get('XT');
    $this_read->{softclipcount} = 0;
    my $type = $const->cigar_types('SOFT_CLIP_CIG');
    if ($cig_str =~ m/$type/){
      $this_read->{softclipcount} = $self->_get_soft_clip_count_from_cigar($algn->cigar_array);
    }
    $this_read->{primaryalnscore} = $a->aux_get('AS');# $algn->get_tag_values('AS');
    $this_read->{qual} = $a->qual;
    $this_read->{start} = $algn->start;
    $this_read->{rdName} = $rdname;
  
    my $this_tum_read_nams_hash = $self->_tum_readnames_hash();
    my $this_tum_readnames_arr = $self->_tum_readnames_arr();
    if(!exists $this_tum_read_nams_hash->{$rdname}){
      push(@$this_tum_readnames_arr, $rdname);
      $this_tum_read_nams_hash->{$rdname} = 0;
    }
    my $this_tum_readnames = $self->_tum_readnames();
    $this_tum_readnames->{$rdname}->{$str} = $this_read;
    $self->_tum_readnames_hash($this_tum_read_nams_hash);
    $self->_tum_readnames_arr($this_tum_readnames_arr);
    $self->_tum_readnames($this_tum_readnames);

  } # End of if this is a covered position, look at deletion event at this location (required for deletion flag)
  return 1;
}

sub populate_muts{
    my ($self, $read) = @_;
    $muts->{'totalTCoverage'} += 1;
    if($read->{str} == 1){
        $muts->{'totalTCoveragePos'} += 1;
    }else{
        $muts->{'totalTCoverageNeg'} += 1;
    }
    if($read->{xt}){
        $muts->{'indelTCount'} += 1;
    }
    push(@{$muts->{'completeMutStrands'}},$read->{str});
    push(@{$muts->{'allTumBases'}},$read->{qbase});
    push(@{$muts->{'allTumBaseQuals'}},$read->{qscore});
    push(@{$muts->{'allTumStrands'}},$read->{str});
    push(@{$muts->{'allTumMapQuals'}},$read->{qual});
    push(@{$muts->{'allMinGapDistances'}},$read->{gapDist});
    push(@{$muts->{'allMutDistPropFromEndOfRead'}},$read->{distFromEnd});
    return if ($self->keepSW() == 0 && defined($read->{xt}) && $read->{xt} eq 'M');

    return if($read->{qscore} < $self->minAnalysedQual());

    $muts->{'tumcvg'} += 1;

    if($read->{str} == 1){
    $muts->{'pcvg'} += 1;
  }else{
    $muts->{'ncvg'} += 1;
    $read->{rdPos} = ($read->{ln} - $read->{rdPos}) + 1;
  }

  return if(uc($read->{qbase}) ne uc($mutBase));

    #Tum quals
    push(@{$muts->{'tqs'}},$read->{qscore});

    #Tum Rd Pos
    push(@{$muts->{'trp'}},$read->{rdPos});

    #Tum rd length
    push(@{$muts->{'trl'}},$read->{ln});

    #Tum XT tags
    push(@{$muts->{'txt'}},$read->{xt});

    #Tum rd start
    push(@{$muts->{'trdst'}},$read->{start});

    #Strands
    push(@{$muts->{'tstr'}},$read->{str});

    #RdNames
    push(@{$muts->{'trn'}},$read->{rdName});

    #Mapping quals
    push(@{$muts->{'tmq'}},$read->{qual});

    #AlnScoresPrm
    push(@{$muts->{'alnp'}},$read->{primaryalnscore});

    #Softclipping
    push(@{$muts->{'sclp'}},$read->{softclipcount});

    return;
}

sub _get_soft_clip_count_from_cigar{
  my ($self, $cig_arr) = @_;
  my $count = 0;
  foreach my $cigentry(@$cig_arr){
    if($cigentry->[0] eq $const->cigar_types('SOFT_CLIP_CIG')){
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
    if($cigSect->[0] eq $const->cigar_types('MATCH_CIG') || $cigSect->[0] eq $const->cigar_types('SKIP_CIG') ||
          $cigSect->[0] eq $const->cigar_types('SOFT_CLIP_CIG')){
      $currentRp+=$cigSect->[1];
    }elsif($cigSect->[0] eq $const->cigar_types('DEL_CIG') || $cigSect->[0] eq $const->cigar_types('INS_CIG')){
      my $dp_start = $currentRp+1;
      my $dp_end = $currentRp+$cigSect->[1];
      my $tmp_dist = 0;
      if($rdPosIndexOfInterest>$dp_end || $rdPosIndexOfInterest<$dp_start){
        $tmp_dist = min(abs($rdPosIndexOfInterest-$dp_start),abs($dp_end-$rdPosIndexOfInterest));
      }
      if($tmp_dist < $min_gap_dist || $min_gap_dist == -1){
        $min_gap_dist = $tmp_dist;
      }
    }
  }
  return $min_gap_dist;
}

sub _getReadPositionFromAlignment{
  my ($self, $currentRefPos, $cigar_array) = @_; # 0-based pos ($a->pos)
  my $rdPosIndexOfInterest = 0;
  foreach my $cigSect(@{$cigar_array}){
    if($cigSect->[0] eq $const->cigar_types('MATCH_CIG')){
      my $op_len = $cigSect->[1];
      if($currentRefPos <= $self->_currentPos() && ($currentRefPos+$op_len) >= $self->_currentPos()){
        for(0..($op_len - 1)) {
          $rdPosIndexOfInterest++;
          $currentRefPos++;
          if($self->_currentPos() == $currentRefPos){
            return ($rdPosIndexOfInterest,$currentRefPos);
          }
        }
      }else{
        $rdPosIndexOfInterest += $op_len;
        $currentRefPos += $op_len;
      }
    }elsif($cigSect->[0] eq $const->cigar_types('DEL_CIG') || $cigSect->[0] eq $const->cigar_types('SKIP_CIG')){
      $currentRefPos += $cigSect->[1];
    }elsif($cigSect->[0] eq $const->cigar_types('INS_CIG') || $cigSect->[0] eq $const->cigar_types('SOFT_CLIP_CIG')){
      $rdPosIndexOfInterest += $cigSect->[1];
    }
  }
}

sub _isCurrentPosCoveredFromAlignment{
  my ($pos, $cigar_array, $pos_of_interest) = @_; # 0-based pos, 1 based current pos
  foreach my $cigSect(@{$cigar_array}){

    if($cigSect->[0] eq $const->cigar_types('MATCH_CIG')){
      if($pos_of_interest >= ($pos + 1) && $pos_of_interest <= ($pos+$cigSect->[1])){
        return 1;
      }
      $pos+= $cigSect->[1];
    }elsif($cigSect->[0] eq $const->cigar_types('DEL_CIG') || $cigSect->[0] eq $const->cigar_types('SKIP_CIG')){
      if($pos_of_interest >= ($pos + 1) && $pos_of_interest <= ($pos+$cigSect->[1])){
        return 0 if ($cigSect->[0] eq $const->cigar_types('SKIP_CIG'));
        return -1 if ($cigSect->[0] eq $const->cigar_types('DEL_CIG'));
      }
      $pos+= $cigSect->[1];
    }
  }
  return 0;
}

sub populate_norms{
    my ($self, $read) = @_;
    if(!defined($muts->{'totalNCoverage'})){
        $muts->{'totalNCoverage'} = 0;
    }
    $muts->{'totalNCoverage'} += 1;

    if(!defined($muts->{'allNormBases'})){
        $muts->{'allNormBases'} = [];
    }
    push(@{$muts->{'allNormBases'}},$read->{qbase});

    if(!defined($muts->{'allNormBaseQuals'})){
        $muts->{'allNormBaseQuals'} = [];
    }
    push(@{$muts->{'allNormBaseQuals'}},$read->{qscore});

    if(!defined($muts->{'allNormStrands'})){
        $muts->{'allNormStrands'} = [];
    }
    push(@{$muts->{'allNormStrands'}},$read->{str});

    return if ($self->keepSW == 0 && defined($read->{xt}) && $read->{xt} eq 'M');

  return if($read->{qscore} < $self->minAnalysedQual());

    if(!defined($muts->{'normcvg'})){
        $muts->{'normcvg'} = 0;
    }
    $muts->{'normcvg'} += 1;

    if($read->{str} == +1){
        $muts->{'npcvg'} += 1;
    }else{
        $muts->{'nncvg'} += 1;
        $read->{rdPos} = ($read->{ln} - $read->{rdPos}) + 1;
    }

  return if(uc($read->{qbase}) ne uc($mutBase));

    #Tum quals
    if(!defined($muts->{'nqs'})){
        my @empty = ();
        $muts->{'nqs'} = \@empty;
    }
    push(@{$muts->{'nqs'}},$read->{qscore});

    #Tum Rd Pos
    if(!defined($muts->{'nrp'})){
        my @empty = ();
        $muts->{'nrp'} = \@empty;
    }
    push(@{$muts->{'nrp'}},$read->{rdPos});

    #Tum rd length
    if(!defined($muts->{'nrl'})){
        my @empty = ();
        $muts->{'nrl'} = \@empty;
    }
    push(@{$muts->{'nrl'}},$read->{ln});

    return;
}

sub _callbackMatchedNormFetch{
  my ($self, $a, $hts) = @_;
  my $flagValue = $a->flag;
  #Check read and mate are mapped.
  return if((int($flagValue) & 2) != 2); # Proper pair check
  return if((int($flagValue) & 3852) != 0);
  # Ensure that we keep
  return if((int($flagValue) & 16) != 0 && (int($flagValue) & 32) != 0);
  return if((int($flagValue) & 16) == 0 && (int($flagValue) & 32) == 0);

  my $algn = Bio::DB::HTS::AlignWrapper->new($a,$hts);
  my $pos = $a->pos;
  my $cigar_array = $algn->cigar_array; # expensive and reused so save to variable
  #Quick check that were covering the base with this read (skips/indels are ignored)
  my $is_covered = _isCurrentPosCoveredFromAlignment($pos, $cigar_array, $self->_currentPos()); #1 is covered, -1 is covered but within indel
  if($is_covered == 1){
    my $this_read;

    my $rdname = $a->qname;
    #Read strand, faster than using $a->strand
    my $str = 1;
    if($algn->reversed){
      $str = -1;
    }
    my $cig_str = $algn->cigar_str; # expensive and reused so save to variable
    #Read population
    $this_read->{str} = $str;

    #Get the correct read position.
    my ($rdPosIndexOfInterest,$currentRefPos) = $self->_getReadPositionFromAlignment($pos, $cigar_array);
    $this_read->{qbase} = substr $a->qseq, $rdPosIndexOfInterest-1, 1;
    $this_read->{qscore} = unpack('C*', substr($a->_qscore, $rdPosIndexOfInterest-1, 1));
    $this_read->{rdPos} = $rdPosIndexOfInterest;

    #Read population
    $this_read->{matchesindel} = ($cig_str =~ m/[ID]/);
    $this_read->{xt} = $a->aux_get('XT');
    $this_read->{ln} = $a->l_qseq;
    $this_read->{softclipcount} = 0;
    my $type = $const->cigar_types('SOFT_CLIP_CIG');
    if ($cig_str =~ m/$type/){
      $this_read->{softclipcount} = $self->_get_soft_clip_count_from_cigar($cigar_array);
    }
    $this_read->{primaryalnscore} = $a->aux_get('AS');# $algn->get_tag_values('AS');
    $this_read->{qual} = $a->qual;
    $this_read->{start} = $algn->start;
    $this_read->{rdName} = $rdname;

    my $this_norm_read_nams_hash = $self->_norm_readnames_hash();
    my $this_norm_readnames_arr = $self->_norm_readnames_arr();
    if(!exists $this_norm_read_nams_hash->{$rdname}){
      push(@$this_norm_readnames_arr, $rdname);
      $this_norm_read_nams_hash->{$rdname} = 0;
    }
    my $this_norm_readnames = $self->_norm_readnames();
    $norm_readnames->{$rdname}->{$str} = $this_read;
    $self->_norm_readnames_hash($this_norm_read_nams_hash);
    $self->_norm_readnames_arr($this_norm_readnames_arr);
    $self->_norm_readnames($this_norm_readnames);
  } # End of if this is a position covered by this alignment
  return 1;
}

sub DESTROY{
  my $self = shift;
  $refBase = undef;
  $mutBase = undef;
  $muts = undef;
  $norms = undef;
  $tum_readnames = undef;
  $norm_readnames = undef;
  $tum_readnames_hash = undef;
  $norm_readnames_hash = undef;
  $tum_readnames_arr = undef;
  $norm_readnames_arr = undef;
  #warn "Base::DESTROY\n";
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
  $self->{'cmnp'} = undef;
  $self->{'gapflg'} = undef;
  return 1;
}

#-----------------------------
#  Getters/setters
#-----------------------------

sub meanMapQualGapFlag{
  my ($self, $p) = @_;
  if(defined($p)){
    $self->{'minMeanMapQualGapFlag'} = $p;
  }else{
    if(!defined($self->{'minMeanMapQualGapFlag'})){
      $self->{'minMeanMapQualGapFlag'} = $const->default_flag_values('MEAN_MAPQ_GAPFLAG');
    }
  }
  return $self->{'minMeanMapQualGapFlag'};
}

sub maxCavemanMatchedNormalProportion{
    my ($self, $val) = @_;
    if(defined($val)){
		 $self->{'cmnmmp'} = $val;
	}else{
		if(!defined($self->{'cmnmmp'})){
			$self->{'cmnmmp'} = $const->default_flag_values('CAVEMAN_MATCHED_NORMAL_MAX_MUT_PROP');
		}
	}
	return $self->{'cmnmmp'};
}

sub maxGapFlagDistFromEndOfReadProp{
    my ($self, $p) = @_;
    if(defined($p)){
         $self->{'maxGapFlagDistFromEndOfReadProp'} = $p;
    }else{
        if(!defined($self->{'maxGapFlagDistFromEndOfReadProp'})){
            $self->{'maxGapFlagDistFromEndOfReadProp'} = $const->default_flag_values('MAX_GAP_DIST_FROM_EOR');
        }
    }
    return $self->{'maxGapFlagDistFromEndOfReadProp'};
}

sub withinXBpOfDeletion{
    my ($self, $p) = @_;
    if(defined($p)){
         $self->{'getWithinXBpOfDeletion'} = $p;
    }else{
        if(!defined($self->{'getWithinXBpOfDeletion'})){
            $self->{'getWithinXBpOfDeletion'} = $const->default_flag_values('WITHIN_XBP_OF_DEL');
        }
    }
    return $self->{'getWithinXBpOfDeletion'};
}

sub minGapFlagDistEndOfReadPercent{
  my ($self, $p) = @_;
    if(defined($p)){
         $self->{'minGapFlagDistEndOfReadPercent'} = $p;
    }else{
        if(!defined($self->{'minGapFlagDistEndOfReadPercent'})){
            $self->{'minGapFlagDistEndOfReadPercent'} = $const->default_flag_values('MIN_GAP_DIST_PCT');
        }
    }
    return $self->{'minGapFlagDistEndOfReadPercent'};
}

sub minGapPresentInPercentReads{
    my ($self, $p) = @_;
    if(defined($p)){
         $self->{'minGapPresentInReads'} = $p;
    }else{
        if(!defined($self->{'minGapPresentInReads'})){
            $self->{'minGapPresentInReads'} = $const->default_flag_values('MIN_GAP_IN_PCT_READS');
        }
    }
    return $self->{'minGapPresentInReads'};
}

sub maxMatchedNormalAlleleProportion{
  my ($self,$p) = @_;
  if(defined($p)){
     $self->{'maxMatchedNormalAlleleProportion'} = $p;
  }else{
    if(!defined($self->{'maxMatchedNormalAlleleProportion'})){
      $self->{'maxMatchedNormalAlleleProportion'} = $const->default_flag_values('MAX_MATCHED_NORM_MUT_ALLELE_PROP');
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
      $self->{'maxPhasingMinorityStrandReadProportion'} = $const->default_flag_values('MAX_PHASING_MINORITY_STRAND_PROP');
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
      $self->{'readPosBeginningOfReadIgnoreProportion'} = $const->default_flag_values('RD_POS_BEGINNING_OF_RD_PROP');
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
      $self->{'readPosTwoThirdsOfReadExtendProportion'} = $const->default_flag_values('RD_POS_END_OF_TWOTHIRDS_EXTEND_PROP');
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
      $self->{'pentamerMinPassAvgQual'} = $const->default_flag_values('MIN_PASS_AVG_QUAL_PENTAMER');
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
      $self->{'pctSamePos'} = $const->default_flag_values('SAME_RD_POS_PERCENT');
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
      $self->{'maxTumIndelProp'} = $const->default_flag_values('MAX_TUM_INDEL_PROP');
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
      $self->{'maxNormIndelProp'} = $const->default_flag_values('MAX_NORM_INDEL_PROP');
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
      $self->{'minAvgMq'} = $const->default_flag_values('MIN_AVG_MAP_QUAL');
    }
  }
  return $self->{'minAvgMq'};
}

sub minPassAvgBaseQualPhasing{
  my ($self,$bq) = @_;
  if(defined($bq)){
    $self->{'phaseQual'} = $bq;
  }elsif(!defined($self->{'phaseQual'})){
    $self->{'phaseQual'} = $const->default_flag_values('MIN_AVG_PHASING_BASE_QUAL');
  }
  return $self->{'phaseQual'};
}

sub depthCutoffProportion{
  my ($self,$dp) = @_;
  if(defined($dp)){
    $self->{'dc'} = $dp;
  }else{
    if(!defined($self->{'dc'})){
      $self->{'dc'} = $const->default_flag_values('DEPTH_CUTOFF_PROP');
    }
  }
  return $self->{'dc'};
}

sub minDepthQual{
  my ($self,$q) = @_;
  if(defined($q)){
    $self->{'d'} = $q;
  }else{
    if(!defined($self->{'d'})){
      $self->{'d'} = $const->default_flag_values('MIN_DEPTH_QUAL');
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
      $self->{'q'} = $const->default_flag_values('MIN_NORM_MUT_ALLELE_BASE_QUAL');
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
      $self->{'minRdPsDpth'} = $const->default_flag_values('MIN_RD_POS_DEPTH');
    }
  }
  return $self->{'minRdPsDpth'};
}

#-----------------------------
#  Post processing flags
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
    $t_bam->header->parse_region($self->_chromosome().":".$self->_currentPos()."-".$self->_currentPos()),
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
    my $readPosOfMut = ($self->_currentPos() - ($rd->pos + 1)) + 1;
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
    if($overCutoff >= ($depth * $self->depthCutoffProportion)){
      return 1;
    }
  }
  return 0;
}

sub getCavemanMatchedNormalResult{
  my ($self, $vcf, $x, $index) = @_;
  my $norm_title = $const->vcf_columns('VCF_COLUMN_NORMAL');
  my $tum_title = $const->vcf_columns('VCF_COLUMN_TUMOUR');
  my $format_title = $const->vcf_columns('VCF_COLUMN_FORMAT');
  my $normal_col = $vcf->get_column($x,$norm_title);
  my $tumour_col = $vcf->get_column($x,$tum_title);
  my $format = $vcf->get_column($x,$format_title);
  if(!defined($self->{'cmnp'})){
    $self->{'cmnp'} = $self->_checkCavemanMatchedNormal($normal_col, $tumour_col, $format);
  }
  return $self->{'cmnp'};
}

sub _checkCavemanMatchedNormal{
  my ($self, $normal_col, $tumour_col, $format) = @_;
  my @splitnorm = split(/:/,$normal_col);
  my @splittum = split(/:/,$tumour_col);
  my @splitformat = split(/:/,$format);
  my $old_format = $const->allele_format('OLD_ALLELE_VCF_FORMAT');
  my $new_format = $const->allele_format('NEW_ALLELE_VCF_FORMAT');
  if($format !~ m/$old_format/ && $format !~ m/$new_format/){
    croak("VCF input format $format for cavemanMatchedNormal doesn't match a known CaVEMan VCF output format");
  }
  $is_stranded_format = 0 if($format =~ m/$old_format/);
  my $total_norm_cvg = 0;
  my $mut_allele_cvg = 0;
  my $total_tumm_cvg = 0;
  my $mut_allele_tum_cvg = 0;
  my $decode = $const->allele_format_idx_old($mutBase);
  if($is_stranded_format==1){
    $decode = $const->allele_format_idx_new($mutBase);
    $total_norm_cvg = sum(@splitnorm[1..8]);
    $total_tumm_cvg = sum(@splittum[1..8]);
  }else{
    $total_norm_cvg = sum(@splitnorm[1..4]);
    $total_tumm_cvg = sum(@splittum[1..4]);
  }
  for my $idx(@{$decode}){
    $mut_allele_cvg += $splitnorm[$idx];
    $mut_allele_tum_cvg += $splittum[$idx];
  }
  my $norm_prop = $mut_allele_cvg/$total_norm_cvg;
  my $tum_prop = $mut_allele_tum_cvg/$total_tumm_cvg;
  #Fail if the difference is less than the given proportion/percentage
  return 0 if($norm_prop > 0 && ($tum_prop - $norm_prop) < $self->maxCavemanMatchedNormalProportion());
  return 1;
}

sub getReadGapFlagResult{
    my ($self) = @_;
    if(!defined($self->{'gapflg'})){
    $self->{'gapflg'} = $self->_checkReadGap();
  }
  return $self->{'gapflg'};
}

sub _checkReadGap{
  my ($self) = @_;
  my $meanMapQ = sum(@{$self->_muts->{'allTumMapQuals'}})/scalar(@{$self->_muts->{'allTumMapQuals'}});
  return 1 if($meanMapQ < $self->meanMapQualGapFlag); #Pass as likely mismapping
  my @tum_base_dist_prop = [];
  my $tum_base_dist_count = 0;  
  my @non_tum_base_dist = [];
  my $norm_base_dist_count = 0;
  foreach (zip($self->_muts->{'allTumBases'},$self->_muts->{'allMinGapDistances'},$self->_muts->{'allMutDistPropFromEndOfRead'})){
    my ($base, $distance, $dist_eor) = @{$_};
    #Pass flag is we find a called variant base within the limits if a deletion 
    if($base eq $mutBase){ 
      return 1 if($distance != -1 && $distance <= $self->withinXBpOfDeletion);
      push(@tum_base_dist_prop, $dist_eor);
      $tum_base_dist_count++;  
    }else{ #Not a variant base - start counting the reads with an indel
      if($distance != -1 && $distance <= $self->withinXBpOfDeletion){
        push(@non_tum_base_dist, $distance);
        $norm_base_dist_count++;
      }
    }
  }
  return 1 if($norm_base_dist_count==0); #Pass if zero reference showing reads with gap
  my $total_reads = scalar(@{$self->_muts->{'allTumMapQuals'}});
  my $percentage_reads_present = ($norm_base_dist_count/$total_reads) * 100;
  my $percent_tum_reads_within_eor_dist = (scalar(grep {$_ != -1 && $_ <= $self->maxGapFlagDistFromEndOfReadProp} @tum_base_dist_prop)/$tum_base_dist_count)*100;
  #Pass if > minGapFlagDistEndOfReadPercent % reads have a distance from end of read > maxGapFlagDistFromEndOfReadProp
  return 1 if($percent_tum_reads_within_eor_dist < $self->minGapFlagDistEndOfReadPercent);
  return 0 if($percentage_reads_present >= $self->minGapPresentInPercentReads);
  return 1;
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
#  DESTROY
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
                                            'depthCutoffProportion' => (1/3),
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
  at least depthCutoffProportion of mutant alleles are of base quality >= minDepthQual
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

=item * maxGapFlagDistFromEndOfReadProp
Sets and returns the maximum gap from end of read allowed in read gap flag.

=item * getWithinXBpOfDeletion
Sets and returns the within X BP of deletion parameter. Used in the readgap flag.

=item * minGapPresentInPercentReads
Sets and returns the minimum gap present in percentage of reads parameter. Used in the readgap flag.

=item * maxGapDistance
Sets and returns the maximum distance from a read gap in the readgap flag.

=back

=head1 AUTHOR

David Jones (drj@sanger.ac.uk)

=head1 COPYRIGHT

=head1 SEE ALSO

perl(1).
Bio::DB::SAM

=cut
