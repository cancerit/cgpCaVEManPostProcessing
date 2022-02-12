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

package Sanger::CGP::CavemanPostProcessor::MNVPostProcessor;

use strict;
use Bio::DB::HTS;
use Bio::DB::HTS::Constants;
use Bio::DB::HTS::Alignment;
use POSIX qw(strftime ceil);
use List::Util qw (sum zip min max);
use Carp;
use Data::Dumper;
use Const::Fast qw(const);

use Sanger::CGP::CavemanPostProcessor;
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

sub _init{
  my ($self,$inputs) = @_;
  if(!defined($inputs->{'tumBam'}) || !defined($inputs->{'normBam'})){
    croak("tumBam and normBam are required for initialisation.\n");
  }
  $self->tumBam($inputs->{'tumBam'}, $inputs->{'ref'});
  $self->normBam($inputs->{'normBam'}, $inputs->{'ref'});
  $self->keepSW($inputs->{'keepSW'}) if exists $inputs->{'keepSW'};
  $self->minAnalysedQual($inputs->{'minAnalysedQual'}) if exists $inputs->{'minAnalysedQual'};
  $self->SUPER::_init($inputs);

  return $self;
}

sub process_hashed_reads{
  my ($self, $code, $hashed_reads, $readname_arr) = @_;

  my %loc_counts = (
    1  => {},
    -1 => {},
  );

  foreach my $readnom(@$readname_arr){
    my $read = $hashed_reads->{$readnom};
    my ($read_to_use, $read_to_use_2);
    if (exists $read->{1} && ! exists $loc_counts{1}{$read->{1}->{qbase}}){
      $loc_counts{1}{$read->{1}->{qbase}} = 0;
    }
    if (exists $read->{-1} && ! exists $loc_counts{-1}{$read->{-1}->{qbase}}){
      $loc_counts{-1}{$read->{-1}->{qbase}} = 0;
    }
    if(exists $read->{1} && exists $read->{-1}) {
      # if change is different have to include both versions
      if( $read->{1}->{qbase} ne $read->{-1}->{qbase}) {
        $read_to_use = $read->{1};
        $read_to_use_2 = $read->{-1};
        $loc_counts{1}{$read->{1}->{qbase}}++;
        $loc_counts{-1}{$read->{-1}->{qbase}}++;
      }
      # score is same check loc counts for that allele to determine where to put
      elsif(sum($read->{1}->{qscore}) == sum($read->{-1}->{qscore})) {
        if($loc_counts{1}{$read->{1}->{qbase}} <= $loc_counts{-1}{$read->{-1}->{qbase}}) {
          $read_to_use = $read->{1};
          $loc_counts{1}{$read->{1}->{qbase}}++;
        }
        else {
          $read_to_use = $read->{-1};
          $loc_counts{-1}{$read->{-1}->{qbase}}++;
        }
      }
      elsif(sum($read->{1}->{qscore}) > sum($read->{-1}->{qscore})) {
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
  #1 is covered, -1 is covered but within indel
  my $is_covered = _isCurrentPosCoveredFromAlignment($pos, 
                                                    $cigar_array, 
                                                    $self->_currentPos(), 
                                                    ($self->_currentPos()+length($self->mutBase())-1));
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
      #Get the correct read position (start of the MNV).
      my ($rdPosIndexOfInterest,$currentRefPos) = $self->_getReadPositionFromAlignment($pos, $cigar_array);
      $this_read->{qbase} = substr $a->qseq, $rdPosIndexOfInterest-1, length($self->mutBase());
      $this_read->{qscore} = [unpack('C*', substr($a->_qscore, $rdPosIndexOfInterest-1, length($self->mutBase())))];
      $this_read->{rdPos} = $rdPosIndexOfInterest;
      $this_read->{gapDist} = _getDistanceFromGapInRead($algn->cigar_array,$rdPosIndexOfInterest, ($rdPosIndexOfInterest+length($self->mutBase()))-1);
      $this_read->{distFromEnd} = min(
                    ($rdPosIndexOfInterest/$this_read->{ln}),
                    (($this_read->{ln}-($rdPosIndexOfInterest+(length($self->mutBase())-1)))/$this_read->{ln}));
    }
    $this_read->{matchesindel} = ($cig_str =~ m/[ID]/);
    $this_read->{xt} = $a->aux_get('XT');
    $this_read->{softclipcount} = 0;
    my $type = $const->cigar_types('SOFT_CLIP_CIG');
    if ($cig_str =~ m/$type/){
      $this_read->{softclipcount} = _get_soft_clip_count_from_cigar($algn->cigar_array);
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
  $self->_muts->{'totalTCoverage'} += 1;
  if($read->{str} == 1){
    $self->_muts->{'totalTCoveragePos'} += 1;
  }else{
    $self->_muts->{'totalTCoverageNeg'} += 1;
  }
  if($read->{xt}){
    $self->_muts->{'indelTCount'} += 1;
  }
  push(@{$self->_muts->{'completeMutStrands'}},$read->{str});
  push(@{$self->_muts->{'allTumBases'}},$read->{qbase});
  push(@{$self->_muts->{'allTumBaseQuals'}},$read->{qscore});
  push(@{$self->_muts->{'allTumStrands'}},$read->{str});
  push(@{$self->_muts->{'allTumMapQuals'}},$read->{qual});
  push(@{$self->_muts->{'allMinGapDistances'}},$read->{gapDist});
  push(@{$self->_muts->{'allMutDistPropFromEndOfRead'}},$read->{distFromEnd});
  return if ($self->keepSW() == 0 && defined($read->{xt}) && $read->{xt} eq 'M');

  return if(min($read->{qscore}) < $self->minAnalysedQual());

  $self->_muts->{'tumcvg'} += 1;

  if($read->{str} == 1){
    $self->_muts->{'pcvg'} += 1;
  }else{
    $self->_muts->{'ncvg'} += 1;
    $read->{rdPos} = ($read->{ln} - $read->{rdPos}) + 1;
  }

  return if(uc($read->{qbase}) ne uc($self->_mutBase()));

  #Tum quals
  push(@{$self->_muts->{'tqs'}},$read->{qscore});

  #Tum Rd Pos
  push(@{$self->_muts->{'trp'}},$read->{rdPos});

  #Tum rd length
  push(@{$self->_muts->{'trl'}},$read->{ln});

  #Tum XT tags
  push(@{$self->_muts->{'txt'}},$read->{xt});

  #Tum rd start
  push(@{$self->_muts->{'trdst'}},$read->{start});

  #Strands
  push(@{$self->_muts->{'tstr'}},$read->{str});

  #RdNames
  push(@{$self->_muts->{'trn'}},$read->{rdName});

  #Mapping quals
  push(@{$self->_muts->{'tmq'}},$read->{qual});

  #AlnScoresPrm
  push(@{$self->_muts->{'alnp'}},$read->{primaryalnscore});

  #Softclipping
  push(@{$self->_muts->{'sclp'}},$read->{softclipcount});

  return;
}

sub _getDistanceFromGapInRead{
  my ($cigar_array,$rdPosIndexOfInterest_start, $rdPosIndexOfInterest_end) = @_;
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
      if($rdPosIndexOfInterest_start>$dp_end || ($rdPosIndexOfInterest_start<$dp_start && $rdPosIndexOfInterest_end < $dp_start)){
        $tmp_dist = min(
          min(abs($rdPosIndexOfInterest_start-$dp_start),abs($dp_end-$rdPosIndexOfInterest_start)),
          min(abs($rdPosIndexOfInterest_end-$dp_start),abs($dp_end-$rdPosIndexOfInterest_end))
        );
      }
      if($tmp_dist < $min_gap_dist || $min_gap_dist == -1){
        $min_gap_dist = $tmp_dist;
      }
    }
  }
  return $min_gap_dist;
}


sub _isCurrentPosCoveredFromAlignment{
  my ($pos, $cigar_array, $pos_of_interest_start, $pos_of_interest_stop) = @_; # 0-based pos, 1 based current pos
  foreach my $cigSect(@{$cigar_array}){

    if($cigSect->[0] eq $const->cigar_types('MATCH_CIG')){
      if($pos_of_interest_start >= ($pos + 1) && $pos_of_interest_stop <= ($pos+$cigSect->[1])){
        return 1;
      }
      $pos+= $cigSect->[1];
    }elsif($cigSect->[0] eq $const->cigar_types('DEL_CIG') || $cigSect->[0] eq $const->cigar_types('SKIP_CIG')){
      if($pos_of_interest_start >= ($pos + 1) && $pos_of_interest_stop <= ($pos+$cigSect->[1]) ||
        $pos_of_interest_start >= ($pos + 1) && $pos_of_interest_start <= ($pos+$cigSect->[1]) ||
         $pos_of_interest_stop >= ($pos + 1) && $pos_of_interest_stop <= ($pos+$cigSect->[1]) ){
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
    if(!defined($self->_muts->{'totalNCoverage'})){
        $self->_muts->{'totalNCoverage'} = 0;
    }
    $self->_muts->{'totalNCoverage'} += 1;

    if(!defined($self->_muts->{'allNormBases'})){
        $self->_muts->{'allNormBases'} = [];
    }
    push(@{$self->_muts->{'allNormBases'}},$read->{qbase});

    if(!defined($self->_muts->{'allNormBaseQuals'})){
        $self->_muts->{'allNormBaseQuals'} = [];
    }
    push(@{$self->_muts->{'allNormBaseQuals'}},$read->{qscore});

    if(!defined($self->_muts->{'allNormStrands'})){
        $self->_muts->{'allNormStrands'} = [];
    }
    push(@{$self->_muts->{'allNormStrands'}},$read->{str});

    return if ($self->keepSW() == 0 && defined($read->{xt}) && $read->{xt} eq 'M');

  return if(min($read->{qscore}) < $self->minAnalysedQual());

    if(!defined($self->_muts->{'normcvg'})){
        $self->_muts->{'normcvg'} = 0;
    }
    $self->_muts->{'normcvg'} += 1;

    if($read->{str} == +1){
        $self->_muts->{'npcvg'} += 1;
    }else{
        $self->_muts->{'nncvg'} += 1;
        $read->{rdPos} = ($read->{ln} - $read->{rdPos}) + 1;
    }

  return if(uc($read->{qbase}) ne uc($self->_mutBase()));

    #Tum quals
    if(!defined($self->_muts->{'nqs'})){
        my @empty = ();
        $self->_muts->{'nqs'} = \@empty;
    }
    push(@{$self->_muts->{'nqs'}},$read->{qscore});

    #Tum Rd Pos
    if(!defined($self->_muts->{'nrp'})){
        my @empty = ();
        $self->_muts->{'nrp'} = \@empty;
    }
    push(@{$self->_muts->{'nrp'}},$read->{rdPos});

    #Tum rd length
    if(!defined($self->_muts->{'nrl'})){
        my @empty = ();
        $self->_muts->{'nrl'} = \@empty;
    }
    push(@{$self->_muts->{'nrl'}},$read->{ln});

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
  my $is_covered = _isCurrentPosCoveredFromAlignment($pos, 
                                            $cigar_array, 
                                            $self->_currentPos(), 
                                            ($self->_currentPos()+length($self->mutBase())-1)
                  );
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
    $this_read->{qscore} = [unpack('C*', substr($a->_qscore, $rdPosIndexOfInterest-1, length($self->mutBase())))];
    $this_read->{rdPos} = $rdPosIndexOfInterest;

    #Read population
    $this_read->{matchesindel} = ($cig_str =~ m/[ID]/);
    $this_read->{xt} = $a->aux_get('XT');
    $this_read->{ln} = $a->l_qseq;
    $this_read->{softclipcount} = 0;
    my $type = $const->cigar_types('SOFT_CLIP_CIG');
    if ($cig_str =~ m/$type/){
      $this_read->{softclipcount} = _get_soft_clip_count_from_cigar($cigar_array);
    }
    $this_read->{primaryalnscore} = $a->aux_get('AS');# $algn->get_tag_values('AS');
    $this_read->{qual} = $a->qual;
    $this_read->{start} = $algn->start;
    $this_read->{rdName} = $rdname;

   my $this_norm_read_names_hash = $self->_norm_readnames_hash();
   my $this_norm_readnames_arr = $self->_norm_readnames_arr();
   if(!exists $this_norm_read_names_hash->{$rdname}){
     push(@$this_norm_readnames_arr, $rdname);
     $this_norm_read_names_hash->{$rdname} = 0;
   }
   my $this_norm_readnames = $self->_norm_readnames();
   $this_norm_readnames->{$rdname}->{$str} = $this_read;
   $self->_norm_readnames_hash($this_norm_read_names_hash);
   $self->_norm_readnames_arr($this_norm_readnames_arr);
   $self->_norm_readnames($this_norm_readnames);

  } # End of if this is a position covered by this alignment
  return 1;
}

sub getDepthResult{
  my ($self) = @_;
  if(!defined($self->{'depth'})){
    $self->{'depth'} = $self->_checkDepth();
  }
  return $self->{'depth'};
}

sub _checkDepth{
  my ($self) = @_;
  my $depth = scalar(@{$self->_muts->{'tqs'}});
  my $overCutoff = 0;
  foreach my $q(@{$self->_muts->{'tqs'}}){
    $q = sum($q) / scalar(@$q); # Use mean of all qualities in the MNV
    if($q >= $self->minDepthQual){
      $overCutoff++;
    }
    if($overCutoff >= ($depth * $self->depthCutoffProportion)){
      return 1;
    }
  }
  return 0;
}

return 1;
