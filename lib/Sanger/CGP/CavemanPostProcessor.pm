##########LICENCE##########
# Copyright (c) 2014-2019 Genome Research Ltd.
#
#Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
#
# This file is part ofcgpCaVEManPostProcessing.
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

package Sanger::CGP::CavemanPostProcessor;

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

our $VERSION = '1.8.9';
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
my $tum_readnames;
my $norm_readnames;
my $tum_readnames_arr;
my $norm_readnames_arr;
my $tum_readnames_hash;
my $norm_readnames_hash;


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
  $self->tumBam($inputs->{'tumBam'}, $inputs->{'ref'});
  $self->normBam($inputs->{'normBam'}, $inputs->{'ref'});
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
    $tum_readnames = undef;
    $tum_readnames_arr = undef;
    $tum_readnames_hash = undef;
  $self->{'tb'}->fetch($chr.':'.$start.'-'.$stop,\&_callbackTumFetch);
    process_hashed_reads(\&populate_muts, $tum_readnames, $tum_readnames_arr);
    $norm_readnames = undef;
    $norm_readnames_arr = undef;
    $norm_readnames_hash = undef;
  $self->{'nb'}->fetch($chr.':'.$start.'-'.$stop,\&_callbackMatchedNormFetch);
    process_hashed_reads(\&populate_norms, $norm_readnames, $norm_readnames_arr);

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
  my ($code, $hashed_reads, $readname_arr) = @_;

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
  my ($algn) = @_;
  my $flagValue = $algn->flag;
  #Check read and mate are mapped. If not return.
    return if((int($flagValue) & 2) != 2); # Proper pair
    return if((int($flagValue) & 3852) != 0);
    # Ensure that we keep
    return if((int($flagValue) & 16) != 0 && (int($flagValue) & 32) != 0);
    return if((int($flagValue) & 16) == 0 && (int($flagValue) & 32) == 0);
    #Calculate other stuff
    my $totalPCovg = 0;
    my $totalNCovg = 0;
    my $indelRdCount = 0;

  if(_isCurrentPosCoveredFromAlignment($algn) == 1){
        return unless ($algn->proper_pair == 1);
        # Ensure that we keep
        return if((int($flagValue) & 16) != 0 && (int($flagValue) & 32) != 0);
        return if((int($flagValue) & 16) == 0 && (int($flagValue) & 32) == 0);
        my $this_read;
        #Get the correct read position.
    my ($rdPosIndexOfInterest,$currentRefPos) = _getReadPositionFromAlignment($algn);
        my @splt = split(//,$algn->qseq);

        my $rdname = $algn->qname;

        #Read strand
    my $str = 1;
    if($algn->reversed){
      $str = -1;
    }

        #Read base
        $this_read->{str} = $str;
        $this_read->{qbase} = $splt[$rdPosIndexOfInterest-1];
        $this_read->{matchesindel} = ($algn->cigar_str =~ m/[ID]/);
        $this_read->{qscore} = $algn->qscore->[$rdPosIndexOfInterest-1];
        $this_read->{xt} = $algn->aux_get('XT');
        $this_read->{ln} = $algn->l_qseq;
        $this_read->{rdPos} = $rdPosIndexOfInterest;
        $this_read->{softclipcount} = 0;
        if ($algn->cigar_str =~ m/$SOFT_CLIP_CIG/){
           $this_read->{softclipcount} = _get_soft_clip_count_from_cigar($algn->cigar_array);
        }
        $this_read->{primaryalnscore} = $algn->get_tag_values('AS');
        $this_read->{qual} = $algn->qual;
        $this_read->{start} = $algn->start;
        $this_read->{rdName} = $rdname;

        if(!exists $tum_readnames_hash->{$rdname}){
            push(@$tum_readnames_arr, $rdname);
            $tum_readnames_hash->{$rdname} = 0;
        }
        $tum_readnames->{$rdname}->{$str} = $this_read;

    } # End of if this is a covered position

  return 1;
}

sub populate_muts{
    my ($read) = @_;
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

    return if ($keepSW == 0 && defined($read->{xt}) && $read->{xt} eq 'M');

  return if($read->{qscore} < $minAnalysedQual);

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

sub populate_norms{
    my ($read) = @_;
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

    return if ($keepSW == 0 && defined($read->{xt}) && $read->{xt} eq 'M');

  return if($read->{qscore} < $minAnalysedQual);

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
  my ($algn) = @_;
  my $flagValue = $algn->flag;
  #Check read and mate are mapped.
    return if((int($flagValue) & 2) != 2); # Proper pair check
  return if((int($flagValue) & 3852) != 0);
    # Ensure that we keep
    return if((int($flagValue) & 16) != 0 && (int($flagValue) & 32) != 0);
    return if((int($flagValue) & 16) == 0 && (int($flagValue) & 32) == 0);
  #Quick check that were covering the base with this read (skips/indels are ignored)

    #Calculate other stuff
    my $totalPCovg = 0;
    my $totalNCovg = 0;
    my $indelRdCount = 0;

    if(_isCurrentPosCoveredFromAlignment($algn) == 1){
        return unless ($algn->proper_pair == 1);
        # Ensure that we keep
        return if((int($flagValue) & 16) != 0 && (int($flagValue) & 32) != 0);
        return if((int($flagValue) & 16) == 0 && (int($flagValue) & 32) == 0);
        my $this_read;
    #Get the correct read position.
    my ($rdPosIndexOfInterest,$currentRefPos) = _getReadPositionFromAlignment($algn,$currentPos);
    my @splt = split(//,$algn->qseq);

        my $rdname = $algn->qname;

        #Read strand
    my $str = 1;
    if($algn->reversed){
      $str = -1;
    }

        #Read population
        $this_read->{str} = $str;
        $this_read->{qbase} = $splt[$rdPosIndexOfInterest-1];
        $this_read->{matchesindel} = ($algn->cigar_str =~ m/[ID]/);
        $this_read->{qscore} = $algn->qscore->[$rdPosIndexOfInterest-1];
        $this_read->{xt} = $algn->aux_get('XT');
        $this_read->{ln} = $algn->l_qseq;
        $this_read->{rdPos} = $rdPosIndexOfInterest;
        $this_read->{softclipcount} = 0;
        if ($algn->cigar_str =~ m/$SOFT_CLIP_CIG/){
           $this_read->{softclipcount} = _get_soft_clip_count_from_cigar($algn->cigar_array);
        }
        $this_read->{primaryalnscore} = $algn->get_tag_values('AS');
        $this_read->{qual} = $algn->qual;
        $this_read->{start} = $algn->start;
        $this_read->{rdName} = $rdname;

        if(!exists $norm_readnames_hash->{$rdname}){
            push(@$norm_readnames_arr, $rdname);
            $norm_readnames_hash->{$rdname} = 0;
        }

        $norm_readnames->{$rdname}->{$str} = $this_read;

  } # End of if this is a position covered by this alignment
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
    $tum_readnames = undef;
    $norm_readnames = undef;
  #warn "Base::DESTROY\n";
}

return 1;
