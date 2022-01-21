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

package Sanger::CGP::CavemanPostProcessor::Constants;

use strict;

use warnings FATAL => 'all';
use autodie qw(:all);
use Const::Fast qw(const);

const my %DEFAULT_FLAG_VALUES => (
  'MATCHED_NORMAL_ALLELE_HICVG_CUTOFF' => 2,
  'MAX_MATCHED_NORMAL_ALLELE_HICVG_PROPORTION' => 0.03,
  'WITHIN_XBP_OF_DEL' => 10,
  'MIN_GAP_IN_PCT_READS' => 30,
  'MEAN_MAPQ_GAPFLAG' => 10,
  'MAX_GAP_DIST_FROM_EOR' => 0.13,
  'MIN_GAP_DIST_PCT' => 75,
  'MIN_SINGLE_END_CVG' => 10,
  'CAVEMAN_MATCHED_NORMAL_MAX_MUT_PROP' => 0.2,
  'MATCHED_NORMAL_MAX_MUT_PROP' => 0.2,
  'MAX_MATCHED_NORM_MUT_ALLELE_PROP' => 0.05,
  'MAX_PHASING_MINORITY_STRAND_PROP' => 0.04,
  'RD_POS_BEGINNING_OF_RD_PROP' => 0.08,
  'RD_POS_END_OF_TWOTHIRDS_EXTEND_PROP' => 0.08,
  'MIN_PASS_AVG_QUAL_PENTAMER' => 20,
  'SAME_RD_POS_PERCENT' => 80,
  'MAX_TUM_INDEL_PROP' => 10,
  'MAX_NORM_INDEL_PROP' => 10,
  'MIN_AVG_MAP_QUAL' => 21,
  'MIN_AVG_PHASING_BASE_QUAL' => 21,
  'MIN_DEPTH_QUAL' => 25,
  'MIN_NORM_MUT_ALLELE_BASE_QUAL' => 15,
  'MIN_RD_POS_DEPTH' => 8,
  'DEPTH_CUTOFF_PROP' => 0.333333,
);

sub default_flag_values {
  my ($class,$item) = @_;
  return %DEFAULT_FLAG_VALUES{$item};
}

const my %CIGAR_TYPES => (
  'MATCH_CIG' => 'M',
  'SKIP_CIG' => 'N',
  'INS_CIG' => 'I',
  'DEL_CIG' => 'D',
  'SOFT_CLIP_CIG' => 'S',
  'HARD_CLIP_CIG' => 'H',
);

sub cigar_types {
  my ($class,$item) = @_;
  return %CIGAR_TYPES{$item};
}

const my %ALLELE_FORMAT => (
  'OLD_ALLELE_VCF_FORMAT' => 'GT:AA:CA:GA:TA:PM',
  'NEW_ALLELE_VCF_FORMAT' => 'GT:FAZ:FCZ:FGZ:FTZ:RAZ:RCZ:RGZ:RTZ:PM',
);

sub allele_format{
  my ($class,$item) = @_;
  return %ALLELE_FORMAT{$item};
}

const my %ALLELE_FORMAT_IDX_OLD => (
  'A'=>[1,5], 
  'C'=>[2,6], 
  'G'=>[3,7], 
  'T'=>[4,8],
);

const my %ALLELE_FORMAT_IDX_NEW => (
  'A'=>[1,5],
  'C' =>[2,6],
  'G'=>[3,7],
  'T'=>[4,8],
);

sub allele_format_idx_old{
  my ($class,$item) = @_;
  return %ALLELE_FORMAT_IDX_OLD{$item};
}

sub allele_format_idx_new{
  my ($class,$item) = @_;
  return %ALLELE_FORMAT_IDX_NEW{$item};
}

1;
