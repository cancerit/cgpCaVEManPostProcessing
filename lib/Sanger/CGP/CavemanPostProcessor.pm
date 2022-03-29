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

package Sanger::CGP::CavemanPostProcessor;

use strict;
use Sanger::CGP::CavemanPostProcessor::Constants;
use POSIX qw(strftime);
use Carp;
use Const::Fast qw(const);
use Attribute::Abstract;
use Data::Dumper;
use base 'Exporter';

our $VERSION = '1.11.0';
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
  $self->_init(\%inputs);
  return $self;
}

sub VERSION{
  return $VERSION;
}

return 1;
