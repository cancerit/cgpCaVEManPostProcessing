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
use List::Util qw (sum zip);
use Carp;
use Const::Fast qw(const);

use Sanger::CGP::CavemanPostProcessor;
our $VERSION = Sanger::CGP::CavemanPostProcessor->VERSION;

sub _init{
  my ($self,$inputs) = @_;

  $self->SUPER::_init($inputs);

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
  $self->{'cmnp'} = undef;
  $self->{'gapflg'} = undef;
  return 1;
}

return 1;
