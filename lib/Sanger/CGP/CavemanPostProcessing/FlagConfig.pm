##########LICENCE##########
# Copyright (c) 2014-2016 Genome Research Ltd.
#
# Author: Cancer Genome Project cgpit@sanger.ac.uk
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
package Sanger::CGP::CavemanPostProcessing::FlagConfig;

use strict;
use warnings FATAL=>'all';
use autodie;

use Sanger::CGP::CavemanPostProcessing;

sub new {
	my ($proto) = shift;
	my $class = ref($proto) || $proto;
 	my $self = {};
  bless($self, $class);
  return $self;
}

sub name{
  my ($self,$nom) = @_;
  if(defined($nom)){
    $self->{'n'} = $nom;
  }
  return $self->{'n'}
}

sub is_intersect{
  my ($self,$val) = @_;
  if(defined($val)){
    $self->{'in'} = $val;
  }
  return $self->{'in'};
}

sub id{
  my ($self,$val) = @_;
  if(defined($val)){
    $self->{'id'} = $val;
  }
  return $self->{'id'};
}

sub description{
  my ($self,$val) = @_;
  if(defined($val)){
    $self->{'d'} = $val;
  }
  return $self->{'d'};
}

sub is_info{
  my ($self,$val) = @_;
  if(defined($val)){
    $self->{'info'} = $val;
  }
  return $self->{'info'};
}

sub value{
  my ($self,$val) = @_;
  if(defined($val)){
    $self->{'v'} = $val;
  }
  return $self->{'v'};
}

sub type{
  my ($self,$val) = @_;
  if(defined($val)){
    $self->{'t'} = $val;
  }
  return $self->{'t'};
}

sub intersect_file{
  my ($self,$val) = @_;
  if(defined($val)){
    $self->{'if'} = $val;
  }
  return $self->{'if'};
}

1;