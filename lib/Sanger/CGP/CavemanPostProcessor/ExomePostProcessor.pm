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

package Sanger::CGP::CavemanPostProcessor::ExomePostProcessor;

use strict;
use Carp;
use Const::Fast qw(const);

use Sanger::CGP::CavemanPostProcessor;
our $VERSION = Sanger::CGP::CavemanPostProcessor->VERSION;

use base qw(Sanger::CGP::CavemanPostProcessor::PostProcessor);

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

#---------------
#	Init methods
#---------------

sub _init{
	my ($self,$inputs) = @_;
	$self->SUPER::_init($inputs);
	return $self;
}

#-----------------------------
#	Post processing filter methods
#-----------------------------


#-----------------
#	Getters/setters
#-----------------


#----------
#	DESTROY!
#----------

sub DESTROY{
	my $self = shift;
	#warn "ExomePostProcessor::DESTROY\n";
	$self->SUPER::DESTROY;
}

return 1;
