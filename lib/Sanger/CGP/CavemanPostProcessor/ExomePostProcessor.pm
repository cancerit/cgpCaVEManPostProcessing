####################################################
# Copyright (c) 2012 Genome Research Ltd.
# Author: Cancer Genome Project, cgpit@sanger.ac.uk
# See LICENCE.TXT for details
####################################################

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
