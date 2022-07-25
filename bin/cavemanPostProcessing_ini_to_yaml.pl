#!/usr/bin/perl
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

use strict;
use warnings FATAL => 'all';

use Carp qw(croak cluck);
use English qw( -no_match_vars );

use Data::Dumper;

use Getopt::Long qw(:config pass_through);
use Pod::Usage;
use Sanger::CGP::CavemanPostProcessor;
use Sanger::CGP::CavemanPostProcessor::ConfigParser;

my $opts = option_builder();
validateInput($opts);
main($opts);

sub main {
  my ($opts) = @_;
  Sanger::CGP::CavemanPostProcessor::ConfigParser::convert_ini_to_yaml($opts->{'f'}, $opts->{'o'});
  return;
}

sub option_builder {
  my ($factory) = @_;

  my %opts = ();

  my $result = &GetOptions (
    'h|help' => \$opts{'h'},
    'i|input=s' => \$opts{'f'},
    'o|outFile=s' => \$opts{'o'},
    'version' => \$opts{'version'},
  );
  return \%opts;
}

sub validateInput {
  my $opts = shift;
  pod2usage(0) if($opts->{'h'});

  if(defined $opts->{'version'}) {
    print sprintf "VERSION: %s\n", Sanger::CGP::CavemanPostProcessor->VERSION;
    exit 0;
  }
  delete $opts->{'version'}; # needs to be deleted or breaks tests

  croak( "Unknown parameter: ".$ARGV[0]) if(scalar(@ARGV) > 0);
  pod2usage("Missing parameter 'input'") if(!defined($opts->{'f'}));
  pod2usage("Missing parameter 'outfile'") if(!defined($opts->{'o'}));
  unless(-e $opts->{'f'} && -r $opts->{'f'}){
    pod2usage("Inpit .ini file does not exist or has incorrect permissions: ".$opts->{'f'}."\n");
  }
  return;
}

__END__

=head1 NAME

cavemanPostProcessing_ini_to_yaml.pl - Utility script to convert old caveman .ini format to .yaml.

=head1 SYNOPSIS

cavemanPostProcessing_ini_to_yaml.pl [-h] -f flag_config.ini -o flag_config.yml 

  General Options:

    --help                 (-h)       Brief documentation

    --version              (-version) Output the version number and exit

    --input                (-i)       The VCF input file to flag.

    --outfile              (-o)       The VCF output file to write.

  Examples:

    cavemanPostProcessing_ini_to_yaml.pl [-h] -f flag_config.ini -o flag_config.yml

=cut
