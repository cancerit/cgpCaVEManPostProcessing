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

use strict;
use Data::Dumper;
use Test::More;
use List::Util qw(first);
use File::Find;
use Cwd;
use Try::Tiny qw(try finally);
use File::Spec;

use FindBin qw($Bin);
my $lib_path = "$Bin/../lib";
my $test_data_path = "$Bin/../testData/";

# Add modules here that cannot be instantiated (should be extended and have no 'new')
# or need a set of inputs - these should be tested in own test script
use constant MODULE_SKIP => qw( Sanger::CGP::CavemanPostProcessor Sanger::CGP::CavemanPostProcessor::PostProcessor Sanger::CGP::CavemanPostProcessing);
use constant VER_MODULE_SKIP => qw ( Sanger::CGP::CavemanPostProcessor::ConfigParser );

my $init_cwd = getcwd;

my @modules;
try {
  chdir($lib_path);
  find({ wanted => \&build_module_set, no_chdir => 1 }, './');
} finally {
  chdir $init_cwd;
  die "The try block died with: @_\n" if(@_);
};

for my $mod(@modules) {
  use_ok($mod) or BAIL_OUT("Unable to 'use' module $mod");
}

for my $mod(@modules) {
  ok($mod->VERSION, "Check version inheritance exists ($mod)") unless (first {$mod eq $_} VER_MODULE_SKIP);
  if($mod->can('new')) { # only try new on things that have new defined
    new_ok($mod => [tumBam=> $test_data_path."test.bam",normBam => $test_data_path."test.bam"]) unless( first {$mod eq $_} MODULE_SKIP );
  }
}

done_testing();

sub build_module_set {
  if($_ =~ m/\.pm$/) {

    my ($dir_str,$file) = (File::Spec->splitpath( $_ ))[1,2];
    $file =~ s/\.pm$//;
    my @dirs = File::Spec->splitdir( $dir_str );
    shift @dirs;
    push @modules, (join '::', @dirs).$file;
  }
}
