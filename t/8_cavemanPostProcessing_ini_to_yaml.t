##########LICENCE##########
# Copyright (c) 2014-2021 Genome Research Ltd.
#
# Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
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

use strict;
use Test::More tests => 2;
use strict;
use warnings FATAL => 'all';
use Carp;
use Data::Dumper;
use Capture::Tiny qw(capture);
use YAML::XS qw(LoadFile);

use FindBin qw($Bin);

my $script_path = "$Bin/../bin/";
my $test_data_path = "$Bin/../testData/";
my $script = $script_path.'cavemanPostProcessing_ini_to_yaml.pl';

my $test_out_file = $test_data_path.'test_config_output.yaml';
my $test_input_file = $test_data_path.'flag.vcf.config_test.ini';
my $expect_output = $test_data_path.'test_config.yaml';

my $perl_exec = "$^X -I ".(join ':', @INC);

main(@ARGV);

sub main{
  run_flag();
  compare($test_out_file, $expect_output, __LINE__);
  unlink($test_out_file);
}

sub run_flag{
  my $cmd = "$perl_exec $script -i $test_input_file -o $test_out_file";
  my ($out, $err, $exit) = capture{ system($cmd) };
  if($exit!=0){
    warn Dumper($err, $out, $exit);
  }
  is($exit, 0,'Conversion from ini to yaml ran correctly');
}

sub compare{
  my ($test_out_file, $expect_output, $line_no) = @_;
  #Read new file
  my $new = LoadFile($test_out_file);
  my $exp = LoadFile($expect_output);
  is_deeply($new,$exp);
}
