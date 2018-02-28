##########LICENCE##########
# Copyright (c) 2014-2018 Genome Research Ltd.
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
use Sanger::CGP::CavemanPostProcessing::FlagConfig;
use Data::Dumper;
use Const::Fast qw(const);

use Test::More tests => 2;

use FindBin qw($Bin);
const my $lib_path => "$Bin/../lib";
const my $test_data_path => "$Bin/../testData/";

subtest 'FlagConfig new' => sub {
  my $cfg = new_ok('Sanger::CGP::CavemanPostProcessing::FlagConfig');
  done_testing();
};

subtest 'FlagConfig methods' => sub {
  my $cfg = new_ok('Sanger::CGP::CavemanPostProcessing::FlagConfig');
  $cfg->name("test");
  ok($cfg->name eq "test","Test name");
  $cfg->description("description");
  ok($cfg->description eq "description","Test description");
  $cfg->type("type");
  ok($cfg->type eq "type","Test type");
  $cfg->type("intersect_file");
  ok($cfg->type eq "intersect_file","Test intersect_file");
  $cfg->type("id");
  ok($cfg->type eq "id","Test id");
  $cfg->value("id");
  ok($cfg->value eq "id","Test value");
  $cfg->is_intersect(1);
  ok($cfg->is_intersect==1,"Test intersect on");
  $cfg->is_intersect(0);
  ok($cfg->is_intersect==0,"Test intersect off");
  $cfg->is_info(1);
  ok($cfg->is_info==1,"Test is_info on");
  $cfg->is_info(0);
  ok($cfg->is_info==0,"Test is_info off");
  done_testing();
};
