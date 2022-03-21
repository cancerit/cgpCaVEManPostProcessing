#!/usr/bin/perl
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

use ExtUtils::MakeMaker;
use File::ShareDir::Install;

install_share dist => 'config';

WriteMakefile(
  NAME          => 'cgpCaVEManPostProcessing',
  VERSION_FROM  => 'lib/Sanger/CGP/CavemanPostProcessor.pm',
  EXE_FILES     => [qw( bin/cgpFlagCaVEMan.pl bin/convertVCFUnmatchedToBed.pl bin/cavemanPostProcessing_ini_to_yaml.pl )],
  PREREQ_PM     => {
                     'Const::Fast' => 0.014,
                     'Attribute::Abstract' => 1.02,
                     'Bio::DB::HTS' => 1.11,
                     'Test::Fatal' => 0.013,
                     'Devel::Cover' => 1.09,
                     'Pod::Coverage' => 0.23,
                     'Config::IniFiles' => 2.83,
                     'File::ShareDir' => 1.102,
                     'LWP::Simple' => 6.00,
                     'Set::IntervalTree' => 0.12,
                     'PerlIO::gzip' => 0.20,
                     'List::Util'=> 1.56,
                     'YAML::XS' => 0.83,
                     'Test::MockObject' => 1.20,
                   }
);

package MY;
use File::ShareDir::Install qw( postamble );
