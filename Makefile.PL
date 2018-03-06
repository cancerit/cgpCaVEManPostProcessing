#!/usr/bin/perl

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


use ExtUtils::MakeMaker;
use File::ShareDir::Install;

install_share dist => 'config';

WriteMakefile(
  NAME          => 'cgpCaVEManPostProcessing',
  VERSION_FROM  => 'lib/Sanger/CGP/CavemanPostProcessing.pm',
  EXE_FILES     => [qw( bin/cgpFlagCaVEMan.pl bin/convertVCFUnmatchedToBed.pl)],
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
                   }
);

package MY;
use File::ShareDir::Install qw( postamble );
