#!/bin/bash
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

SOURCE_BEDTOOLS="https://github.com/arq5x/bedtools2/releases/download/v2.21.0/bedtools-2.21.0.tar.gz"

EXP_BTV="2.21.0"

version_gt () {
  test $(printf '%s\n' $@ | sort -V | head -n 1) == "$1";
}


get_distro () {
  if hash curl 2>/dev/null; then
    curl -sSL -o $1.tar.gz -C - --retry 10 $2
  else
    wget -nv -O $1.tar.gz $2
  fi
  mkdir -p $1
  tar --strip-components 1 -C $1 -zxf $1.tar.gz
}

if [ "$#" -ne "1" ] ; then
  echo "Please provide an installation path  such as /opt/ICGC"
  exit 0
fi

INST_PATH=$1

# get current directory
INIT_DIR=`pwd`

# cleanup inst_path
mkdir -p $INST_PATH/bin
cd $INST_PATH
INST_PATH=`pwd`
cd $INIT_DIR

# make sure that build is self contained
unset PERL5LIB
ARCHNAME=`perl -e 'use Config; print $Config{archname};'`
PERLROOT=$INST_PATH/lib/perl5
export PERL5LIB="$PERLROOT"
export PATH="$INST_PATH/bin:$PATH"

#create a location to build dependencies
SETUP_DIR=$INIT_DIR/install_tmp
mkdir -p $SETUP_DIR

# re-initialise log file
echo > $INIT_DIR/setup.log

# log information about this system
(
    echo '============== System information ===='
    set -x
    lsb_release -a
    uname -a
    sw_vers
    system_profiler
    grep MemTotal /proc/meminfo
    set +x
    echo
) >>$INIT_DIR/setup.log 2>&1

CGPVCF=`perl -le 'eval "require $ARGV[0]" and print $ARGV[0]->VERSION' Sanger::CGP::Vcf`
if [[ "x$CGPVCF" == "x" ]] ; then
  echo "PREREQUISITE: Please install cgpVcf before proceeding:"
  echo "  https://github.com/cancerit/cgpVcf/releases"
  exit 1;
fi

BIODBHTS=`perl -le 'eval "require $ARGV[0]" and print $ARGV[0]->VERSION' Bio::DB::HTS`
if [[ "x$BIODBHTS" == "x" ]] ; then
  echo "PREREQUISITE: Please install Bio::DB::HTS before proceeding"
  exit 1;
fi

CPU=`grep -c ^processor /proc/cpuinfo`
if [ $? -eq 0 ]; then
  if [ "$CPU" -gt "6" ]; then
    CPU=6
  fi
else
  CPU=1
fi
echo "Max compilation CPUs set to $CPU"

set -e

#add bin path for install tests
export PATH="$INST_PATH/bin:$PATH"

BTV=`bedtools --version | cut -d v -f 2`
echo -n "Building bedtools ..."
if [ -e $SETUP_DIR/bedtools.success ]; then
  echo -n " previously installed ...";
else
  if [[ "x$BTV" == "x" ]] ; then
    echo "PREREQUISITE: bedtools version >= $EXP_BTV. Installing... "
    set -x
    cd $SETUP_DIR
    if [ ! -e bedtools ]; then
      get_distro "bedtools2" $SOURCE_BEDTOOLS
    fi
    make -C bedtools2 -j$CPU
    cp bedtools2/bin/* $INST_PATH/bin/.
    set +x
    touch $SETUP_DIR/bedtools.success
  else
    if version_gt $BTV $EXP_BTV; then
      echo "  bedtools version is good ($BTV)"
      touch $SETUP_DIR/bedtools.success
    else
      echo "PREREQUISITE: bedtools version >= $EXP_BTV but found version $BTV). Installing ..."
      set -x
      cd $SETUP_DIR
      if [ ! -e bedtools ]; then
        get_distro "bedtools2" $SOURCE_BEDTOOLS
      fi
      make -C bedtools2 -j$CPU
      cp bedtools2/bin/* $INST_PATH/bin/.
      set +x
      touch $SETUP_DIR/bedtools.success
    fi
  fi
fi

cd $INIT_DIR

echo -n "Installing Perl prerequisites ..."
set -x
perl $INST_PATH/bin/cpanm -v --mirror http://cpan.metacpan.org -l $INST_PATH/ --installdeps .
set +x

echo -n "Installing cgpCaVEManPostProcessing ..."
set -xe
perl Makefile.PL INSTALL_BASE=$INST_PATH
make
make test
make install
set +x

# cleanup all junk
rm -rf $SETUP_DIR

echo
echo
echo "Please add the following to beginning of path:"
echo "  $INST_PATH/bin"
echo "Please add the following to beginning of PERL5LIB:"
echo "  $PERLROOT"
echo

exit 0
