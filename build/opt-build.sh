#! /bin/bash
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
set -xe

if [[ -z "${TMPDIR}" ]]; then
  TMPDIR=/tmp
fi

set -u

if [ "$#" -lt "1" ] ; then
  echo "Please provide an installation path such as /opt/ICGC"
  exit 1
fi

# get path to this script
SCRIPT_PATH=`dirname $0`;
SCRIPT_PATH=`(cd $SCRIPT_PATH && pwd)`

# get the location to install to
INST_PATH=$1
mkdir -p $1
INST_PATH=`(cd $1 && pwd)`
echo $INST_PATH

# get current directory
INIT_DIR=`pwd`

CPU=`grep -c ^processor /proc/cpuinfo`
if [ $? -eq 0 ]; then
  if [ "$CPU" -gt "6" ]; then
    CPU=6
  fi
else
  CPU=1
fi
echo "Max compilation CPUs set to $CPU"

SETUP_DIR=$INIT_DIR/install_tmp
mkdir -p $SETUP_DIR/distro # don't delete the actual distro directory until the very end
mkdir -p $INST_PATH/bin
cd $SETUP_DIR

# make sure tools installed can see the install loc of libraries
set +u
export LD_LIBRARY_PATH=`echo $INST_PATH/lib:$LD_LIBRARY_PATH | perl -pe 's/:\$//;'`
export PATH=`echo $INST_PATH/bin:$PATH | perl -pe 's/:\$//;'`
export MANPATH=`echo $INST_PATH/man:$INST_PATH/share/man:$MANPATH | perl -pe 's/:\$//;'`
export PERL5LIB=`echo $INST_PATH/lib/perl5:$PERL5LIB | perl -pe 's/:\$//;'`
export C_INCLUDE_PATH=`echo $INST_PATH/include/ | perl -pe 's/:\$//;'`
set -u

## vcftools
if [ ! -e $SETUP_DIR/vcftools.success ]; then
  curl -sSL --retry 10 https://github.com/vcftools/vcftools/releases/download/v${VER_VCFTOOLS}/vcftools-${VER_VCFTOOLS}.tar.gz > distro.tar.gz
  rm -rf distro/*
  tar --strip-components 2 -C distro -xzf distro.tar.gz
  cd distro
  ./configure --prefix=$INST_PATH --with-pmdir=lib/perl5
  make -j$CPU
  make install
  cd $SETUP_DIR
  rm -rf distro.* distro/*
  touch $SETUP_DIR/vcftools.success
fi

### cgpVcf
if [ ! -e $SETUP_DIR/cgpVcf.success ]; then
  curl -sSL --retry 10 https://github.com/cancerit/cgpVcf/archive/${VER_CGPVCF}.tar.gz > distro.tar.gz
  rm -rf distro/*
  tar --strip-components 1 -C distro -xzf distro.tar.gz
  cd distro
  cpanm --no-interactive --notest --mirror http://cpan.metacpan.org --notest -l $INST_PATH --installdeps .
  cpanm -v --no-interactive --mirror http://cpan.metacpan.org -l $INST_PATH .
  cd $SETUP_DIR
  rm -rf distro.* distro/*
  touch $SETUP_DIR/cgpVcf.success
fi

### bedtools for cgpCaVEManPostProcessing
if [ ! -e $SETUP_DIR/bedtools.success ]; then
  curl -sSL --retry 10 https://github.com/arq5x/bedtools2/releases/download/v${VER_BEDTOOLS}/bedtools-${VER_BEDTOOLS}.tar.gz > distro.tar.gz
  rm -rf distro/*
  tar --strip-components 1 -C distro -xzf distro.tar.gz
  make -C distro -j$CPU
  cp distro/bin/* $INST_PATH/bin/.
  cd $SETUP_DIR
  rm -rf distro.* distro/*
  touch $SETUP_DIR/bedtools.success
fi

### cgpCaVEManPostProcessing
if [ ! -e $SETUP_DIR/cgpCaVEManPostProcessing.success ]; then
  cpanm --no-interactive --notest --mirror http://cpan.metacpan.org --notest -l $INST_PATH File::ShareDir::Install
  curl -sSL --retry 10 https://github.com/cancerit/cgpCaVEManPostProcessing/archive/${VER_CGPCAVEPOSTPROC}.tar.gz > distro.tar.gz
  rm -rf distro/*
  tar --strip-components 1 -C distro -xzf distro.tar.gz
  cd distro
  cpanm --no-interactive --notest --mirror http://cpan.metacpan.org --notest -l $INST_PATH --installdeps .
  cpanm -v --no-interactive --mirror http://cpan.metacpan.org -l $INST_PATH .
  cd $SETUP_DIR
  rm -rf distro.* distro/*
  touch $SETUP_DIR/cgpCaVEManPostProcessing.success
fi
