#!/bin/bash

##########LICENCE##########
# Copyright (c) 2014 Genome Research Ltd.
#
# Author: Cancer Genome Project cgpit@sanger.ac.uk
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

SOURCE_BEDTOOLS="https://github.com/arq5x/bedtools2/releases/download/v2.21.0/bedtools-2.21.0.tar.gz"

done_message () {
    if [ $? -eq 0 ]; then
        echo " done."
        if [ "x$1" != "x" ]; then
            echo $1
        fi
    else
        echo " failed.  See setup.log file for error messages." $2
        echo "    Please check INSTALL file for items that should be installed by a package manager"
        exit 1
    fi
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
PERLARCH=$PERLROOT/$ARCHNAME
export PERL5LIB="$PERLROOT:$PERLARCH"

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

BIODBSAM=`perl -le 'eval "require $ARGV[0]" and print $ARGV[0]->VERSION' Bio::DB::Sam`
if [[ "x$BIODBSAM" == "x" ]] ; then
  echo "PREREQUISITE: Please install Bio::DB::Sam before proceeding"
  exit 1;
fi

#add bin path for install tests
export PATH="$INST_PATH/bin:$PATH"

echo -n "Building bedtools ..."
if [ -e $SETUP_DIR/bedtools.success ]; then
  echo -n " previously installed ...";
else
  cd $SETUP_DIR
  (
  set -x
  if [ ! -e bedtools ]; then
    get_distro "bedtools2" $SOURCE_BEDTOOLS
  fi
  make -C bedtools2 -j$CPU
  cp bedtools2/bin/* $INST_PATH/bin/.
  touch $SETUP_DIR/bedtools.success
  )>>$INIT_DIR/setup.log 2>&1
fi
done_message "" "Failed to build bedtools."

cd $INIT_DIR

echo -n "Installing Perl prerequisites ..."
if ! ( perl -MExtUtils::MakeMaker -e 1 >/dev/null 2>&1); then
    echo
    echo "WARNING: Your Perl installation does not seem to include a complete set of core modules.  Attempting to cope with this, but if installation fails please make sure that at least ExtUtils::MakeMaker is installed.  For most users, the best way to do this is to use your system's package manager: apt, yum, fink, homebrew, or similar."
fi
(
  set -x
  $INIT_DIR/bin/cpanm -v --mirror http://cpan.metacpan.org -notest -l $INST_PATH/ --installdeps . < /dev/null
  set +x
) >>$INIT_DIR/setup.log 2>&1
done_message "" "Failed during installation of core dependencies."

echo -n "Installing cgpCaVEManPostProcessing ..."
(
  perl Makefile.PL INSTALL_BASE=$INST_PATH
  make
  make test
  make install
) >>$INIT_DIR/setup.log 2>&1
done_message "" "cgpCaVEManPostProcessing install failed."

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
