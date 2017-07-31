#!/usr/bin/env bash
#
# Teem: Tools to process and visualize scientific data and images             .
# Copyright (C) 2013, 2012, 2011, 2010, 2009  University of Chicago
# Copyright (C) 2008, 2007, 2006, 2005  Gordon Kindlmann
# Copyright (C) 2004, 2003, 2002, 2001, 2000, 1999, 1998  University of Utah
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# (LGPL) as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# The terms of redistributing and/or modifying this software also
# include exceptions to the LGPL that facilitate static linking.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; if not, write to Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#

set -o nounset

# The various directories used for this process
# directory used for making source modifications
TEEM_SRC=~/teem
# clean check-out; no other modifications
TEEM_SVN=~/teem-svn
# CMake build directory for teem-svn
TEEM_SVN_BUILD=~/teem-svn-build
# target of "make install" from teem-svn-build
TEEM_SVN_INSTALL=~/teem-svn-install

if [[ ! (-d "$TEEM_SRC" &&
         -d "$TEEM_SVN" &&
         -d "$TEEM_SVN_BUILD" &&
         -d "$TEEM_SVN_INSTALL") ]]; then
  echo "Not all of following directories exist:"
  echo "  TEEM_SRC=$TEEM_SRC"
  echo "  TEEM_SVN=$TEEM_SVN"
  echo "  TEEM_SVN_BUILD=$TEEM_SVN_BUILD"
  echo "  TEEM_SVN_INSTALL=$TEEM_SVN_INSTALL"
  echo "Sorry for the inconvenience.  This script is used mainly by GLK"
  echo "to periodically refresh teem.py, especially prior to releases,"
  echo "and it isn't more generally useful (yet)."
  exit 1
fi

function doo {
  echo "==== $1"
  eval $1
  ret=$?
  if [ $ret != 0 ]; then
    echo "==== ERROR (status $ret)"
    exit $ret
  fi
}

doo "cd $TEEM_SVN"
doo "svn update"
doo "cd $TEEM_SVN_BUILD"
doo "make"
doo "make install"
doo "cd $TEEM_SRC/python/ctypes"
doo "python gen-teem.py ctypeslib-gccxml-0.9 $TEEM_SVN_INSTALL"
dfile=$(mktemp /tmp/svndiff.XXXXXXXXX)
doo "svn diff teem.py | tee $dfile"
if [[ -s $dfile ]]; then
  echo "===="
  echo "==== NOTE: There were new differences; consider \"svn commit teem.py\""
  echo "===="
fi
rm -f $dfile
