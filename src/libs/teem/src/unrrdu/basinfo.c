/*
  Teem: Tools to process and visualize scientific data and images             .
  Copyright (C) 2013, 2012, 2011, 2010, 2009  University of Chicago
  Copyright (C) 2008, 2007, 2006, 2005  Gordon Kindlmann
  Copyright (C) 2004, 2003, 2002, 2001, 2000, 1999, 1998  University of Utah

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.
  The terms of redistributing and/or modifying this software also
  include exceptions to the LGPL that facilitate static linking.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with this library; if not, write to Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include "unrrdu.h"
#include "privateUnrrdu.h"

#define INFO "Modify whole-array attributes (not per-axis)"
static const char *_unrrdu_basinfoInfoL =
(INFO
 ", which is called \"basic info\" in Nrrd terminology. "
 "The only attributes which are set are those for which command-line "
 "options are given.\n "
 "* Uses no particular function; just sets fields in NrrdAxisInfo");

int
unrrdu_basinfoMain(int argc, const char **argv, const char *me,
                   hestParm *hparm) {
  /* these are stock for unrrdu */
  hestOpt *opt = NULL;
  airArray *mop;
  int pret;
  char *err;
  /* these are stock for things using the usual -i and -o */
  char *out;
  Nrrd *nin, *nout;
  /* these are specific to this command */
  NrrdIoState *nio;
  char *spcStr, *_origStr, *origStr;
  int space;
  unsigned int spaceDim;

  hestOptAdd(&opt, "spc,space", "space", airTypeString, 1, 1, &spcStr, "",
             "identify the space (e.g. \"RAS\", \"LPS\") in which the array "
             "conceptually lives, from the nrrdSpace airEnum, which in turn "
             "determines the dimension of the space.  Or, use an integer>0 to"
             "give the dimension of a space that nrrdSpace doesn't know about. "
             "By default (not using this option), the enclosing space is "
             "set as unknown.");
  hestOptAdd(&opt, "orig,origin", "origin", airTypeString, 1, 1, &_origStr, "",
             "(NOTE: must quote vector) the origin in space of the array: "
             "the location of the center "
             "of the first sample, of the form \"(x,y,z)\" (or however "
             "many coefficients are needed for the chosen space). Quoting the "
             "vector is needed to stop interpretation from the shell");
  OPT_ADD_NIN(nin, "input nrrd");
  OPT_ADD_NOUT(out, "output nrrd");

  mop = airMopNew();
  airMopAdd(mop, opt, (airMopper)hestOptFree, airMopAlways);
  nio = nrrdIoStateNew();
  airMopAdd(mop, nio, (airMopper)nrrdIoStateNix, airMopAlways);

  USAGE(_unrrdu_basinfoInfoL);
  PARSE();
  airMopAdd(mop, opt, (airMopper)hestParseFree, airMopAlways);

  nout = nrrdNew();
  airMopAdd(mop, nout, (airMopper)nrrdNuke, airMopAlways);
  if (nrrdCopy(nout, nin)) {
    airMopAdd(mop, err = biffGetDone(NRRD), airFree, airMopAlways);
    fprintf(stderr, "%s: error copying input:\n%s", me, err);
    airMopError(mop);
    return 1;
  }

  /* HEY: copy and paste from unrrdu/make.c */
  if (airStrlen(spcStr)) {
    space = airEnumVal(nrrdSpace, spcStr);
    if (!space) {
      /* couldn't parse it as space, perhaps its a uint */
      if (1 != sscanf(spcStr, "%u", &spaceDim)) {
        fprintf(stderr, "%s: couldn't parse \"%s\" as a nrrdSpace "
                "or as a uint", me, spcStr);
        airMopError(mop); return 1;
      }
      /* else we did parse it as a uint */
      nout->space = nrrdSpaceUnknown;
      nout->spaceDim = spaceDim;
    } else {
      /* we did parse a known space */
      nrrdSpaceSet(nout, space);
    }
  }

  /* HEY: copy and paste from unrrdu/make.c */
  if (airStrlen(_origStr)) {
    /* why this is necessary is a bit confusing to me, both the check for
       enclosing quotes, and the need to use to a separate variable (isn't
       hest doing memory management of addresses, not variables?) */
    if ('\"' == _origStr[0] && '\"' == _origStr[strlen(_origStr)-1]) {
      _origStr[strlen(_origStr)-1] = 0;
      origStr = _origStr + 1;
    } else {
      origStr = _origStr;
    }
    /* same hack about using NrrdIoState->line as basis for parsing */
    nio->line = origStr;
    nio->pos = 0;
    if (nrrdFieldInfoParse[nrrdField_space_origin](NULL, nout,
                                                   nio, AIR_TRUE)) {
      airMopAdd(mop, err = biffGetDone(NRRD), airFree, airMopAlways);
      fprintf(stderr, "%s: trouble with origin \"%s\":\n%s",
              me, origStr, err);
      nio->line = NULL; airMopError(mop); return 1;
    }
    nio->line = NULL;
  }

  SAVE(out, nout, NULL);

  airMopOkay(mop);
  return 0;
}

UNRRDU_CMD(basinfo, INFO);
