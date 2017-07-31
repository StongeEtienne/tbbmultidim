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

#include "teem/nrrd.h"

/*
** Tests:
** nrrdBoundarySpecNew
** nrrdBoundarySpecNix
** nrrdBoundarySpecParse
** nrrdBoundarySpecSprint
** nrrdBoundarySpecCheck
** nrrdHestBoundarySpec
*/

void
checkFailOrDie(const NrrdBoundarySpec *bspec, airArray *mop) {
  static const char me[]="checkFailOrDie";
  char *err;

  if (!nrrdBoundarySpecCheck(bspec)) {
    /* what, no error? */
    fprintf(stderr, "%s: did not get expected check failure!\n", me);
    airMopError(mop);
    exit(1);
  }
  /* good, there's an error */
  err = biffGetDone(NRRD);
  airMopAdd(mop, err, airFree, airMopAlways);
  fprintf(stderr, "(expecting error): %s", err);
  return;
}

void
parseFailOrDie(NrrdBoundarySpec *bspec, const char *str, airArray *mop) {
  static const char me[]="parseFailOrDie";
  char *err;

  if (!nrrdBoundarySpecParse(bspec, str)) {
    /* what, no error? */
    fprintf(stderr, "%s: did not get expected parse fail on \"%s\"\n", me, str);
    airMopError(mop);
    exit(1);
  }
  /* good, there's an error */
  err = biffGetDone(NRRD);
  airMopAdd(mop, err, airFree, airMopAlways);
  fprintf(stderr, "(expecting error): %s", err);
  return;
}

/*
** yes, the equality check should really be on boundary specs,
** not on the strings, but that isn't implemented yet
*/
void
psLoopOrDie(NrrdBoundarySpec *bsp, const char *str,
                     airArray *mop) {
  static const char me[]="psLoopOrDie";
  char *err;
  char buff[AIR_STRLEN_LARGE];

  if (nrrdBoundarySpecParse(bsp, str)
      || nrrdBoundarySpecSprint(buff, bsp)) {
    err = biffGetDone(NRRD);
    airMopAdd(mop, err, airFree, airMopAlways);
    fprintf(stderr, "%s: error: %s", me, err);
    airMopError(mop);
    exit(1);
  }
  if (strcmp(str, buff)) {
    fprintf(stderr, "%s: parse->sprint->\"%s\" != given \"%s\"\n", me,
            buff, str);
    airMopError(mop);
    exit(1);
  }
  fprintf(stderr, "(looped okay): %s\n", str);
  return;
}

static const char *tbspecInfo =
  "for testing handling of boundary specifications";

int
main(int argc, const char *argv[]) {
  /* stock variables */
  const char *me;
  hestOpt *hopt=NULL;
  hestParm *hparm;
  airArray *mop;
  /* variables specific to this program */
  char *err;
  NrrdBoundarySpec *bspec, **bsv;
  unsigned int bsNum, bsIdx;
  char buff[AIR_STRLEN_LARGE];

  me = argv[0];
  mop = airMopNew();
  hparm = hestParmNew();
  airMopAdd(mop, hparm, (airMopper)hestParmFree, airMopAlways);
  hestOptAdd(&hopt, "bs", "bspec0", airTypeOther, 1, -1, &bsv, NULL,
             "bspecs", &bsNum, NULL, nrrdHestBoundarySpec);
  hestParseOrDie(hopt, argc-1, argv+1, hparm, me, tbspecInfo,
                 AIR_TRUE, AIR_TRUE, AIR_TRUE);
  airMopAdd(mop, hopt, (airMopper)hestOptFree, airMopAlways);
  airMopAdd(mop, hopt, (airMopper)hestParseFree, airMopAlways);

  printf("command line options:\n");
  for (bsIdx=0; bsIdx<bsNum; bsIdx++) {
    if (nrrdBoundarySpecSprint(buff, bsv[bsIdx])) {
      airMopAdd(mop, err = biffGetDone(NRRD), airFree, airMopAlways);
      fprintf(stderr, "%s: problem with bsv[%u]:%s", me, bsIdx, err);
      airMopError(mop); return 1;
    }
    printf("%s: bspec[%u] = %s\n", me, bsIdx, buff);
  }

  bspec = nrrdBoundarySpecNew();
  airMopAdd(mop, bspec, (airMopper)nrrdBoundarySpecNix, airMopAlways);

  fprintf(stderr, "%s: bogus value\n", me);
  bspec->boundary = -10;
  checkFailOrDie(bspec, mop);

  fprintf(stderr, "%s: pad with bogus padValue\n", me);
  bspec->boundary = nrrdBoundaryPad;
  bspec->padValue = AIR_POS_INF;
  checkFailOrDie(bspec, mop);

  fprintf(stderr, "%s: checking parse failures\n", me);
  parseFailOrDie(bspec, "bingo", mop);
  parseFailOrDie(bspec, "wrap:10", mop);
  parseFailOrDie(bspec, "pad", mop);
  parseFailOrDie(bspec, "pad:nan", mop);
  parseFailOrDie(bspec, "pad:bob", mop);

  fprintf(stderr, "%s: checking string->bpsec->string\n", me);
  psLoopOrDie(bspec, airEnumStr(nrrdBoundary, nrrdBoundaryBleed), mop);
  psLoopOrDie(bspec, airEnumStr(nrrdBoundary, nrrdBoundaryWrap), mop);
  psLoopOrDie(bspec, airEnumStr(nrrdBoundary, nrrdBoundaryWeight), mop);
  psLoopOrDie(bspec, airEnumStr(nrrdBoundary, nrrdBoundaryMirror), mop);
  psLoopOrDie(bspec, "pad:0", mop);
  psLoopOrDie(bspec, "pad:3", mop);

  airMopOkay(mop);
  return 0;
}
