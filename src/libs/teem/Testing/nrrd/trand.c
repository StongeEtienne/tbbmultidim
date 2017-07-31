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
#include <testDataPath.h>

/*
** Tests:
** airSrandMT
** airNormalRand
** nrrdNew
** nrrdAlloc_va
** nrrdHisto
** nrrdHistoDraw
** nrrdSave (to .pgm file)
** nrrdNuke
*/

#define BINS 1000
#define HGHT 1000

int
main(int argc, const char *argv[]) {
  const char *me;
  size_t vi, ii, qvalLen;
  Nrrd *nval, *nhist, *nimg, *nread, *ncorr, *ninmem[3];
  double aa, bb, *val;
  airArray *mop;
  char *name, explain[AIR_STRLEN_LARGE];
#define VALS 0
#define HIST 1
#define IMAG 2
  /* PGM image since this Teem build might not support PNG */
  char *mine[3] = { "vals.nrrd",
                    "histo.nrrd",
                    "histo.pgm" };
  char *corr[3] = { "test/trandvals.nrrd",
                    "test/trandhisto.nrrd",
                    "test/trandhisto.pgm"};
  char *what[3] = { "value",
                    "histogram",
                    "histogram image" };
  int differ, wi;

  AIR_UNUSED(argc);
  me = argv[0];
  mop = airMopNew();

  qvalLen = 10*BINS;
  nrrdAlloc_va(nval=nrrdNew(), nrrdTypeDouble, 1, 4*qvalLen);
  airMopAdd(mop, nval, (airMopper)nrrdNuke, airMopAlways);
  val = AIR_CAST(double*, nval->data);

  nhist=nrrdNew();
  airMopAdd(mop, nhist, (airMopper)nrrdNuke, airMopAlways);
  nimg=nrrdNew();
  airMopAdd(mop, nimg, (airMopper)nrrdNuke, airMopAlways);
  nread = nrrdNew();
  airMopAdd(mop, nread, (airMopper)nrrdNuke, airMopAlways);
  ncorr = nrrdNew();
  airMopAdd(mop, ncorr, (airMopper)nrrdNuke, airMopAlways);

  airSrandMT(999);
  vi = 0;
  /* without first casting to float, the platform-dependent
     differences in the values from airNormalRand() would
     lead to testing errors, e.g.:
     correct (test/trandvals.nrrd) and generated values differ:
     valA[0]=0.36654774192269141 < valB[0]=0.36654774192269146 by 5.55112e-17
     Would be nice to figure out exactly what the origin of that is ... */
  for (ii=0; ii<qvalLen; ii++) {
    airNormalRand(&aa, NULL);
    val[vi++] = AIR_CAST(float, aa);
  }
  for (ii=0; ii<qvalLen; ii++) {
    airNormalRand(NULL, &bb);
    val[vi++] = AIR_CAST(float, bb);
  }
  for (ii=0; ii<qvalLen; ii++) {
    airNormalRand(&aa, &bb);
    val[vi++] = AIR_CAST(float, aa);
    val[vi++] = AIR_CAST(float, bb);
  }

  if (nrrdSave(mine[VALS], nval, NULL)
      || nrrdHisto(nhist, nval, NULL, NULL, BINS, nrrdTypeInt)
      || nrrdSave(mine[HIST], nhist, NULL)
      || nrrdHistoDraw(nimg, nhist, HGHT, AIR_TRUE, 0.0)
      || nrrdSave(mine[IMAG], nimg, NULL)) {
    char *err;
    airMopAdd(mop, err = biffGetDone(NRRD), airFree, airMopAlways);
    fprintf(stderr, "%s: trouble:\n%s", me, err);
    airMopError(mop); return 1;
  }

  ninmem[VALS] = nval;
  ninmem[HIST] = nhist;
  ninmem[IMAG] = nimg;
  for (wi=0; wi<3; wi++) {
    name = testDataPathPrefix(corr[wi]);
    airMopAdd(mop, name, airFree, airMopAlways);
    if (nrrdLoad(ncorr, name, NULL)
        || nrrdLoad(nread, mine[wi], NULL)) {
      char *err;
      airMopAdd(mop, err = biffGetDone(NRRD), airFree, airMopAlways);
      fprintf(stderr, "%s: trouble reading %s:\n%s", me, err, what[wi]);
      airMopError(mop); return 1;
    }
    if (nrrdCompare(ninmem[wi], nread, AIR_FALSE /* onlyData */,
                    0.0 /* epsilon */, &differ, explain)) {
      char *err;
      airMopAdd(mop, err = biffGetDone(NRRD), airFree, airMopAlways);
      fprintf(stderr, "%s: trouble comparing in-mem and from-disk %s:\n%s",
              me, what[wi], err);
      airMopError(mop); return 1;
    }
    if (differ) {
      fprintf(stderr, "%s: in-mem and from-disk (%s) %ss differ: %s\n",
              me, mine[wi], what[wi], explain);
      airMopError(mop); return 1;
    } else {
      printf("%s: good: in-mem and from-disk %ss same\n", me, what[wi]);
    }
    if (nrrdCompare(ncorr, nread, AIR_FALSE /* onlyData */,
                    0.0 /* epsilon */, &differ, explain)) {
      char *err;
      airMopAdd(mop, err = biffGetDone(NRRD), airFree, airMopAlways);
      fprintf(stderr, "%s: trouble comparing correct and generated %s:\n%s",
              me, what[wi], err);
      airMopError(mop); return 1;
    }
    if (differ) {
      fprintf(stderr, "%s: correct (%s) and generated %ss differ: %s\n",
              me, corr[wi], what[wi], explain);
      airMopError(mop); return 1;
    } else {
      printf("%s: good: correct and generated %ss same\n", me, what[wi]);
    }
  }

  airMopOkay(mop);
  return 0;
}
