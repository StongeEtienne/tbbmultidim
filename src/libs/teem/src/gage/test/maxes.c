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

#include "../gage.h"

void
maxes1(Nrrd *nout, const Nrrd *nin) {
  unsigned int sx, xi, xd;
  int xo, ismax;
  double val, (*lup)(const void *, size_t),
    (*ins)(void *, size_t, double);

  lup = nrrdDLookup[nin->type];
  ins = nrrdDInsert[nin->type];
  sx = AIR_CAST(unsigned int, nin->axis[0].size);
  for (xi=0; xi<sx; xi++) {
    ismax = AIR_TRUE;
    val = lup(nin->data, xi);
    for (xo=-1; xo<=1; xo+=2) {
      xd = !xi && xo < 0 ? sx-1 : AIR_MOD(xi+xo, sx);
      ismax &= lup(nin->data, xd) < val;
      if (!ismax) break;
    }
    ins(nout->data, xi, ismax);
  }
  return;
}

void
maxes2(Nrrd *nout, const Nrrd *nin) {
  unsigned int sx, sy, xi, yi, xd, yd;
  int xo, yo, ismax;
  double val, (*lup)(const void *, size_t),
    (*ins)(void *, size_t, double);

  lup = nrrdDLookup[nin->type];
  ins = nrrdDInsert[nin->type];
  sx = AIR_CAST(unsigned int, nin->axis[0].size);
  sy = AIR_CAST(unsigned int, nin->axis[1].size);
  for (yi=0; yi<sy; yi++) {
    for (xi=0; xi<sx; xi++) {
      ismax = AIR_TRUE;
      val = lup(nin->data, xi + sx*yi);
      for (yo=-1; yo<=1; yo++) {
        yd = !yi && yo < 0 ? sy-1 : AIR_MOD(yi+yo, sy);
        for (xo=-1; xo<=1; xo++) {
          if (!xo && !yo) continue;
          xd = !xi && xo < 0 ? sx-1 : AIR_MOD(xi+xo, sx);
          ismax &= lup(nin->data, xd + sx*yd) < val;
        }
      }
      ins(nout->data, xi + sx*yi, ismax);
    }
  }
  return;
}

void
maxes3(Nrrd *nout, const Nrrd *nin) {
  unsigned int sx, sy, sz, xi, yi, zi, xd, yd, zd;
  int xo, yo, zo, ismax;
  double val, (*lup)(const void *, size_t),
    (*ins)(void *, size_t, double);

  lup = nrrdDLookup[nin->type];
  ins = nrrdDInsert[nin->type];
  sx = AIR_CAST(unsigned int, nin->axis[0].size);
  sy = AIR_CAST(unsigned int, nin->axis[1].size);
  sz = AIR_CAST(unsigned int, nin->axis[2].size);
  for (zi=0; zi<sz; zi++) {
    for (yi=0; yi<sy; yi++) {
      for (xi=0; xi<sx; xi++) {
        ismax = AIR_TRUE;
        val = lup(nin->data, xi + sx*(yi + sy*zi));
        for (zo=-1; zo<=1; zo++) {
          zd = !zi && zo < 0 ? sz-1 : AIR_MOD(zi+zo, sz);
          for (yo=-1; yo<=1; yo++) {
            yd = !yi && yo < 0 ? sy-1 : AIR_MOD(yi+yo, sy);
            for (xo=-1; xo<=1; xo++) {
              if (!xo && !yo && !zo) continue;
              xd = !xi && xo < 0 ? sx-1 : AIR_MOD(xi+xo, sx);
              ismax &= lup(nin->data, xd + sx*(yd + sy*zd)) < val;
            }
          }
        }
        ins(nout->data, xi + sx*(yi + sy*zi), ismax);
      }
    }
  }
  return;
}

char *maxesInfo = ("sets maxima to 1, else 0");

int
main(int argc, const char *argv[]) {
  const char *me;
  char *outS;
  hestOpt *hopt;
  hestParm *hparm;
  airArray *mop;

  char *err;
  Nrrd *nin, *nout;

  me = argv[0];
  mop = airMopNew();
  hparm = hestParmNew();
  hopt = NULL;
  airMopAdd(mop, hparm, (airMopper)hestParmFree, airMopAlways);
  hestOptAdd(&hopt, "i", "nin", airTypeOther, 1, 1, &nin, NULL,
             "input image", NULL, NULL, nrrdHestNrrd);
  hestOptAdd(&hopt, "o", "filename", airTypeString, 1, 1, &outS, NULL,
             "output volume");
  hestParseOrDie(hopt, argc-1, argv+1, hparm,
                 me, maxesInfo, AIR_TRUE, AIR_TRUE, AIR_TRUE);
  airMopAdd(mop, hopt, (airMopper)hestOptFree, airMopAlways);
  airMopAdd(mop, hopt, (airMopper)hestParseFree, airMopAlways);

  if (!AIR_IN_CL(1, nin->dim, 3)) {
    airMopAdd(mop, err = biffGetDone(NRRD), airFree, airMopAlways);
    fprintf(stderr, "%s: nin->dim %u not 1, 2, or 3", me, nin->dim);
    airMopError(mop); return 1;
  }

  nout = nrrdNew();
  airMopAdd(mop, nout, (airMopper)nrrdNuke, airMopAlways);
  if (nrrdCopy(nout, nin)) {
    airMopAdd(mop, err = biffGetDone(NRRD), airFree, airMopAlways);
    fprintf(stderr, "%s: couldn't allocate output:\n%s", me, err);
    airMopError(mop); return 1;
  }

  if (1 == nin->dim) {
    maxes1(nout, nin);
  } else if (2 == nin->dim) {
    maxes2(nout, nin);
  } else {
    maxes3(nout, nin);
  }

  if (nrrdSave(outS, nout, NULL)) {
    airMopAdd(mop, err = biffGetDone(NRRD), airFree, airMopAlways);
    fprintf(stderr, "%s: couldn't save output:\n%s", me, err);
    airMopError(mop); return 1;
  }

  airMopOkay(mop);
  exit(0);
}
