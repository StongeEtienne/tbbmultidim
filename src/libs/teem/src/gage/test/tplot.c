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

char *tplotInfo = ("for plotting tau");

int
main(int argc, const char *argv[]) {
  const char *me;
  hestOpt *hopt;
  hestParm *hparm;
  airArray *mop;

  double dommm[2];
  unsigned int num;
  char *err;
  char *outS;

  Nrrd *nout;
  double *out, dd;
  unsigned int ii;
  int tee;

  me = argv[0];
  mop = airMopNew();
  hparm = hestParmNew();
  hopt = NULL;
  airMopAdd(mop, hparm, (airMopper)hestParmFree, airMopAlways);
  hestOptAdd(&hopt, "tee", NULL, airTypeInt, 0, 0, &tee, NULL,
             "domain of plot is t, not sigma");
  hestOptAdd(&hopt, "mm", "min max", airTypeDouble, 2, 2, dommm, NULL,
             "range of domain to plot");
  hestOptAdd(&hopt, "n", "# samp", airTypeUInt, 1, 1, &num, NULL,
             "number of samples");
  hestOptAdd(&hopt, "o", "filename", airTypeString, 1, 1, &outS, NULL,
             "output volume");
  hestParseOrDie(hopt, argc-1, argv+1, hparm,
                 me, tplotInfo, AIR_TRUE, AIR_TRUE, AIR_TRUE);
  airMopAdd(mop, hopt, (airMopper)hestOptFree, airMopAlways);
  airMopAdd(mop, hopt, (airMopper)hestParseFree, airMopAlways);

  nout = nrrdNew();
  airMopAdd(mop, nout, (airMopper)nrrdNuke, airMopAlways);
  if (nrrdMaybeAlloc_va(nout, nrrdTypeDouble, 1,
                        AIR_CAST(size_t, num))) {
    airMopAdd(mop, err = biffGetDone(NRRD), airFree, airMopAlways);
    fprintf(stderr, "%s: couldn't allocate output:\n%s", me, err);
    airMopError(mop); return 1;
  }
  out = AIR_CAST(double *, nout->data);
  nout->axis[0].min = dommm[0];
  nout->axis[0].max = dommm[1];
  nout->axis[0].center = nrrdCenterNode;
  for (ii=0; ii<num; ii++) {
    dd = AIR_AFFINE(0, ii, num-1, dommm[0], dommm[1]);
    out[ii] = tee ? gageTauOfTee(dd) : gageTauOfSig(dd);
  }

  if (nrrdSave(outS, nout, NULL)) {
    airMopAdd(mop, err = biffGetDone(NRRD), airFree, airMopAlways);
    fprintf(stderr, "%s: couldn't save output:\n%s", me, err);
    airMopError(mop); return 1;
  }

  airMopOkay(mop);
  exit(0);
}
