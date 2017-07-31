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

char *optsigInfo = ("Computes tables of optimal sigmas for hermite-spline "
                    "reconstruction of scale space");

int
main(int argc, const char *argv[]) {
  const char *me;
  hestOpt *hopt;
  hestParm *hparm;
  airArray *mop;

  char *err, *outS;
  double sigma[2], convEps, cutoff;
  int measr[2], tentRecon;
  unsigned int sampleNum[2], dim, measrSampleNum, maxIter, num, ii;
  gageOptimSigContext *osctx;
  double *scalePos, *out, info[512];
  Nrrd *nout;
  NrrdKernelSpec *kss;
  double *plotPos; unsigned int plotPosNum;

  me = argv[0];
  mop = airMopNew();
  hparm = hestParmNew();
  hopt = NULL;
  airMopAdd(mop, hparm, (airMopper)hestParmFree, airMopAlways);
  hestOptAdd(&hopt, "num", "# samp", airTypeUInt, 2, 2, sampleNum, NULL,
             "min and max number of samples to optimize");
  hestOptAdd(&hopt, "dim", "dimension", airTypeUInt, 1, 1, &dim, "3",
             "dimension of image to work with");
  hestOptAdd(&hopt, "sig", "min max", airTypeDouble, 2, 2, sigma, NULL,
             "sigmas to use for bottom and top sample; normal use "
             "should have 0 for bottom scale");
  hestOptAdd(&hopt, "cut", "cut", airTypeDouble, 1, 1, &cutoff, "4",
             "at how many sigmas to cut-off discrete gaussian");
  hestOptAdd(&hopt, "maxi", "max", airTypeUInt, 1, 1, &maxIter, "1000",
             "maximum # iterations");
  hestOptAdd(&hopt, "N", "# samp", airTypeUInt, 1, 1, &measrSampleNum, "300",
             "number of samples in the measurement of error across scales");
  hestOptAdd(&hopt, "eps", "eps", airTypeDouble, 1, 1, &convEps, "0.0001",
             "convergence threshold for optimization");
  hestOptAdd(&hopt, "m", "m1 m2", airTypeEnum, 2, 2, measr, "l2 l2",
             "how to measure error across image, and across scales",
             NULL, nrrdMeasure);
  hestOptAdd(&hopt, "p", "s0 s1", airTypeDouble, 1, -1, &plotPos, "nan nan",
             "hack: don't do optimization; just plot the recon error given "
             "these (two or more) samples along scale.  OR, hackier hack: "
             "if there is only one value given here, do a different "
             "plot: of the recon error within a small window (with the width "
             "given here) in rho, as rho moves through its range",
             &plotPosNum);
  hestOptAdd(&hopt, "kss", "kern", airTypeOther, 1, 1, &kss, "hermite",
             "kernel for gageKernelStack", NULL, NULL, nrrdHestKernelSpec);
  hestOptAdd(&hopt, "tent", NULL, airTypeInt, 0, 0, &tentRecon, NULL,
             "same hack: plot error with tent recon, not hermite");
  hestOptAdd(&hopt, "o", "nout", airTypeString, 1, 1, &outS, NULL,
             "output array");
  hestParseOrDie(hopt, argc-1, argv+1, hparm,
                 me, optsigInfo, AIR_TRUE, AIR_TRUE, AIR_TRUE);
  airMopAdd(mop, hopt, (airMopper)hestOptFree, airMopAlways);
  airMopAdd(mop, hopt, (airMopper)hestParseFree, airMopAlways);

  nout = nrrdNew();
  airMopAdd(mop, nout, AIR_CAST(airMopper, nrrdNuke), airMopAlways);

  osctx = gageOptimSigContextNew(dim, sampleNum[1],
                                 measrSampleNum,
                                 sigma[0], sigma[1],
                                 cutoff);
  if (!osctx) {
    airMopAdd(mop, err = biffGetDone(GAGE), airFree, airMopAlways);
    fprintf(stderr, "%s: trouble:\n%s", me, err);
    airMopError(mop); return 1;
  }
  airMopAdd(mop, osctx, AIR_CAST(airMopper, gageOptimSigContextNix),
            airMopAlways);

  scalePos = AIR_CALLOC(sampleNum[1], double);
  airMopAdd(mop, scalePos, airFree, airMopAlways);

  if (1 == plotPosNum && AIR_EXISTS(plotPos[0])) {
    /* hackity hack: a different kind of plotting requested */
    if (gageOptimSigErrorPlotSliding(osctx, nout,
                                     plotPos[0],
                                     measrSampleNum,
                                     kss, measr[0])) {
      airMopAdd(mop, err = biffGetDone(GAGE), airFree, airMopAlways);
      fprintf(stderr, "%s: trouble:\n%s", me, err);
      airMopError(mop); return 1;
    }
  } else if (AIR_EXISTS(plotPos[0]) && AIR_EXISTS(plotPos[1])) {
    /* hack: plotting requested */
    if (gageOptimSigErrorPlot(osctx, nout,
                              plotPos, plotPosNum,
                              kss,
                              measr[0])) {
      airMopAdd(mop, err = biffGetDone(GAGE), airFree, airMopAlways);
      fprintf(stderr, "%s: trouble:\n%s", me, err);
      airMopError(mop); return 1;
    }
  } else {
    /* do sample position optimization */
    if (nrrdMaybeAlloc_va(nout, nrrdTypeDouble, 2,
                          AIR_CAST(size_t, sampleNum[1]+1),
                          AIR_CAST(size_t, sampleNum[1]+1))) {
      airMopAdd(mop, err = biffGetDone(NRRD), airFree, airMopAlways);
      fprintf(stderr, "%s: trouble allocating output:\n%s", me, err);
      airMopError(mop); return 1;
    }
    out = AIR_CAST(double *, nout->data);
    /* hacky way of saving some of the computation information */
    /* HEY: this is what KVPs are for */
    info[0] = cutoff;
    info[1] = measrSampleNum;
    info[2] = measr[0];
    info[3] = measr[1];
    info[4] = convEps;
    info[5] = maxIter;
    for (ii=0; ii<sampleNum[1]+1; ii++) {
      out[ii] = info[ii];
    }
    for (num=sampleNum[0]; num<=sampleNum[1]; num++) {
      printf("\n%s: ======= optimizing %u/%u samples (sigma %g--%g) \n\n",
             me, num, sampleNum[1]+1, sigma[0], sigma[1]);
      if (gageOptimSigCalculate(osctx, scalePos, num,
                                kss, measr[0], measr[1],
                                maxIter, convEps)) {
        airMopAdd(mop, err = biffGetDone(GAGE), airFree, airMopAlways);
        fprintf(stderr, "%s: trouble:\n%s", me, err);
        airMopError(mop); return 1;
      }
      for (ii=0; ii<num; ii++) {
        out[ii + (sampleNum[1]+1)*num] = scalePos[ii];
      }
      out[sampleNum[1] + (sampleNum[1]+1)*num] = osctx->finalErr;
    }
  }
  if (nrrdSave(outS, nout, NULL)) {
    airMopAdd(mop, err = biffGetDone(NRRD), airFree, airMopAlways);
    fprintf(stderr, "%s: trouble saving output:\n%s\n", me, err);
    airMopError(mop);
    return 1;
  }

  airMopOkay(mop);
  exit(0);
}
