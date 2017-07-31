/*
  Teem: Tools to process and visualize scientific data and images             .
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


#include "../dye.h"

/*
** hist[0]: hue vs sat
** hist[1]: hue vs val
*/
void
imageProc(Nrrd *nhproj[3], Nrrd *nhist[2], unsigned int sH,
          float *rgb, unsigned int size0, unsigned int sXY,
          unsigned int overSampleNum, float overSampleScale) {
  unsigned int xyi, hi, si, vi, oi;
  float rr, gg, bb, hh, ss, vv, *hist[2];
  double rndA, rndB;

  nrrdZeroSet(nhist[0]);
  nrrdZeroSet(nhist[1]);
  hist[0] = AIR_CAST(float *, nhist[0]->data);
  hist[1] = AIR_CAST(float *, nhist[1]->data);
  for (xyi=0; xyi<sXY; xyi++) {
    rr = AIR_CLAMP(0, rgb[0], 255);
    gg = AIR_CLAMP(0, rgb[1], 255);
    bb = AIR_CLAMP(0, rgb[2], 255);
    rr = AIR_AFFINE(-1, rr, 256, 0, 1);
    gg = AIR_AFFINE(-1, gg, 256, 0, 1);
    bb = AIR_AFFINE(-1, bb, 256, 0, 1);
    dyeRGBtoHSV(&hh, &ss, &vv, rr, gg, bb);
    si = airIndexClamp(0, ss, 1, sH);
    vi = airIndexClamp(0, vv, 1, sH);

#define UPDATE_HIST(rnd)                                                \
    hi = airIndexClamp(0, hh + overSampleScale*(1-ss)*(rnd), 1, sH);    \
    hist[0][hi + sH*si] += 1.0/overSampleNum;                           \
    hist[1][hi + sH*vi] += 1.0/overSampleNum

    if (overSampleNum % 2 == 1) {
      airNormalRand(&rndA, NULL);
      UPDATE_HIST(rndA);
      overSampleNum -= 1;
    }
    for (oi=0; oi<overSampleNum; oi+=2) {
      airNormalRand(&rndA, &rndB);
      UPDATE_HIST(rndA);
      UPDATE_HIST(rndB);
    }
    rgb += size0;
  }
  nrrdProject(nhproj[0], nhist[0], 1, nrrdMeasureHistoMean, nrrdTypeFloat);
  nrrdProject(nhproj[1], nhist[1], 1, nrrdMeasureHistoMean, nrrdTypeFloat);
  nrrdProject(nhproj[2], nhist[1], 1, nrrdMeasureSum, nrrdTypeFloat);
}

const char *mchistInfo = ("Makes a color histogram of frames. ");

int
main(int argc, const char *argv[]) {
  const char *me;
  hestOpt *hopt;
  hestParm *hparm;
  airArray *mop;

  char **ninStr, *err, *outS, doneStr[13];
  Nrrd *nin0, *nin, *nrgb, *nout, *nhist[2], *npreout, *nhproj[3];
  float *rgb;
  float *out, *preout, *hist[2], maxSum,
    upSample, overSampleScale;
  unsigned int size0, sX, sY, sH, ninLen, ti, overSampleNum;
  NrrdResampleContext *rsmc;
  NrrdKernelSpec *ksp;

  me = argv[0];
  mop = airMopNew();
  hopt = NULL;
  hparm = hestParmNew();
  airMopAdd(mop, hparm, (airMopper)hestParmFree, airMopAlways);
  hparm->respFileEnable = AIR_TRUE;
  hestOptAdd(&hopt, "i", "images", airTypeString, 1, -1, &ninStr, NULL,
             "input image sequence", &ninLen, NULL, NULL);
  hestOptAdd(&hopt, "sh", "histo size", airTypeUInt, 1, 1, &sH, "500",
             "histogram size");
  hestOptAdd(&hopt, "k", "kern", airTypeOther, 1, 1, &ksp,
             "tent", "kernel for upsampling images",
             NULL, NULL, nrrdHestKernelSpec);
  hestOptAdd(&hopt, "us", "upsampling", airTypeFloat, 1, 1, &upSample,
             "1", "amount of upsampling of image");
  hestOptAdd(&hopt, "osn", "# oversmp", airTypeUInt, 1, 1, &overSampleNum,
             "1", "number of sample per (upsampled) pixel");
  hestOptAdd(&hopt, "osc", "scaling", airTypeFloat, 1, 1, &overSampleScale,
             "1", "scaling with oversampling");
  hestOptAdd(&hopt, "ms", "max sum", airTypeFloat, 1, 1, &maxSum,
             "10", "per-hue histogram summation is non-linearly and "
             "asymptotically clamped to this maximum");
  hestOptAdd(&hopt, "o", "nout", airTypeString, 1, 1, &outS, "-",
             "output filename", NULL);

  hestParseOrDie(hopt, argc-1, argv+1, hparm,
                 me, mchistInfo, AIR_TRUE, AIR_TRUE, AIR_TRUE);
  airMopAdd(mop, hopt, (airMopper)hestOptFree, airMopAlways);
  airMopAdd(mop, hopt, (airMopper)hestParseFree, airMopAlways);

  if (0 == overSampleNum) {
    fprintf(stderr, "%s: overSampleNum must be > 0\n", me);
    airMopError(mop);
    return 1;
  }
  nin0 = nrrdNew();
  airMopAdd(mop, nin0, (airMopper)nrrdNuke, airMopAlways);
  if (nrrdLoad(nin0, ninStr[0], NULL)) {
    airMopAdd(mop, err=biffGetDone(NRRD), airFree, airMopAlways);
    fprintf(stderr, "%s: couldn't load first image:\n%s", me, err);
    airMopError(mop);
    return 1;
  }
  if (!( (3 == nin0->axis[0].size || 4 == nin0->axis[0].size)
         && 3 == nin0->dim
         && nrrdTypeUChar == nin0->type )) {
    fprintf(stderr, "%s: 1st image not 3D (3-or-4)-by-X-by-Y %s array "
            "(got %u-D %s array)\n", me,
            airEnumStr(nrrdType, nrrdTypeUChar),
            nin0->dim,
            airEnumStr(nrrdType, nin0->type));
    airMopError(mop);
    return 1;
  }
  rsmc = nrrdResampleContextNew();
  airMopAdd(mop, rsmc, (airMopper)nrrdResampleContextNix, airMopAlways);
  size0 = AIR_CAST(unsigned int, nin0->axis[0].size);
  sX = AIR_CAST(unsigned int, upSample*nin0->axis[1].size);
  sY = AIR_CAST(unsigned int, upSample*nin0->axis[2].size);
  nrgb = nrrdNew();
  airMopAdd(mop, nrgb, (airMopper)nrrdNuke, airMopAlways);
  if (nrrdResampleDefaultCenterSet(rsmc, nrrdCenterCell)
      || nrrdResampleInputSet(rsmc, nin0)
      || nrrdResampleKernelSet(rsmc, 1, ksp->kernel, ksp->parm)
      || nrrdResampleKernelSet(rsmc, 2, ksp->kernel, ksp->parm)
      || nrrdResampleSamplesSet(rsmc, 1, sX)
      || nrrdResampleSamplesSet(rsmc, 2, sY)
      || nrrdResampleRangeFullSet(rsmc, 1)
      || nrrdResampleRangeFullSet(rsmc, 2)
      || nrrdResampleTypeOutSet(rsmc, nrrdTypeFloat)
      || nrrdResampleRenormalizeSet(rsmc, AIR_TRUE)
      || nrrdResampleExecute(rsmc, nrgb)) {
    airMopAdd(mop, err=biffGetDone(NRRD), airFree, airMopAlways);
    fprintf(stderr, "%s: error resampling slice:\n%s", me, err);
    airMopError(mop);
    return 1;
  }

  nhist[0] = nrrdNew();
  airMopAdd(mop, nhist[0], (airMopper)nrrdNuke, airMopAlways);
  nhist[1] = nrrdNew();
  airMopAdd(mop, nhist[1], (airMopper)nrrdNuke, airMopAlways);
  if (nrrdMaybeAlloc_va(nhist[0], nrrdTypeFloat, 2,
                        AIR_CAST(size_t, sH),
                        AIR_CAST(size_t, sH))
      || nrrdMaybeAlloc_va(nhist[1], nrrdTypeFloat, 2,
                           AIR_CAST(size_t, sH),
                           AIR_CAST(size_t, sH))) {
    airMopAdd(mop, err=biffGetDone(NRRD), airFree, airMopAlways);
    fprintf(stderr, "%s: error allocating histos:\n%s", me, err);
    airMopError(mop);
    return 1;
  }
  nhist[0]->axis[0].min = nhist[0]->axis[1].min = 0.0;
  nhist[0]->axis[0].max = nhist[0]->axis[1].max = 1.0;
  nhist[1]->axis[0].min = nhist[1]->axis[1].min = 0.0;
  nhist[1]->axis[0].max = nhist[1]->axis[1].max = 1.0;
  nhproj[0] = nrrdNew();
  airMopAdd(mop, nhproj[0], (airMopper)nrrdNix, airMopAlways);
  nhproj[1] = nrrdNew();
  airMopAdd(mop, nhproj[1], (airMopper)nrrdNix, airMopAlways);
  nhproj[2] = nrrdNew();
  airMopAdd(mop, nhproj[2], (airMopper)nrrdNix, airMopAlways);

  printf("working ...       ");
  hist[0] = AIR_CAST(float *, nhist[0]->data);
  hist[1] = AIR_CAST(float *, nhist[1]->data);
  nin = nrrdNew();
  airMopAdd(mop, nin, (airMopper)nrrdNuke, airMopAlways);
  npreout = nrrdNew();
  airMopAdd(mop, npreout, (airMopper)nrrdNuke, airMopAlways);
  if (nrrdMaybeAlloc_va(npreout, nrrdTypeFloat, 3,
                        AIR_CAST(size_t, 3),
                        AIR_CAST(size_t, sH),
                        AIR_CAST(size_t, ninLen))) {
    airMopAdd(mop, err=biffGetDone(NRRD), airFree, airMopAlways);
    fprintf(stderr, "%s: error allocating pre-output:\n%s", me, err);
    airMopError(mop);
    return 1;
  }
  preout = AIR_CAST(float *, npreout->data);
  for (ti=0; ti<ninLen; ti++) {
    printf("%s", airDoneStr(0, ti, ninLen, doneStr)); fflush(stdout);
    if (nrrdLoad(nin, ninStr[ti], NULL)) {
      airMopAdd(mop, err=biffGetDone(NRRD), airFree, airMopAlways);
      fprintf(stderr, "%s: couldn't load image[%u]:\n%s", me, ti, err);
      airMopError(mop);
      return 1;
    }
    if (!nrrdSameSize(nin0, nin, AIR_TRUE)) {
      airMopAdd(mop, err=biffGetDone(NRRD), airFree, airMopAlways);
      fprintf(stderr, "%s: nin[%u] not like nin[0]:\n%s", me, ti, err);
      airMopError(mop);
      return 1;
    }
    if (nrrdResampleInputSet(rsmc, nin)
        || nrrdResampleExecute(rsmc, nrgb)) {
      airMopAdd(mop, err=biffGetDone(NRRD), airFree, airMopAlways);
      fprintf(stderr, "%s: trouble resampling nin[%u]:\n%s", me, ti, err);
      airMopError(mop);
      return 1;
    }
    if (nrrdWrap_va(nhproj[0], preout + 0*sH, nrrdTypeFloat, 1, AIR_CAST(size_t, sH)) ||
        nrrdWrap_va(nhproj[1], preout + 1*sH, nrrdTypeFloat, 1, AIR_CAST(size_t, sH)) ||
        nrrdWrap_va(nhproj[2], preout + 2*sH, nrrdTypeFloat, 1, AIR_CAST(size_t, sH))) {
      airMopAdd(mop, err=biffGetDone(NRRD), airFree, airMopAlways);
      fprintf(stderr, "%s: trouble wrapping output[%u]:\n%s", me, ti, err);
      airMopError(mop);
      return 1;
    }
    rgb = AIR_CAST(float *, nrgb->data);
    imageProc(nhproj, nhist, sH,
              rgb, size0, sX*sY,
              overSampleNum, overSampleScale);
    preout += 3*sH;
  }
  printf("%s\n", airDoneStr(0, ti, ninLen, doneStr)); fflush(stdout);

  nout = nrrdNew();
  airMopAdd(mop, nout, (airMopper)nrrdNuke, airMopAlways);
  if (nrrdMaybeAlloc_va(nout, nrrdTypeFloat, 3,
                        AIR_CAST(size_t, 3),
                        AIR_CAST(size_t, sH),
                        AIR_CAST(size_t, ninLen))) {
    airMopAdd(mop, err=biffGetDone(NRRD), airFree, airMopAlways);
    fprintf(stderr, "%s: error allocating output:\n%s", me, err);
    airMopError(mop);
    return 1;
  }
  out = AIR_CAST(float *, nout->data);
  preout = AIR_CAST(float *, npreout->data);
  for (ti=0; ti<ninLen; ti++) {
    unsigned int hi;
    float hh, vv, ss, scl;
    for (hi=0; hi<sH; hi++) {
      hh = AIR_AFFINE(0, hi, sH, 0, 1);
      if (!preout[hi + 2*sH]) {
        ELL_3V_SET(out + 3*hi, 0, 0, 0);
      } else {
        ss = preout[hi + 2*sH];
        scl = ss/(maxSum + ss);
        vv = scl*preout[hi + 1*sH];
        dyeHSVtoRGB(out + 0 + 3*hi, out + 1 + 3*hi, out + 2 + 3*hi,
                    hh, preout[hi + 0*sH], vv);
      }
    }
    out += 3*sH;
    preout += 3*sH;
  }

  if (nrrdSave(outS, nout, NULL)) {
    airMopAdd(mop, err=biffGetDone(NRRD), airFree, airMopAlways);
    fprintf(stderr, "%s: error saving output:\n%s", me, err);
    airMopError(mop);
    return 1;
  }

  airMopOkay(mop);
  return 0;
}
