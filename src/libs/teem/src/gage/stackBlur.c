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

#include "gage.h"
#include "privateGage.h"

const char *
_gageSigmaSamplingStr[] = {
  "(unknown_sampling)",
  "unisig", /* "uniform-sigma", */
  "unitau", /* "uniform-tau", */
  "optil2" /* "optimal-3d-l2l2" */
};

const char *
_gageSigmaSamplingDesc[] = {
  "unknown sampling",
  "uniform samples along sigma",
  "uniform samples along Lindeberg's tau",
  "optimal sampling (3D L2 image error and L2 error across scales)"
};

const char *
_gageSigmaSamplingStrEqv[] = {
  "uniform-sigma", "unisigma", "unisig",
  "uniform-tau", "unitau",
  "optimal-3d-l2l2", "optimal-l2l2", "optil2",
  ""
};

const int
_gageSigmaSamplingValEqv[] = {
  gageSigmaSamplingUniformSigma, gageSigmaSamplingUniformSigma,
  /* */ gageSigmaSamplingUniformSigma,
  gageSigmaSamplingUniformTau, gageSigmaSamplingUniformTau,
  gageSigmaSamplingOptimal3DL2L2, gageSigmaSamplingOptimal3DL2L2,
  /* */ gageSigmaSamplingOptimal3DL2L2
};

const airEnum
_gageSigmaSampling_enum = {
  "sigma sampling strategy",
  GAGE_SIGMA_SAMPLING_MAX,
  _gageSigmaSamplingStr, NULL,
  _gageSigmaSamplingDesc,
  _gageSigmaSamplingStrEqv, _gageSigmaSamplingValEqv,
  AIR_FALSE
};
const airEnum *const
gageSigmaSampling = &_gageSigmaSampling_enum;


void
gageStackBlurParmInit(gageStackBlurParm *parm) {

  if (parm) {
    parm->num = 0;
    parm->sigmaRange[0] = AIR_NAN;
    parm->sigmaRange[1] = AIR_NAN;
    parm->sigmaSampling = gageSigmaSamplingUnknown;
    parm->sigma = airFree(parm->sigma);
    parm->kspec = nrrdKernelSpecNix(parm->kspec);
    /* this will be effectively moot when nrrdKernelDiscreteGaussian is used
       with a bit cut-off, and will only help with smaller cut-offs and with
       any other kernel, and will be moot for FFT-based blurring */
    parm->renormalize = AIR_TRUE;
    parm->bspec = nrrdBoundarySpecNix(parm->bspec);
    parm->oneDim = AIR_FALSE;
    /* the cautious application of the FFT--based blurring justifies enables
       it by default */
    parm->needSpatialBlur = AIR_FALSE;
    parm->verbose = 1; /* HEY: this may be revisited */
    parm->dgGoodSigmaMax = nrrdKernelDiscreteGaussianGoodSigmaMax;
  }
  return;
}

/*
** does not use biff
*/
gageStackBlurParm *
gageStackBlurParmNew() {
  gageStackBlurParm *parm;

  parm = AIR_CALLOC(1, gageStackBlurParm);
  gageStackBlurParmInit(parm);
  return parm;
}

gageStackBlurParm *
gageStackBlurParmNix(gageStackBlurParm *sbp) {

  if (sbp) {
    airFree(sbp->sigma);
    nrrdKernelSpecNix(sbp->kspec);
    nrrdBoundarySpecNix(sbp->bspec);
    free(sbp);
  }
  return NULL;
}

/*
** *differ is set to 0 or 1; not useful for sorting
*/
int
gageStackBlurParmCompare(const gageStackBlurParm *aa, const char *_nameA,
                         const gageStackBlurParm *bb, const char *_nameB,
                         int *differ, char explain[AIR_STRLEN_LARGE]) {
  static const char me[]="gageStackBlurParmCompare",
    baseA[]="A", baseB[]="B";
  const char *nameA, *nameB;
  unsigned int si, warnLen = AIR_STRLEN_LARGE/4;
  char stmp[2][AIR_STRLEN_LARGE], subexplain[AIR_STRLEN_LARGE];

  if (!(aa && bb && differ)) {
    biffAddf(GAGE, "%s: got NULL pointer (%p %p %p)", me,
             AIR_VOIDP(aa), AIR_VOIDP(bb), AIR_VOIDP(differ));
    return 1;
  }
  nameA = _nameA ? _nameA : baseA;
  nameB = _nameB ? _nameB : baseB;
  if (strlen(nameA) + strlen(nameB) > warnLen) {
    biffAddf(GAGE, "%s: names (len %s, %s) might lead to overflow", me,
             airSprintSize_t(stmp[0], strlen(nameA)),
             airSprintSize_t(stmp[1], strlen(nameB)));
    return 1;
  }
  /*
  ** HEY: really ambivalent about not doing this check:
  ** its unusual in Teem to not take an opportunity to do this kind
  ** of sanity check when its available, but we don't really know the
  ** circumstances of when this will be called, and if that includes
  ** some interaction with hest, there may not yet have been the chance
  ** to complete the sbp.
  if (gageStackBlurParmCheck(aa)) {
    biffAddf(GAGE, "%s: problem with sbp %s", me, nameA);
    return 1;
  }
  if (gageStackBlurParmCheck(bb)) {
    biffAddf(GAGE, "%s: problem with sbp %s", me, nameB);
    return 1;
  }
  */
#define CHECK(VAR, FMT)                                                 \
  if (aa->VAR != bb->VAR) {                                             \
    if (explain) {                                                      \
      sprintf(explain, "%s->" #VAR "=" #FMT " != %s->" #VAR "=" #FMT,   \
              nameA, aa->VAR, nameB,  bb->VAR);                         \
    }                                                                   \
    *differ = 1;                                                        \
    return 0;                                                           \
  }
  CHECK(num, %u);
  CHECK(sigmaRange[0], %.17g);
  CHECK(sigmaRange[1], %.17g);
  CHECK(renormalize, %d);
  CHECK(oneDim, %d);
  CHECK(needSpatialBlur, %d);
  /* This is sketchy: the apparent point of the function is to see if two
     sbp's are different.  But a big role of the function is to enable
     leeching in meet.  And for leeching, a difference in verbose is moot */
  /* CHECK(verbose, %d); */
  CHECK(dgGoodSigmaMax, %.17g);
#undef CHECK
  if (aa->sigmaSampling != bb->sigmaSampling) {
    if (explain) {
      sprintf(explain, "%s->sigmaSampling=%s != %s->sigmaSampling=%s",
              nameA, airEnumStr(gageSigmaSampling, aa->sigmaSampling),
              nameB, airEnumStr(gageSigmaSampling, bb->sigmaSampling));
    }
    *differ = 1; return 0;
  }
  for (si=0; si<aa->num; si++) {
    if (aa->sigma[si] != bb->sigma[si]) {
      if (explain) {
        sprintf(explain, "%s->sigma[%u]=%.17g != %s->sigma[%u]=%.17g",
                nameA, si, aa->sigma[si], nameB, si, bb->sigma[si]);
      }
      *differ = 1; return 0;
    }
  }
  if (nrrdKernelSpecCompare(aa->kspec, bb->kspec,
                            differ, subexplain)) {
    biffMovef(GAGE, NRRD, "%s: trouble comparing kernel specs", me);
    return 1;
  }
  if (*differ) {
    if (explain) {
      sprintf(explain, "kernel specs different: %s", subexplain);
    }
    *differ = 1; return 0;
  }
  if (nrrdBoundarySpecCompare(aa->bspec, bb->bspec,
                              differ, subexplain)) {
    biffMovef(GAGE, NRRD, "%s: trouble comparing boundary specs", me);
    return 1;
  }
  if (*differ) {
    if (explain) {
      sprintf(explain, "boundary specs different: %s", subexplain);
    }
    *differ = 1; return 0;
  }
  /* no differences so far */
  *differ = 0;
  return 0;
}

int
gageStackBlurParmCopy(gageStackBlurParm *dst,
                      const gageStackBlurParm *src) {
  static const char me[]="gageStackBlurParmCopy";
  int differ;
  char explain[AIR_STRLEN_LARGE];

  if (!(dst && src)) {
    biffAddf(GAGE, "%s: got NULL pointer", me);
    return 1;
  }
  if (gageStackBlurParmCheck(src)) {
    biffAddf(GAGE, "%s: given src parm has problems", me);
    return 1;
  }
  if (gageStackBlurParmSigmaSet(dst, src->num,
                                src->sigmaRange[0], src->sigmaRange[1],
                                src->sigmaSampling)
      || gageStackBlurParmKernelSet(dst, src->kspec)
      || gageStackBlurParmRenormalizeSet(dst, src->renormalize)
      || gageStackBlurParmDgGoodSigmaMaxSet(dst, src->dgGoodSigmaMax)
      || gageStackBlurParmBoundarySpecSet(dst, src->bspec)
      || gageStackBlurParmNeedSpatialBlurSet(dst, src->needSpatialBlur)
      || gageStackBlurParmVerboseSet(dst, src->verbose)
      || gageStackBlurParmOneDimSet(dst, src->oneDim)) {
    biffAddf(GAGE, "%s: problem setting dst parm", me);
    return 1;
  }
  if (gageStackBlurParmCompare(dst, "copy", src, "original",
                               &differ, explain)) {
    biffAddf(GAGE, "%s: trouble assessing correctness of copy", me);
    return 1;
  }
  if (differ) {
    biffAddf(GAGE, "%s: problem: copy not equal: %s", me, explain);
    return 1;
  }
  return 0;
}

int
gageStackBlurParmSigmaSet(gageStackBlurParm *sbp, unsigned int num,
                          double sigmaMin, double sigmaMax,
                          int sigmaSampling) {
  static const char me[]="gageStackBlurParmSigmaSet";
  unsigned int ii;

  if (!( sbp )) {
    biffAddf(GAGE, "%s: got NULL pointer", me);
    return 1;
  }
  airFree(sbp->sigma);
  sbp->sigma = NULL;
  if (!( 0 <= sigmaMin )) {
    biffAddf(GAGE, "%s: need sigmaMin >= 0 (not %g)", me, sigmaMin);
    return 1;
  }
  if (!( sigmaMin < sigmaMax )) {
    biffAddf(GAGE, "%s: need sigmaMax %g > sigmaMin %g",
             me, sigmaMax, sigmaMin);
    return 1;
  }
  if (airEnumValCheck(gageSigmaSampling, sigmaSampling)) {
    biffAddf(GAGE, "%s: %d is not a valid %s", me,
             sigmaSampling, gageSigmaSampling->name);
    return 1;
  }
  if (!( num >= 2 )) {
    biffAddf(GAGE, "%s: need # scale samples >= 2 (not %u)", me, num);
    return 1;
  }
  sbp->sigma = AIR_CALLOC(num, double);
  if (!sbp->sigma) {
    biffAddf(GAGE, "%s: couldn't alloc scale for %u", me, num);
    return 1;
  }
  sbp->num = num;
  sbp->sigmaRange[0] = sigmaMin;
  sbp->sigmaRange[1] = sigmaMax;
  sbp->sigmaSampling = sigmaSampling;

  switch (sigmaSampling) {
    double tau0, tau1, tau;
    unsigned int sigmax;
  case gageSigmaSamplingUniformSigma:
    for (ii=0; ii<num; ii++) {
      sbp->sigma[ii] = AIR_AFFINE(0, ii, num-1, sigmaMin, sigmaMax);
    }
    break;
  case gageSigmaSamplingUniformTau:
    tau0 = gageTauOfSig(sigmaMin);
    tau1 = gageTauOfSig(sigmaMax);
    for (ii=0; ii<num; ii++) {
      tau = AIR_AFFINE(0, ii, num-1, tau0, tau1);
      sbp->sigma[ii] = gageSigOfTau(tau);
    }
    break;
  case gageSigmaSamplingOptimal3DL2L2:
    sigmax = AIR_CAST(unsigned int, sigmaMax);
    if (0 != sigmaMin) {
      biffAddf(GAGE, "%s: sigmaMin %g != 0", me, sigmaMin);
      return 1;
    }
    if (sigmax != sigmaMax) {
      biffAddf(GAGE, "%s: sigmaMax %g not an integer", me, sigmaMax);
      return 1;
    }
    if (gageOptimSigSet(sbp->sigma, num, sigmax)) {
      biffAddf(GAGE, "%s: trouble setting optimal sigmas", me);
      return 1;
    }
    break;
  default:
    biffAddf(GAGE, "%s: sorry, sigmaSampling %s (%d) not implemented", me,
             airEnumStr(gageSigmaSampling, sigmaSampling), sigmaSampling);
    return 1;
  }
  if (sbp->verbose > 1) {
    fprintf(stderr, "%s: %u samples in [%g,%g] via %s:\n", me,
            num, sigmaMin, sigmaMax,
            airEnumStr(gageSigmaSampling, sigmaSampling));
    for (ii=0; ii<num; ii++) {
      if (ii) {
        fprintf(stderr, "%s:           "
                "| deltas: %g\t               %g\n", me,
                sbp->sigma[ii] - sbp->sigma[ii-1],
                gageTauOfSig(sbp->sigma[ii])
                - gageTauOfSig(sbp->sigma[ii-1]));
      }
      fprintf(stderr, "%s: sigma[%02u]=%g%s\t         tau=%g\n", me, ii,
              sbp->sigma[ii], !sbp->sigma[ii] ? "     " : "",
              gageTauOfSig(sbp->sigma[ii]));
    }
  }

  return 0;
}

int
gageStackBlurParmScaleSet(gageStackBlurParm *sbp,
                          unsigned int num,
                          double smin, double smax,
                          int uniform, int optimal) {
  static const char me[]="gageStackBlurParmScaleSet";
  int sampling;

  fprintf(stderr, "\n%s: !!! This function is deprecated; use "
          "gageStackBlurParmSigmaSet instead !!!\n\n", me);
  if (uniform && optimal) {
    biffAddf(GAGE, "%s: can't have both uniform and optimal sigma sampling",
             me);
    return 1;
  }
  sampling = (uniform
              ? gageSigmaSamplingUniformSigma
              : (optimal
                 ? gageSigmaSamplingOptimal3DL2L2
                 : gageSigmaSamplingUniformTau));
  if (gageStackBlurParmSigmaSet(sbp, num, smin, smax, sampling)) {
    biffAddf(GAGE, "%s: trouble", me);
    return 1;
  }
  return 0;
}

int
gageStackBlurParmKernelSet(gageStackBlurParm *sbp,
                           const NrrdKernelSpec *kspec) {
  static const char me[]="gageStackBlurParmKernelSet";

  if (!( sbp && kspec )) {
    biffAddf(GAGE, "%s: got NULL pointer", me);
    return 1;
  }
  nrrdKernelSpecNix(sbp->kspec);
  sbp->kspec = nrrdKernelSpecCopy(kspec);
  return 0;
}

int
gageStackBlurParmRenormalizeSet(gageStackBlurParm *sbp,
                                int renormalize) {
  static const char me[]="gageStackBlurParmRenormalizeSet";

  if (!sbp) {
    biffAddf(GAGE, "%s: got NULL pointer", me);
    return 1;
  }
  sbp->renormalize = renormalize;
  return 0;
}

int
gageStackBlurParmBoundarySet(gageStackBlurParm *sbp,
                             int boundary, double padValue) {
  static const char me[]="gageStackBlurParmBoundarySet";

  if (!sbp) {
    biffAddf(GAGE, "%s: got NULL pointer", me);
    return 1;
  }
  nrrdBoundarySpecNix(sbp->bspec);
  sbp->bspec = nrrdBoundarySpecNew();
  sbp->bspec->boundary = boundary;
  sbp->bspec->padValue = padValue;
  if (nrrdBoundarySpecCheck(sbp->bspec)) {
    biffMovef(GAGE, NRRD, "%s: problem", me);
    return 1;
  }
  return 0;
}

int
gageStackBlurParmBoundarySpecSet(gageStackBlurParm *sbp,
                                 const NrrdBoundarySpec *bspec) {
  static const char me[]="gageStackBlurParmBoundarySet";

  if (!sbp) {
    biffAddf(GAGE, "%s: got NULL pointer", me);
    return 1;
  }
  nrrdBoundarySpecNix(sbp->bspec);
  sbp->bspec = nrrdBoundarySpecCopy(bspec);
  if (nrrdBoundarySpecCheck(sbp->bspec)) {
    biffMovef(GAGE, NRRD, "%s: problem", me);
    return 1;
  }
  return 0;
}

int
gageStackBlurParmOneDimSet(gageStackBlurParm *sbp, int oneDim) {
  static const char me[]="gageStackBlurParmOneDimSet";

  if (!sbp) {
    biffAddf(GAGE, "%s: got NULL pointer", me);
    return 1;
  }
  sbp->oneDim = oneDim;
  return 0;
}

int
gageStackBlurParmNeedSpatialBlurSet(gageStackBlurParm *sbp,
                                    int needSpatialBlur) {
  static const char me[]="gageStackBlurParmNeedSpatialBlurSet";

  if (!sbp) {
    biffAddf(GAGE, "%s: got NULL pointer", me);
    return 1;
  }
  sbp->needSpatialBlur = needSpatialBlur;
  return 0;
}

int
gageStackBlurParmVerboseSet(gageStackBlurParm *sbp, int verbose) {
  static const char me[]="gageStackBlurParmVerboseSet";

  if (!sbp) {
    biffAddf(GAGE, "%s: got NULL pointer", me);
    return 1;
  }
  sbp->verbose = verbose;
  return 0;
}

int
gageStackBlurParmDgGoodSigmaMaxSet(gageStackBlurParm *sbp,
                                 double dgGoodSigmaMax) {
  static const char me[]="gageStackBlurParmDgGoodSigmaMaxSet";

  if (!sbp) {
    biffAddf(GAGE, "%s: got NULL pointer", me);
    return 1;
  }
  if (!dgGoodSigmaMax > 0) {
    biffAddf(GAGE, "%s: given dgGoodSigmaMax %g not > 0", me, dgGoodSigmaMax);
    return 1;
  }
  sbp->dgGoodSigmaMax = dgGoodSigmaMax;
  return 0;
}

int
gageStackBlurParmCheck(const gageStackBlurParm *sbp) {
  static const char me[]="gageStackBlurParmCheck";
  unsigned int ii;

  if (!sbp) {
    biffAddf(GAGE, "%s: got NULL pointer", me);
    return 1;
  }
  if (!( sbp->num >= 2)) {
    biffAddf(GAGE, "%s: need num >= 2, not %u", me, sbp->num);
    return 1;
  }
  if (!sbp->sigma) {
    biffAddf(GAGE, "%s: sigma vector not allocated", me);
    return 1;
  }
  if (!sbp->kspec) {
    biffAddf(GAGE, "%s: blurring kernel not set", me);
    return 1;
  }
  if (!sbp->bspec) {
    biffAddf(GAGE, "%s: boundary specification not set", me);
    return 1;
  }
  for (ii=0; ii<sbp->num; ii++) {
    if (!AIR_EXISTS(sbp->sigma[ii])) {
      biffAddf(GAGE, "%s: sigma[%u] = %g doesn't exist", me, ii,
               sbp->sigma[ii]);
      return 1;
    }
    if (ii) {
      if (!( sbp->sigma[ii-1] < sbp->sigma[ii] )) {
        biffAddf(GAGE, "%s: sigma[%u] = %g not < sigma[%u] = %g", me,
                 ii, sbp->sigma[ii-1], ii+1, sbp->sigma[ii]);
        return 1;
      }
    }
  }
  /* HEY: no sanity check on kernel because there is no
     nrrdKernelSpecCheck(), but there should be! */
  if (nrrdBoundarySpecCheck(sbp->bspec)) {
    biffMovef(GAGE, NRRD, "%s: problem with boundary", me);
    return 1;
  }
  return 0;
}

int
gageStackBlurParmParse(gageStackBlurParm *sbp,
                       int extraFlags[256],
                       char **extraParmsP,
                       const char *_str) {
  static const char me[]="gageStackBlurParmParse";
  char *str, *mnmfS, *stok, *slast=NULL, *parmS, *eps;
  int flagSeen[256];
  double sigmaMin, sigmaMax, dggsm;
  unsigned int sigmaNum, parmNum;
  int haveFlags, verbose, verboseGot=AIR_FALSE, dggsmGot=AIR_FALSE,
    sampling, samplingGot=AIR_FALSE, E;
  airArray *mop, *epsArr;
  NrrdKernelSpec *kspec=NULL;
  NrrdBoundarySpec *bspec=NULL;

  if (!( sbp && _str )) {
    biffAddf(GAGE, "%s: got NULL pointer", me);
    return 1;
  }
  if (!( str = airStrdup(_str) )) {
    biffAddf(GAGE, "%s: couldn't copy input", me);
    return 1;
  }
  mop = airMopNew();
  airMopAdd(mop, str, airFree, airMopAlways);
  if (extraParmsP) {
    /* start with empty string */
    epsArr = airArrayNew(AIR_CAST(void **, &eps), NULL, sizeof(char), 42);
    airMopAdd(mop, epsArr, (airMopper)airArrayNuke, airMopAlways);
    airArrayLenIncr(epsArr, 1);
    *eps = '\0';
  } else {
    epsArr = NULL;
  }

  /* working with assumption that '/' does not appear
     in mnmfS <minScl>-<#smp>-<maxScl>[-<flags>] */
  if ( (parmS = strchr(str, '/')) ) {
    /* there are in fact parms */
    *parmS = '\0';
    parmS++;
  } else {
    parmS = NULL;
  }
  mnmfS = str;
  if (!( 3 == airStrntok(mnmfS, "-") || 4 == airStrntok(mnmfS, "-") )) {
    biffAddf(GAGE, "%s: didn't get 3 or 4 \"-\"-separated tokens in \"%s\"",
             me, mnmfS);
    airMopError(mop); return 1;
  }
  haveFlags = (4 == airStrntok(mnmfS, "-"));
  stok = airStrtok(mnmfS, "-", &slast);
  if (1 != sscanf(stok, "%lg", &sigmaMin)) {
    biffAddf(GAGE, "%s: couldn't parse \"%s\" as max sigma", me, stok);
    airMopError(mop); return 1;
  }
  stok = airStrtok(NULL, "-", &slast);
  if (1 != sscanf(stok, "%u", &sigmaNum)) {
    biffAddf(GAGE, "%s: couldn't parse \"%s\" as # scale samples", me, stok);
    airMopError(mop); return 1;
  }
  stok = airStrtok(NULL, "-", &slast);
  if (1 != sscanf(stok, "%lg", &sigmaMax)) {
    biffAddf(GAGE, "%s: couldn't parse \"%s\" as max scale", me, stok);
    airMopError(mop); return 1;
  }
  bzero(flagSeen, sizeof(flagSeen));
  if (extraFlags) {
    /* not sizeof(extraFlags) == sizeof(int*) */
    bzero(extraFlags, sizeof(flagSeen));
  }
  if (haveFlags) {
    char *flags, *ff;
    /* look for various things in flags */
    flags = airToLower(airStrdup(airStrtok(NULL, "-", &slast)));
    airMopAdd(mop, flags, airFree, airMopAlways);
    ff = flags;
    while (*ff && '+' != *ff) {
      /* '1': oneDim
         'r': turn OFF spatial kernel renormalize
         'u': uniform (in sigma) sampling
         'o': optimized (3d l2l2) sampling
         'p': need spatial blur
      */
      if (strchr("1ruop", *ff)) {
        flagSeen[AIR_CAST(unsigned char, *ff)] = AIR_TRUE;
      } else {
        if (extraFlags) {
          extraFlags[AIR_CAST(unsigned char, *ff)] = AIR_TRUE;
        } else {
          biffAddf(GAGE, "%s: got extra flag '%c' but NULL extraFlag",
                   me, *ff);
          airMopError(mop); return 1;
        }
      }
      ff++;
    }
    if (flagSeen['u'] && flagSeen['o']) {
      biffAddf(GAGE, "%s: can't have both optimal ('o') and uniform ('u') "
               "flags set in \"%s\"", me, flags);
      airMopError(mop); return 1;
    }
    if (ff && '+' == *ff) {
      biffAddf(GAGE, "%s: sorry, can no longer indicate a derivative "
               "normalization bias via '+' in \"%s\" in flags \"%s\"; "
               "use \"dnbias=\" parm instead", me, ff, flags);
      airMopError(mop); return 1;
    }
  }
  if (parmS) {
    unsigned int parmIdx;
    char *pval, xeq[AIR_STRLEN_SMALL];
    parmNum = airStrntok(parmS, "/");
    for (parmIdx=0; parmIdx<parmNum; parmIdx++) {
      if (!parmIdx) {
        stok = airStrtok(parmS, "/", &slast);
      } else {
        stok = airStrtok(NULL, "/", &slast);
      }
      if (strcpy(xeq, "k=") && stok == strstr(stok, xeq)) {
        pval = stok + strlen(xeq);
        kspec = nrrdKernelSpecNew();
        airMopAdd(mop, kspec, (airMopper)nrrdKernelSpecNix, airMopAlways);
        if (nrrdKernelSpecParse(kspec, pval)) {
          biffMovef(GAGE, NRRD, "%s: couldn't parse \"%s\" as blurring kernel",
                    me, pval);
          airMopError(mop); return 1;
        }
      } else if (strcpy(xeq, "b=") && strstr(stok, xeq) == stok) {
        pval = stok + strlen(xeq);
        bspec = nrrdBoundarySpecNew();
        airMopAdd(mop, bspec, (airMopper)nrrdBoundarySpecNix, airMopAlways);
        if (nrrdBoundarySpecParse(bspec, pval)) {
          biffMovef(GAGE, NRRD, "%s: couldn't parse \"%s\" as boundary",
                    me, pval);
          airMopError(mop); return 1;
        }
      } else if (strcpy(xeq, "v=") && strstr(stok, xeq) == stok) {
        pval = stok + strlen(xeq);
        if (1 != sscanf(pval, "%d", &verbose)) {
          biffAddf(GAGE, "%s: couldn't parse \"%s\" as verbose int", me, pval);
          airMopError(mop); return 1;
        }
        verboseGot = AIR_TRUE;
      } else if (strcpy(xeq, "s=") && strstr(stok, xeq) == stok) {
        pval = stok + strlen(xeq);
        sampling = airEnumVal(gageSigmaSampling, pval);
        if (gageSigmaSamplingUnknown == sampling) {
          biffAddf(GAGE, "%s: couldn't parse \"%s\" as %s", me, pval,
                   gageSigmaSampling->name);
          airMopError(mop); return 1;
        }
        samplingGot = AIR_TRUE;
      } else if (strcpy(xeq, "dggsm=") && strstr(stok, xeq) == stok) {
        pval = stok + strlen(xeq);
        if (1 != sscanf(pval, "%lg", &dggsm)) {
          biffAddf(GAGE, "%s: couldn't parse \"%s\" as dgGoodSigmaMax double",
                   me, pval);
          airMopError(mop); return 1;
        }
        dggsmGot = AIR_TRUE;
      } else {
        /* doesn't match any of the parms we know how to parse */
        if (extraParmsP) {
          airArrayLenIncr(epsArr, AIR_CAST(int, 2 + strlen(stok)));
          if (strlen(eps)) {
            strcat(eps, "/");
          }
          strcat(eps, stok);
        } else {
          biffAddf(GAGE, "%s: got extra parm \"%s\" but NULL extraParmsP",
                   me, stok);
          airMopError(mop); return 1;
        }
      }
    }
  }
  /* have parsed everything, now error checking and making sense */
  if (flagSeen['u'] && flagSeen['o']) {
    biffAddf(GAGE, "%s: can't use flags 'u' and 'o' at same time", me);
    airMopError(mop); return 1;
  }
  if ((flagSeen['u'] || flagSeen['o']) && samplingGot) {
    biffAddf(GAGE, "%s: can't use both 'u','o' flags and parms to "
             "specify sigma sampling", me);
    airMopError(mop); return 1;
  }
  if (!samplingGot) {
    /* have to set sampling from flags */
    if (flagSeen['u']) {
      sampling = gageSigmaSamplingUniformSigma;
    } else if (flagSeen['o']) {
      sampling = gageSigmaSamplingOptimal3DL2L2;
    } else {
      sampling = gageSigmaSamplingUniformTau;
    }
  }
  /* setting sbp fields */
  E = 0;
  if (!E) E |= gageStackBlurParmSigmaSet(sbp, sigmaNum,
                                         sigmaMin, sigmaMax, sampling);
  if (kspec) {
    if (!E) E |= gageStackBlurParmKernelSet(sbp, kspec);
  }
  if (flagSeen['r']) {
    if (!E) E |= gageStackBlurParmRenormalizeSet(sbp, AIR_FALSE);
  }
  if (dggsmGot) {
    if (!E) E |= gageStackBlurParmDgGoodSigmaMaxSet(sbp, dggsm);
  }
  if (bspec) {
    if (!E) E |= gageStackBlurParmBoundarySpecSet(sbp, bspec);
  }
  if (flagSeen['p']) {
    if (!E) E |= gageStackBlurParmNeedSpatialBlurSet(sbp, AIR_TRUE);
  }
  if (verboseGot) {
    if (!E) E |= gageStackBlurParmVerboseSet(sbp, verbose);
  }
  if (flagSeen['1']) {
    if (!E) E |= gageStackBlurParmOneDimSet(sbp, AIR_TRUE);
  }
  /* NOT doing the final check, because if this is being called from
     hest, the caller won't have had time to set the default info in
     the sbp (like the default kernel), so it will probably look
     incomplete.
     if (!E) E |= gageStackBlurParmCheck(sbp); */
  if (E) {
    biffAddf(GAGE, "%s: problem with blur parm specification", me);
    airMopError(mop); return 1;
  }
  if (extraParmsP) {
    if (airStrlen(eps)) {
      *extraParmsP = airStrdup(eps);
    } else {
      *extraParmsP = NULL;
    }
  }
  airMopOkay(mop);
  return 0;
}

int
gageStackBlurParmSprint(char str[AIR_STRLEN_LARGE],
                        const gageStackBlurParm *sbp,
                        int extraFlag[256],
                        char *extraParm) {
  static const char me[]="gageStackBlurParmSprint";
  char *out, stmp[AIR_STRLEN_LARGE];
  int needFlags, hef;
  unsigned int fi;

  if (!(str && sbp)) {
    biffAddf(GAGE, "%s: got NULL pointer", me);
    return 1;
  }

  out = str;
  sprintf(out, "%.17g-%u-%.17g",
          sbp->sigmaRange[0], sbp->num, sbp->sigmaRange[1]);
  out += strlen(out);
  hef = AIR_FALSE;
  if (extraFlag) {
    for (fi=0; fi<256; fi++) {
      hef |= extraFlag[fi];
    }
  }
  needFlags = (sbp->oneDim
               || sbp->renormalize
               || sbp->needSpatialBlur
               || hef);
  if (needFlags) {
    strcat(out, "-");
    if (sbp->oneDim)          { strcat(out, "1"); }
    if (sbp->renormalize)     { strcat(out, "r"); }
    if (sbp->needSpatialBlur) { strcat(out, "p"); }
    if (hef) {
      for (fi=0; fi<256; fi++) {
        if (extraFlag[fi]) {
          sprintf(stmp, "%c", AIR_CAST(char, fi));
          strcat(out, stmp);
        }
      }
    }
  }

  if (sbp->kspec) {
    strcat(out, "/");
    if (nrrdKernelSpecSprint(stmp, sbp->kspec)) {
      biffMovef(GAGE, NRRD, "%s: problem with kernel", me);
      return 1;
    }
    strcat(out, "k="); strcat(out, stmp);
  }

  if (sbp->bspec) {
    strcat(out, "/");
    if (nrrdBoundarySpecSprint(stmp, sbp->bspec)) {
      biffMovef(GAGE, NRRD, "%s: problem with boundary", me);
      return 1;
    }
    strcat(out, "b="); strcat(out, stmp);
  }

  if (!airEnumValCheck(gageSigmaSampling, sbp->sigmaSampling)) {
    strcat(out, "/s=");
    strcat(out, airEnumStr(gageSigmaSampling, sbp->sigmaSampling));
  }

  if (sbp->verbose) {
    sprintf(stmp, "/v=%d", sbp->verbose);
    strcat(out, stmp);
  }

  if (sbp->kspec
      && nrrdKernelDiscreteGaussian == sbp->kspec->kernel
      && nrrdKernelDiscreteGaussianGoodSigmaMax != sbp->dgGoodSigmaMax) {
    sprintf(stmp, "/dggsm=%.17g", sbp->dgGoodSigmaMax);
    strcat(out, stmp);
  }

  if (extraParm) {
    strcat(out, "/");
    strcat(out, extraParm);
  }

  return 0;
}

int
_gageHestStackBlurParmParse(void *ptr, char *str, char err[AIR_STRLEN_HUGE]) {
  gageStackBlurParm **sbp;
  char me[]="_gageHestStackBlurParmParse", *nerr;

  if (!(ptr && str)) {
    sprintf(err, "%s: got NULL pointer", me);
    return 1;
  }
  sbp = (gageStackBlurParm **)ptr;
  if (!strlen(str)) {
    /* got an empty string; we trying to emulate an "optional"
       command-line option, in a program that likely still has
       the older -ssn, -ssr, -kssb, and which may or may not
       need any scale-space functionality */
    *sbp = NULL;
  } else {
    *sbp = gageStackBlurParmNew();
    /* NOTE: no way to retrieve extraFlags or extraParms from hest */
    if (gageStackBlurParmParse(*sbp, NULL, NULL, str)) {
      nerr = biffGetDone(GAGE);
      airStrcpy(err, AIR_STRLEN_HUGE, nerr);
      gageStackBlurParmNix(*sbp);
      free(nerr);
      return 1;
    }
  }
  return 0;
}

hestCB
_gageHestStackBlurParm = {
  sizeof(gageStackBlurParm*),
  "stack blur specification",
  _gageHestStackBlurParmParse,
  (airMopper)gageStackBlurParmNix
};

hestCB *
gageHestStackBlurParm = &_gageHestStackBlurParm;

static int
_checkNrrd(Nrrd *const nblur[], const Nrrd *const ncheck[],
           unsigned int blNum, int checking,
           const Nrrd *nin, const gageKind *kind) {
  static const char me[]="_checkNrrd";
  unsigned int blIdx;

  for (blIdx=0; blIdx<blNum; blIdx++) {
    if (checking) {
      if (nrrdCheck(ncheck[blIdx])) {
        biffMovef(GAGE, NRRD, "%s: bad ncheck[%u]", me, blIdx);
        return 1;
      }
    } else {
      if (!nblur[blIdx]) {
        biffAddf(GAGE, "%s: NULL blur[%u]", me, blIdx);
        return 1;
      }
    }
  }
  if (3 + kind->baseDim != nin->dim) {
    biffAddf(GAGE, "%s: need nin->dim %u (not %u) with baseDim %u", me,
             3 + kind->baseDim, nin->dim, kind->baseDim);
    return 1;
  }
  return 0;
}

#define KVP_NUM 9

static const char
_blurKey[KVP_NUM][AIR_STRLEN_LARGE] = {/*  0  */ "gageStackBlur",
                                       /*  1 */  "cksum",
                                       /*  2  */ "scale",
                                       /*  3  */ "kernel",
                                       /*  4  */ "renormalize",
                                       /*  5  */ "boundary",
                                       /*  6  */ "onedim",
                                       /*  7  */ "spatialblurred",
#define KVP_SBLUR_IDX                      7
                                       /*  8  */ "dgGoodSigmaMax"
#define KVP_DGGSM_IDX                      8
                                       /* (9 == KVP_NUM, above) */};

typedef struct {
  char val[KVP_NUM][AIR_STRLEN_LARGE];
} blurVal_t;

static blurVal_t *
_blurValAlloc(airArray *mop, gageStackBlurParm *sbp, NrrdKernelSpec *kssb,
              const Nrrd *nin, int spatialBlurred) {
  static const char me[]="_blurValAlloc";
  blurVal_t *blurVal;
  unsigned int blIdx, cksum;

  blurVal = AIR_CAST(blurVal_t *, calloc(sbp->num, sizeof(blurVal_t)));
  if (!blurVal) {
    biffAddf(GAGE, "%s: couldn't alloc blurVal for %u", me, sbp->num);
    return NULL;
  }
  cksum = nrrdCRC32(nin, airEndianLittle);

  for (blIdx=0; blIdx<sbp->num; blIdx++) {
    kssb->parm[0] = sbp->sigma[blIdx];
    sprintf(blurVal[blIdx].val[0], "true");
    sprintf(blurVal[blIdx].val[1], "%u", cksum);
    sprintf(blurVal[blIdx].val[2], "%.17g", sbp->sigma[blIdx]);
    nrrdKernelSpecSprint(blurVal[blIdx].val[3], kssb);
    sprintf(blurVal[blIdx].val[4], "%s", sbp->renormalize ? "true" : "false");
    nrrdBoundarySpecSprint(blurVal[blIdx].val[5], sbp->bspec);
    sprintf(blurVal[blIdx].val[6], "%s",
            sbp->oneDim ? "true" : "false");
    sprintf(blurVal[blIdx].val[7], "%s",
            spatialBlurred ? "true" : "false");
    sprintf(blurVal[blIdx].val[8], "%.17g", sbp->dgGoodSigmaMax);
  }
  airMopAdd(mop, blurVal, airFree, airMopAlways);
  return blurVal;
}

/*
** some spot checks suggest that where the PSF of lindeberg-gaussian blurring
** should be significantly non-zero, this is more accurate than the current
** "discrete gauss" kernel, but for small values near zero, the spatial
** blurring is more accurate.  This is due to how with limited numerical
** precision, the FFT can produce very low amplitude noise.
*/
static int
_stackBlurDiscreteGaussFFT(Nrrd *const nblur[], gageStackBlurParm *sbp,
                           const Nrrd *nin, const gageKind *kind) {
  static const char me[]="_stackBlurDiscreteGaussFFT";
  size_t sizeAll[NRRD_DIM_MAX], *size, ii, xi, yi, zi, nn;
  Nrrd *ninC, /* complex version of input, same as input type */
    *ninFT,   /* FT of input, type double */
    *noutFT,  /* FT of output, values set manually from ninFT, as double */
    *noutCd,  /* complex version of output, still as double */
    *noutC;   /* complex version of output, as input type;
                 will convert/clamp this to get nblur[i] */
  double (*lup)(const void *, size_t), (*ins)(void *, size_t, double),
    *ww[3], tblur, theta, *inFT, *outFT;
  unsigned int blIdx, axi, ftaxes[3] = {1,2,3};
  airArray *mop;
  int axmap[NRRD_DIM_MAX];

  mop = airMopNew();
  ninC = nrrdNew();
  airMopAdd(mop, ninC, (airMopper)nrrdNuke, airMopAlways);
  ninFT = nrrdNew();
  airMopAdd(mop, ninFT, (airMopper)nrrdNuke, airMopAlways);
  noutFT = nrrdNew();
  airMopAdd(mop, noutFT, (airMopper)nrrdNuke, airMopAlways);
  noutCd = nrrdNew();
  airMopAdd(mop, noutCd, (airMopper)nrrdNuke, airMopAlways);
  noutC = nrrdNew();
  airMopAdd(mop, noutC, (airMopper)nrrdNuke, airMopAlways);

  if (gageKindScl != kind) {
    biffAddf(GAGE, "%s: sorry, non-scalar kind not yet implemented", me);
    /* but it really wouldn't be that hard ... */
    airMopError(mop); return 1;
  }
  if (3 != nin->dim) {
    biffAddf(GAGE, "%s: sanity check fail: nin->dim %u != 3", me, nin->dim);
    airMopError(mop); return 1;
  }
  lup = nrrdDLookup[nin->type];
  ins = nrrdDInsert[nin->type];
  /* unurrdu/fft.c handles real input by doing an axis insert and then
     padding, but that is overkill; this is a more direct way, which perhaps
     should be migrated to unurrdu/fft.c */
  sizeAll[0] = 2;
  size = sizeAll + 1;
  nrrdAxisInfoGet_nva(nin, nrrdAxisInfoSize, size);
  if (nrrdMaybeAlloc_nva(ninC, nin->type, nin->dim+1, sizeAll)) {
    biffMovef(GAGE, NRRD, "%s: couldn't allocate complex-valued input", me);
    airMopError(mop); return 1;
  }
  for (axi=0; axi<3; axi++) {
    if (!( ww[axi] = AIR_CALLOC(size[axi], double) )) {
      biffAddf(GAGE, "%s: couldn't allocate axis %u buffer", me, axi);
      airMopError(mop); return 1;
    }
    airMopAdd(mop, ww[axi], airFree, airMopAlways);
  }
  nn = size[0]*size[1]*size[2];
  for (ii=0; ii<nn; ii++) {
    ins(ninC->data, 0 + 2*ii, lup(nin->data, ii));
    ins(ninC->data, 1 + 2*ii, 0.0);
  }
  for (axi=0; axi<4; axi++) {
    if (!axi) {
      axmap[axi] = -1;
    } else {
      axmap[axi] = axi-1;
    }
  }
  if (nrrdAxisInfoCopy(ninC, nin, axmap, NRRD_AXIS_INFO_NONE)
      || nrrdBasicInfoCopy(ninC, nin,
                           (NRRD_BASIC_INFO_DATA_BIT |
                            NRRD_BASIC_INFO_DIMENSION_BIT |
                            NRRD_BASIC_INFO_CONTENT_BIT |
                            NRRD_BASIC_INFO_COMMENTS_BIT |
                            NRRD_BASIC_INFO_KEYVALUEPAIRS_BIT))) {
    biffMovef(GAGE, NRRD, "%s: couldn't set complex-valued axinfo", me);
    airMopError(mop); return 1;
  }
  ninC->axis[0].kind = nrrdKindComplex; /* should use API */
  /*
  nrrdSave("ninC.nrrd", ninC, NULL);
  */
  /* the copy to noutFT is just to allocate it; the values
     there will be re-written over and over in the loop below */
  if (nrrdFFT(ninFT, ninC, ftaxes, 3,
              +1 /* forward */,
              AIR_TRUE /* rescale */,
              nrrdFFTWPlanRigorEstimate /* should generalize! */)
      || nrrdCopy(noutFT, ninFT)) {
    biffMovef(GAGE, NRRD, "%s: trouble with initial transforms", me);
    airMopError(mop); return 1;
  }
  /*
  nrrdSave("ninFT.nrrd", ninFT, NULL);
  */
  inFT = AIR_CAST(double *, ninFT->data);
  outFT = AIR_CAST(double *, noutFT->data);
  for (blIdx=0; blIdx<sbp->num; blIdx++) {
    if (sbp->verbose) {
      fprintf(stderr, "%s: . . . %u/%u (scale %g, tau %g) . . . ", me,
              blIdx, sbp->num, sbp->sigma[blIdx],
              gageTauOfSig(sbp->sigma[blIdx]));
      fflush(stderr);
    }
    tblur = sbp->sigma[blIdx]*sbp->sigma[blIdx];
    for (axi=0; axi<3; axi++) {
      for (ii=0; ii<size[axi]; ii++) {
        theta = AIR_AFFINE(0, ii, size[axi], 0.0, 2*AIR_PI);
        /* from eq (22) of T. Lindeberg "Scale-Space for Discrete
           Signals", IEEE PAMI 12(234-254); 1990 */
        ww[axi][ii] = exp(tblur*(cos(theta)-1.0));
        /*
        fprintf(stderr, "!%s: ww[%u][%u] = %g\n", me, axi,
                AIR_CAST(unsigned int, ii), ww[axi][ii]);
        */
      }
    }
    ii=0;
    /*
    for (axi=0; axi<3; axi++) {
      fprintf(stderr, "!%s: size[%u] = %u\n", me, axi,
              AIR_CAST(unsigned int, size[axi]));
    }
    */
    for (zi=0; zi<size[2]; zi++) {
      for (yi=0; yi<size[1]; yi++) {
        for (xi=0; xi<size[0]; xi++) {
          double wght;
          wght = sbp->oneDim ? 1.0 : ww[1][yi]*ww[2][zi];
          wght *= ww[0][xi];
          outFT[0 + 2*ii] = wght*inFT[0 + 2*ii];
          outFT[1 + 2*ii] = wght*inFT[1 + 2*ii];
          /*
          fprintf(stderr, "!%s: out[%u] = (%g,%g) = %g * (%g,%g)\n", me,
                  AIR_CAST(unsigned int, ii),
                  outFT[0 + 2*ii], outFT[1 + 2*ii],
                  wght,
                  inFT[0 + 2*ii], inFT[1 + 2*ii]);
          */
          ii++;
        }
      }
    }
    if (nrrdFFT(noutCd, noutFT, ftaxes, 3,
                -1 /* backward */,
                AIR_TRUE /* rescale */,
                nrrdFFTWPlanRigorEstimate /* should generalize! */)
        || (nrrdTypeDouble == nin->type
            ? nrrdCopy(noutC, noutCd)
            : nrrdCastClampRound(noutC, noutCd, nin->type,
                                 AIR_TRUE /* clamp */,
                                 +1 /* roundDir, when needed */))
        || nrrdSlice(nblur[blIdx], noutC, 0, 0)
        || nrrdContentSet_va(nblur[blIdx], "blur", nin, "")) {
      biffMovef(GAGE, NRRD, "%s: trouble with back transform %u", me, blIdx);
      airMopError(mop); return 1;
    }
    if (sbp->verbose) {
      fprintf(stderr, "done\n");
    }
    /*
    if (0) {
      char fname[AIR_STRLEN_SMALL];
      sprintf(fname, "noutFT-%03u.nrrd", blIdx);
      nrrdSave(fname, noutFT, NULL);
      sprintf(fname, "noutCd-%03u.nrrd", blIdx);
      nrrdSave(fname, noutCd, NULL);
      sprintf(fname, "noutC-%03u.nrrd", blIdx);
      nrrdSave(fname, noutC, NULL);
      sprintf(fname, "nblur-%03u.nrrd", blIdx);
      nrrdSave(fname, nblur[blIdx], NULL);
    }
    */
  }

  airMopOkay(mop);
  return 0;
}

static int
_stackBlurSpatial(Nrrd *const nblur[], gageStackBlurParm *sbp,
                  NrrdKernelSpec *kssb,
                  const Nrrd *nin, const gageKind *kind) {
  static const char me[]="_stackBlurSpatial";
  NrrdResampleContext *rsmc;
  Nrrd *niter;
  unsigned int axi, blIdx;
  int E, iterative, rsmpType;
  double timeStepMax, /* max length of diffusion time allowed per blur,
                         as determined by sbp->dgGoodSigmaMax */
    timeDone,         /* amount of diffusion time just applied */
    timeLeft;         /* amount of diffusion time left to do */
  airArray *mop;

  mop = airMopNew();
  rsmc = nrrdResampleContextNew();
  airMopAdd(mop, rsmc, (airMopper)nrrdResampleContextNix, airMopAlways);
  if (nrrdKernelDiscreteGaussian == kssb->kernel) {
    iterative = AIR_TRUE;
    /* we don't want to lose precision when iterating */
    rsmpType = nrrdResample_nt;
    /* may be used with iterative diffusion */
    niter = nrrdNew();
    airMopAdd(mop, niter, (airMopper)nrrdNuke, airMopAlways);
  } else {
    iterative = AIR_FALSE;
    rsmpType = nrrdTypeDefault;
    niter = NULL;
  }

  E = 0;
  if (!E) E |= nrrdResampleDefaultCenterSet(rsmc, nrrdDefaultCenter);
  /* the input for the first scale is indeed nin, regardless
     of iterative */
  if (!E) E |= nrrdResampleInputSet(rsmc, nin);
  if (kind->baseDim) {
    unsigned int bai;
    for (bai=0; bai<kind->baseDim; bai++) {
      if (!E) E |= nrrdResampleKernelSet(rsmc, bai, NULL, NULL);
    }
  }
  for (axi=0; axi<3; axi++) {
    if (!E) E |= nrrdResampleSamplesSet(rsmc, kind->baseDim + axi,
                                        nin->axis[kind->baseDim + axi].size);
    if (!E) E |= nrrdResampleRangeFullSet(rsmc, kind->baseDim + axi);
  }
  if (!E) E |= nrrdResampleBoundarySpecSet(rsmc, sbp->bspec);
  if (!E) E |= nrrdResampleTypeOutSet(rsmc, rsmpType);
  if (!E) E |= nrrdResampleClampSet(rsmc, AIR_TRUE); /* probably moot */
  if (!E) E |= nrrdResampleRenormalizeSet(rsmc, sbp->renormalize);
  if (E) {
    biffAddf(GAGE, "%s: trouble setting up resampling", me);
    airMopError(mop); return 1;
  }

  timeDone = 0;
  timeStepMax = (sbp->dgGoodSigmaMax)*(sbp->dgGoodSigmaMax);
  for (blIdx=0; blIdx<sbp->num; blIdx++) {
    if (sbp->verbose) {
      fprintf(stderr, "%s: . . . blurring %u / %u (scale %g) . . . ",
              me, blIdx, sbp->num, sbp->sigma[blIdx]);
      fflush(stderr);
    }
    if (iterative) {
      double timeNow = sbp->sigma[blIdx]*sbp->sigma[blIdx];
      unsigned int passIdx = 0;
      timeLeft = timeNow - timeDone;
      if (sbp->verbose) {
        fprintf(stderr, "\n");
        fprintf(stderr, "%s: scale %g == time %g (tau %g);\n"
                "               timeLeft %g = %g - %g\n",
                me, sbp->sigma[blIdx], timeNow, gageTauOfTee(timeNow),
                timeLeft, timeNow, timeDone);
        if (timeLeft > timeStepMax) {
          fprintf(stderr, "%s: diffusing for time %g in steps of %g\n", me,
                  timeLeft, timeStepMax);
        }
        fflush(stderr);
      }
      do {
        double timeDo;
        if (blIdx || passIdx) {
          /* either we're past the first scale (blIdx >= 1), or
             (unlikely) we're on the first scale but after the first
             pass of a multi-pass blurring, so we have to feed the
             previous result back in as input.
             AND: the way that niter is being used is very sneaky,
             and probably too clever: the resampling happens in
             multiple passes, among buffers internal to nrrdResample;
             so its okay have the output and input nrrds be the same:
             they're never used at the same time. */
          if (!E) E |= nrrdResampleInputSet(rsmc, niter);
        }
        timeDo = (timeLeft > timeStepMax
                  ? timeStepMax
                  : timeLeft);
        /* it is the repeated re-setting of this parm[0] which motivated
           copying to our own kernel spec, so that the given one in the
           gageStackBlurParm can stay untouched */
        kssb->parm[0] = sqrt(timeDo);
        for (axi=0; axi<3; axi++) {
          if (!sbp->oneDim || !axi) {
            /* we set the blurring kernel on this axis if
               we are NOT doing oneDim, or, we are,
               but this is axi == 0 */
            if (!E) E |= nrrdResampleKernelSet(rsmc, kind->baseDim + axi,
                                               kssb->kernel,
                                               kssb->parm);
          } else {
            /* what to do with oneDom on axi 1, 2 */
            /* you might think that we should just do no resampling at all
               on this axis, but that would undermine the in==out==niter
               trick described above; and produce the mysterious behavior
               that the second scale-space volume is all 0.0 */
            double boxparm[NRRD_KERNEL_PARMS_NUM] = {1.0};
            if (!E) E |= nrrdResampleKernelSet(rsmc, kind->baseDim + axi,
                                               nrrdKernelBox, boxparm);
          }
        }
        if (sbp->verbose) {
          fprintf(stderr, "  pass %u (timeLeft=%g => "
                  "time=%g, sigma=%g) ...\n",
                  passIdx, timeLeft, timeDo, kssb->parm[0]);
        }
        if (!E) E |= nrrdResampleExecute(rsmc, niter);
        /* for debugging problem with getting zero output
        if (!E) {
          NrrdRange *nrange;
          nrange = nrrdRangeNewSet(niter, AIR_FALSE);
          fprintf(stderr, "%s: min/max = %g/%g\n", me,
                  nrange->min, nrange->max);
          if (!nrange->min || !nrange->max) {
            fprintf(stderr, "%s: what? zero zero\n", me);
            biffAddf(GAGE, "%s: no good", me);
            airMopError(mop); return 1;
          }
        }
        */
        timeLeft -= timeDo;
        passIdx++;
      } while (!E && timeLeft > 0.0);
      /* at this point we have to replicate the behavior of the
         last stage of resampling (e.g. _nrrdResampleOutputUpdate
         in nrrd/resampleContext.c), since we've gently hijacked
         the resampling to access the nrrdResample_t blurring
         result (for further blurring) */
      if (!E) E |= nrrdCastClampRound(nblur[blIdx], niter, nin->type,
                                      AIR_TRUE,
                                      nrrdTypeIsIntegral[nin->type]);
      if (!E) E |= nrrdContentSet_va(nblur[blIdx], "blur", nin, "");
      timeDone = timeNow;
    } else { /* do blurring in one shot */
      kssb->parm[0] = sbp->sigma[blIdx];
      for (axi=0; axi<(sbp->oneDim ? 1 : 3); axi++) {
        if (!E) E |= nrrdResampleKernelSet(rsmc, kind->baseDim + axi,
                                           kssb->kernel,
                                           kssb->parm);
      }
      if (!E) E |= nrrdResampleExecute(rsmc, nblur[blIdx]);
    }
    if (E) {
      if (sbp->verbose) {
        fprintf(stderr, "problem!\n");
      }
      biffMovef(GAGE, NRRD, "%s: trouble w/ %u of %u (scale %g)",
                me, blIdx, sbp->num, sbp->sigma[blIdx]);
      airMopError(mop); return 1;
    }
    if (sbp->verbose) {
      fprintf(stderr, "  done.\n");
    }
  } /* for blIdx */

  airMopOkay(mop);
  return 0;
}

/*
** little helper function to do pre-blurring of a given nrrd
** of the sort that might be useful for scale-space gage use
**
** nblur has to already be allocated for "blNum" Nrrd*s, AND, they all
** have to point to valid (possibly empty) Nrrds, so they can hold the
** results of blurring
*/
int
gageStackBlur(Nrrd *const nblur[], gageStackBlurParm *sbp,
              const Nrrd *nin, const gageKind *kind) {
  static const char me[]="gageStackBlur";
  unsigned int blIdx, kvpIdx;
  NrrdKernelSpec *kssb;
  blurVal_t *blurVal;
  airArray *mop;
  int E, fftable, spatialBlurred;

  if (!(nblur && sbp && nin && kind)) {
    biffAddf(GAGE, "%s: got NULL pointer", me);
    return 1;
  }
  if (gageStackBlurParmCheck(sbp)) {
    biffAddf(GAGE, "%s: problem with parms", me);
    return 1;
  }
  if (_checkNrrd(nblur, NULL, sbp->num, AIR_FALSE, nin, kind)) {
    biffAddf(GAGE, "%s: problem with input ", me);
    return 1;
  }
  mop = airMopNew();
  kssb = nrrdKernelSpecCopy(sbp->kspec);
  airMopAdd(mop, kssb, (airMopper)nrrdKernelSpecNix, airMopAlways);
  /* see if we can use FFT-based implementation */
  fftable = (!sbp->needSpatialBlur
             && nrrdBoundaryWrap == sbp->bspec->boundary
             && nrrdKernelDiscreteGaussian == sbp->kspec->kernel);
  if (fftable && nrrdFFTWEnabled) {
    /* go directly to FFT-based blurring */
    if (_stackBlurDiscreteGaussFFT(nblur, sbp, nin, kind)) {
      biffAddf(GAGE, "%s: trouble with frequency-space blurring", me);
      airMopError(mop); return 1;
    }
    spatialBlurred = AIR_FALSE;
  } else { /* else either not fft-able, or not it was, but not available;
              in either case we have to do spatial blurring */
    if (fftable && !nrrdFFTWEnabled) {
      if (sbp->verbose) {
        fprintf(stderr, "%s: NOTE: FFT-based blurring applicable but not "
                "available in this Teem build (not built with FFTW)\n", me);
      }
    } else {
      if (sbp->verbose) {
        char kstr[AIR_STRLEN_LARGE], bstr[AIR_STRLEN_LARGE];
        nrrdKernelSpecSprint(kstr, kssb);
        nrrdBoundarySpecSprint(bstr, sbp->bspec);
        fprintf(stderr, "%s: (FFT-based blurring not applicable: "
                "need spatial blur=%s, boundary=%s, kernel=%s)\n", me,
                sbp->needSpatialBlur ? "yes" : "no", bstr, kstr);
      }
    }
    if (_stackBlurSpatial(nblur, sbp, kssb, nin, kind)) {
      biffAddf(GAGE, "%s: trouble with spatial-domain blurring", me);
      airMopError(mop); return 1;
    }
    spatialBlurred = AIR_TRUE;
  }
  /* add the KVPs to document how these were blurred */
  if (!( blurVal = _blurValAlloc(mop, sbp, kssb, nin, spatialBlurred) )) {
    biffAddf(GAGE, "%s: problem getting KVP buffer", me);
    airMopError(mop); return 1;
  }
  E = 0;
  for (blIdx=0; blIdx<sbp->num; blIdx++) {
    for (kvpIdx=0; kvpIdx<KVP_NUM; kvpIdx++) {
      if (KVP_DGGSM_IDX != kvpIdx) {
        if (!E) E |= nrrdKeyValueAdd(nblur[blIdx], _blurKey[kvpIdx],
                                     blurVal[blIdx].val[kvpIdx]);
      } else {
        /* only need to save dgGoodSigmaMax if it was spatially blurred
           with the discrete gaussian kernel */
        if (spatialBlurred
            && nrrdKernelDiscreteGaussian == kssb->kernel) {
          if (!E) E |= nrrdKeyValueAdd(nblur[blIdx], _blurKey[kvpIdx],
                                       blurVal[blIdx].val[kvpIdx]);
        }
      }
    }
  }
  if (E) {
    biffMovef(GAGE, NRRD, "%s: trouble adding KVPs", me);
    airMopError(mop); return 1;
  }

  airMopOkay(mop);
  return 0;
}

/*
******** gageStackBlurCheck
**
** (docs)
**
*/
int
gageStackBlurCheck(const Nrrd *const nblur[],
                   gageStackBlurParm *sbp,
                   const Nrrd *nin, const gageKind *kind) {
  static const char me[]="gageStackBlurCheck";
  gageShape *shapeOld, *shapeNew;
  blurVal_t *blurVal;
  airArray *mop;
  unsigned int blIdx, kvpIdx;
  NrrdKernelSpec *kssb;

  if (!(nblur && sbp && nin && kind)) {
    biffAddf(GAGE, "%s: got NULL pointer", me);
    return 1;
  }
  mop = airMopNew();
  kssb = nrrdKernelSpecCopy(sbp->kspec);
  airMopAdd(mop, kssb, (airMopper)nrrdKernelSpecNix, airMopAlways);
  if (gageStackBlurParmCheck(sbp)
      || _checkNrrd(NULL, nblur, sbp->num, AIR_TRUE, nin, kind)
      || (!( blurVal = _blurValAlloc(mop, sbp, kssb, nin,
                                     (sbp->needSpatialBlur
                                      ? AIR_TRUE
                                      : AIR_FALSE)) )) ) {
    biffAddf(GAGE, "%s: problem", me);
    airMopError(mop); return 1;
  }

  shapeNew = gageShapeNew();
  airMopAdd(mop, shapeNew, (airMopper)gageShapeNix, airMopAlways);
  if (gageShapeSet(shapeNew, nin, kind->baseDim)) {
    biffAddf(GAGE, "%s: trouble setting up reference shape", me);
    airMopError(mop); return 1;
  }
  shapeOld = gageShapeNew();
  airMopAdd(mop, shapeOld, (airMopper)gageShapeNix, airMopAlways);

  for (blIdx=0; blIdx<sbp->num; blIdx++) {
    if (nin->type != nblur[blIdx]->type) {
      biffAddf(GAGE, "%s: nblur[%u]->type %s != nin type %s\n", me,
               blIdx, airEnumStr(nrrdType, nblur[blIdx]->type),
               airEnumStr(nrrdType, nin->type));
      airMopError(mop); return 1;
    }
    /* check to see if nblur[blIdx] is as expected */
    if (gageShapeSet(shapeOld, nblur[blIdx], kind->baseDim)
        || !gageShapeEqual(shapeOld, "nblur", shapeNew, "nin")) {
      biffAddf(GAGE, "%s: trouble, or nblur[%u] shape != nin shape",
               me, blIdx);
      airMopError(mop); return 1;
    }
    /* see if the KVPs are there with the required values */
    for (kvpIdx=0; kvpIdx<KVP_NUM; kvpIdx++) {
      char *tmpval;
      tmpval = nrrdKeyValueGet(nblur[blIdx], _blurKey[kvpIdx]);
      airMopAdd(mop, tmpval, airFree, airMopAlways);
      if (KVP_DGGSM_IDX != kvpIdx) {
        if (!tmpval) {
          biffAddf(GAGE, "%s: didn't see key \"%s\" in nblur[%u]", me,
                   _blurKey[kvpIdx], blIdx);
          airMopError(mop); return 1;
        }
        if (KVP_SBLUR_IDX == kvpIdx) {
          /* this KVP is handled differently */
          if (!sbp->needSpatialBlur) {
            /* we don't care if it was frequency-domain or spatial-domain
               blurring, so there's no right answer */
            continue;
          }
          /* else we do need spatial domain blurring; so do the check below */
        }
        if (strcmp(tmpval, blurVal[blIdx].val[kvpIdx])) {
          biffAddf(GAGE, "%s: found key[%s] \"%s\" != wanted \"%s\"", me,
                   _blurKey[kvpIdx], tmpval, blurVal[blIdx].val[kvpIdx]);
          airMopError(mop); return 1;
        }
      } else {
        /* KVP_DGGSM_IDX == kvpIdx; handled differently since KVP isn't saved
           when not needed */
        if (!sbp->needSpatialBlur) {
          /* if you don't care about needing spatial blurring, you
             lose the right to care about a difference in dggsm */
          continue;
        }
        if (tmpval) {
          /* therefore it was spatially blurred, so there's a saved DGGSM,
             and we do need spatial blurring, so we compare them */
          if (strcmp(tmpval, blurVal[blIdx].val[kvpIdx])) {
            biffAddf(GAGE, "%s: found key[%s] \"%s\" != wanted \"%s\"", me,
                     _blurKey[kvpIdx], tmpval, blurVal[blIdx].val[kvpIdx]);
            /* HEY: a change in the discrete Gaussian cut-off will result
               in a recomputation of the blurrings, even with FFT-based
               blurring, though cut-off is completely moot then */
            airMopError(mop); return 1;
          }
        }
        continue;
      }
    }
  }

  airMopOkay(mop);
  return 0;
}

int
gageStackBlurGet(Nrrd *const nblur[], int *recomputedP,
                 gageStackBlurParm *sbp,
                 const char *format,
                 const Nrrd *nin, const gageKind *kind) {
  static const char me[]="gageStackBlurGet";
  airArray *mop;
  int recompute;
  unsigned int ii;


  if (!( nblur && sbp && nin && kind )) {
    biffAddf(GAGE, "%s: got NULL pointer", me);
    return 1;
  }
  for (ii=0; ii<sbp->num; ii++) {
    if (!nblur[ii]) {
      biffAddf(GAGE, "%s: nblur[%u] NULL", me, ii);
      return 1;
    }
  }
  if (gageStackBlurParmCheck(sbp)) {
    biffAddf(GAGE, "%s: trouble with blur parms", me);
    return 1;
  }
  mop = airMopNew();

  /* set recompute flag */
  if (!airStrlen(format)) {
    /* no info about files to load, obviously have to recompute */
    if (sbp->verbose) {
      fprintf(stderr, "%s: no file info, must recompute blurrings\n", me);
    }
    recompute = AIR_TRUE;
  } else {
    char *fname, *suberr;
    int firstExists;
    FILE *file;
    /* do have info about files to load, but may fail in many ways */
    fname = AIR_CALLOC(strlen(format) + AIR_STRLEN_SMALL, char);
    if (!fname) {
      biffAddf(GAGE, "%s: couldn't allocate fname", me);
      airMopError(mop); return 1;
    }
    airMopAdd(mop, fname, airFree, airMopAlways);
    /* see if we can get the first file (number 0) */
    sprintf(fname, format, 0);
    firstExists = !!(file = fopen(fname, "r"));
    airFclose(file);
    if (!firstExists) {
      if (sbp->verbose) {
        fprintf(stderr, "%s: no file \"%s\"; will recompute blurrings\n",
                me, fname);
      }
      recompute = AIR_TRUE;
    } else if (nrrdLoadMulti(nblur, sbp->num, format, 0, NULL)) {
      airMopAdd(mop, suberr = biffGetDone(NRRD), airFree, airMopAlways);
      if (sbp->verbose) {
        fprintf(stderr, "%s: will recompute blurrings that couldn't be "
                "read:\n%s\n", me, suberr);
      }
      recompute = AIR_TRUE;
    } else if (gageStackBlurCheck(AIR_CAST(const Nrrd*const*, nblur),
                                  sbp, nin, kind)) {
      airMopAdd(mop, suberr = biffGetDone(GAGE), airFree, airMopAlways);
      if (sbp->verbose) {
        fprintf(stderr, "%s: will recompute blurrings (from \"%s\") "
                "that don't match:\n%s\n", me, format, suberr);
      }
      recompute = AIR_TRUE;
    } else {
      /* else precomputed blurrings could all be read, and did match */
      if (sbp->verbose) {
        fprintf(stderr, "%s: will reuse %u %s pre-blurrings.\n", me,
                sbp->num, format);
      }
      recompute = AIR_FALSE;
    }
  }
  if (recompute) {
    if (gageStackBlur(nblur, sbp, nin, kind)) {
      biffAddf(GAGE, "%s: trouble computing blurrings", me);
      airMopError(mop); return 1;
    }
  }
  if (recomputedP) {
    *recomputedP = recompute;
  }

  airMopOkay(mop);
  return 0;
}

/*
******** gageStackBlurManage
**
** does the work of gageStackBlurGet and then some:
** allocates the array of Nrrds, allocates an array of doubles for scale,
** and saves output if recomputed
*/
int
gageStackBlurManage(Nrrd ***nblurP, int *recomputedP,
                    gageStackBlurParm *sbp,
                    const char *format,
                    int saveIfComputed, NrrdEncoding *enc,
                    const Nrrd *nin, const gageKind *kind) {
  static const char me[]="gageStackBlurManage";
  Nrrd **nblur;
  unsigned int ii;
  airArray *mop;
  int recomputed;

  if (!( nblurP && sbp && nin && kind )) {
    biffAddf(GAGE, "%s: got NULL pointer", me);
    return 1;
  }
  nblur = *nblurP = AIR_CALLOC(sbp->num, Nrrd *);
  if (!nblur) {
    biffAddf(GAGE, "%s: couldn't alloc %u Nrrd*s", me, sbp->num);
    return 1;
  }

  mop = airMopNew();
  airMopAdd(mop, nblurP, (airMopper)airSetNull, airMopOnError);
  airMopAdd(mop, *nblurP, airFree, airMopOnError);
  for (ii=0; ii<sbp->num; ii++) {
    nblur[ii] = nrrdNew();
    airMopAdd(mop, nblur[ii], (airMopper)nrrdNuke, airMopOnError);
  }
  if (gageStackBlurGet(nblur, &recomputed, sbp, format, nin, kind)) {
    biffAddf(GAGE, "%s: trouble getting nblur", me);
    airMopError(mop); return 1;
  }
  if (recomputedP) {
    *recomputedP = recomputed;
  }
  if (recomputed && format && saveIfComputed) {
    NrrdIoState *nio;
    int E;
    E = 0;
    if (enc) {
      if (!enc->available()) {
        biffAddf(GAGE, "%s: requested %s encoding which is not "
                 "available in this build", me, enc->name);
        airMopError(mop); return 1;
      }
      nio = nrrdIoStateNew();
      airMopAdd(mop, nio, (airMopper)nrrdIoStateNix, airMopAlways);
      if (!E) E |= nrrdIoStateEncodingSet(nio, nrrdEncodingGzip);
    } else {
      nio = NULL;
    }
    if (!E) E |= nrrdSaveMulti(format, AIR_CAST(const Nrrd *const *,
                                                nblur),
                               sbp->num, 0, nio);
    if (E) {
      biffMovef(GAGE, NRRD, "%s: trouble saving blurrings", me);
      airMopError(mop); return 1;
    }
  }

  airMopOkay(mop);
  return 0;
}
