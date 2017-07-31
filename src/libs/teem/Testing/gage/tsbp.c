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

#include "teem/gage.h"

/*
** Tests:
** gageStackBlurParmParse
** gageStackBlurParmSprint
** gageHestStackBlurParm
** gageStackBlurParmCompare
*/

/* experimenting with globals for the sake of simplifying tests */
int extraFlag[256];
char *extraParm;
gageStackBlurParm *sbp, *sbq;
airArray *mop;
char *err;
char buff[AIR_STRLEN_LARGE];

void
parseFailOrDie(const char *str) {
  static const char me[]="parseFailOrDie";

  fprintf(stderr, "%s(\"%s\"): ", me, str);
  if (!gageStackBlurParmParse(sbp, extraFlag, &extraParm, str)) {
    fprintf(stderr, "problem did NOT fail!\n");
    airMopError(mop); exit(1);
  }
  /* else good, we got an error */
  airMopAdd(mop, err = biffGetDone(GAGE), airFree, airMopAlways);
  fprintf(stderr, "did get error: %s", err);
  return;
}

void
parseOrDie(const char *str) {
  static const char me[]="parseOrDie";
  int differ;
  char explain[AIR_STRLEN_LARGE]="!explain_not_set!";

  if (gageStackBlurParmParse(sbp, extraFlag, &extraParm, str)
      || gageStackBlurParmSprint(buff, sbp, extraFlag, extraParm)) {
    airMopAdd(mop, err = biffGetDone(GAGE), airFree, airMopAlways);
    fprintf(stderr, "%s: failed to parse \"%s\" or then sprint: %s\n",
            me, str, err);
    airMopError(mop); exit(1);
  }
  if (gageStackBlurParmParse(sbq, extraFlag, &extraParm, str)) {
    airMopAdd(mop, err = biffGetDone(GAGE), airFree, airMopAlways);
    fprintf(stderr, "%s: failed to parse \"%s\": %s\n", me, str, err);
    airMopError(mop); exit(1);
  }
  /* can tweak sbq here to make sure compare is working */
  /*
  sbp->kspec = nrrdKernelSpecNew();
  nrrdKernelSpecParse(sbp->kspec, "cubic:1,1,0");
  sbq->kspec = nrrdKernelSpecNew();
  nrrdKernelSpecParse(sbq->kspec, "gauss:1,5");
  */
  /*
  sbp->bspec = nrrdBoundarySpecNew();
  sbq->bspec = nrrdBoundarySpecNew();
  nrrdBoundarySpecParse(sbp->bspec, "pad:3");
  nrrdBoundarySpecParse(sbq->bspec, "pad:2");
  */
  if (gageStackBlurParmCompare(sbp, "first", sbq, "second",
                               &differ, explain)) {
    airMopAdd(mop, err = biffGetDone(GAGE), airFree, airMopAlways);
    fprintf(stderr, "%s: failed to compare: %s", me, err);
    airMopError(mop); exit(1);
  }
  if (differ) {
    fprintf(stderr, "%s: two sbps from same string differ: %s\n",
            me, explain);
    airMopError(mop); exit(1);
  }
  printf("%s: \"%s\" -> \"%s\"\n", me, str, buff);
  gageStackBlurParmInit(sbp);
  gageStackBlurParmInit(sbq);
  return;
}

static const char *sbpInfo =
  "for testing handling of stack blur parms";
int
main(int argc, const char **argv) {
  /* stock variables */
  const char *me;
  hestOpt *hopt=NULL;
  hestParm *hparm;
  /* variables specific to this program */
  gageStackBlurParm **sbpv;
  unsigned int sbpNum, sbpIdx;

  me = argv[0];
  mop = airMopNew();
  hparm = hestParmNew();
  airMopAdd(mop, hparm, (airMopper)hestParmFree, airMopAlways);
  hestOptAdd(&hopt, "sbp", "sbp0", airTypeOther, 1, -1, &sbpv, NULL,
             "stack blur parms", &sbpNum, NULL, gageHestStackBlurParm);
  hestParseOrDie(hopt, argc-1, argv+1, hparm, me, sbpInfo,
                 AIR_TRUE, AIR_TRUE, AIR_TRUE);
  airMopAdd(mop, hopt, (airMopper)hestOptFree, airMopAlways);
  airMopAdd(mop, hopt, (airMopper)hestParseFree, airMopAlways);

  printf("%s: command line options: ---------- \n", me);
  for (sbpIdx=0; sbpIdx<sbpNum; sbpIdx++) {
    if (gageStackBlurParmSprint(buff, sbpv[sbpIdx], NULL, NULL)) {
      airMopAdd(mop, err = biffGetDone(GAGE), airFree, airMopAlways);
      fprintf(stderr, "%s: problem with sbpv[%u]:%s", me, sbpIdx, err);
      airMopError(mop); return 1;
    }
    printf("%s: sbp[%u] = \"%s\"\n", me, sbpIdx, buff);
  }
  printf("\n");

  sbp = gageStackBlurParmNew();
  airMopAdd(mop, sbp, (airMopper)gageStackBlurParmNix, airMopAlways);
  sbq = gageStackBlurParmNew();
  airMopAdd(mop, sbq, (airMopper)gageStackBlurParmNix, airMopAlways);

  printf("%s: testing various broken strings ---------- \n", me);
  parseFailOrDie("0-8");
  parseFailOrDie("0-8-joe");
  parseFailOrDie("0-n-8.3");
  parseFailOrDie("0-n-8.3-r");
  parseFailOrDie("0-4-8.3-uo");
  parseFailOrDie("0-4-8.3-u+4");
  parseFailOrDie("0-4-8.3-u/k=bingo");
  parseFailOrDie("0-4-8.3-u/k=dg:1,5/b=bingo");
  parseFailOrDie("0-4-8.3-u/k=dg:1,5/b=pad:joe");
  parseFailOrDie("0-4-8.3-u/k=dg:1,5/b=pad:0/v=n");
  parseFailOrDie("0-4-8.3-u/k=dg:1,5/b=pad:0/v=1/s=optiL2");
  parseFailOrDie("0-4-8.3/k=dg:1,5/b=pad:0/v=1/s=optiL2/dggsm=bingo");
  printf("\n");

  printf("%s: testing various okay strings ---------- \n", me);
  parseOrDie("0-4-8.3");
  parseOrDie("0-4-8-o");
  parseOrDie("0-4-8.3-u");
  parseOrDie("0-4-8.3-u1rpn/k=dg:1,5");
  parseOrDie("0-4-8.3-u1rpn/k=dg:1,5/b=pad:42/v=1/dggsm=8");
  printf("\n");

  airMopOkay(mop);
  printf("%s: all okay!\n", me);
  return 0;
}
