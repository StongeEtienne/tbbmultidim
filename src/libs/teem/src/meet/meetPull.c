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

#include "meet.h"

meetPullVol *
meetPullVolNew(void) {
  meetPullVol *ret;

  ret = AIR_CALLOC(1, meetPullVol);
  if (ret) {
    ret->kind = NULL;
    ret->fileName = ret->volName = NULL;
    ret->sbp = NULL;
    ret->leeching = AIR_FALSE;
    ret->derivNormSS = AIR_FALSE;
    ret->recomputedSS = AIR_FALSE;
    ret->derivNormBiasSS = 0.0;
    ret->nin = NULL;
    ret->ninSS = NULL;
  }
  return ret;
}

/*
** DOES use biff
*/
meetPullVol *
meetPullVolCopy(const meetPullVol *mpv) {
  static const char me[]="meetPullVolCopy";
  meetPullVol *ret;
  unsigned int si;
  airArray *mop;

  mop = airMopNew();
  ret = meetPullVolNew();
  airMopAdd(mop, ret, (airMopper)meetPullVolNix, airMopOnError);
  /* HEY: hope this is okay for dynamic kinds */
  ret->kind = mpv->kind;
  ret->fileName = airStrdup(mpv->fileName);
  ret->volName = airStrdup(mpv->volName);
  if (mpv->sbp) {
    ret->sbp = gageStackBlurParmNew();
    if (gageStackBlurParmCopy(ret->sbp, mpv->sbp)) {
      biffMovef(MEET, GAGE, "%s: problem", me);
      airMopError(mop); return NULL;
    }
  }
  ret->leeching = AIR_FALSE;
  ret->derivNormSS = mpv->derivNormSS;
  ret->recomputedSS = AIR_FALSE;
  ret->derivNormBiasSS = mpv->derivNormBiasSS;
  if (mpv->sbp) {
    ret->nin = NULL;
    ret->ninSS = AIR_CALLOC(ret->sbp->num, Nrrd *);
    for (si=0; si<mpv->sbp->num; si++) {
      ret->ninSS[si] = nrrdNew();
      if (nrrdCopy(ret->ninSS[si], mpv->ninSS[si])) {
        biffMovef(MEET, NRRD, "%s: problem with ninSS[%u]", me, si);
        airMopError(mop); return NULL;
      }
    }
  } else {
    ret->nin = nrrdNew();
    if (nrrdCopy(ret->nin, mpv->nin)) {
      biffMovef(MEET, NRRD, "%s: problem with nin", me);
      airMopError(mop); return NULL;
    }
    ret->ninSS = NULL;
  }
  airMopOkay(mop);
  return ret;
}

/*
******** meetPullVolParse
**
** parses a string to extract all the information necessary to create
** the pullVolume (this function originated in Deft/src/main-pull.c)
*/
int
meetPullVolParse(meetPullVol *mpv, const char *_str) {
  static const char me[]="meetPullVolParse";
#define VFMT_SHRT "<fileName>:<kind>:<volName>"
  /* there are other flags and parms but these are the main ones */
#define SFMT "<minScl>-<#smp>-<maxScl>[-n][/k=kss][/b=bspec][/s=smpling]"
#define VFMT_LONG "<fileName>:<kind>:" SFMT ":<volName>"
  char *str, *ctok, *clast=NULL;
  airArray *mop;
  int wantSS;
  unsigned int ctokn;

  if (!(mpv && _str)) {
    biffAddf(MEET, "%s: got NULL pointer", me);
    return 1;
  }
  if (!( str = airStrdup(_str) )) {
    biffAddf(MEET, "%s: couldn't strdup input", me);
    return 1;
  }

  mop = airMopNew();
  airMopAdd(mop, str, airFree, airMopAlways);
  ctokn = airStrntok(str, ":");
  /* An annoying side-effect of putting the blurring kernel specification and
     boundary specification inside the string representation of the
     gageStackBlurParm, is that the colon in dg:1,6" or "pad:10" is now
     confused as a delimiter in the (e.g.) "vol.nrrd:scalar:0-8-5.5:V" string
     representation of meetPullVol, as in
     "vol.nrrd:scalar:0-8-5.5/k=dg:1,6/b=pad:1:V". So we have to be more
     permissive in the number of tokens (hacky) */
  if (!( ctokn >= 3 )) {
    biffAddf(MEET, "%s: didn't get at least 3 \":\"-separated tokens in "
             "\"%s\"; not of form " VFMT_SHRT " or " VFMT_LONG , me, _str);
    airMopError(mop); return 1;
  }
  wantSS = (ctokn > 3);

  ctok = airStrtok(str, ":", &clast);
  if (!(mpv->fileName = airStrdup(ctok))) {
    biffAddf(MEET, "%s: couldn't strdup fileName", me);
    airMopError(mop); return 1;
  }
  airMopAdd(mop, &(mpv->fileName), (airMopper)airSetNull, airMopOnError);
  airMopAdd(mop, mpv->fileName, airFree, airMopOnError);
  ctok = airStrtok(NULL, ":", &clast);
  if (!(mpv->kind = meetConstGageKindParse(ctok))) {
    biffAddf(MEET, "%s: couldn't parse \"%s\" as kind", me, ctok);
    airMopError(mop); return 1;
  }
  if (wantSS) {
    int extraFlag[256]; char *extraParm=NULL, *ptok, *plast;
    unsigned int efi, cti;
    char *sbps;
    /* the hack to make the ":" inside a blurring kernel specification or
       boundary specification be unlike the ":" that delimits the real
       meetPullVol fields */
    sbps = airStrdup(_str); /* to have a buffer big enough */
    airMopAdd(mop, sbps, airFree, airMopAlways);
    strcpy(sbps, "");
    for (cti=0; cti<ctokn-3; cti++) {
      if (cti) {
        strcat(sbps, ":");
      }
      ctok = airStrtok(NULL, ":", &clast);
      strcat(sbps, ctok);
    }
    mpv->sbp = gageStackBlurParmNix(mpv->sbp);
    mpv->sbp = gageStackBlurParmNew();
    if (gageStackBlurParmParse(mpv->sbp, extraFlag, &extraParm, sbps)) {
      biffMovef(MEET, GAGE, "%s: problem parsing sbp from \"%s\"", me, sbps);
      airMopError(mop); return 1;
    }
    mpv->derivNormSS = !!extraFlag['n'];
    extraFlag['n'] = AIR_FALSE;
    for (efi=0; efi<256; efi++) {
      if (extraFlag[AIR_CAST(unsigned char, efi)]) {
        biffAddf(MEET, "%s: got unknown extra flag '%c' in \"%s\"", me,
                 AIR_CAST(char, efi), sbps);
        airMopError(mop); return 1;
      }
    }
    if (extraParm) {
      unsigned int pmi, pmn;
      static const char dnbiase[]="dnbias=";
      airMopAdd(mop, extraParm, airFree, airMopAlways);
      pmn = airStrntok(extraParm, "/");
      for (pmi=0; pmi<pmn; pmi++) {
        ptok = airStrtok(!pmi ? extraParm : NULL, "/", &plast);
        if (strstr(ptok, dnbiase) == ptok) {
          if (1 != sscanf(ptok + strlen(dnbiase), "%lg",
                          &(mpv->derivNormBiasSS))) {
            biffAddf(MEET, "%s: couldn't parse \"%s\" as double in \"%s\"",
                     me, ptok + strlen(dnbiase), ptok);
            airMopError(mop); return 1;
          }
        } else {
          biffAddf(MEET, "%s: got unknown extra parm %s in \"%s\"",
                   me, ptok, extraParm);
          airMopError(mop); return 1;
        }
      }
    }
  } else {
    /* no scale-space stuff wanted */
    mpv->sbp = NULL;
    mpv->ninSS = NULL;
  }
  ctok = airStrtok(NULL, ":", &clast);
  if (!(mpv->volName = airStrdup(ctok))) {
    biffAddf(MEET, "%s: couldn't strdup volName", me);
    airMopError(mop); return 1;
  }
  airMopAdd(mop, &(mpv->volName), (airMopper)airSetNull, airMopOnError);
  airMopAdd(mop, mpv->volName, airFree, airMopOnError);

  if (strchr(ctok, '-')) {
    biffAddf(MEET, "%s: you probably don't want \"-\" in volume name \"%s\"; "
             "forgot last \":\" in scale sampling specification?", me, ctok);
    airMopError(mop); return 1;
  }

  airMopOkay(mop);
  return 0;
}

int
meetHestPullVolParse(void *ptr, char *str, char err[AIR_STRLEN_HUGE]) {
  static const char me[]="meetHestPullVolParse";
  meetPullVol *mpv, **mpvP;
  airArray *mop;

  if (!(ptr && str)) {
    sprintf(err, "%s: got NULL pointer", me);
    return 1;
  }
  mop = airMopNew();
  mpvP = AIR_CAST(meetPullVol **, ptr);
  *mpvP = mpv = meetPullVolNew();
  airMopAdd(mop, mpvP, (airMopper)airSetNull, airMopOnError);
  airMopAdd(mop, mpv, (airMopper)meetPullVolNix, airMopOnError);
  if (meetPullVolParse(mpv, str)) {
    char *ler;
    airMopAdd(mop, ler = biffGetDone(MEET), airFree, airMopOnError);
    airStrcpy(err, AIR_STRLEN_HUGE, ler);
    airMopError(mop);
    return 1;
  }
  airMopOkay(mop);
  return 0;
}

/*
******** meetPullVolNix
**
** this frees stuff allocated meetPullVolParse and meetPullVolLoadMulti
*/
meetPullVol *
meetPullVolNix(meetPullVol *mpv) {

  if (mpv) {
    if (!mpv->leeching && mpv->nin) {
      nrrdNuke(mpv->nin);
    }
    if (mpv->sbp) {
      unsigned int ssi;
      if (mpv->ninSS) {
        /* need this check because the mpv may not have benefitted
           from meetPullVolLoadMulti, so it might be incomplete */
        for (ssi=0; ssi<mpv->sbp->num; ssi++) {
          if (!mpv->leeching) {
            nrrdNuke(mpv->ninSS[ssi]);
          }
        }
        airFree(mpv->ninSS);
      }
      gageStackBlurParmNix(mpv->sbp);
    }
    airFree(mpv->fileName);
    airFree(mpv->volName);
    airFree(mpv);
  }
  return NULL;
}

hestCB
_meetHestPullVol = {
  sizeof(meetPullVol *),
  "meetPullVol",
  meetHestPullVolParse,
  (airMopper)meetPullVolNix,
};

hestCB *
meetHestPullVol = &_meetHestPullVol;

/*
******** meetPullVolLeechable
**
** indicates whether lchr can leech from orig (saved in *can), and if not,
** explanation is saved in explain (if non-NULL)
**
** always uses biff
*/
int
meetPullVolLeechable(const meetPullVol *lchr,
                     const meetPullVol *orig,
                     int *can, char explain[AIR_STRLEN_LARGE]) {
  static const char me[]="meetPullVolLeechable";
  char subexplain[AIR_STRLEN_LARGE];

  if (!( lchr && orig && can )) {
    biffAddf(MEET, "%s: got NULL pointer (%p %p %p)", me, lchr, orig, can);
    return 1;
  }
  /* can leech, if not reading from stdin */
  *can = !!strcmp(orig->fileName, "-");
  if (!*can) {
    if (explain) {
      sprintf(explain, "original loaded from stdin");
    }
    return 0;
  }
  /* can, if coming from same file */
  *can = !strcmp(orig->fileName, lchr->fileName);
  if (!*can) {
    if (explain) {
      sprintf(explain, "not from same file (\"%s\" vs \"%s\")\n",
              lchr->fileName, orig->fileName);
    }
    return 0;
  }
  /* can, if same kind */
  *can = (orig->kind == lchr->kind);
  if (!*can) {
    if (explain) {
      sprintf(explain, "not same kind (%s vs %s)\n",
              lchr->kind->name, orig->kind->name);
    }
    return 0;
  }
  /* can, if both using or both not using scale-space */
  *can = (!!lchr->sbp == !!orig->sbp);
  if (!*can) {
    if (explain) {
      sprintf(explain, "not agreeing on use of scale-space (%s vs %s)\n",
              lchr->sbp ? "yes" : "no", orig->sbp ? "yes" : "no");
    }
    return 0;
  }
  if (orig->sbp) {
    int differ;
    if (gageStackBlurParmCompare(lchr->sbp, "potential leecher",
                                 orig->sbp, "original",
                                 &differ, subexplain)) {
      biffAddf(MEET, "%s: problem comparing sbps", me);
      return 1;
    }
    if (differ) {
      if (explain) {
        sprintf(explain, "different uses of scale-space: %s", subexplain);
      }
      *can = AIR_FALSE;
      return 0;
    }
  }
  /* DO allow difference in derivNormSS (the main reason for leeching),
     as well as derivNormBiasSS */
  /* no differences so far */
  *can = AIR_TRUE;
  return 0;
}

void
meetPullVolLeech(meetPullVol *vol,
                 const meetPullVol *volPrev) {

  if (vol && volPrev) {
    vol->nin = volPrev->nin;
    if (vol->sbp) {
      unsigned int ni;
      /* have to allocate ninSS here; in volPrev it was probably allocated
         by gageStackBlurManage */
      vol->ninSS = AIR_CALLOC(vol->sbp->num, Nrrd *);
      for (ni=0; ni<vol->sbp->num; ni++) {
        vol->ninSS[ni] = volPrev->ninSS[ni];
      }
    }
    vol->leeching = AIR_TRUE;
  }
  return;
}

/*
** This is kind of a sad function. The big re-write of gageStackBlurParm in
** late August 2013 was motivated by the frustration of how there was no
** centralized and consistent way of representing (by text or by command-line
** options) all the things that determine scale-space "stack" creation.
** Having re-organized gageStackBlurParm, the meetPullVol was re-organized to
** include one inside, which is clearly better than the previous reduplication
** of the stack blur parms inside the meetPullVol. Parsing the meetPullVol
** from the command-line (as in done in puller) highlights the annoying fact
** that hest wants to be the origin of information: you can't have hest
** supplement existing information with whatever it learns from the
** command-line, especially when hest is parsing 1 or more of something, and
** especially when the existing information would be coming from other
** command-line arguments.
**
** So, this sad function says, "ok all you meetPullVol parsed from the
** command-line: if you don't already have a boundary or a kernel set, here's
** one to use". What makes it sad is how the whole point of the
** gageStackBlurParm re-org was that knowledge about the internals of the
** gageStackBlurParm was now going to be entirely localized in gage. But here
** we are listing off two of its fields as parameters to this function, which
** means its API might change the next time the gageStackBlurParm is updated.
**
** To help keep track of what info was actually used, *kssSetP and *bspSetP
** (if non-NULL) are set to the number of kernel and boundary specs that are
** "finished" in this way.
*/
int
meetPullVolStackBlurParmFinishMulti(meetPullVol **mpv, unsigned int mpvNum,
                                    unsigned int *kssSetP,
                                    unsigned int *bspSetP,
                                    const NrrdKernelSpec *kssblur,
                                    const NrrdBoundarySpec *bspec) {
  static const char me[]="meetPullVolStackBlurParmFinishMulti";
  unsigned int ii, kssSet, bspSet;

  if (!mpv || !mpvNum) {
    biffAddf(MEET, "%s: got NULL mpv (%p) or 0 mpvNum (%u)",
             me, AIR_VOIDP(mpv), mpvNum);
    return 1;
  }
  kssSet = bspSet = 0;
  for (ii=0; ii<mpvNum; ii++) {
    if (kssblur && mpv[ii]->sbp && !(mpv[ii]->sbp->kspec)) {
      if (gageStackBlurParmKernelSet(mpv[ii]->sbp, kssblur)) {
        biffMovef(MEET, GAGE, "%s: trouble w/ kernel on mpv[%u]", me, ii);
        return 1;
      }
      kssSet++;
    }
    if (bspec && mpv[ii]->sbp && !(mpv[ii]->sbp->bspec)) {
      if (gageStackBlurParmBoundarySpecSet(mpv[ii]->sbp, bspec)) {
        biffMovef(MEET, GAGE, "%s: trouble w/ boundary on mpv[%u]", me, ii);
        return 1;
      }
      bspSet++;
    }
  }
  if (kssSetP) {
    *kssSetP = kssSet;
  }
  if (bspSetP) {
    *bspSetP = bspSet;
  }
  return 0;
}

/*
******** meetPullVolLoadMulti
**
** at this point the only per-pullVolume information required for
** loading/creating the volumes, which isn't already in the
** meetPullVol, is the cachePath, so that is passed explicitly.
*/
int
meetPullVolLoadMulti(meetPullVol **mpv, unsigned int mpvNum,
                     char *cachePath, int verbose) {
  static const char me[]="meetPullVolLoadMulti";
  unsigned int mpvIdx;
  airArray *mop;
  meetPullVol *vol;

  if (!( mpv && cachePath)) {
    biffAddf(MEET, "%s: got NULL pointer", me);
    return 1;
  }
  mop = airMopNew();
  for (mpvIdx=0; mpvIdx<mpvNum; mpvIdx++) {
    unsigned int pvi;
    int leechable;
    char explain[AIR_STRLEN_LARGE];
    vol = mpv[mpvIdx];
    for (pvi=0; pvi<mpvIdx; pvi++) {
      if (meetPullVolLeechable(vol, mpv[pvi], &leechable, explain)) {
        biffAddf(MEET, "%s: problem testing leechable(v[%u]->v[%u])",
                 me, mpvIdx, pvi);
        return 1;
      }
      if (leechable) {
        meetPullVolLeech(vol, mpv[pvi]);
        break; /* prevent a chain of leeching */
      } else {
        if (verbose) {
          fprintf(stderr, "%s: mpv[%u] cannot leech mpv[%u]: %s\n", me,
                  mpvIdx, pvi, explain);
        }
      }
    }
    if (pvi < mpvIdx) {
      /* we leeched this one, move on */
      if (verbose) {
        fprintf(stderr, "%s: vspec[%u] (%s) leeching off vspec[%u] (%s)\n",
                me, mpvIdx, vol->volName, pvi, mpv[pvi]->volName);
      }
      continue;
    }
    /* else we're not leeching */
    vol->leeching = AIR_FALSE;
    vol->nin = nrrdNew();
    airMopAdd(mop, &(vol->nin), (airMopper)airSetNull, airMopOnError);
    airMopAdd(mop, vol->nin, (airMopper)nrrdNuke, airMopOnError);
    if (nrrdLoad(vol->nin, vol->fileName, NULL)) {
      biffMovef(MEET, NRRD, "%s: trouble loading mpv[%u]->nin (\"%s\")",
                me, mpvIdx, vol->volName);
      airMopError(mop); return 1;
    }
    if (vol->sbp) {
      char formatSS[AIR_STRLEN_LARGE];
      sprintf(formatSS, "%s/%s-%%03u-%03u.nrrd",
              cachePath, vol->volName, vol->sbp->num);
      if (verbose) {
        fprintf(stderr, "%s: managing %s ... \n", me, formatSS);
      }
      if (gageStackBlurManage(&(vol->ninSS), &(vol->recomputedSS), vol->sbp,
                              formatSS, AIR_TRUE, NULL,
                              vol->nin, vol->kind)) {
        biffMovef(MEET, GAGE, "%s: trouble getting volume stack (\"%s\")",
                  me, formatSS);
        airMopError(mop); return 1;
      }
      if (verbose) {
        fprintf(stderr, "%s: ... done\n", me);
      }
    }
  }
  airMopOkay(mop);
  return 0;
}

/*
******** meetPullVolAddMulti
**
** the spatial (k00, k11, k22) and scale (kSSrecon) reconstruction
** kernels are not (yet) part of the meetPullVol, so have to be passed
** in here
*/
int
meetPullVolAddMulti(pullContext *pctx,
                    meetPullVol **mpv, unsigned int mpvNum,
                    const NrrdKernelSpec *k00,
                    const NrrdKernelSpec *k11,
                    const NrrdKernelSpec *k22,
                    const NrrdKernelSpec *kSSrecon) {
  static const char me[]="meetPullVolAddMulti";
  unsigned int mpvIdx;

  if (!( pctx && mpv )) {
    biffAddf(MEET, "%s: got NULL pointer", me);
    return 1;
  }
  for (mpvIdx=0; mpvIdx<mpvNum; mpvIdx++) {
    meetPullVol *vol;
    int E;
    vol = mpv[mpvIdx];
    if (!vol->sbp) {
      E = pullVolumeSingleAdd(pctx, vol->kind, vol->volName,
                              vol->nin, k00, k11, k22);
    } else {
      E = pullVolumeStackAdd(pctx, vol->kind, vol->volName, vol->nin,
                             AIR_CAST(const Nrrd *const *,
                                      vol->ninSS),
                             vol->sbp->sigma, vol->sbp->num,
                             vol->derivNormSS, vol->derivNormBiasSS,
                             k00, k11, k22, kSSrecon);
    }
    if (E) {
      biffMovef(MEET, PULL, "%s: trouble adding volume %u/%u (\"%s\")",
                me, mpvIdx, mpvNum, vol->volName);
      return 1;
    }
  }

  return 0;
}

meetPullInfo *
meetPullInfoNew(void) {
  meetPullInfo *ret;

  ret = AIR_CALLOC(1, meetPullInfo);
  if (ret) {
    ret->info = 0;
    ret->source = pullSourceUnknown;
    ret->prop = pullPropUnknown;
    ret->constraint = AIR_FALSE;
    ret->volName = ret->itemStr = NULL;
    ret->zero = ret->scale = AIR_NAN;
  }
  return ret;
}

meetPullInfo *
meetPullInfoNix(meetPullInfo *minf) {

  if (minf) {
    airFree(minf->volName);
    airFree(minf->itemStr);
    free(minf);
  }
  return NULL;
}

static int
zeroScaleSet(meetPullInfo *minf, int haveZS, char **lastP) {
  static const char me[]="_zeroScaleSet";
  char *tok;

  if (haveZS) {
    tok = airStrtok(NULL, ":", lastP);
    if (1 != sscanf(tok, "%lf", &(minf->zero))) {
      biffAddf(MEET, "%s: couldn't parse %s as zero (double)", me, tok);
      return 1;
    }
    tok = airStrtok(NULL, ":", lastP);
    if (1 != sscanf(tok, "%lf", &(minf->scale))) {
      biffAddf(MEET, "%s: couldn't parse %s as scale (double)", me, tok);
      return 1;
    }
  } else {
    minf->zero = minf->scale = AIR_NAN;
  }
  return 0;
}

int
meetPullInfoParse(meetPullInfo *minf, const char *_str) {
  static const char me[]="meetPullInfoParse";
#define IFMT_GAGE "<info>[-c]:<volname>:<item>[:<zero>:<scale>]"
#define IFMT_PROP "<info>:prop=<prop>[:<zero>:<scale>]"
#define PROP_PREFIX "prop="   /* has to end with = */
  char *str, *tok, *last=NULL, *iflags;
  airArray *mop;
  int haveZS, source;

  if (!(minf && _str)) {
    biffAddf(MEET, "%s: got NULL pointer", me);
    return 1;
  }
  if ( (3 == airStrntok(_str, ":") || 5 == airStrntok(_str, ":"))
       && 1 == airStrntok(_str, "=") ) {
    source = pullSourceGage;
    haveZS = (5 == airStrntok(_str, ":"));
  } else if ( (2 == airStrntok(_str, ":") || 4 == airStrntok(_str, ":"))
              && 2 == airStrntok(_str, "=") ) {
    source = pullSourceProp;
    haveZS = (4 == airStrntok(_str, ":"));
  } else {
    biffAddf(MEET, "%s: \"%s\" not of form " IFMT_GAGE " or " IFMT_PROP,
             me, _str);
    return 1;
  }

  mop = airMopNew();
  if (!( str = airStrdup(_str) )) {
    biffAddf(MEET, "%s: couldn't strdup input", me);
    return 1;
  }
  airMopAdd(mop, str, airFree, airMopAlways);

  minf->source = source;
  if (pullSourceGage == source) {
    tok = airStrtok(str, ":", &last);
    iflags = strchr(tok, '-');
    if (iflags) {
      *iflags = '\0';
      iflags++;
    }
    if (!(minf->info = airEnumVal(pullInfo, tok))) {
      biffAddf(MEET, "%s: couldn't parse \"%s\" as %s",
               me, tok, pullInfo->name);
      airMopError(mop); return 1;
    }
    if (iflags) {
      if (strchr(iflags, 'c')) {
        minf->constraint = AIR_TRUE;
      }
    }
    tok = airStrtok(NULL, ":", &last);
    airFree(minf->volName);
    minf->volName = airStrdup(tok);
    airMopAdd(mop, minf->volName, airFree, airMopOnError);
    tok = airStrtok(NULL, ":", &last);
    airFree(minf->itemStr);
    minf->itemStr = airStrdup(tok);
    airMopAdd(mop, minf->itemStr, airFree, airMopOnError);
    if (zeroScaleSet(minf, haveZS, &last)) {
      biffAddf(MEET, "%s: couldn't parse zero or scale",  me);
      airMopError(mop); return 1;
    }
  } else if (pullSourceProp == source) {
    /* "<info>:prop=<prop>[:<zero>:<scale>]" */
    tok = airStrtok(str, ":", &last);
    if (!(minf->info = airEnumVal(pullInfo, tok))) {
      biffAddf(MEET, "%s: couldn't parse \"%s\" as %s",
               me, tok, pullInfo->name);
      airMopError(mop); return 1;
    }
    tok = airStrtok(NULL, ":", &last);
    if (strncmp(PROP_PREFIX, tok, strlen(PROP_PREFIX))) {
      biffAddf(MEET, "%s: property info didn't start with %s",
               me, PROP_PREFIX);
    }
    tok += strlen(PROP_PREFIX);
    if (!(minf->prop = airEnumVal(pullProp, tok))) {
      biffAddf(MEET, "%s: couldn't parse \"%s\" as %s",
               me, tok, pullProp->name);
      airMopError(mop); return 1;
    }
    if (zeroScaleSet(minf, haveZS, &last)) {
      biffAddf(MEET, "%s: couldn't parse zero or scale",  me);
      airMopError(mop); return 1;
    }
  } else {
    biffAddf(MEET, "%s: sorry, source %s not handled",
             me, airEnumStr(pullSource, source));
    airMopError(mop); return 1;
  }

  airMopOkay(mop);
  return 0;
}

int
meetHestPullInfoParse(void *ptr, char *str, char err[AIR_STRLEN_HUGE]) {
  static const char me[]="meetHestPullInfoParse";
  airArray *mop;
  meetPullInfo **minfP, *minf;

  if (!(ptr && str)) {
    sprintf(err, "%s: got NULL pointer", me);
    return 1;
  }
  mop = airMopNew();
  minfP = AIR_CAST(meetPullInfo **, ptr);
  *minfP = minf = meetPullInfoNew();
  airMopAdd(mop, minfP, (airMopper)airSetNull, airMopOnError);
  airMopAdd(mop, minf, (airMopper)meetPullInfoNix, airMopOnError);
  if (meetPullInfoParse(minf, str)) {
    char *ler;
    airMopAdd(mop, ler = biffGetDone(MEET), airFree, airMopOnError);
    airStrcpy(err, AIR_STRLEN_HUGE, ler);
    airMopError(mop);
    return 1;
  }
  airMopOkay(mop);
  return 0;
}

hestCB
_meetHestPullInfo = {
  sizeof(meetPullInfo *),
  "meetPullInfo",
  meetHestPullInfoParse,
  (airMopper)meetPullInfoNix
};

hestCB *
meetHestPullInfo = &_meetHestPullInfo;

int
meetPullInfoAddMulti(pullContext *pctx,
                     meetPullInfo **minf, unsigned int minfNum) {
  static const char me[]="meetPullInfoAddMulti";
  const pullVolume *vol;
  unsigned int ii;
  airArray *mop;

  if (!( pctx && minf )) {
    biffAddf(MEET, "%s: got NULL pointer", me);
    return 1;
  }

  mop = airMopNew();
  for (ii=0; ii<minfNum; ii++) {
    pullInfoSpec *ispec;
    ispec = pullInfoSpecNew();
    airMopAdd(mop, ispec, (airMopper)pullInfoSpecNix, airMopOnError);
    ispec->volName = airStrdup(minf[ii]->volName);
    ispec->source = minf[ii]->source;
    ispec->info = minf[ii]->info;
    ispec->prop = minf[ii]->prop;
    ispec->zero = minf[ii]->zero;
    ispec->scale = minf[ii]->scale;
    ispec->constraint = minf[ii]->constraint;
    /* the item is the one thing that takes some work to recover;
       we need to find the volume and find the item from its kind->enm */
    if (pullSourceGage == ispec->source) {
      if (!( vol = pullVolumeLookup(pctx, minf[ii]->volName) )) {
        biffMovef(MEET, PULL, "%s: can't find volName \"%s\" for minf[%u]",
                  me, minf[ii]->volName, ii);
        airMopError(mop); return 1;
      }
      if (!( ispec->item = airEnumVal(vol->kind->enm, minf[ii]->itemStr))) {
        biffAddf(MEET, "%s: can't parse \"%s\" as item of %s kind (minf[%u])\n",
                 me, minf[ii]->itemStr, vol->kind->name, ii);
        airMopError(mop); return 1;
      }
    }
    if (pullInfoSpecAdd(pctx, ispec)) {
      biffMovef(MEET, PULL, "%s: trouble adding ispec from minf[%u]", me, ii);
      airMopError(mop); return 1;
    }
    /* else added the ispec okay. If we have an error with a different
       ispec later in this loop, who's job is it to free up the ispecs
       that have been added successfully?  In teem/src/bin/puller, that
       is done by pullContextNix. So we now extricate ourself from the
       business of freeing this ispec in case of error; one of the few
       times that airMopSub is really needed */
    airMopSub(mop, ispec, (airMopper)pullInfoSpecNix);
  }

  airMopOkay(mop);
  return 0;
}
