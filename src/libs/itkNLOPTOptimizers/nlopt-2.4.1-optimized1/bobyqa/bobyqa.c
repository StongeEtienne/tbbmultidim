/*
   Original Fortran code by Powell (2009).  Converted via v2c,
   cleaned up, and incorporated into NLopt by S. G. Johnson (2009).
   See README. */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "bobyqa.h"
/***************************************************************
/* CPB: two-dimensioan array allocator and memory free
/*      could be put in bobyqa.h
/*      set to type double for use here
/***************************************************************/
size_t cache_line_size() {
    FILE * p = 0;

static size_t cachelinesize = 0;
/*  
 *  The check for cache line size may be dependent on Linux vendor.
 */
	if ( cachelinesize==0 )
	{
	    p = fopen("/sys/devices/system/cpu/cpu0/cache/index0/coherency_line_size", "r");
	    unsigned int i = 0;
	    if (p) {
		fscanf(p, "%d", &i);
		fclose(p);
	    }
	    cachelinesize=(size_t)i;
	}
    return cachelinesize;
}

double **Allocate2DArray( int nRows, int nCols)
{
    //(step 0) Determine cacheline size and adjust nCols to use full lines
    size_t line_size = cache_line_size();
    int doubles_per_line = line_size / sizeof(double);
    int rem = nCols % doubles_per_line;
    if (rem != 0) nCols += (doubles_per_line - rem);
    
    //(step 1) allocate memory for array of elements of column
    double **ppi = malloc(sizeof(double *)*nRows);
    
    //(step 2) allocate memory for array of elements of each row
    double *curPtr = _mm_malloc(sizeof(double) * nRows * nCols, 64);
    
    // Now point the pointers in the right place
    int i;
    for( i = 0; i < nRows; ++i)
    {
        *(ppi + i) = curPtr;
        curPtr += nCols;
    }
    return ppi;
}

void Free2DArray(double** Array)
{
    _mm_free(*Array);
    free(Array);
}
/***************************************************************/

typedef double (*bobyqa_func)(int n, const double *x, void *func_data);

#define MIN2(a,b) ((a) <= (b) ? (a) : (b))
#define MAX2(a,b) ((a) >= (b) ? (a) : (b))
#define IABS(x) ((x) < 0 ? -(x) : (x))

#define U(n) ((unsigned) (n))

static void update_(int n, int npt, int ndim, double **bmat, 
    double **zmat, double *vlag, double beta, 
    double denom, int knew)
{
    /* System generated locals */
    double d__1, d__2, d__3;

    /* Local variables */
    int i, j, k, jl, jp;
    double one, tau, temp;
    int nptm;
    double zero, alpha, tempa, tempb, ztest;
    double *w = _mm_malloc(sizeof(double) * (ndim + 8), 64);


/*     The arrays BMAT and ZMAT are updated, as required by the new position */
/*     of the interpolation point that has the index KNEW. The vector VLAG has */
/*     N+NPT components, set on entry to the first NPT and last N components */
/*     of the product Hw in equation (4.11) of the Powell (2006) paper on */
/*     NEWUOA. Further, BETA is set on entry to the value of the parameter */
/*     with that name, and DENOM is set to the denominator of the updating */
/*     formula. Elements of ZMAT may be treated as zero if their moduli are */
/*     at most ZTEST. The first NDIM elements of W are used for working space. */

__assume_aligned(vlag, 64);

    /* Function Body */
    one = 1.;
    zero = 0.;
    nptm = npt - n - 1;
    ztest = zero;
    
   
    for (j = 0; j < nptm; ++j) {
        for (k = 0; k < npt; ++k) {
/* L10: */
            ztest = MAX2(ztest,fabs(zmat[j][k]));
        }
    }
    ztest *= 1e-20;

    /*     Apply the rotations that put zeros in the KNEW-th row of ZMAT. */

    jl = 1;
  
    for (j = 1; j < nptm; ++j) {	/* CPB: original was DO J=2,NPTM */
        if (fabs(zmat[j][knew]) > ztest) {
            d__1 = zmat[0][knew];
            d__2 = zmat[j][knew];
            temp = sqrt(d__1 * d__1 + d__2 * d__2);
            tempa = zmat[0][knew] / temp;
            tempb = zmat[j][knew] / temp;
 
            for (i = 0; i < npt; ++i) {
                temp = tempa * zmat[0][i] + tempb * zmat[j][i];
                zmat[j][i]  = tempa * zmat[j][i] - tempb * zmat[0][i];
/* L20: */
                zmat[0][i] = temp;
            }
        }
        zmat[j][knew] = zero;
/* L30: */
    }

/*     Put the first NPT components of the KNEW-th column of HLAG into W, */
/*     and calculate the parameters of the updating formula. */


    for (i = 0; i < npt; ++i) {
        w[i] = zmat[0][knew] * zmat[0][i];
/* L40: */
    }
    alpha = w[knew];
    tau = vlag[knew];
    vlag[knew] -= one;

/*     Complete the updating of ZMAT. */

    temp = sqrt(denom);
    tempb = zmat[0][knew] / temp;
    tempa = tau / temp;

    for (i = 0; i < npt; ++i) {
/* L50: */
        zmat[0][i] = tempa * zmat[0][i] - tempb * vlag[i];
    }

/*     Finally, update the matrix BMAT. */


    for (j = 0; j < n; ++j) {
        jp = npt + j;
        w[jp] = bmat[j][knew];
        tempa = (alpha * vlag[jp] - tau * w[jp]) / denom;
        tempb = (-(beta) * w[jp] - tau * vlag[jp]) / denom;

        for (i = 0; i <= jp; ++i) {
            bmat[j][i] = bmat[j][i] + tempa * vlag[i] + tempb * w[i];
            if (i >= npt) {
                bmat[(i - npt)][jp] = bmat[j][i];
            }
/* L60: */
        }
    }
    _mm_free(w);
} /* update_ */

static nlopt_result rescue_(int n, int npt, const double *xl, const double *xu, 
            /* int *maxfun */
            nlopt_stopping *stop,
            bobyqa_func calfun, void *calfun_data,
    double *xbase, double **xpt, double *fval, double *xopt, double *gopt,
    double *hq, double *pq, int ndim, double **bmat, 
    double **zmat, double *sl, double *su, /* int nf,  */
    double delta, int *kopt, double *vlag, double **ptsaux, 
    double *ptsid)
{
    /* System generated locals */
    double d__1, d__2, d__3, d__4;

    /* Local variables */
    double f;
    int i, j, k, ih, jp, ip, iq, np, iw;
    double xp, xq, den;
    int ihp;
    double one;
    int ihq, jpn, kpt;
    double sum, diff, half, beta;
    int kold;
    double winc;
    int nrem, knew;
    double temp, bsum;
    int nptm;
    double zero, hdiag, fbase, sfrac, denom, vquad, sumpq;
    double dsqmin, distsq, vlmxsq;
    double *V_distsq, *V_sum, *V_hdiag;
    double *w = _mm_malloc(sizeof(double) * (ndim+npt+8), 64); 
    nlopt_result rCode = NLOPT_SUCCESS;


/*     The arguments N, NPT, XL, XU, MAXFUN, XBASE, XPT, FVAL, XOPT, */
/*       GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU have the same meanings as */
/*       the corresponding arguments of BOBYQB on the entry to RESCUE. */
/*     NF is maintained as the number of calls of CALFUN so far, except that */
/*       NF is set to -1 if the value of MAXFUN prevents further progress. */
/*     KOPT is maintained so that FVAL(KOPT) is the least calculated function */
/*       value. Its correct value must be given on entry. It is updated if a */
/*       new least function value is found, but the corresponding changes to */
/*       XOPT and GOPT have to be made later by the calling program. */
/*     DELTA is the current trust region radius. */
/*     VLAG is a working space vector that will be used for the values of the */
/*       provisional Lagrange functions at each of the interpolation points. */
/*       They are part of a product that requires VLAG to be of length NDIM. */
/*     PTSAUX is also a working space array. For J=1,2,...,N, PTSAUX(1,J) and */
/*       PTSAUX(2,J) specify the two positions of provisional interpolation */
/*       points when a nonzero step is taken along e_J (the J-th coordinate */
/*       direction) through XBASE+XOPT, as specified below. Usually these */
/*       steps have length DELTA, but other lengths are chosen if necessary */
/*       in order to satisfy the given bounds on the variables. */
/*     PTSID is also a working space array. It has NPT components that denote */
/*       provisional new positions of the original interpolation points, in */
/*       case changes are needed to restore the linear independence of the */
/*       interpolation conditions. The K-th point is a candidate for change */
/*       if and only if PTSID(K) is nonzero. In this case let p and q be the */
/*       int parts of PTSID(K) and (PTSID(K)-p) multiplied by N+1. If p */
/*       and q are both positive, the step from XBASE+XOPT to the new K-th */
/*       interpolation point is PTSAUX(1,p)*e_p + PTSAUX(1,q)*e_q. Otherwise */
/*       the step is PTSAUX(1,p)*e_p or PTSAUX(2,q)*e_q in the cases q=0 or */
/*       p=0, respectively. */
/*     The first NDIM+NPT elements of the array W are used for working space. */
/*     The final elements of BMAT and ZMAT are set in a well-conditioned way */
/*       to the values that are appropriate for the new interpolation points. */
/*     The elements of GOPT, HQ and PQ are also revised to the values that are */
/*       appropriate to the final quadratic model. */

__assume_aligned(xbase, 64);
__assume_aligned(fval, 64);
__assume_aligned(xopt, 64);
__assume_aligned(gopt, 64);
__assume_aligned(hq, 64);
__assume_aligned(pq, 64);
__assume_aligned(sl, 64);
__assume_aligned(su, 64);
__assume_aligned(vlag, 64);
__assume_aligned(ptsid, 64);

/*     Set some constants. */

    /* Function Body */
    half = .5;
    one = 1.;
    zero = 0.;
    np = n + 1;
    sfrac = half / (double) np;
    nptm = npt - np;
    V_distsq = _mm_malloc(sizeof(double) * (npt + 8), 64);
    V_sum = _mm_malloc(sizeof(double) * (npt + 8), 64);
    V_hdiag = _mm_malloc(sizeof(double) * (npt + 8), 64);

/*     Shift the interpolation points so that XOPT becomes the origin, and set */
/*     the elements of ZMAT to zero. The value of SUMPQ is required in the */
/*     updating of HQ below. The squares of the distances from XOPT to the */
/*     other interpolation points are set at the end of W. Increments of WINC */
/*     may be added later to these squares to balance the consideration of */
/*     the choice of point that is going to become current. */

    sumpq = zero;
    winc = zero;
    V_distsq[0:npt] = zero;

    for (j = 0; j < n; ++j) {
        for (k = 0; k < npt; ++k) {
            xpt[j][k] -= xopt[j];
/* L10: */
            d__1 = xpt[j][k];
            V_distsq[k] += d__1 * d__1;
        }
    }
    for (k = 0; k < npt; ++k) {
        sumpq += pq[k];
        w[ndim + k] = V_distsq[k];
        winc = MAX2(winc,V_distsq[k]);
/* L20: */
    }

//    for (j = 0; j < nptm; ++j) {
//        for (k = 0; k < npt; ++k)
//            zmat[j][k] = zero;
//    }

    zmat[0:nptm][0:npt] = zero;

/*     Update HQ so that HQ and PQ define the second derivatives of the model */
/*     after XBASE has been shifted to the trust region centre. */

    ih = 0;
    for (j = 0; j < n; ++j) {
        w[j] = half * sumpq * xopt[j];

        for (k = 0; k < npt; ++k) {
/* L30: */
            w[j] += pq[k] * xpt[j][k];
        }

        for (i = 0; i <= j; ++i) {
            ++ih;
/* L40: */
            hq[ih-1] = hq[ih-1] + w[i] * xopt[j] + w[j] * xopt[i];
        }
    }

/*     Shift XBASE, SL, SU and XOPT. Set the elements of BMAT to zero, and */
/*     also set the elements of PTSAUX. */

    for (j = 0; j < n; ++j) {
        xbase[j] += xopt[j];
        sl[j] -= xopt[j];
        su[j] -= xopt[j];
        xopt[j] = zero;
        ptsaux[0][j] = MIN2(delta, su[j]);
        ptsaux[1][j] = MAX2(-(delta), sl[j]);
        if (ptsaux[0][j] + ptsaux[1][j] < zero) {
            temp = ptsaux[0][j];
            ptsaux[0][j] = ptsaux[1][j];
            ptsaux[1][j] = temp;
        }
        if (fabs(ptsaux[1][j]) < half * (fabs(ptsaux[0][j]))) {
            ptsaux[1][j] = half * ptsaux[0][j];
        }

        for (i = 0; i < ndim; ++i) {
/* L50: */
            bmat[j][i] = zero;
        }
    }
    fbase = fval[*kopt];

/*     Set the identifiers of the artificial interpolation points that are */
/*     along a coordinate direction from XOPT, and set the corresponding */
/*     nonzero elements of BMAT and ZMAT. */

    ptsid[0] = sfrac;

    for (j = 0; j < n; ++j) {
        jp = j + 1;
        jpn = jp + n;
        ptsid[jp] = (double) (j+1) + sfrac; /* CPB: j is counter here */
        if (jpn <= (npt - 1)) {
            ptsid[jpn] = (double) (j+1) / (double) np + sfrac;
            temp = one / (ptsaux[0][j] - ptsaux[1][j]);
            bmat[j][jp] = -temp + one / ptsaux[0][j];
            bmat[j][jpn] = temp + one / ptsaux[1][j];
            bmat[j][0] = -bmat[j][jp] - bmat[j][jpn];
            zmat[j][0] = sqrt(2.) / fabs(ptsaux[0][j] * ptsaux[1][j]);
            zmat[j][jp] = zmat[j][0] * ptsaux[1][j] * temp;
            zmat[j][jpn] = -zmat[j][0] * ptsaux[0][j] * temp;
        } 
        else {
            bmat[j][0] = -one / ptsaux[0][j];
            bmat[j][jp] = one / ptsaux[0][j];
            bmat[j][j + npt] = -half * (ptsaux[0][j] * ptsaux[0][j]);
        }
/* L60: */
    }

/*     Set any remaining identifiers with their nonzero elements of ZMAT. */

    if (npt >= n + np) {

        for (k = 2 * np; k < npt; ++k) {	
            iw = (int) (((double) (k - np) - half) / (double) (n));
            ip = k - np - iw *(n);
            iq = ip + iw;
            if (iq > n) {
                iq -= n;
            }
/* CPB: offset counters by -1 when used as indices on zero-based index */
            ptsid[k-1] = (double) ip + (double) iq / (double) np + sfrac;
            temp = one / (ptsaux[0][ip-1] * ptsaux[0][iq-1]);
            zmat[0][k-np-1] = temp;
            zmat[k-np-1][ip] = -temp;
            zmat[k-np-1][iq] = -temp;
/* L70: */
            zmat[k-np-1][k-1] = temp;
        }
    }
    nrem = npt;
    kold = 0;
    knew = *kopt;

/*     Reorder the provisional points in the way that exchanges PTSID(KOLD) */
/*     with PTSID(KNEW). */

L80:

    for (j = 0; j < n; ++j) {
        temp = bmat[j][kold];
        bmat[j][kold] = bmat[j][knew];
/* L90: */
        bmat[j][knew] = temp;
    }

    for (j = 0; j < nptm; ++j) {
        temp = zmat[j][kold];
        zmat[j][kold] = zmat[j][knew];
/* L100: */
        zmat[j][knew ] = temp;
    }
    ptsid[kold] = ptsid[knew];
    ptsid[knew] = zero;
    w[ndim + knew] = zero;
    --nrem;
    if (knew != *kopt) {
        temp = vlag[kold];
        vlag[kold] = vlag[knew];
        vlag[knew] = temp;

/*     Update the BMAT and ZMAT matrices so that the status of the KNEW-th */
/*     interpolation point can be changed from provisional to original. The */
/*     branch to label 350 occurs if all the original points are reinstated. */
/*     The nonnegative values of W(NDIM+K) are required in the search below. */

        update_(n, npt, ndim, bmat, zmat, vlag, beta, denom, knew);  
        if (nrem == 0) {
            goto L350;
        }

        for (k = 0; k < npt; ++k) {
/* L110: */
            w[ndim + k] = fabs(w[ndim + k]);
        }
    }

/*     Pick the index KNEW of an original interpolation point that has not */
/*     yet replaced one of the provisional interpolation points, giving */
/*     attention to the closeness to XOPT and to previous tries with KNEW. */

L120:
    dsqmin = zero;

    for (k = 0; k < npt; ++k) {
        if (w[ndim + k] > zero) {
            if (dsqmin == zero || w[ndim + k] < dsqmin) {
                knew = k;
                dsqmin = w[ndim + k];
            }
        }
/* L130: */
    }
    if (dsqmin == zero) {
        goto L260;
    }

/*     Form the W-vector of the chosen original interpolation point. */


    for (j = 0; j < n; ++j) {
/* L140: */
        w[npt + j] = xpt[j][knew];
    }

    for (k = 0; k < npt; ++k) {
        sum = zero;
        if (k == *kopt) {
        } 
        else if (ptsid[k] == zero) {

            for (j = 0; j < n; ++j) {
/* L150: */
                sum += w[npt + j] * xpt[j][k];
            }
        } 
        else {
            ip = (int) ptsid[k];
            if (ip > 0) {
                sum = w[npt + ip-1] * ptsaux[0][ip-1];
            }
            iq = (int) ((double) np * ptsid[k] - (double) (ip * np));
            if (iq > 0) {
                iw = 0;
                if (ip == 0) {
                    iw = 1;
                }
                sum += w[npt + iq-1] * ptsaux[iw][iq-1];
            }
        }
/* L160: */
        w[k] = half * sum * sum;
    }

/*     Calculate VLAG and BETA for the required updating of the H matrix if */
/*     XPT(KNEW,.) is reinstated in the set of interpolation points. */


    V_sum[0:npt] = zero;

    for (j = 0; j < n; ++j) {
       for (k = 0; k < npt; ++k) {
/* L170: */
            V_sum[k] += bmat[j][k] * w[npt + j];
        }
    }
/* L180: */
    vlag[0:npt] = V_sum[0:npt];

    beta = zero;
    
    for (j = 0; j < nptm; ++j) {
        sum = zero;
   
        for (k = 0; k < npt; ++k) {
/* L190: */
            sum += zmat[j][k] * w[k];
        }
        beta -= sum * sum;
  
        for (k = 0; k < npt; ++k) {
/* L200: */
            vlag[k] += sum * zmat[j][k];
        }
    }
    bsum = zero;
    distsq = zero;
 
    for (j = 0; j < n; ++j) {
        sum = zero;

        for (k = 0; k < npt; ++k) {
/* L210: */
            sum += bmat[j][k] * w[k];
        }
        jp = j + npt;
        bsum += sum * w[jp];

        for (ip = npt; ip < ndim; ++ip) {
/* L220: */
            sum += bmat[j][ip] * w[ip];
        }
        bsum += sum * w[jp];
        vlag[jp] = sum;
/* L230: */
        distsq += xpt[j][knew] * xpt[j][knew];
    }
    beta = half * distsq * distsq + beta - bsum;
    vlag[*kopt] += one;

/*     KOLD is set to the index of the provisional interpolation point that is */
/*     going to be deleted to make way for the KNEW-th original interpolation */
/*     point. The choice of KOLD is governed by the avoidance of a small value */
/*     of the denominator in the updating calculation of UPDATE. */

    denom = zero;
    vlmxsq = zero;
    V_hdiag[0:npt] = zero;

    for (j = 0; j < nptm; ++j) {
       for (k = 0; k < npt; ++k) 
/* L240: */
          V_hdiag[k] += zmat[j][k] * zmat[j][k];
    }
    for (k = 0; k < npt; ++k) {
        if (ptsid[k] != zero) {
            den = beta * V_hdiag[k] + vlag[k] * vlag[k];
            if (den > denom) {
                kold = k;
                denom = den;
            }
        }
    }

    for (k = 0; k < npt; ++k) {
/* L250: */
        vlmxsq = MAX2(vlmxsq, vlag[k] * vlag[k]);
    }

    if (denom <= vlmxsq * .01) {
        w[ndim + knew] = -w[ndim + knew] - winc;
        goto L120;
    }
    goto L80;

/*     When label 260 is reached, all the final positions of the interpolation */
/*     points have been chosen although any changes have not been included yet */
/*     in XPT. Also the final BMAT and ZMAT matrices are complete, but, apart */
/*     from the shift of XBASE, the updating of the quadratic model remains to */
/*     be done. The following cycle through the new interpolation points begins */
/*     by putting the new point in XPT(KPT,.) and by setting PQ(KPT) to zero, */
/*     except that a RETURN occurs if MAXFUN prohibits another value of F. */

L260:
   
    for (kpt = 0; kpt < npt; ++kpt) {
        if (ptsid[kpt] == zero) continue;
            //goto L340;

        if (nlopt_stop_forced(stop)) {rCode = NLOPT_FORCED_STOP; goto L350;}
        else if (nlopt_stop_evals(stop)) {rCode = NLOPT_MAXEVAL_REACHED; goto L350;}
        else if (nlopt_stop_time(stop)) {rCode = NLOPT_MAXTIME_REACHED; goto L350;}

        ih = 0;
  
        for (j = 0; j < n; ++j) {
            w[j] = xpt[j][kpt];
            xpt[j][kpt] = zero;
            temp = pq[kpt] * w[j];
  
            for (i = 0; i <= j; ++i) {
                ++ih;
/* L270: */
                hq[ih-1] += temp * w[i];
            }    
        }
        pq[kpt] = zero;
        ip = (int) ptsid[kpt];
        iq = (int) ((double) np * ptsid[kpt] - (double) (ip * np));
        if (ip > 0) {
            xp = ptsaux[0][ip-1];
            xpt[ip-1][kpt] = xp;
        }
        if (iq > 0) {
            xq = ptsaux[0][iq-1];
            if (ip == 0) {
                xq = ptsaux[1][iq-1];
            }
            xpt[iq-1][kpt] = xq;
        }

/*     Set VQUAD to the value of the current model at the new point. */

        vquad = fbase;
        if (ip > 0) {
            ihp = (ip + ip * ip) / 2;
            vquad += xp * (gopt[ip-1] + half * xp * hq[ihp-1]);
        }
        if (iq > 0) {
            ihq = (iq + iq * iq) / 2;
            vquad += xq * (gopt[iq-1] + half * xq * hq[ihq-1]);
            if (ip > 0) {
                iw = MAX2(ihp,ihq) - IABS(ip - iq);
                vquad += xp * xq * hq[iw-1];
            }
        }
 
        for (k = 0; k < npt; ++k) {
            temp = zero;
            if (ip > 0) {
                temp += xp * xpt[ip-1][k];
            }
            if (iq > 0) {
                temp += xq * xpt[iq-1][k];
            }
/* L280: */
            vquad += half * pq[k] * temp * temp;
        }

/*     Calculate F at the new interpolation point, and set DIFF to the factor */
/*     that is going to multiply the KPT-th Lagrange function when the model */
/*     is updated to provide interpolation to the new function value. */

        for (i = 0; i < n; ++i) {
            d__4 = xbase[i] + xpt[i][kpt];
            d__1 = MAX2(xl[i],d__4);
            w[i] = MIN2(d__1,xu[i]);
            if (xpt[i][kpt] == sl[i]) {
                w[i] = xl[i];
            }
            if (xpt[i][kpt] == su[i]) {
                w[i] = xu[i];
            }
/* L290: */
        }

        stop->nevals++;
        f = calfun(n, &w[0], calfun_data);
        fval[kpt] = f;
        if (f < fval[*kopt]) {
            *kopt = kpt;
        }
        if (nlopt_stop_forced(stop)) {rCode = NLOPT_FORCED_STOP; goto L350;}
        else if (f < stop->minf_max) {rCode = NLOPT_MINF_MAX_REACHED; goto L350;}
        else if (nlopt_stop_evals(stop)) {rCode = NLOPT_MAXEVAL_REACHED; goto L350;}
        else if (nlopt_stop_time(stop)) {rCode = NLOPT_MAXTIME_REACHED; goto L350;}

        diff = f - vquad;

/*     Update the quadratic model. The RETURN from the subroutine occurs when */
/*     all the new interpolation points are included in the model. */

        for (i = 0; i < n; ++i) {
/* L310: */
            gopt[i] += diff * bmat[i][kpt];
        }
        
        for (k = 0; k < npt; ++k) {
            sum = zero;
       
            for (j = 0; j < nptm; ++j) {
/* L320: */
                sum += zmat[j][k] * zmat[j][kpt];
            }
            temp = diff * sum;
            if (ptsid[k] == zero) {
                pq[k] += temp;
            } 
            else {
                ip = (int) ptsid[k];
                iq = (int) ((double) np * ptsid[k] - (double) (ip * np));
                ihq = (iq * iq + iq) / 2;
                if (ip == 0) {
                    hq[ihq-1] += temp * (ptsaux[1][iq-1] * ptsaux[1][iq-1]);
                } 
                else {
                    ihp = (ip * ip + ip) / 2;
                    hq[ihp-1] += temp * (ptsaux[0][ip-1] * ptsaux[0][ip-1]);
                    if (iq > 0) {
                        hq[ihq-1] += temp * (ptsaux[0][iq-1] * ptsaux[0][iq-1]);
                        iw = MAX2(ihp,ihq) - IABS(iq - ip);
                        hq[iw-1] += temp * ptsaux[0][ip-1] * ptsaux[0][iq-1];
                    }
                }
            }
/* L330: */
        }
        ptsid[kpt] = zero;
/* L340:
        ; */
    }
L350:
    _mm_free(w);
    _mm_free(V_distsq);
    _mm_free(V_sum);
    _mm_free(V_hdiag);
    return rCode;
} /* rescue_ */


static void altmov_(int n, int npt, double **xpt,
    double *xopt, int ndim, double **bmat, double **zmat, 
    double *sl, double *su, int *kopt, int *knew, 
    double *adelt, double *xnew, double *xalt, double *alpha,
    double *cauchy, double *glag, double *hcol )
{
    /* System generated locals */
    int xpt_dim1, xpt_offset, bmat_dim1, bmat_offset, zmat_dim1, 
        zmat_offset, i__1, i__2;
    double d__1, d__2, d__3, d__4;

    /* Local variables */
    int i, j, k;
    double ha, gw, one, diff, half;
    int ilbd, isbd;
    double slbd;
    int iubd;
    double vlag, subd, temp;
    int ksav;
    double step, zero, curv;
    int iflag;
    double scale, csave, tempa, tempb, tempd, OneSQRT2, sumin, ggfree;
    int ibdsav;
    double dderiv, bigstp, predsq, presav, distsq, stpsav, wfixsq, wsqsav;
    double *V_temp;
    double *w = _mm_malloc(sizeof(double) * (2*n+8), 64); 
    /* CPB: values to denote when i == 0 for zero-based index value into XNEW */
    int suZERO = n+1, slZERO = n+1;


/*     The arguments N, NPT, XPT, XOPT, BMAT, ZMAT, NDIM, SL and SU all have */
/*       the same meanings as the corresponding arguments of BOBYQB. */
/*     KOPT is the index of the optimal interpolation point. */
/*     KNEW is the index of the interpolation point that is going to be moved. */
/*     ADELT is the current trust region bound. */
/*     XNEW will be set to a suitable new position for the interpolation point */
/*       XPT(KNEW,.). Specifically, it satisfies the SL, SU and trust region */
/*       bounds and it should provide a large denominator in the next call of */
/*       UPDATE. The step XNEW-XOPT from XOPT is restricted to moves along the */
/*       straight lines through XOPT and another interpolation point. */
/*     XALT also provides a large value of the modulus of the KNEW-th Lagrange */
/*       function subject to the constraints that have been mentioned, its main */
/*       difference from XNEW being that XALT-XOPT is a constrained version of */
/*       the Cauchy step within the trust region. An exception is that XALT is */
/*       not calculated if all components of GLAG (see below) are zero. */
/*     ALPHA will be set to the KNEW-th diagonal element of the H matrix. */
/*     CAUCHY will be set to the square of the KNEW-th Lagrange function at */
/*       the step XALT-XOPT from XOPT for the vector XALT that is returned, */
/*       except that CAUCHY is set to zero if XALT is not calculated. */
/*     GLAG is a working space vector of length N for the gradient of the */
/*       KNEW-th Lagrange function at XOPT. */
/*     HCOL is a working space vector of length NPT for the second derivative */
/*       coefficients of the KNEW-th Lagrange function. */
/*     W is a working space vector of length 2N that is going to hold the */
/*       constrained Cauchy step from XOPT of the Lagrange function, followed */
/*       by the downhill version of XALT when the uphill step is calculated. */

/*     Set the first NPT components of W to the leading elements of the */
/*     KNEW-th column of the H matrix. */


    /* Function Body */
    half = .5;
    one = 1.;
    zero = 0.;
    OneSQRT2 = one + sqrt(2.);
    V_temp = _mm_malloc(sizeof(double) * (npt + 8), 64);

    for (k = 0; k < npt; ++k) {
/* L10: */
        hcol[k] = zero;
    }

    for (j = 0; j < npt - n - 1; ++j) {
        temp = zmat[j][*knew];

        for (k = 0; k < npt; ++k) {
/* L20: */
            hcol[k] += temp * zmat[j][k];
        }
    }
    *alpha = hcol[*knew];
    ha = half * (*alpha);

/*     Calculate the gradient of the KNEW-th Lagrange function at XOPT. */

    for (i = 0; i < n; ++i) {
/* L30: */
        glag[i] = bmat[i][*knew];
    }
    
    V_temp[0:npt] = zero;
   
    for (j = 0; j < n; ++j) {
        temp = xopt[j];
        for (k = 0; k < npt; ++k) {
/* L40: */
            V_temp[k] += xpt[j][k] * temp;
        }
    }
    for (k = 0; k < npt; ++k) 
        V_temp[k] = hcol[k] * V_temp[k];
  
    for (i = 0; i < n; ++i) {
        for (k = 0; k < npt; ++k) {
/* L50: */
            glag[i] += V_temp[k] * xpt[i][k];
        }
    }

/*     Search for a large denominator along the straight lines through XOPT */
/*     and another interpolation point. SLBD and SUBD will be lower and upper */
/*     bounds on the step along each of these lines in turn. PREDSQ will be */
/*     set to the square of the predicted denominator for each line. PRESAV */
/*     will be set to the largest admissible value of PREDSQ that occurs. */

    presav = zero;
    
    for (k = 0; k < npt; ++k) {
        if (k == *kopt) continue;
            //goto L80;
        dderiv = zero;
        distsq = zero;
        
        for (i = 0; i < n; ++i) {
            temp = xpt[i][k] - xopt[i];
            dderiv += glag[i] * temp;
/* L60: */
            distsq += temp * temp;
        }
        subd = *adelt / sqrt(distsq);
        slbd = -subd;
        ilbd = 0;
        iubd = 0;
        sumin = MIN2(one,subd);

/*     Revise SLBD and SUBD if necessary because of the bounds in SL and SU. */

        for (i = 0; i < n; ++i) {
            temp = xpt[i][k] - xopt[i];
            if (temp > zero) {
                if (slbd * temp < sl[i] - xopt[i]) {
                    slbd = (sl[i] - xopt[i]) / temp;
                    ilbd = (i==0) ? -slZERO : -i;
                }
                if (subd * temp > su[i] - xopt[i]) {
                    subd = MAX2(sumin, (su[i] - xopt[i]) / temp);
                    iubd = (i==0) ? suZERO : i;
                }
            } 
            else if (temp < zero) {
                if (slbd * temp > su[i] - xopt[i]) {
                    slbd = (su[i] - xopt[i]) / temp;
                    ilbd = (i==0) ? slZERO : i;
                }
                if (subd * temp < sl[i] - xopt[i]) {
                    subd = MAX2(sumin,(sl[i] - xopt[i]) / temp);
                    iubd = (i==0) ? -suZERO : -i;
                }
            }
/* L70: */
        }

/*     Seek a large modulus of the KNEW-th Lagrange function when the index */
/*     of the other interpolation point on the line through XOPT is KNEW. */

        if (k == *knew) {
            diff = dderiv - one;
            step = slbd;
            vlag = slbd * (dderiv - slbd * diff);
            isbd = ilbd;
            temp = subd * (dderiv - subd * diff);
            if (fabs(temp) > fabs(vlag)) {
                step = subd;
                vlag = temp;
                isbd = iubd;
            }
            tempd = half * dderiv;
            tempa = tempd - diff * slbd;
            tempb = tempd - diff * subd;
            if (tempa * tempb < zero) {
                temp = tempd * tempd / diff;
                if (fabs(temp) > fabs(vlag)) {
                step = tempd / diff;
                vlag = temp;
                isbd = 0;
                }
            }

/*     Search along each of the other lines through XOPT and another point. */

        } 
        else {
            step = slbd;
            vlag = slbd * (one - slbd);
            isbd = ilbd;
            temp = subd * (one - subd);
            if (fabs(temp) > fabs(vlag)) {
                step = subd;
                vlag = temp;
                isbd = iubd;
            }
            if (subd > half) {
                if (fabs(vlag) < .25) {
                    step = half;
                    vlag = .25;
                    isbd = 0;
                }
            }
            vlag *= dderiv;
        }

/*     Calculate PREDSQ for the current line search and maintain PRESAV. */

        temp = step * (one - step) * distsq;
        predsq = vlag * vlag * (vlag * vlag + ha * temp * temp);
        if (predsq > presav) {
            presav = predsq;
            ksav = k;
            stpsav = step;
            ibdsav = isbd;
        }
/*L80:
    ; */
    }

/*     Construct XNEW in a way that satisfies the bound constraints exactly. */

    for (i = 0; i < n; ++i) {
        temp = xopt[i] + stpsav * (xpt[i][ksav] - xopt[i]);
/* L90: */
        d__2 = MIN2(su[i],temp);
        xnew[i] = MAX2(sl[i],d__2);
    }
/* CPB: Use variables slZERO and suZERO to denote that the [0] element 
 * was used to set ilbd and iubd between L70 and L80. This leaves the zero 
 * value (in isbd) as the indication to not update xnew after loop L90. 
 * Test for the saved (ibdsav) index modified to test for s_ZERO value to 
 * signify the zero result stored in isbd. */
            
    if (ibdsav < 0) {
        if (ibdsav == -slZERO || ibdsav == -suZERO) 
            xnew[0] = sl[0];
        else    
            xnew[-ibdsav] = sl[-ibdsav];
    }
    if (ibdsav > 0) {
        if (ibdsav == slZERO || ibdsav == suZERO) 
            xnew[0] = su[0];
        else    
            xnew[ibdsav] = su[ibdsav];
    }

/*     Prepare for the iterative method that assembles the constrained Cauchy */
/*     step in W. The sum of squares of the fixed components of W is formed in */
/*     WFIXSQ, and the free components of W are set to BIGSTP. */

    bigstp = *adelt + *adelt;
    iflag = 0;
L100:
    wfixsq = zero;
    ggfree = zero;
   
    for (i = 0; i < n; ++i) {
        w[i] = zero;
        tempa = MIN2(xopt[i] - sl[i], glag[i]);
        tempb = MAX2(xopt[i] - su[i], glag[i]);
        if (tempa > zero || tempb < zero) {
            w[i] = bigstp;
            ggfree += glag[i] * glag[i];
        }
/* L110: */
    }
    if (ggfree == zero) {
        *cauchy = zero;
        goto L200;  /* RETURN */
    }

/*     Investigate whether more components of W can be fixed. */

L120:
    temp = *adelt * (*adelt) - wfixsq;
    if (temp > zero) {
        wsqsav = wfixsq;
        step = sqrt(temp / ggfree);
        ggfree = zero;
       
        for (i = 0; i < n; ++i) {
            if (w[i] == bigstp) {
                temp = xopt[i] - step * glag[i];
                if (temp <= sl[i]) {
                    w[i] = sl[i] - xopt[i];
                    wfixsq += w[i] * w[i];
                } 
                else if (temp >= su[i]) {
                    w[i] = su[i] - xopt[i];
                    wfixsq += w[i] * w[i];
                } 
                else {
                    ggfree += glag[i] * glag[i];
                }
            }
/* L130: */
        }
        if (wfixsq > wsqsav && ggfree > zero) { /* repeat - until */
            goto L120;
        }
    }

/*     Set the remaining free components of W and all components of XALT, */
/*     except that W may be scaled later. */

    gw = zero;

    for (i = 0; i < n; ++i) {
        if (w[i] == bigstp) {
            w[i] = -step * glag[i];
            d__2 = MIN2(su[i],xopt[i] + w[i]);
            xalt[i] = MAX2(sl[i],d__2);
        } 
        else if (w[i] == zero) {
           xalt[i] = xopt[i];
        } 
        else if (glag[i] > zero) {
           xalt[i] = sl[i];
        } 
        else {
           xalt[i] = su[i];
        }
/* L140: */
        gw += glag[i] * w[i];
    }

/*     Set CURV to the curvature of the KNEW-th Lagrange function along W. */
/*     Scale W by a factor less than one if that can reduce the modulus of */
/*     the Lagrange function at XOPT+W. Set CAUCHY to the final value of */
/*     the square of this function. */

    curv = zero;

    V_temp[0:npt] = zero;

    for (j = 0; j < n; ++j) {
       temp = w[j];
       for (k = 0; k < npt; ++k) 
/* L150: */
            V_temp[k] += xpt[j][k] * temp;
    }

    for (k = 0; k < npt; ++k) {
/* L160: */
        curv += hcol[k] * V_temp[k] * V_temp[k];
    }

    if (iflag == 1) {
        curv = -curv;
    }
    if (curv > -gw && curv < -OneSQRT2 * gw) {
        scale = -gw / curv;

        for (i = 0; i < n; ++i) {
            temp = xopt[i] + scale * w[i];
/* L170: */
            d__2 = MIN2(su[i],temp);
            xalt[i] = MAX2(sl[i],d__2);
        }
        *cauchy = (half * gw * scale) * (half * gw * scale);
    } 
    else {
        *cauchy = (gw + half * curv) * (gw + half * curv);
    }

/*     If IFLAG is zero, then XALT is calculated as before after reversing */
/*     the sign of GLAG. Thus two XALT vectors become available. The one that */
/*     is chosen is the one that gives the larger value of CAUCHY. */

    if (iflag == 0) {

        for (i = 0; i < n; ++i) {
            glag[i] = -glag[i];
/* L180: */
            w[n + i] = xalt[i];
        }
        csave = *cauchy;
        iflag = 1;
        goto L100;
    }
    if (csave > *cauchy) {

        for (i = 0; i < n; ++i) {
/* L190: */
            xalt[i] = w[n + i];
        }
        *cauchy = csave;
    }
L200:
    _mm_free(V_temp);
    _mm_free(w);
    return;
} /* altmov_ */

static void trsbox_(int n, int npt, double **xpt,
    double *xopt, double *gopt, double *hq, double *pq, 
    double *sl, double *su, double *delta, double *xnew, 
    double *d, double *gnew, double *xbdi, double *s, 
    double *hs, double *hred, double *dsq, double *crvmin)
{
    /* System generated locals */
    int xpt_dim1, xpt_offset, i__1, i__2;
    double d__1, d__2, d__3, d__4;

    /* Local variables */
    int i, j, k, ih;
    double ds;
    int iu;
    double dhd, dhs, cth, one, shs, sth, ssq, half, beta, sdec, blen;
    int iact, nact;
    double angt, qred;
    int isav;
    double temp, zero, xsav, xsum, angbd, dredg, sredg;
    int iterc;
    double resid, delsq, ggsav, tempa, tempb, ratio, sqstp, redmax, 
        dredsq, redsav, onemin, gredsq, rednew;
    int itcsav;
    double rdprev, rdnext, stplen, stepsq;
    int itermax;
    double *V_temp;


/*     The arguments N, NPT, XPT, XOPT, GOPT, HQ, PQ, SL and SU have the same */
/*       meanings as the corresponding arguments of BOBYQB. */
/*     DELTA is the trust region radius for the present calculation, which */
/*       seeks a small value of the quadratic model within distance DELTA of */
/*       XOPT subject to the bounds on the variables. */
/*     XNEW will be set to a new vector of variables that is approximately */
/*       the one that minimizes the quadratic model within the trust region */
/*       subject to the SL and SU constraints on the variables. It satisfies */
/*       as equations the bounds that become active during the calculation. */
/*     D is the calculated trial step from XOPT, generated iteratively from an */
/*       initial value of zero. Thus XNEW is XOPT+D after the final iteration. */
/*     GNEW holds the gradient of the quadratic model at XOPT+D. It is updated */
/*       when D is updated. */
/*     XBDI is a working space vector. For I=1,2,...,N, the element XBDI(I) is */
/*       set to -1.0, 0.0, or 1.0, the value being nonzero if and only if the */
/*       I-th variable has become fixed at a bound, the bound being SL(I) or */
/*       SU(I) in the case XBDI(I)=-1.0 or XBDI(I)=1.0, respectively. This */
/*       information is accumulated during the construction of XNEW. */
/*     The arrays S, HS and HRED are also used for working space. They hold the */
/*       current search direction, and the changes in the gradient of Q along S */
/*       and the reduced D, respectively, where the reduced D is the same as D, */
/*       except that the components of the fixed variables are zero. */
/*     DSQ will be set to the square of the length of XNEW-XOPT. */
/*     CRVMIN is set to zero if D reaches the trust region boundary. Otherwise */
/*       it is set to the least curvature of H that occurs in the conjugate */
/*       gradient searches that are not restricted by any constraints. The */
/*       value CRVMIN=-1.0D0 is set, however, if all of these searches are */
/*       constrained. */

/*     A version of the truncated conjugate gradient is applied. If a line */
/*     search is restricted by a constraint, then the procedure is restarted, */
/*     the values of the variables that are at their bounds being fixed. If */
/*     the trust region boundary is reached, then further changes may be made */
/*     to D, each one being in the two dimensional space that is spanned */
/*     by the current D and the gradient of Q at XOPT+D, staying on the trust */
/*     region boundary. Termination occurs when the reduction in Q seems to */
/*     be close to the greatest reduction that can be achieved. */

/*     Set some constants. */

    /* Function Body */
    half = .5;
    one = 1.;
    onemin = -1.;
    zero = 0.;
    V_temp = _mm_malloc(sizeof(double) * (npt + 8), 64);

/*     The sign of GOPT(I) gives the sign of the change to the I-th variable */
/*     that will reduce Q from its value at XOPT. Thus XBDI(I) shows whether */
/*     or not to fix the I-th variable at one of its bounds initially, with */
/*     NACT being set to the number of fixed variables. D and GNEW are also */
/*     set for the first iteration. DELSQ is the upper bound on the sum of */
/*     squares of the free variables. QRED is the reduction in Q so far. */

    iterc = 0;
    nact = 0;
    sqstp = zero;
    
    for (i = 0; i < n; ++i) {
        xbdi[i] = zero;
        if (xopt[i] <= sl[i]) {
            if (gopt[i] >= zero) {
                xbdi[i] = onemin;
            }
        } 
        else if (xopt[i] >= su[i]) {
            if (gopt[i] <= zero) {
                xbdi[i] = one;
            }
        }
        if (xbdi[i] != zero) {
            ++nact;
        }
        d[i] = zero;
/* L10: */
        gnew[i] = gopt[i];
    }
    delsq = *delta * *delta;
    qred = zero;
    *crvmin = onemin;

/*     Set the next search direction of the conjugate gradient method. It is */
/*     the steepest descent direction initially and when the iterations are */
/*     restarted because a variable has just been fixed by a bound, and of */
/*     course the components of the fixed variables are zero. ITERMAX is an */
/*     upper bound on the indices of the conjugate gradient iterations. */

L20:
    beta = zero;
L30:
    stepsq = zero;

    for (i = 0; i < n; ++i) {
        if (xbdi[i] != zero) {
            s[i] = zero;
        } 
        else if (beta == zero) {
            s[i] = -gnew[i];
        } 
        else {
            s[i] = beta * s[i] - gnew[i];
        }
/* L40: */
        stepsq += s[i] * s[i];
    }
    if (stepsq == zero) {
        goto L190;
    }
    if (beta == zero) {
        gredsq = stepsq;
        itermax = iterc + n - nact;
    }
    if (gredsq * delsq <= qred * 1e-4 * qred) {
        goto L190;
    }

/*     Multiply the search direction by the second derivative matrix of Q and */
/*     calculate some scalars for the choice of steplength. Then set BLEN to */
/*     the length of the the step to the trust region boundary and STPLEN to */
/*     the steplength, ignoring the simple bounds. */

    goto L210;
L50:
    resid = delsq;
    ds = zero;
    shs = zero;

    for (i = 0; i < n; ++i) {
        if (xbdi[i] == zero) {
            resid -= d[i] * d[i];
            ds += s[i] * d[i];
            shs += s[i] * hs[i];
        }
/* L60: */
    }
    if (resid <= zero) {
        goto L90;
    }
    temp = sqrt(stepsq * resid + ds * ds);
    if (ds < zero) {
        blen = (temp - ds) / stepsq;
    } 
    else {
        blen = resid / (temp + ds);
    }
    stplen = blen;
    if (shs > zero) {
        stplen = MIN2(blen, gredsq / shs);
    }

/*     Reduce STPLEN if necessary in order to preserve the simple bounds, */
/*     letting IACT be the index of the new constrained variable. */

    iact = -1; /* CPB: Use -1 as marker that iact not set */

    for (i = 0; i < n; ++i) {
        if (s[i] != zero) {
            xsum = xopt[i] + d[i];
            if (s[i] > zero) {
                temp = (su[i] - xsum) / s[i];
            } 
            else {
                temp = (sl[i] - xsum) / s[i];
            }
            if (temp < stplen) {
                stplen = temp;
                iact = i;
            }
        }
/* L70: */
    }

/*     Update CRVMIN, GNEW and D. Set SDEC to the decrease that occurs in Q. */

    sdec = zero;
    if (stplen > zero) {
        ++iterc;
        temp = shs / stepsq;
        if (iact == -1 && temp > zero) {
            *crvmin = MIN2(*crvmin,temp);
            if (*crvmin == onemin) {
                *crvmin = temp;
            }
        }
        ggsav = gredsq;
        gredsq = zero;

        for (i = 0; i < n; ++i) {
            gnew[i] += stplen * hs[i];
            if (xbdi[i] == zero) {
                gredsq += gnew[i] * gnew[i];
            }
/* L80: */
            d[i] += stplen * s[i];
        }
        d__1 = stplen * (ggsav - half * stplen * shs);
        sdec = MAX2(d__1,zero);
        qred += sdec;
    }

/*     Restart the conjugate gradient method if it has hit a new bound. */

    if (iact >= 0) {
        ++nact;
        xbdi[iact] = one;
        if (s[iact] < zero) {
            xbdi[iact] = onemin;
        }
        delsq -= d[iact] * d[iact];
        if (delsq <= zero) {
            goto L90;
        }
        goto L20;
    }

/*     If STPLEN is less than BLEN, then either apply another conjugate */
/*     gradient iteration or RETURN. */

    if (stplen < blen) {
        if (iterc == itermax) {
            goto L190;
        }
        if (sdec <= qred * .01) {
            goto L190;
        }
        beta = gredsq / ggsav;
        goto L30;
    }
L90:
    *crvmin = zero;

/*     Prepare for the alternative iteration by calculating some scalars */
/*     and by multiplying the reduced D by the second derivative matrix of */
/*     Q, where S holds the reduced D in the call of GGMULT. */

L100:
    if (nact >= n - 1) {
        goto L190;
    }
    dredsq = zero;
    dredg = zero;
    gredsq = zero;

    for (i = 0; i < n; ++i) {
        if (xbdi[i] == zero) {
            dredsq += d[i] * d[i];
            dredg += d[i] * gnew[i];
            gredsq += gnew[i] * gnew[i];
            s[i] = d[i];
        } 
        else {
            s[i] = zero;
        }
/* L110: */
    }
    itcsav = iterc;
    goto L210;

/*     Let the search direction S be a linear combination of the reduced D */
/*     and the reduced G that is orthogonal to the reduced D. */

L120:
    ++iterc;
    temp = gredsq * dredsq - dredg * dredg;
    if (temp <= qred * 1e-4 * qred) {
        goto L190;
    }
    temp = sqrt(temp);
 
    for (i = 0; i < n; ++i) {
        if (xbdi[i] == zero) {
            s[i] = (dredg * d[i] - dredsq * gnew[i]) / temp;
        } 
        else {
            s[i] = zero;
        }
/* L130: */
    }
    sredg = -temp;

/*     By considering the simple bounds on the variables, calculate an upper */
/*     bound on the tangent of half the angle of the alternative iteration, */
/*     namely ANGBD, except that, if already a free variable has reached a */
/*     bound, there is a branch back to label 100 after fixing that variable. */

    angbd = one;
    iact = -1;  /* CPB: Use -1 as marker that iact not set */

    for (i = 0; i < n; ++i) {
        if (xbdi[i] == zero) {
            tempa = xopt[i] + d[i] - sl[i];
            tempb = su[i] - xopt[i] - d[i];
            if (tempa <= zero) {
                ++nact;
                xbdi[i] = onemin;
                goto L100;
            } 
            else if (tempb <= zero) {
                ++nact;
                xbdi[i] = one;
                goto L100;
            }
            ratio = one;
            ssq = d[i] * d[i] + s[i] * s[i];
            temp = ssq - (xopt[i] - sl[i]) * (xopt[i] - sl[i]);
            if (temp > zero) {
                temp = sqrt(temp) - s[i];
                if (angbd * temp > tempa) {
                    angbd = tempa / temp;
                    iact = i;
                    xsav = onemin;
                }
            }
            temp = ssq - (su[i] - xopt[i]) * (su[i] - xopt[i]);
            if (temp > zero) {
                temp = sqrt(temp) + s[i];
                if (angbd * temp > tempb) {
                    angbd = tempb / temp;
                    iact = i;
                    xsav = one;
                }
            }
        }    
/* L140: */
    }

/*     Calculate HHD and some curvatures for the alternative iteration. */

    goto L210;
L150:
    shs = zero;
    dhs = zero;
    dhd = zero;
   
    for (i = 0; i < n; ++i) {
        if (xbdi[i] == zero) {
            shs += s[i] * hs[i];
            dhs += d[i] * hs[i];
            dhd += d[i] * hred[i];
        }
/* L160: */
    }

/*     Seek the greatest reduction in Q for a range of equally spaced values */
/*     of ANGT in [0,ANGBD], where ANGT is the tangent of half the angle of */
/*     the alternative iteration. */

    redmax = zero;
    isav = -1; /* Use -1 as marker that isav not set */
    redsav = zero;
    iu = (int) (angbd * 17. + 3.1);

    for (i = 0; i < iu; ++i) {
        angt = angbd * (double) (i+1) / (double) iu;
        sth = (angt + angt) / (one + angt * angt);
        temp = shs + angt * (angt * dhd - dhs - dhs);
        rednew = sth * (angt * dredg - sredg - half * sth * temp);
        if (rednew > redmax) {
            redmax = rednew;
            isav = i;
            rdprev = redsav;
        } 
        else if (i == isav + 1) {
            rdnext = rednew;
        }
/* L170: */
        redsav = rednew;
    }

/*     Return if the reduction is zero. Otherwise, set the sine and cosine */
/*     of the angle of the alternative iteration, and calculate SDEC. */

    if (isav < 0) {
        goto L190;
    }

/* CPB: incremented isav here since it is off by 1 (from Fortran 
 * version) and is compared to iu (but not used as index) */
    ++isav;

    if (isav < iu) {
        temp = (rdnext - rdprev) / (redmax + redmax - rdprev - rdnext);
        angt = angbd * ((double) isav + half * temp) / (double) iu;
    }
    cth = (one - angt * angt) / (one + angt * angt);
    sth = (angt + angt) / (one + angt * angt);
    temp = shs + angt * (angt * dhd - dhs - dhs);
    sdec = sth * (angt * dredg - sredg - half * sth * temp);
    if (sdec <= zero) {
        goto L190;
    }

/*     Update GNEW, D and HRED. If the angle of the alternative iteration */
/*     is restricted by a bound on a free variable, that variable is fixed */
/*     at the bound. */

    dredg = zero;
    gredsq = zero;

    for (i = 0; i < n; ++i) {
        gnew[i] = gnew[i] + (cth - one) * hred[i] + sth * hs[i];
        if (xbdi[i] == zero) {
            d[i] = cth * d[i] + sth * s[i];
            dredg += d[i] * gnew[i];
            gredsq += gnew[i] * gnew[i];
        }
/* L180: */
        hred[i] = cth * hred[i] + sth * hs[i];
    }
    qred += sdec;
    if (iact >= 0 && isav == iu) {
        ++nact;
        xbdi[iact] = xsav;
        goto L100;
    }

/*     If SDEC is sufficiently small, then RETURN after setting XNEW to */
/*     XOPT+D, giving careful attention to the bounds. */

    if (sdec > qred * .01) {
        goto L120;
    }
L190:
    *dsq = zero;

    for (i = 0; i < n; ++i) {
        d__1 = MIN2((xopt[i] + d[i]),su[i]);
        xnew[i] = MAX2(d__1,sl[i]);
        if (xbdi[i] == onemin) {
            xnew[i] = sl[i];
        }
        if (xbdi[i] == one) {
            xnew[i] = su[i];
        }
        d[i] = xnew[i] - xopt[i];
/* L200: */
        *dsq += d[i] * d[i];
    }
    //printf("trsbox OUT:\ngnew, gopt: %f %f %f %f\n", gnew[0], gnew[1], gopt[0], gopt[1]);
    //printf("s, hs: %f %f %f %f\n", s[0], s[1], hs[0], hs[1]);
    //printf("hq: %f %f %f\n\n", hq[0], hq[1], hq[2]);

    _mm_free(V_temp);
    return;
/*     The following instructions multiply the current S-vector by the second */
/*     derivative matrix of the quadratic model, putting the product in HS. */
/*     They are reached from three different parts of the software above and */
/*     they can be regarded as an external subroutine. */

L210:
    ih = 0;

    for (j = 0; j < n; ++j) {
        hs[j] = zero;

        for (i = 0; i <= j; ++i) {
            ++ih;
            if (i < j) {
                hs[j] += hq[ih-1] * s[i];
            }
/* L220: */
            hs[i] += hq[ih-1] * s[j];
        }
    }

/* CPB: Restructuring of loops to access xpt in row-major order will compute */
/* excess steps, with some elements of V_temp == zero when pq[k] == zero */
/* Done to enhance utilization of vectorization.
 */

    V_temp[0:npt] = zero;

    for (j = 0; j < n; ++j) {
        for (k = 0; k < npt; ++k) {
/* L230: */
            V_temp[k] += xpt[j][k] * s[j];
        }
    }

    for (k = 0; k < npt; ++k) 
        V_temp[k] *= pq[k];

    for (i = 0; i < n; ++i) {
        for (k = 0; k < npt; ++k) {
/* L240: */
            hs[i] += V_temp[k] * xpt[i][k];
        }
/* L250: */
    }
    if (*crvmin != zero) {
        goto L50;
    }
    if (iterc > itcsav) {
        goto L150;
    }
    //i__2 = n;
    for (i = 0; i < n; ++i) {
/* L260: */
        hred[i] = hs[i];
    }
    goto L120;
} /* trsbox_ */

static nlopt_result prelim_(int n, int npt, double *x, 
    const double *xl, const double *xu, double *rhobeg, 
            nlopt_stopping *stop,
            bobyqa_func calfun, void *calfun_data,
     double *xbase, double **xpt, double *fval,
     double *gopt, double *hq, double *pq, int ndim, double **bmat, 
    double **zmat, double *sl, double *su, 
            int *kopt)
{
    /* System generated locals */
    double d__1, d__2, d__3, d__4;

    /* Local variables */
    double f;
    int i, j, k, ih, np, nfm;
    double one;
    int nfx, ipt, jpt;
    double two, fbeg, diff, half, temp, zero, recip, stepa, stepb;
    int itemp;
    double rhosq;

    int nf;
/*
 * CPB: index versions of counter variables for zero-based indexing
 */
    int nf_idx, nfm_idx, nfx_idx;

/*     The arguments N, NPT, X, XL, XU, RHOBEG, and MAXFUN are the */
/*       same as the corresponding arguments in SUBROUTINE BOBYQA. */
/*     The arguments XBASE, XPT, FVAL, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU */
/*       are the same as the corresponding arguments in BOBYQB, the elements */
/*       of SL and SU being set in BOBYQA. */
/*     GOPT is usually the gradient of the quadratic model at XOPT+XBASE, but */
/*       it is set by PRELIM to the gradient of the quadratic model at XBASE. */
/*       If XOPT is nonzero, BOBYQB will change it to its usual value later. */
/*     NF is maintaned as the number of calls of CALFUN so far. */
/*     KOPT will be such that the least calculated value of F so far is at */
/*       the point XPT(KOPT,.)+XBASE in the space of the variables. */

/*     SUBROUTINE PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ, */
/*     BMAT and ZMAT for the first iteration, and it maintains the values of */
/*     NF and KOPT. The vector X is also changed by PRELIM. */

/*     Set some constants. */

    /* Function Body */
    half = .5;
    one = 1.;
    two = 2.;
    zero = 0.;
    rhosq = *rhobeg * *rhobeg;
    recip = one / rhosq;
    np = n + 1;

/*     Set XBASE to the initial vector of variables, and set the initial */
/*     elements of XPT, BMAT, HQ, PQ and ZMAT to zero. */

    /* CPB: initializations can be done in single statement with CEAN */

    for (j = 0; j < n; ++j) {
        xbase[j] = x[j];
    
        for (k = 0; k < npt; ++k) {
/* L10: */
            xpt[j][k] = zero;
        }
   
        for (i = 0; i < ndim; ++i) {
/* L20: */
            bmat[j][i] = zero;
        }
    }
  
    for (ih = 0; ih < (n * np) / 2; ++ih) {
/* L30: */
        hq[ih] = zero;
    }
 
    for (k = 0; k < npt; ++k) {
        pq[k] = zero;

        for (j = 0; j < npt - np; ++j) {
/* L40: */
            zmat[j][k] = zero;
        }
    }

/*    Begin the initialization procedure. NF becomes one more than the number */
/*    of function values so far. The coordinates of the displacement of the */
/*    next initial interpolation point from XBASE are set in XPT(NF+1,.). */

    nf = 0;
    nf_idx = -1;
L50:
    nfm = nf;
    nfm_idx = nf_idx;
    nfx = nf - n;
    nfx_idx = nfx - 1;
    ++(nf);
    ++(nf_idx);
    if (nfm <= n * 2) {        
        if (nfm >= 1 && nfm <= n) {
            stepa = *rhobeg;
            if (su[nfm_idx] == zero) {
                stepa = -stepa;
            }
            xpt[nfm_idx][nf_idx] = stepa;
        } 
        else if (nfm > n) {        
            stepa = xpt[nfx_idx][nf_idx - n]; 
            stepb = -(*rhobeg);
            if (sl[nfx_idx] == zero) {
                stepb = MIN2(two * *rhobeg,su[nfx_idx]);
            }
            if (su[nfx_idx] == zero) {
                stepb = MAX2(-two * *rhobeg,sl[nfx_idx]);
            }
            xpt[nfx_idx][nf_idx] = stepb;
        }
    } 
    else {
        itemp = (nfm - np) / n;
        jpt = nfm - itemp * n - n;
        ipt = jpt + itemp;
        if (ipt > n) {
            itemp = jpt;
            jpt = ipt - n;
            ipt = itemp;
        }
        /* CPB: used -1 for 0-based indexing */
        xpt[ipt-1][nf_idx] = xpt[ipt-1][ipt]; 
        xpt[jpt-1][nf_idx] = xpt[jpt-1][jpt];
    }
/*     Calculate the next value of F. The least function value so far and */
/*     its index are required. */


    for (j = 0; j < n; ++j) {
        d__1 = MAX2(xl[j],xbase[j] + xpt[j][nf_idx]);
        x[j] = MIN2(d__1,xu[j]);
        if (xpt[j][nf_idx] == sl[j]) {
            x[j] = xl[j];
        }
        if (xpt[j][nf_idx] == su[j]) {
            x[j] = xu[j];
        }
/* L60: */
    }
    stop->nevals++;
    f = calfun(n, x, calfun_data);
    fval[nf_idx] = f;
    if (nf == 1) {
        fbeg = f;
        *kopt = 0;
    } 
    else if (f < fval[*kopt]) {
        *kopt = nf_idx;
    }

/*     Set the nonzero initial elements of BMAT and the quadratic model in the */
/*     cases when NF is at most 2*N+1. If NF exceeds N+1, then the positions */
/*     of the NF-th and (NF-N)-th interpolation points may be switched, in */
/*     order that the function value at the first of them contributes to the */
/*     off-diagonal second derivative terms of the initial quadratic model. */

    if (nf <= (n * 2) + 1) {
        if (nf >= 2 && nf <= n + 1) {
            gopt[nfm_idx] = (f - fbeg) / stepa;
            if (npt < nf + n) {
                bmat[0][nfm_idx] = -one / stepa;
                bmat[nfm_idx][nf_idx] = one / stepa;
                bmat[nfm_idx][npt + nfm_idx] = -half * rhosq; /* WATCH */
            }
        } 
        else if (nf >= n + 2) {        /* CPB: was (nf >= n + 2) */
            ih = nfx * (nfx + 1) / 2;
            temp = (f - fbeg) / stepb;
            diff = stepb - stepa;
            hq[ih-1] = two * (temp - gopt[nfx_idx]) / diff;
            gopt[nfx_idx] = (gopt[nfx_idx] * stepb - temp * stepa) / diff;
            if (stepa * stepb < zero) {
                if (f < fval[nf_idx - n]) {
                    fval[nf_idx] = fval[nf_idx - n];
                    fval[nf_idx - n] = f;
                    if (*kopt == nf_idx) {
                        *kopt = nf_idx - n;
                    }
                    xpt[nfx_idx][nf_idx - n] = stepb;
                    xpt[nfx_idx][nf_idx] = stepa;
                }
            }
            bmat[nfx_idx][0] = -(stepa + stepb) / (stepa * stepb);
            bmat[nfx_idx][nf_idx] = -half / xpt[nfx_idx][nf_idx - n];
            bmat[nfx_idx][nf_idx - n] = -bmat[nfx_idx][0] - bmat[nfx_idx][nf_idx];
            zmat[nfx_idx][0] = sqrt(two) / (stepa * stepb);
            zmat[nfx_idx][nf_idx] = sqrt(half) / rhosq;
            zmat[nfx_idx][nf_idx - n] = -zmat[nfx_idx][0] - zmat[nfx_idx][nf_idx];
        }

/*     Set the off-diagonal second derivatives of the Lagrange functions and */
/*     the initial quadratic model. */

    } 
    else { /* CPB Adjusted ipt and jpt values by -1 when used as index values */
        ih = ipt * (ipt - 1) / 2 + jpt;
        zmat[nfx_idx][0] = recip;
        zmat[nfx_idx][nf_idx] = recip;
        zmat[nfx_idx][ipt] = -recip; /* CPB: Removed +1 on ipt */
        zmat[nfx_idx][jpt] = -recip; /* CPB: Removed +1 on jpt */
        temp = xpt[ipt-1][nf_idx] * xpt[jpt-1][nf_idx];
        hq[ih-1] = (fbeg - fval[ipt] - fval[jpt] + f) / temp; 
        /* CPB: Removed +1 on ipt and jpt in line above */
    }
    if (nlopt_stop_forced(stop)) return NLOPT_FORCED_STOP;
    else if (f < stop->minf_max) return NLOPT_MINF_MAX_REACHED;
    else if (nlopt_stop_evals(stop)) return NLOPT_MAXEVAL_REACHED;
    else if (nlopt_stop_time(stop))  return NLOPT_MAXTIME_REACHED;
    
    if (nf < npt) {
        goto L50;
    }
    return NLOPT_SUCCESS;
} /* prelim_ */

static nlopt_result bobyqb_(int n, int npt, double *x, 
    const double *xl, const double *xu, double *rhobeg, double *
    rhoend, 
                nlopt_stopping *stop,
                bobyqa_func calfun, void *calfun_data,
                double *minf,
        double *xbase, 
    double **xpt, double *fval, double *xopt, double *gopt,
     double *hq, double *pq, int ndim, double **bmat, 
    double **zmat, 
    double *sl, double *su, double *xnew, 
    double *xalt, double *d, double *vlag)
{
    /* System generated locals */
    double d__1, d__2;

    /* Local variables */
    double f;
    int i, j, k, ih, jj, nh, ip, jp;
    double dx;
    int np;
    double den, one, ten, dsq, rho, sum, two, diff, half, beta, gisq;
    int knew;
    double temp, suma, sumb, bsum, fopt;
    int kopt, nptm;
    double zero, curv;
    int ksav;
    double gqsq, dist, sumw, sumz, diffa, diffb, diffc, hdiag;
    int kbase;
    double alpha, delta, adelt, denom, fsave, bdtol, delsq, inv_delsq;
    int nresc, nfsav;
    double ratio, dnorm, vquad, pqold, tenth;
    int itest;
    double sumpq, scaden;
    double errbig, cauchy, fracsq, biglsq, densav;
    double bdtest;
    double crvmin, frhosq;
    double distsq;
    int ntrits;
    double xoptsq;
/*
 * CPB: try _mm_malloc(p, 64) and __assume_aligned(V_a, 64) before vector loops
 *
 */
    double *V_hdiag, *V_den, *V_distsq, *V_d__2, *V_temp;
    double *V_sum, *V_suma, *V_sumb;
    double *w_tp;
    double *w = _mm_malloc(sizeof(double) * (3*ndim + 8), 64);
    /* for rescue_ */
    double **ptsaux, *ptsid;
    /* for altmov_ */
    double *glag, *hcol;
    /* for trsbox_ */
    double *gnew, *xbdi, *s, *hs, *hred;

    nlopt_result rc = NLOPT_SUCCESS, rc2;

/*     The arguments N, NPT, X, XL, XU, RHOBEG, RHOEND, and MAXFUN */
/*       are identical to the corresponding arguments in SUBROUTINE BOBYQA. */
/*     XBASE holds a shift of origin that should reduce the contributions */
/*       from rounding errors to values of the model and Lagrange functions. */
/*     XPT is a two-dimensional array that holds the coordinates of the */
/*       interpolation points relative to XBASE. */
/*     FVAL holds the values of F at the interpolation points. */
/*     XOPT is set to the displacement from XBASE of the trust region centre. */
/*     GOPT holds the gradient of the quadratic model at XBASE+XOPT. */
/*     HQ holds the explicit second derivatives of the quadratic model. */
/*     PQ contains the parameters of the implicit second derivatives of the */
/*       quadratic model. */
/*     BMAT holds the last N columns of H. */
/*     ZMAT holds the factorization of the leading NPT by NPT submatrix of H, */
/*       this factorization being ZMAT times ZMAT^T, which provides both the */
/*       correct rank and positive semi-definiteness. */
/*     NDIM is the first dimension of BMAT and has the value NPT+N. */
/*     SL and SU hold the differences XL-XBASE and XU-XBASE, respectively. */
/*       All the components of every XOPT are going to satisfy the bounds */
/*       SL(I) .LEQ. XOPT(I) .LEQ. SU(I), with appropriate equalities when */
/*       XOPT is on a constraint boundary. */
/*     XNEW is chosen by SUBROUTINE TRSBOX or ALTMOV. Usually XBASE+XNEW is the */
/*       vector of variables for the next call of CALFUN. XNEW also satisfies */
/*       the SL and SU constraints in the way that has just been mentioned. */
/*     XALT is an alternative to XNEW, chosen by ALTMOV, that may replace XNEW */
/*       in order to increase the denominator in the updating of UPDATE. */
/*     D is reserved for a trial step from XOPT, which is usually XNEW-XOPT. */
/*     VLAG contains the values of the Lagrange functions at a new point X. */
/*       They are part of a product that requires VLAG to be of length NDIM. */

__assume_aligned(xbase, 64);
__assume_aligned(fval, 64);
__assume_aligned(xopt, 64);
__assume_aligned(gopt, 64);
__assume_aligned(hq, 64);
__assume_aligned(pq, 64);
__assume_aligned(sl, 64);
__assume_aligned(su, 64);
__assume_aligned(xnew, 64);
__assume_aligned(xalt, 64);
__assume_aligned(d, 64);
__assume_aligned(vlag, 64);


    V_hdiag = _mm_malloc(sizeof(double)*(npt+8), 64);
    V_den = _mm_malloc(sizeof(double)*(npt+8), 64);
    V_distsq = _mm_malloc(sizeof(double)*(npt+8), 64);
    V_d__2 = _mm_malloc(sizeof(double)*(npt+8), 64);
    V_temp = _mm_malloc(sizeof(double)*(npt+8), 64);
    V_sum = _mm_malloc(sizeof(double)*(npt+8), 64);
    V_suma = _mm_malloc(sizeof(double)*(npt+8), 64);
    V_sumb = _mm_malloc(sizeof(double)*(npt+8), 64);
    
    /* Function Body */
    half = .5;
    one = 1.;
    ten = 10.;
    tenth = .1;
    two = 2.;
    zero = 0.;
    np = n + 1;
    nptm = npt - np;
    nh = n * np / 2;

#define ZERO 0.

    glag =  _mm_malloc(sizeof(double) * (U(np)+8), 64);       if (!glag) { rc = NLOPT_OUT_OF_MEMORY; goto done; }
    hcol =  _mm_malloc(sizeof(double) * (U(ndim)+8), 64);    if (!hcol) { rc = NLOPT_OUT_OF_MEMORY; goto done; }
    gnew =  _mm_malloc(sizeof(double) * (U(np)+8), 64);       if (!gnew) { rc = NLOPT_OUT_OF_MEMORY; goto done; }
    xbdi =  _mm_malloc(sizeof(double) * (U(n)+8), 64);       if (!xbdi) { rc = NLOPT_OUT_OF_MEMORY; goto done; }
    s    =  _mm_malloc(sizeof(double) * (U(n)+8), 64);       if (!s)    { rc = NLOPT_OUT_OF_MEMORY; goto done; }
    hs   =  _mm_malloc(sizeof(double) * (U(n)+8), 64);       if (!hs)   { rc = NLOPT_OUT_OF_MEMORY; goto done; } 
    hred =  _mm_malloc(sizeof(double) * (U(n)+8), 64);       if (!hred) { rc = NLOPT_OUT_OF_MEMORY; goto done; }
    ptsaux = Allocate2DArray(2, n);                             if (!ptsaux) { rc = NLOPT_OUT_OF_MEMORY; goto done; }
    ptsid =  _mm_malloc(sizeof(double) * (U(npt)+8), 64);    if (!ptsid) { rc = NLOPT_OUT_OF_MEMORY; goto done; }

    
/*     The call of PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ, */
/*     BMAT and ZMAT for the first iteration, with the corresponding values of */
/*     of NF and KOPT, which are the number of calls of CALFUN so far and the */
/*     index of the interpolation point at the trust region centre. Then the */
/*     initial XOPT is set too. The branch to label 720 occurs if MAXFUN is */
/*     less than NPT. GOPT will be updated if KOPT is different from KBASE. */

    rc2 = prelim_(n, npt, x, xl, xu, rhobeg, stop, calfun, calfun_data,
        xbase, xpt, fval, gopt, hq, pq, ndim, bmat, zmat, sl, su, &kopt);
    xoptsq = zero;
    for (i = 0; i < n; ++i) {
        xopt[i] = xpt[i][kopt];
/* L10: */
        xoptsq += xopt[i] * xopt[i];
    }
    fsave = fval[0];
    if (rc2 != NLOPT_SUCCESS) {
      rc = rc2;
      goto L720;
    }
    kbase = 0;

/*     Complete the settings that are required for the iterative procedure. */

    rho = *rhobeg;
    delta = rho;
    nresc = stop->nevals;
    ntrits = 0;
    diffa = zero;
    diffb = zero;
    itest = 0;
    nfsav = stop->nevals;

/*     Update GOPT if necessary before the first iteration and after each */
/*     call of RESCUE that makes a call of CALFUN. */

L20:
    if (kopt != kbase) {
        ih = 0;
        for (j = 0; j < n; ++j) {
            //i__2 = j;
            for (i = 0; i <= j; ++i) {
                ++ih;
                if (i < j) {
                    gopt[j] += hq[ih-1] * xopt[i];
                }
/* L30: */
                gopt[i] += hq[ih-1] * xopt[j];
            }
        }
        if (stop->nevals > npt) {
            //i__2 = npt;
            for (k = 0; k < npt; ++k) {
                temp = zero;
                for (j = 0; j < n; ++j) {
/* L40: */
                    temp += xpt[j][k] * xopt[j];
                }
                temp = pq[k] * temp;
                for (i = 0; i < n; ++i) {
/* L50: */
                    gopt[i] += temp * xpt[i][k];
                }
            }
        }
    }

/*     Generate the next point in the trust region that provides a small value */
/*     of the quadratic model subject to the constraints on the variables. */
/*     The int NTRITS is set to the number "trust region" iterations that */
/*     have occurred since the last "alternative" iteration. If the length */
/*     of XNEW-XOPT is less than HALF*RHO, however, then there is a branch to */
/*     label 650 or 680 with NTRITS=-1, instead of calculating F at XNEW. */

L60:
    trsbox_(n, npt, xpt, xopt, gopt, hq, pq, sl, su, &delta, xnew, d, 
        gnew, xbdi, s, hs, hred, &dsq, &crvmin);

    dnorm = MIN2(delta,sqrt(dsq));
    if (dnorm < half * rho) {
        ntrits = -1;
        distsq = (ten * rho) * (ten * rho);
        if (stop->nevals <= nfsav + 2) {
            goto L650;
        }

/*     The following choice between labels 650 and 680 depends on whether or */
/*     not our work with the current RHO seems to be complete. Either RHO is */
/*     decreased or termination occurs if the errors in the quadratic model at */
/*     the last three interpolation points compare favourably with predictions */
/*     of likely improvements to the model within distance HALF*RHO of XOPT. */

        d__1 = MAX2(diffa,diffb);
        errbig = MAX2(d__1,diffc);
        frhosq = rho * .125 * rho;
        if (crvmin > zero && errbig > frhosq * crvmin) {
            goto L650;
        }
        bdtol = errbig / rho;
        for (j = 0; j < n; ++j) {
            bdtest = bdtol;
            if (xnew[j] == sl[j]) {
                bdtest = w[j];
            }
            if (xnew[j] == su[j]) {
                bdtest = -w[j];
            }
            if (bdtest < bdtol) {
                curv = hq[(j + j * j) / 2];
                //i__2 = npt;
                for (k = 0; k < npt; ++k) {
/* L70: */
                    curv += pq[k] * (xpt[j][k] * xpt[j][k]);
                }
                bdtest += half * curv * rho;
                if (bdtest < bdtol) {
                    goto L650;
                }
            }
/* L80: */
        }    
        goto L680;
    }
    ++ntrits;

/*     Severe cancellation is likely to occur if XOPT is too far from XBASE. */
/*     If the following test holds, then XBASE is shifted so that XOPT becomes */
/*     zero. The appropriate changes are made to BMAT and to the second */
/*     derivatives of the current model, beginning with the changes to BMAT */
/*     that do not depend on ZMAT. VLAG is used temporarily for working space. */

L90:
    if (dsq <= xoptsq * .001) {
        fracsq = xoptsq * .25;
        sumpq = zero;
        for (k = 0; k < npt; ++k) 
            sumpq += pq[k];

        V_sum[0:npt] = -half * xoptsq;
        for (i = 0; i < n; ++i) {
            for (k = 0; k < npt; ++k) 
/* L100: */
               V_sum[k] += xpt[i][k] * xopt[i];
        }
        w_tp = &w[npt];
        w_tp[0:npt] = V_sum[0:npt];

        V_temp[0:npt] = fracsq - half * V_sum[0:npt];
        for (k = 0; k < npt; ++k) {
            for (i = 0; i < n; ++i) {
                w[i] = bmat[i][k];
                vlag[i] = V_sum[k] * xpt[i][k] + V_temp[k] * xopt[i];
                ip = npt + i;

                for (j = 0; j <= i; ++j) {
/* L110: */
                    bmat[j][ip] += w[i] * vlag[j] + vlag[i] * w[j];
                }
            }
        }

/*     Then the revisions of BMAT that depend on ZMAT are calculated. */

        w_tp = &w[npt];
        for (jj = 0; jj < nptm; ++jj) {
            sumz = zero;
            sumw = zero;
        
            for (k = 0; k < npt; ++k) {
                sumz += zmat[jj][k];
                vlag[k] = w_tp[k] * zmat[jj][k];
/* L120: */
                sumw += vlag[k];
            }
       
            for (j = 0; j < n; ++j) {
                sum = (fracsq * sumz - half * sumw) * xopt[j];
                
                for (k = 0; k < npt; ++k) {
/* L130: */
                    sum += vlag[k] * xpt[j][k];
                }
                w[j] = sum;
                
                for (k = 0; k < npt; ++k) {
/* L140: */
                    bmat[j][k] += sum * zmat[jj][k];
                }
            }
            
            for (i = 0; i < n; ++i) {
                ip = i + npt;
                temp = w[i];
          
                for (j = 0; j <= i; ++j) {
/* L150: */
                    bmat[j][ip] += temp * w[j];
                }
            }
        }

/*     The following instructions complete the shift, including the changes */
/*     to the second derivative parameters of the quadratic model. */

        ih = 0;

        for (j = 0; j < n; ++j) {
            w[j] = -half * sumpq * xopt[j];

            for (k = 0; k < npt; ++k) {
                w[j] += pq[k] * xpt[j][k];
/* L160: */
                xpt[j][k] -= xopt[j];
            }

            for (i = 0; i <= j; ++i) {
                ++ih;
                hq[ih-1] = hq[ih-1] + w[i] * xopt[j] + xopt[i] * w[j];
/* L170: */
                bmat[j][npt + i] = bmat[i][npt + j];
            }
        }

        for (i = 0; i < n; ++i) {
            xbase[i] += xopt[i];
            xnew[i] -= xopt[i];
            sl[i] -= xopt[i];
            su[i] -= xopt[i];
/* L180: */
            xopt[i] = zero;
        }
        xoptsq = zero;
    }
    if (ntrits == 0) {
        goto L210;
    }
    goto L230;

/*     XBASE is also moved to XOPT by a call of RESCUE. This calculation is */
/*     more expensive than the previous shift, because new matrices BMAT and */
/*     ZMAT are generated from scratch, which may include the replacement of */
/*     interpolation points whose positions seem to be causing near linear */
/*     dependence in the interpolation conditions. Therefore RESCUE is called */
/*     only if rounding errors have reduced by at least a factor of two the */
/*     denominator of the formula for updating the H matrix. It provides a */
/*     useful safeguard, but is not invoked in most applications of BOBYQA. */


L190:
    nfsav = stop->nevals;
    kbase = kopt;
    rc2 = rescue_(n, npt, xl, xu, 
          stop, calfun, calfun_data,
          xbase, xpt, fval, xopt, gopt,
          hq, pq, ndim, bmat, zmat,
          sl, su, delta, &kopt, vlag,
          ptsaux, ptsid);

/*     XOPT is updated now in case the branch below to label 720 is taken. */
/*     Any updating of GOPT occurs after the branch below to label 20, which */
/*     leads to a trust region iteration as does the branch to label 60. */

    xoptsq = zero;
    if (kopt != kbase) {
        for (i = 0; i < n; ++i) {
            xopt[i] = xpt[i][kopt];
/* L200: */
            xoptsq += xopt[i] * xopt[i];
        }
    }
    if (rc2 != NLOPT_SUCCESS) { 
        rc = rc2;
        goto L720; 
    }
    nresc = stop->nevals;
    if (nfsav < stop->nevals) {
        nfsav = stop->nevals;
        goto L20;
    }
    if (ntrits > 0) {
        goto L60;
    }

/*     Pick two alternative vectors of variables, relative to XBASE, that */
/*     are suitable as new positions of the KNEW-th interpolation point. */
/*     Firstly, XNEW is set to the point on a line through XOPT and another */
/*     interpolation point that minimizes the predicted value of the next */
/*     denominator, subject to ||XNEW - XOPT|| .LEQ. ADELT and to the SL */
/*     and SU bounds. Secondly, XALT is set to the best feasible point on */
/*     a constrained version of the Cauchy step of the KNEW-th Lagrange */
/*     function, the corresponding value of the square of this function */
/*     being returned in CAUCHY. The choice between these alternatives is */
/*     going to be made when the denominator is calculated. */

L210:
    altmov_(n, npt, xpt, xopt, ndim, bmat, zmat,
        sl, su, &kopt, &knew, &adelt, xnew,
        xalt, &alpha, &cauchy, glag, hcol);
  
    for (i = 0; i < n; ++i) {
/* L220: */
        d[i] = xnew[i] - xopt[i];
    }

// CPB
//printf("bobyqb_: n, npt, ndim, knew %d %d %d %d\n", n, npt, ndim, knew);

/*     Calculate VLAG and BETA for the current choice of D. The scalar */
/*     product of D with XPT(K,.) is going to be held in W(NPT+K) for */
/*     use when VQUAD is calculated. */

L230:
    V_suma[0:npt] = zero;
    V_sumb[0:npt] = zero;
    V_sum[0:npt] = zero;

    for (j = 0; j < n; ++j) {
        for (k = 0; k < npt; ++k) {
            V_suma[k] += xpt[j][k] * d[j];
            V_sumb[k] += xpt[j][k] * xopt[j];
/* L240: */
            V_sum[k] += bmat[j][k] * d[j];
        }
    }
    w_tp = &w[npt];

    for (k = 0; k < npt; ++k) {
        w[k] = V_suma[k] * (half * V_suma[k] + V_sumb[k]);
        vlag[k] = V_sum[k];
/* L250: */
        w_tp[k] = V_suma[k];
    }
    beta = zero;
    for (jj = 0; jj < nptm; ++jj) {
        sum = zero;
        for (k = 0; k < npt; ++k) {
/* L260: */
            sum += zmat[jj][k] * w[k];
        }
        beta -= sum * sum;
        
        for (k = 0; k < npt; ++k) {
/* L270: */
            vlag[k] += sum * zmat[jj][k];
        }
    }
    dsq = zero;
    bsum = zero;
    dx = zero;

    for (j = 0; j < n; ++j) {
        dsq += d[j] * d[j];
        sum = zero;
        
        for (k = 0; k < npt; ++k) {
/* L280: */
            sum += w[k] * bmat[j][k];
        }
        bsum += sum * d[j];
        jp = npt + j;

        for (i = 0; i < n; ++i) {
/* L290: */
            sum += bmat[i][jp] * d[i];
        }
        vlag[jp] = sum;
        bsum += sum * d[j];
/* L300: */
        dx += d[j] * xopt[j];
    }
    beta = dx * dx + dsq * (xoptsq + dx + dx + half * dsq) + beta - bsum;
    vlag[kopt] += one;

/*     If NTRITS is zero, the denominator may be increased by replacing */
/*     the step D of ALTMOV by a Cauchy step. Then RESCUE may be called if */
/*     rounding errors have damaged the chosen denominator. */

    double *xpt_row;
    if (ntrits == 0) {
        denom = vlag[knew] * vlag[knew] + alpha * beta;
        if (denom < cauchy && cauchy > zero) {
        
            for (i = 0; i < n; ++i) {
                xnew[i] = xalt[i];
/* L310: */
                d[i] = xnew[i] - xopt[i];
            }
            cauchy = zero;
            goto L230;
        }
        if (denom <= half * (vlag[knew] * vlag[knew])) {
            if (stop->nevals > nresc) {
                goto L190;
            }
        /* Return from BOBYQA because of much cancellation in a
           denominator. */
            rc = NLOPT_ROUNDOFF_LIMITED;
            goto L720;
        }

/*     Alternatively, if NTRITS is positive, then set KNEW to the index of */
/*     the next interpolation point to be deleted to make room for a trust */
/*     region step. Again RESCUE may be called if rounding errors have damaged */
/*     the chosen denominator, which is the reason for attempting to select */
/*     KNEW before calculating the next value of the objective function. */

    } 
    else {
        delsq = delta * delta;
        inv_delsq = 1.0 / delsq;
        scaden = zero;
        biglsq = zero;
        knew = 0;
        V_hdiag[0:npt] = zero;
       
        for (jj = 0; jj < nptm; ++jj) {
/* L330: */
            for (k = 0; k < npt; ++k) {
                V_hdiag[k] += zmat[jj][k] * zmat[jj][k];
            }
        }
        for (k = 0; k < npt; ++k) {
            V_den[k] = beta * V_hdiag[k] + vlag[k] * vlag[k];
        }
        
// CPB: Single division operation to change original divides to multiply op

        V_distsq[0:npt] = zero;
        for (j = 0; j < n; ++j) {
/* L340: */
            for (k = 0; k < npt; ++k) {
                temp = (xpt[j][k] - xopt[j]);
                V_distsq[k] += temp * temp;
                V_d__2[k] = (V_distsq[k] * inv_delsq) * (V_distsq[k] * inv_delsq);
                V_temp[k] = MAX2(one,V_d__2[k]);
            }
        }
        for (k = 0; k < npt; ++k) {
              if (k != kopt) { 
                 if (V_temp[k] * V_den[k] > scaden) {
                    scaden = V_temp[k] * V_den[k];
                    knew = k;
                    denom = V_den[k];
                 }
                 V_d__2[k] = V_temp[k] * (vlag[k] * vlag[k]);
                 biglsq = MAX2(biglsq,V_d__2[k]);
              }
/*
L350:
        ;
*/
        }
        if (scaden <= half * biglsq) {
            if (stop->nevals > nresc) {
                goto L190;
            }
        /* Return from BOBYQA because of much cancellation in a
           denominator. */
            rc = NLOPT_ROUNDOFF_LIMITED;
            goto L720;
        }
    }

/*     Put the variables for the next calculation of the objective function */
/*       in XNEW, with any adjustments for the bounds. */


/*     Calculate the value of the objective function at XBASE+XNEW, unless */
/*       the limit on the number of calculations of F has been reached. */

L360:

    for (i = 0; i < n; ++i) {
        d__1 = MAX2(xl[i],xbase[i] + xnew[i]);
        x[i] = MIN2(d__1,xu[i]);
        if (xnew[i] == sl[i]) {
            x[i] = xl[i];
        }
        if (xnew[i] == su[i]) {
            x[i] = xu[i];
        }
/* L380: */
    }

    if (nlopt_stop_forced(stop)) rc = NLOPT_FORCED_STOP;
    else if (nlopt_stop_evals(stop)) rc = NLOPT_MAXEVAL_REACHED;
    else if (nlopt_stop_time(stop)) rc = NLOPT_MAXTIME_REACHED;
    if (rc != NLOPT_SUCCESS) goto L720;

    stop->nevals++;
    f = calfun(n, x, calfun_data);
    if (ntrits == -1) {
        fsave = f;
        rc = NLOPT_XTOL_REACHED;
        if (fsave < fval[kopt]) { *minf = f; return rc; }
        goto L720;
    }

    if (f < stop->minf_max) {
      *minf = f;
      return NLOPT_MINF_MAX_REACHED;
    }

/*     Use the quadratic model to predict the change in F due to the step D, */
/*       and set DIFF to the error of this prediction. */

    fopt = fval[kopt];
    vquad = zero;
    ih = 0;

    for (j = 0; j < n; ++j) {
        vquad += d[j] * gopt[j];

        for (i = 0; i <= j; ++i) {
            ++ih;
            temp = d[i] * d[j];
            if (i == j) {
                temp = half * temp;
            }
/* L410: */
            vquad += hq[ih-1] * temp;
        }
    }

    w_tp = &w[npt];
    for (k = 0; k < npt; ++k) {
/* L420: */
        vquad += half * pq[k] * (w_tp[k] * w_tp[k]);
    }
    diff = f - fopt - vquad;
    diffc = diffb;
    diffb = diffa;
    diffa = fabs(diff);
    if (dnorm > rho) {
        nfsav = stop->nevals;
    }

/*     Pick the next value of DELTA after a trust region step. */

    if (ntrits > 0) {
        if (vquad >= zero) {
      /* Return from BOBYQA because a trust region step has failed
         to reduce Q. */
        rc = NLOPT_ROUNDOFF_LIMITED; /* or FTOL_REACHED? */
        goto L720;
        }
        ratio = (f - fopt) / vquad;
        if (ratio <= tenth) {
            delta = MIN2(half * delta,dnorm);
        } 
        else if (ratio <= .7) {
            delta = MAX2(half * delta,dnorm);
        } 
        else {
            delta = MAX2(half * delta,dnorm + dnorm);
        }
        if (delta <= rho * 1.5) {
            delta = rho;
        }

/*     Recalculate KNEW and DENOM if the new F is less than FOPT. */

        if (f < fopt) {
            ksav = knew;
            densav = denom;
            delsq = delta * delta;
            inv_delsq = 1.0 / delsq;
            scaden = zero;
            biglsq = zero;
            knew = 0;
            V_hdiag[0:npt] = zero;
            V_distsq[0:npt] = zero;
        
            for (jj = 0; jj < nptm; ++jj) {
                for (k = 0; k < npt; ++k) 
/* L440: */
                    V_hdiag[k] += zmat[jj][k] * zmat[jj][k];
            }

            for (k = 0; k < npt; ++k) 
                V_den[k] = beta * V_hdiag[k] + vlag[k] * vlag[k];
                
            for (j = 0; j < n; ++j) {
                for (k = 0; k < npt; ++k) {
/* L450: */
                    temp = (xpt[j][k] - xnew[j]);
                    V_distsq[k] += temp*temp;
                }
            }
            for (k = 0; k < npt; ++k) {
                V_d__2[k] = (V_distsq[k] * inv_delsq) * (V_distsq[k] * inv_delsq);
                V_temp[k] = MAX2(one,V_d__2[k]);
/* L460: */
                d__2 = V_temp[k] * (vlag[k] * vlag[k]);
                biglsq = MAX2(biglsq,d__2);
            }
            for (k = 0; k < npt; ++k) {
                if (V_temp[k] * V_den[k] > scaden) {
                    scaden = V_temp[k] * V_den[k];
                    knew = k;
                    denom = V_den[k];
                }
            }

            if (scaden <= half * biglsq) {
                knew = ksav;
                denom = densav;
            }
        }
    }

/*     Update BMAT and ZMAT, so that the KNEW-th interpolation point can be */
/*     moved. Also update the second derivative terms of the model. */

    update_(n, npt, ndim, bmat, zmat, vlag, beta, denom, knew);
    ih = 0;
    pqold = pq[knew];
    pq[knew] = zero;

    for (i = 0; i < n; ++i) {
        temp = pqold * xpt[i][knew];

        for (j = 0; j <= i; ++j) {
            ++ih;
/* L470: */
            hq[ih-1] += temp * xpt[j][knew];
        }
    }

    double *zmat_row;
    for (jj = 0; jj < nptm; ++jj) {
        temp = diff * zmat[jj][knew];
        zmat_row = zmat[jj];
__assume_aligned(zmat_row, 64);

        for (k = 0; k < npt; ++k) {
/* L480: */
            pq[k] += temp * zmat_row[k];
        }
    }

/*     Include the new interpolation point, and make the changes to GOPT at */
/*     the old XOPT that are caused by the updating of the quadratic model. */

    fval[knew] = f;

    for (i = 0; i < n; ++i) {
        xpt[i][knew] = xnew[i];
/* L490: */
        w[i] = bmat[i][knew];
    }

    V_suma[0:npt] = zero;
    for (jj = 0; jj < nptm; ++jj) {
        temp = zmat[jj][knew];
        zmat_row = zmat[jj];
__assume_aligned(zmat_row, 64);
        for (k = 0; k < npt; ++k) 
/* L500: */
            V_suma[k] += temp * zmat_row[k];
    }

    for (k = 0; k < npt; ++k) {
        if (nlopt_isinf(V_suma[k])) {
      /* SGJ: detect singularity here (happend if we run
         for too many iterations) ... is there another way to recover? */
            rc = NLOPT_ROUNDOFF_LIMITED;
            goto L720;
        }
    }

    V_sumb[0:npt] = zero;
    for (j = 0; j < n; ++j) {
        temp = xopt[j];
        xpt_row = xpt[j];
__assume_aligned(xpt_row, 64);
        for (k = 0; k < npt; ++k) 
/* L510: */
            V_sumb[k] += temp * xpt_row[k];
    }
    for (k = 0; k < npt; ++k) {
        V_temp[k] = V_suma[k] * V_sumb[k];
    }

    for (i = 0; i < n; ++i) {
        xpt_row = xpt[i];
__assume_aligned(xpt_row, 64);
        for (k = 0; k < npt; ++k) 
/* L520: */
            w[i] += V_temp[k] * xpt_row[k];
    }
    
    for (i = 0; i < n; ++i) {
/* L530: */
        gopt[i] += diff * w[i];
    }

/*     Update XOPT, GOPT and KOPT if the new calculated F is less than FOPT. */

    if (f < fopt) {
        kopt = knew;
        xoptsq = zero;
        ih = 0;

        for (j = 0; j < n; ++j) {
            xopt[j] = xnew[j];
            xoptsq += xopt[j] * xopt[j];
            
            for (i = 0; i <= j; ++i) {
                ++ih;
                if (i < j) {
                    gopt[j] += hq[ih-1] * d[i];
                }
/* L540: */
                gopt[i] += hq[ih-1] * d[j];
            }
        }

        V_temp[0:npt] = zero;
        for (j = 0; j < n; ++j) {
            for (k = 0; k < npt; ++k)
/* L550: */
                V_temp[k] += xpt[j][k] * d[j];
        }
        for (k = 0; k < npt; ++k)
            V_temp[k] *= pq[k];
  
//          V_temp[k] = pq[k] * V_temp[k];
       
        for (i = 0; i < n; ++i) {
            for (k = 0; k < npt; ++k)
/* L560: */
                gopt[i] += V_temp[k] * xpt[i][k];
        }
        if (nlopt_stop_ftol(stop, f, fopt)) {
            rc = NLOPT_FTOL_REACHED;
            goto L720;
        }
    }

/*     Calculate the parameters of the least Frobenius norm interpolant to */
/*     the current data, the gradient of this interpolant at XOPT being put */
/*     into VLAG(NPT+I), I=1,2,...,N. */

    if (ntrits > 0) {
 
        for (k = 0; k < npt; ++k) {
            vlag[k] = fval[k] - fval[kopt];
/* L570: */
            w[k] = zero;
        }

        for (j = 0; j < nptm; ++j) {
            sum = zero;

            for (k = 0; k < npt; ++k) {
/* L580: */
                sum += zmat[j][k] * vlag[k];
            }

            for (k = 0; k < npt; ++k) {
/* L590: */
                w[k] += sum * zmat[j][k];
            }
        }

        V_sum[0:npt] = 0;
        for (j = 0; j <n; ++j) {
            for (k = 0; k < npt; ++k)
/* L600: */
                V_sum[k] += xpt[j][k] * xopt[j];
        }
        w_tp = &w[npt];

        for (k = 0; k < npt; ++k) {
            w_tp[k] = w[k];
/* L610: */
            w[k] = sum * w[k];
        }
        gqsq = zero;
        gisq = zero;

double *bmat_row;

        for (i = 0; i < n; ++i) {
            sum = zero;

            bmat_row = bmat[i];
            xpt_row = xpt[i];
__assume_aligned(bmat_row, 64);
__assume_aligned(xpt_row, 64);

            for (k = 0; k < npt; ++k) {
/* L620: */
                //sum += bmat[i][k] * vlag[k] + xpt[i][k] * w[k];
                sum += bmat_row[k] * vlag[k] + xpt_row[k] * w[k];
            }
            if (xopt[i] == sl[i]) {
                d__1 = MIN2(zero,gopt[i]);
                gqsq += d__1 * d__1;
                d__1 = MIN2(zero,sum);
                gisq += d__1 * d__1;
            } 
            else if (xopt[i] == su[i]) {
                d__1 = MAX2(zero,gopt[i]);
                gqsq += d__1 * d__1;
                d__1 = MAX2(zero,sum);
                gisq += d__1 * d__1;
            } 
            else {
                gqsq += gopt[i] * gopt[i];
                gisq += sum * sum;
            }
/* L630: */
            vlag[npt + i] = sum;
        }

/*     Test whether to replace the new quadratic model by the least Frobenius */
/*     norm interpolant, making the replacement if the test is satisfied. */

        ++itest;
        if (gqsq < ten * gisq) {
            itest = 0;
        }
        if (itest >= 3) {
            d__1 = MAX2(npt,nh);
            for (i = 0; i < d__1; ++i) {
                if (i < n) {        /* CPB: was (i <= n) */
                    gopt[i] = vlag[npt + i];
                }
                if (i < npt) {        /* CPB: was (i <= npt) */
                    pq[i] = w[npt + i];
                }
                if (i < nh) {        /* CPB: was (i <= nh) */
                    hq[i] = zero;
                }
/* L640: */
            }
            itest = 0;         /* CPB: factored out of loop */
        }
    }

/*     If a trust region step has provided a sufficient decrease in F, then */
/*     branch for another trust region calculation. The case NTRITS=0 occurs */
/*     when the new interpolation point was reached by an alternative step. */

    if (ntrits == 0) {
        goto L60;
    }
    if (f <= fopt + tenth * vquad) {
        goto L60;
    }

/*     Alternatively, find out if the interpolation points are close enough */
/*       to the best point so far. */

    d__1 = (two * delta) * (two * delta);
    d__2 = (ten * rho) * (ten * rho);
    distsq = MAX2(d__1,d__2);
L650:
    knew = 0;

    V_sum[0:npt] = zero;

    for (j = 0; j < n; ++j) {
        for (k = 0; k < npt; ++k) {
/* L660: */
            temp = (xpt[j][k] - xopt[j]);
            V_sum[k] += temp * temp;
        }
    }
    for (k = 0; k < npt; ++k) {
/* L670: */
        if (V_sum[k] > distsq) {
            knew = k;
            distsq = V_sum[k];
        }
    }

/*     If KNEW is positive, then ALTMOV finds alternative new positions for */
/*     the KNEW-th interpolation point within distance ADELT of XOPT. It is */
/*     reached via label 90. Otherwise, there is a branch to label 60 for */
/*     another trust region iteration, unless the calculations with the */
/*     current RHO are complete. */

    if (knew > 0) {
        dist = sqrt(distsq);
        if (ntrits == -1) {
            delta = MIN2(tenth * delta,half * dist);
            if (delta <= rho * 1.5) {
                delta = rho;
            }
        }
        ntrits = 0;
        d__1 = MIN2(tenth * dist,delta);
        adelt = MAX2(d__1,rho);
        dsq = adelt * adelt;
        goto L90;
    }
    if (ntrits == -1) {
        goto L680;
    }
    if (ratio > zero) {
        goto L60;
    }
    if (MAX2(delta,dnorm) > rho) {
        goto L60;
    }

/*     The calculations with the current value of RHO are complete. Pick the */
/*       next values of RHO and DELTA. */

L680:
    if (rho > *rhoend) {
        delta = half * rho;
        ratio = rho / *rhoend;
        if (ratio <= 16.) {
            rho = *rhoend;
        }     
        else if (ratio <= 250.) {
            rho = sqrt(ratio) * *rhoend;
        } 
        else {
            rho = tenth * rho;
        }
        delta = MAX2(delta,rho);
        ntrits = 0;
        nfsav = stop->nevals;
        goto L60;
    }

/*     Return from the calculation, after another Newton-Raphson step, if */
/*       it is too short to have been tried before. */

    if (ntrits == -1) {
        goto L360;
    }
L720:
    /* originally: if (fval[kopt] <= fsave) -- changed by SGJ, since
       this seems like a slight optimization to not update x[]
       unnecessarily, at the expense of possibly not returning the
       best x[] found so far if the algorithm is stopped suddenly
       (e.g. runs out of time) ... it seems safer to execute this
       unconditionally, and the efficiency loss seems negligible. */
    
    for (i = 0; i < n; ++i) {
        d__1 = MAX2(xl[i],xbase[i] + xopt[i]);
        x[i] = MIN2(d__1,xu[i]);
        if (xopt[i] == sl[i]) {
            x[i] = xl[i];
        }
        if (xopt[i] == su[i]) {
            x[i] = xu[i];
        }
/* L730: */
    }
    f = fval[kopt];
    *minf = f;
done:
    _mm_free(glag);
    _mm_free(hcol);
    _mm_free(gnew);
    _mm_free(xbdi);
    _mm_free(s);
    _mm_free(hs);
    _mm_free(hred);
    Free2DArray(ptsaux);
    _mm_free(ptsid);

    _mm_free(w);
    _mm_free(V_hdiag);
    _mm_free(V_den);
    _mm_free(V_distsq);
    _mm_free(V_d__2);
    _mm_free(V_temp);
    _mm_free(V_sum);
    _mm_free(V_suma);
    _mm_free(V_sumb);
    return rc;
} /* bobyqb_ */

/**************************************************************************/


typedef struct {
     double *s, *xs;
     nlopt_func f; void *f_data;
} rescale_fun_data;

static double rescale_fun(int n, const double *x, void *d_)
{
     rescale_fun_data *d = (rescale_fun_data*) d_;
     nlopt_unscale(U(n), d->s, x, d->xs);
     return d->f(U(n), d->xs, NULL, d->f_data);
}

nlopt_result bobyqa(int n, int npt, double *x, 
            const double *xl, const double *xu, 
            const double *dx,
            nlopt_stopping *stop, double *minf,
            nlopt_func f, void *f_data)
{
    /* System generated locals */
    int i__1;
    double d__1, d__2;

    /* Local variables */
    int j, id, np, iw, igo, ihq, ixb, ixa, ifv, isl, jsl, ipq, ivl, ixn, 
        ixo, ixp, isu, jsu, ndim;
    double temp, zero;
    int ibmat, izmat;

    double rhobeg, rhoend;
    double *w0 = NULL;
    nlopt_result ret;
    double *s = NULL, *sxl = NULL, *sxu = NULL, *xs = NULL;
    rescale_fun_data calfun_data;
/* CPB: Add arrays for bobyqb parameters */
    double *xbase, *fval, *xopt, *gopt, *hq, *pq;
    /*double *xbase, **xpt, *fval, *xopt, *gopt, *hq, *pq, **bmat, **zmat;*/
    double *sl, *su, *xnew, *xalt, *d, *vlag;
    
    double **xpt;
    double **bmat;
    double **zmat;

    ndim = npt + n;
    np = n + 1;



    /* SGJ 2010: rescale parameters to make the initial step sizes dx
                 equal in all directions */
    s = nlopt_compute_rescaling(U(n), dx);
    if (!s) return NLOPT_OUT_OF_MEMORY;

    /* this statement must go before goto done, so that --x occurs */
    nlopt_rescale(U(n), s, x, x); /*--x;*/

    xs =  malloc(sizeof(double) * (U(n)));

    sxl = nlopt_new_rescaled(U(n), s, xl);
    if (!sxl) { ret = NLOPT_OUT_OF_MEMORY; goto done; }
    xl = sxl;
    sxu = nlopt_new_rescaled(U(n), s, xu);
    if (!sxu) { ret = NLOPT_OUT_OF_MEMORY; goto done; }
    xu = sxu;
    
	// ADDED nlopt 2.4.1
	nlopt_reorder_bounds(n, sxl, sxu);

/*
    xbase =  malloc(sizeof(double) * (U(n)));      if (!xbase) { ret = NLOPT_OUT_OF_MEMORY; goto done; }
    xpt   = Allocate2DArray(n, npt);                         if (!xpt) { ret = NLOPT_OUT_OF_MEMORY; goto done; }
    fval  = malloc(sizeof(double) * (U(npt)));    if (!fval) { ret = NLOPT_OUT_OF_MEMORY; goto done; }
    xopt  = malloc(sizeof(double) * (U(n)));      if (!xopt) { ret = NLOPT_OUT_OF_MEMORY; goto done; }
    gopt  = malloc(sizeof(double) * (U(n)));      if (!gopt) { ret = NLOPT_OUT_OF_MEMORY; goto done; }
    hq    = malloc(sizeof(double) * (U(n*np/2))); if (!hq) { ret = NLOPT_OUT_OF_MEMORY; goto done; }
    pq    = malloc(sizeof(double) * (U(npt)));    if (!pq) { ret = NLOPT_OUT_OF_MEMORY; goto done; }
    bmat  = Allocate2DArray(n, ndim);                        if (!bmat) { ret = NLOPT_OUT_OF_MEMORY; goto done; }
    zmat  = Allocate2DArray(npt-np, npt);                    if (!zmat) { ret = NLOPT_OUT_OF_MEMORY; goto done; }
    sl    = malloc(sizeof(double) * (U(n)));      if (!sl) { ret = NLOPT_OUT_OF_MEMORY; goto done; }
    su    = malloc(sizeof(double) * (U(n)));      if (!su) { ret = NLOPT_OUT_OF_MEMORY; goto done; }
    xnew  = malloc(sizeof(double) * (U(n)));      if (!xnew) { ret = NLOPT_OUT_OF_MEMORY; goto done; }
    xalt  = malloc(sizeof(double) * (U(n)));      if (!xalt) { ret = NLOPT_OUT_OF_MEMORY; goto done; }
    d     = malloc(sizeof(double) * (U(n)));      if (!d) { ret = NLOPT_OUT_OF_MEMORY; goto done; }
    vlag  = malloc(sizeof(double) * (U(ndim)));   if (!vlag) { ret = NLOPT_OUT_OF_MEMORY; goto done; }
*/
    
    xbase =  _mm_malloc(sizeof(double) * (U(n))+8, 64);      if (!xbase) { ret = NLOPT_OUT_OF_MEMORY; goto done; }
    xpt   = Allocate2DArray(n, npt);                         if (!xpt) { ret = NLOPT_OUT_OF_MEMORY; goto done; }
    fval  = _mm_malloc(sizeof(double) * (U(npt))+8, 64);    if (!fval) { ret = NLOPT_OUT_OF_MEMORY; goto done; }
    xopt  = _mm_malloc(sizeof(double) * (U(n))+8, 64);      if (!xopt) { ret = NLOPT_OUT_OF_MEMORY; goto done; }
    gopt  = _mm_malloc(sizeof(double) * (U(n))+8, 64);      if (!gopt) { ret = NLOPT_OUT_OF_MEMORY; goto done; }
    hq    = _mm_malloc(sizeof(double) * (U(n*np/2))+8, 64); if (!hq) { ret = NLOPT_OUT_OF_MEMORY; goto done; }
    pq    = _mm_malloc(sizeof(double) * (U(npt))+8, 64);    if (!pq) { ret = NLOPT_OUT_OF_MEMORY; goto done; }
    bmat  = Allocate2DArray(n, ndim);                        if (!bmat) { ret = NLOPT_OUT_OF_MEMORY; goto done; }
    zmat  = Allocate2DArray(npt-np, npt);                    if (!zmat) { ret = NLOPT_OUT_OF_MEMORY; goto done; }
    sl    = _mm_malloc(sizeof(double) * (U(n))+8, 64);      if (!sl) { ret = NLOPT_OUT_OF_MEMORY; goto done; }
    su    = _mm_malloc(sizeof(double) * (U(n))+8, 64);      if (!su) { ret = NLOPT_OUT_OF_MEMORY; goto done; }
    xnew  = _mm_malloc(sizeof(double) * (U(n))+8, 64);      if (!xnew) { ret = NLOPT_OUT_OF_MEMORY; goto done; }
    xalt  = _mm_malloc(sizeof(double) * (U(n))+8, 64);      if (!xalt) { ret = NLOPT_OUT_OF_MEMORY; goto done; }
    d     = _mm_malloc(sizeof(double) * (U(n))+8, 64);      if (!d) { ret = NLOPT_OUT_OF_MEMORY; goto done; }
    vlag  = _mm_malloc(sizeof(double) * (U(ndim))+8, 64);   if (!vlag) { ret = NLOPT_OUT_OF_MEMORY; goto done; }
    
	// NLOPT 2.4.1     rhobeg = dx[0] / s[0]; /* equals all other dx[i] after rescaling */
	rhobeg = fabs(dx[0] / s[0]);

    calfun_data.s = s;
    calfun_data.xs = xs;
    calfun_data.f = f;
    calfun_data.f_data = f_data;

    /* SGJ, 2009: compute rhoend from NLopt stop info */
    rhoend = stop->xtol_rel * (rhobeg);

	// NLOPT 2.4.1
	//for (j = 0; j < n; ++j)
    // if (rhoend < stop->xtol_abs[j] / s[j])
    //      rhoend = stop->xtol_abs[j] / s[j];
	// replaced by:
    for (j = 0; j < n; ++j)
	 if (rhoend < stop->xtol_abs[j] / fabs(s[j]))
	      rhoend = stop->xtol_abs[j] / fabs(s[j]);


/*     This subroutine seeks the least value of a function of many variables, */
/*     by applying a trust region method that forms quadratic models by */
/*     interpolation. There is usually some freedom in the interpolation */
/*     conditions, which is taken up by minimizing the Frobenius norm of */
/*     the change to the second derivative of the model, beginning with the */
/*     zero matrix. The values of the variables are constrained by upper and */
/*     lower bounds. The arguments of the subroutine are as follows. */

/*     N must be set to the number of variables and must be at least two. */
/*     NPT is the number of interpolation conditions. Its value must be in */
/*       the interval [N+2,(N+1)(N+2)/2]. Choices that exceed 2*N+1 are not */
/*       recommended. */
/*     Initial values of the variables must be set in X(1),X(2),...,X(N). They */
/*       will be changed to the values that give the least calculated F. */
/*     For I=1,2,...,N, XL(I) and XU(I) must provide the lower and upper */
/*       bounds, respectively, on X(I). The construction of quadratic models */
/*       requires XL(I) to be strictly less than XU(I) for each I. Further, */
/*       the contribution to a model from changes to the I-th variable is */
/*       damaged severely by rounding errors if XU(I)-XL(I) is too small. */
/*     RHOBEG and RHOEND must be set to the initial and final values of a trust */
/*       region radius, so both must be positive with RHOEND no greater than */
/*       RHOBEG. Typically, RHOBEG should be about one tenth of the greatest */
/*       expected change to a variable, while RHOEND should indicate the */
/*       accuracy that is required in the final values of the variables. An */
/*       error return occurs if any of the differences XU(I)-XL(I), I=1,...,N, */
/*       is less than 2*RHOBEG. */
/*     MAXFUN must be set to an upper bound on the number of calls of CALFUN. */
/*     The array W will be used for working space. Its length must be at least */
/*       (NPT+5)*(NPT+N)+3*N*(N+5)/2. */

/*     SUBROUTINE CALFUN (N,X,F) has to be provided by the user. It must set */
/*     F to the value of the objective function for the current values of the */
/*     variables X(1),X(2),...,X(N), which are generated automatically in a */
/*     way that satisfies the bounds given in XL and XU. */

/*     Return if the value of NPT is unacceptable. */

    /* Parameter adjustments */
    /* CPB: Not needed if using 0-based addressing for C arrays */
//    --xu;
//    --xl;

    /* Function Body */
    if (npt < n + 2 || npt > (n + 2) * np / 2) {
      /* Return from BOBYQA because NPT is not in the required interval */
//printf("NPT is not in the required interval %d %d %d \n" , npt , n , np);
      ret = NLOPT_INVALID_ARGS;
      goto done;
    }

/*     Partition the working space array, so that different parts of it can */
/*     be treated separately during the calculation of BOBYQB. The partition */
/*     requires the first (NPT+2)*(NPT+N)+3*N*(N+5)/2 elements of W plus the */
/*     space that is taken by the last array in the argument list of BOBYQB. */

//    w0 = _mm_malloc(sizeof(double) * U((npt+5)*(npt+n)+3*n*(n+5)/2), 64);
//    w0 = malloc(sizeof(double) * U((npt+5)*(npt+n)+3*n*(n+5)/2));
//    if (!w0) { ret = NLOPT_OUT_OF_MEMORY; goto done; }

/*   Return if there is insufficient space between the bounds. Modify the */
/*   initial X if necessary in order to avoid conflicts between the bounds */
/*   and the construction of the first quadratic model. The lower and upper */
/*   bounds on moves from the updated X are set now, in the ISL and ISU */
/*   partitions of W, in order to provide useful and exact information about */
/*   components of X that become within distance RHOBEG from their bounds. */

    zero = 0.;
    for (j = 0; j < n; ++j) {
        temp = xu[j] - xl[j];
        if (temp < rhobeg + rhobeg) {
      /* Return from BOBYQA because one of the differences
         XU(I)-XL(I)s is less than 2*RHOBEG. */
            ret = NLOPT_INVALID_ARGS;
            goto done;
        }
//        jsl = isl + j - 1;  // not needed since separate arrays for sl and su are used
//        jsu = jsl + n;
//        w[jsl] = xl[j] - x[j];  /* substitute sl[j] for w[jsl]  */
//        w[jsu] = xu[j] - x[j];  /* substitute su[j] for w[jsu]    */
        sl[j] = xl[j] - x[j];
        su[j] = xu[j] - x[j];        
        if (sl[j] >= -(rhobeg)) {
            if (sl[j] >= zero) {
                x[j] = xl[j];
                sl[j] = zero;
                su[j] = temp;
            } 
            else {
                x[j] = xl[j] + rhobeg;
                sl[j] = -(rhobeg);
                d__1 = xu[j] - x[j];
                su[j] = MAX2(d__1,rhobeg);
            }
        } 
        else if (su[j] <= rhobeg) {
            if (su[j] <= zero) {
                x[j] = xu[j];
                sl[j] = -temp;
                su[j] = zero;
            } 
            else {
                x[j] = xu[j] - rhobeg;
                d__1 = xl[j] - x[j], d__2 = -(rhobeg);
                sl[j] = MIN2(d__1,d__2);
                su[j] = rhobeg;
            }
        }
/* L30: */
    }

/*     Make the call of BOBYQB. */

    ret = bobyqb_(n, npt, x, xl, xu, &rhobeg, &rhoend,
          stop, rescale_fun, &calfun_data, minf,
          xbase, xpt, fval, xopt, gopt, hq, pq, 
          ndim, bmat, zmat, sl, su, xnew, xalt,
          d, vlag);
          
/* CPB: Is the final workspace parameter needed? Can we bury that in bobyqb_ ? */          

done:
    //if (w0) _mm_free(w0);
    if (sxl) free(sxl);
    if (sxu) free(sxu);
    if (xs) free(xs);
    //++x; 
    nlopt_unscale(U(n), s, x, x);
    if (s) free(s);
    if (xbase) _mm_free(xbase);
    if (xpt) Free2DArray(xpt);
    if (fval) _mm_free(fval);
    if (xopt) _mm_free(xopt);
    if (gopt) _mm_free(gopt);
    if (hq) _mm_free(hq);
    if (pq) _mm_free(pq);
    if (bmat) Free2DArray(bmat);
    if (zmat) Free2DArray(zmat);
    if (sl) _mm_free(sl);
    if (su) _mm_free(su);
    if (xnew) _mm_free(xnew);
    if (xalt) _mm_free(xalt);
    if (d) _mm_free(d);
    if (vlag) _mm_free(vlag); 
    
    return ret;
} /* bobyqa_ */


