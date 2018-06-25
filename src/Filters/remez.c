/* SIGTOOLS module by Travis Oliphant

Copyright 2005 Travis Oliphant
Permission to use, copy, modify, and distribute this software without fee
is granted under the SciPy License.
*/

#include "remez.h"
#include <setjmp.h>
#include <stdlib.h>

#include <stdio.h>
#include <math.h>


jmp_buf MALLOC_FAIL;

char *check_malloc(size_t size)
{
    char *the_block = malloc(size);
    if (the_block == NULL) {
        printf("\nERROR: unable to allocate %zu bytes!\n", size);
        longjmp(MALLOC_FAIL,-1);
    }
    return the_block;
}


/************************************************************************
 * Start of portable, non-python specific routines.                     *
 ************************************************************************/

/* Some core routines are written
in a portable way so that they could be used in other applications.  The 
order filtering, however uses python-specific constructs in its guts 
and is therefore Python dependent.  This could be changed in a 
straightforward way but I haven't done it for lack of time.*/

static int index_out_of_bounds(intp *indices, intp *max_indices, int ndims) {
  int bad_index = 0, k = 0;

  while (!bad_index && (k++ < ndims)) {
    bad_index = ((*(indices) >= *(max_indices++)) || (*(indices) < 0));
    indices++;
  }
  return bad_index;
}

/* This maybe could be redone with stride information so it could be 
 * called with non-contiguous arrays:  I think offsets is related to 
 * the difference between the strides.  I'm not sure about init_offset 
 * just yet.  I think it needs to be calculated because of mode_dep
 * but probably with dim1 being the size of the "original, unsliced" array
 */

static intp compute_offsets (uintp *offsets, intp *offsets2, intp *dim1,
                             intp *dim2, intp *dim3, intp *mode_dep,
                             int nd) {
  int k,i;
  intp init_offset = 0;

  for (k = 0; k < nd - 1; k++) 
    {
      init_offset += mode_dep[k];
      init_offset *= dim1[k+1];
    }
  init_offset += mode_dep[k] - 2;
  
  k = nd;
  while(k--) {
    offsets[k] = 0;
    offsets2[k] = 0;
    for (i = k + 1; i < nd - 1; i++) {
      offsets[k] += dim1[i] - dim2[i];
      offsets[k] *= dim1[i+1];

      offsets2[k] += dim1[i] - dim3[i];
      offsets2[k] *= dim1[i+1];
    }

    if (k < nd - 1) {
      offsets[k] += dim1[i] - dim2[i];
      offsets2[k] += dim1[i] - dim3[i];
    }
    offsets[k] += 1;
    offsets2[k] += 1;
  }
  return init_offset;
}

/* increment by 1 the index into an N-D array, doing the necessary
   carrying when the index reaches the dimension along that axis */ 
static int increment(intp *ret_ind, int nd, intp *max_ind) {    
    int k, incr = 1;
    
    k = nd - 1;
    if (++ret_ind[k] >= max_ind[k]) {
      while (k >= 0 && (ret_ind[k] >= max_ind[k]-1)) {
	incr++;
	ret_ind[k--] = 0;
      }
      if (k >= 0) ret_ind[k]++;
    }
    return incr;
}

/********************************************************
 *
 *  Code taken from remez.c by Erik Kvaleberg which was 
 *    converted from an original FORTRAN by
 *
 * AUTHORS: JAMES H. MCCLELLAN
 *
 *         DEPARTMENT OF ELECTRICAL ENGINEERING AND COMPUTER SCIENCE
 *         MASSACHUSETTS INSTITUTE OF TECHNOLOGY
 *         CAMBRIDGE, MASS. 02139
 *
 *         THOMAS W. PARKS
 *         DEPARTMENT OF ELECTRICAL ENGINEERING
 *         RICE UNIVERSITY
 *         HOUSTON, TEXAS 77001
 *
 *         LAWRENCE R. RABINER
 *         BELL LABORATORIES
 *         MURRAY HILL, NEW JERSEY 07974
 *
 *  
 *  Adaptation to C by 
 *      egil kvaleberg
 *      husebybakken 14a
 *      0379 oslo, norway
 *  Email:
 *      egil@kvaleberg.no
 *  Web:
 *      http://www.kvaleberg.com/
 * 
 *
 *********************************************************/


#define BANDPASS       1
#define DIFFERENTIATOR 2
#define HILBERT        3

#define GOBACK goto
#define DOloop(a,from,to) for ( (a) = (from); (a) <= (to); ++(a))
#define PI    3.14159265358979323846
#define TWOPI (PI+PI)

/*
 *-----------------------------------------------------------------------
 * FUNCTION: lagrange_interp (d)
 *  FUNCTION TO CALCULATE THE LAGRANGE INTERPOLATION
 *  COEFFICIENTS FOR USE IN THE FUNCTION gee.
 *-----------------------------------------------------------------------
 */
static double lagrange_interp(int k, int n, int m, double *x)
{
    int j, l;
    double q, retval;

    retval = 1.0;
    q = x[k];
    DOloop(l,1,m) {
	for (j = l; j <= n; j += m) {
	    if (j != k)
		retval *= 2.0 * (q - x[j]);
	}
    }
    return 1.0 / retval;
}

/*
 *-----------------------------------------------------------------------
 * FUNCTION: freq_eval (gee)
 *  FUNCTION TO EVALUATE THE FREQUENCY RESPONSE USING THE
 *  LAGRANGE INTERPOLATION FORMULA IN THE BARYCENTRIC FORM
 *-----------------------------------------------------------------------
 */
static double freq_eval(int k, int n, double *grid, double *x, double *y, double *ad)
{
    int j;
    double p,c,d,xf;

    d = 0.0;
    p = 0.0;
    xf = cos(TWOPI * grid[k]);

    DOloop(j,1,n) {
	c = ad[j] / (xf - x[j]);
	d += c;
	p += c * y[j];
    }

    return p/d;
}


/*
 *-----------------------------------------------------------------------
 * SUBROUTINE: remez
 *  THIS SUBROUTINE IMPLEMENTS THE REMEZ EXCHANGE ALGORITHM
 *  FOR THE WEIGHTED CHEBYSHEV APPROXIMATION OF A CONTINUOUS
 *  FUNCTION WITH A SUM OF COSINES.  INPUTS TO THE SUBROUTINE
 *  ARE A DENSE GRID WHICH REPLACES THE FREQUENCY AXIS, THE
 *  DESIRED FUNCTION ON THIS GRID, THE WEIGHT FUNCTION ON THE
 *  GRID, THE NUMBER OF COSINES, AND AN INITIAL GUESS OF THE
 *  EXTREMAL FREQUENCIES.  THE PROGRAM MINIMIZES THE CHEBYSHEV
 *  ERROR BY DETERMINING THE BSMINEST LOCATION OF THE EXTREMAL
 *  FREQUENCIES (POINTS OF MAXIMUM ERROR) AND THEN CALCULATES
 *  THE COEFFICIENTS OF THE BEST APPROXIMATION.
 *-----------------------------------------------------------------------
 */
static int remez(double *dev, double des[], double grid[], double edge[],  
	   double wt[], int ngrid, int nbands, int iext[], double alpha[],
	   int nfcns, int itrmax, double *work, int dimsize)
		/* dev, iext, alpha                         are output types */
		/* des, grid, edge, wt, ngrid, nbands, nfcns are input types */
{
    int k, k1, kkk, kn, knz, klow, kup, nz, nzz, nm1;
    int cn;
    int j, jchnge, jet, jm1, jp1;
    int l, luck=0, nu, nut, nut1=0, niter;

    double ynz=0.0, comp=0.0, devl, gtemp, fsh, y1=0.0, err, dtemp, delf, dnum, dden;
    double aa=0.0, bb=0.0, ft, xe, xt;

    static double *a, *p, *q;
    static double *ad, *x, *y;

    a = work; p = a + dimsize+1; q = p + dimsize+1; 
    ad = q + dimsize+1; x = ad + dimsize+1; y = x + dimsize+1;
    devl = -1.0;
    nz  = nfcns+1;
    nzz = nfcns+2;
    niter = 0;

    do {
    L100:
	iext[nzz] = ngrid + 1;
	++niter;

	if (niter > itrmax) break;

	printf("ITERATION %2d: ",niter);

	DOloop(j,1,nz) {
	    x[j] = cos(grid[iext[j]]*TWOPI);
	    printf("  j=%02d: iext[j]=%d grid[iext[j]]=%g x[j]=%g\n", j, iext[j], grid[iext[j]], x[j]);
	}
	jet = (nfcns-1) / 15 + 1;

	DOloop(j,1,nz) {
	    ad[j] = lagrange_interp(j,nz,jet,x);
	}

	dnum = 0.0;
	dden = 0.0;
	k = 1;

	DOloop(j,1,nz) {
	    l = iext[j];
	    dnum += ad[j] * des[l];
	    dden += (double)k * ad[j] / wt[l];
	    k = -k;
	}
	*dev = dnum / dden;

	printf("DEVIATION = %1.16f\n",*dev);

	nu = 1;
	if ( (*dev) > 0.0 ) nu = -1;
	(*dev) = -(double)nu * (*dev);
	k = nu;
	DOloop(j,1,nz) {
	    l = iext[j];
	    y[j] = des[l] + (double)k * (*dev) / wt[l];
	    k = -k;
	}
	if ( (*dev) <= devl ) {
	    /* finished */
	    return -1;
	}
	devl = (*dev);
	jchnge = 0;
	k1 = iext[1];
	knz = iext[nz];
	klow = 0;
	nut = -nu;
	j = 1;

    /*
     * SEARCH FOR THE EXTREMAL FREQUENCIES OF THE BEST APPROXIMATION
     */

    L200:
	if (j == nzz) ynz = comp;
	if (j >= nzz) goto L300;
	kup = iext[j+1];
	l = iext[j]+1;
	nut = -nut;
	if (j == 2) y1 = comp;
	comp = (*dev);
	if (l >= kup) goto L220;
	err = (freq_eval(l,nz,grid,x,y,ad)-des[l]) * wt[l];
	if (((double)nut*err-comp) <= 0.0) goto L220;
	comp = (double)nut * err;
    L210:
	if (++l >= kup) goto L215;
	err = (freq_eval(l,nz,grid,x,y,ad)-des[l]) * wt[l];
	if (((double)nut*err-comp) <= 0.0) goto L215;
	comp = (double)nut * err;
	GOBACK L210;

    L215:
	iext[j++] = l - 1;
	klow = l - 1;
	++jchnge;
	GOBACK L200;

    L220:
	--l;
    L225:
	if (--l <= klow) goto L250;
	err = (freq_eval(l,nz,grid,x,y,ad)-des[l]) * wt[l];
	if (((double)nut*err-comp) > 0.0) goto L230;
	if (jchnge <= 0) goto L225;
	goto L260;

    L230:
	comp = (double)nut * err;
    L235:
	if (--l <= klow) goto L240;
	err = (freq_eval(l,nz,grid,x,y,ad)-des[l]) * wt[l];
	if (((double)nut*err-comp) <= 0.0) goto L240;
	comp = (double)nut * err;
	GOBACK L235;
    L240:
	klow = iext[j];
	iext[j] = l+1;
	++j;
	++jchnge;
	GOBACK L200;

    L250:
	l = iext[j]+1;
	if (jchnge > 0) GOBACK L215;

    L255:
	if (++l >= kup) goto L260;
	err = (freq_eval(l,nz,grid,x,y,ad)-des[l]) * wt[l];
	if (((double)nut*err-comp) <= 0.0) GOBACK L255;
	comp = (double)nut * err;

	GOBACK L210;
    L260:
	klow = iext[j++];
	GOBACK L200;

    L300:
	if (j > nzz) goto L320;
	if (k1 > iext[1] ) k1 = iext[1];
	if (knz < iext[nz]) knz = iext[nz];
	nut1 = nut;
	nut = -nu;
	l = 0;
	kup = k1;
	comp = ynz*(1.00001);
	luck = 1;
    L310:
	if (++l >= kup) goto L315;
	err = (freq_eval(l,nz,grid,x,y,ad)-des[l]) * wt[l];
	if (((double)nut*err-comp) <= 0.0) GOBACK L310;
	comp = (double) nut * err;
	j = nzz;
	GOBACK L210;

    L315:
	luck = 6;
	goto L325;

    L320:
	if (luck > 9) goto L350;
	if (comp > y1) y1 = comp;
	k1 = iext[nzz];
    L325:
	l = ngrid+1;
	klow = knz;
	nut = -nut1;
	comp = y1*(1.00001);
    L330:
	if (--l <= klow) goto L340;
	err = (freq_eval(l,nz,grid,x,y,ad)-des[l]) * wt[l];
	if (((double)nut*err-comp) <= 0.0) GOBACK L330;
	j = nzz;
	comp = (double) nut * err;
	luck = luck + 10;
	GOBACK L235;
    L340:
	if (luck == 6) goto L370;
	DOloop(j,1,nfcns) {
	    iext[nzz-j] = iext[nz-j];
	}
	iext[1] = k1;
	GOBACK L100;
    L350:
	kn = iext[nzz];
	DOloop(j,1,nfcns) iext[j] = iext[j+1];
	iext[nz] = kn;

	GOBACK L100;
    L370:
	;
    } while (jchnge > 0);

/*
 *    CALCULATION OF THE COEFFICIENTS OF THE BEST APPROXIMATION
 *    USING THE INVERSE DISCRETE FOURIER TRANSFORM
 */
    nm1 = nfcns - 1;
    fsh = 1.0e-06;
    gtemp = grid[1];
    x[nzz] = -2.0;
    cn  = 2*nfcns - 1;
    delf = 1.0/cn;
    l = 1;
    kkk = 0;

    if (edge[1] == 0.0 && edge[2*nbands] == 0.5) kkk = 1;

    if (nfcns <= 3) kkk = 1;
    if (kkk !=     1) {
	dtemp = cos(TWOPI*grid[1]);
	dnum  = cos(TWOPI*grid[ngrid]);
	aa    = 2.0/(dtemp-dnum);
	bb    = -(dtemp+dnum)/(dtemp-dnum);
    }

    DOloop(j,1,nfcns) {
	ft = (j - 1) * delf;
	xt = cos(TWOPI*ft);
	if (kkk != 1) {
	    xt = (xt-bb)/aa;
#if 0
	    /*XX* ckeck up !! */
	    xt1 = sqrt(1.0-xt*xt);
	    ft = atan2(xt1,xt)/TWOPI;
#else
	    ft = acos(xt)/TWOPI;
#endif
	}
L410:
	xe = x[l];
	if (xt > xe) goto L420;
	if ((xe-xt) < fsh) goto L415;
	++l;
	GOBACK L410;
L415:
	a[j] = y[l];
	goto L425;
L420:
	if ((xt-xe) < fsh) GOBACK L415;
	grid[1] = ft;
	a[j] = freq_eval(1,nz,grid,x,y,ad);
L425:
	if (l > 1) l = l-1;
    }

    grid[1] = gtemp;
    dden = TWOPI / cn;
    DOloop (j,1,nfcns) {
	dtemp = 0.0;
	dnum = (j-1) * dden;
	if (nm1 >= 1) {
	    DOloop(k,1,nm1) {
		dtemp += a[k+1] * cos(dnum*k);
	    }
	}
	alpha[j] = 2.0 * dtemp + a[1];
    }

    DOloop(j,2,nfcns) alpha[j] *= 2.0 / cn;
    alpha[1] /= cn;

    if (kkk != 1) {
	p[1] = 2.0*alpha[nfcns]*bb+alpha[nm1];
	p[2] = 2.0*aa*alpha[nfcns];
	q[1] = alpha[nfcns-2]-alpha[nfcns];
	DOloop(j,2,nm1) {
	    if (j >= nm1) {
		aa *= 0.5;
		bb *= 0.5;
	    }
	    p[j+1] = 0.0;
	    DOloop(k,1,j) {
		a[k] = p[k];
		p[k] = 2.0 * bb * a[k];
	    }
	    p[2] += a[1] * 2.0 *aa;
	    jm1 = j - 1;
	    DOloop(k,1,jm1) p[k] += q[k] + aa * a[k+1];
	    jp1 = j + 1;
	    DOloop(k,3,jp1) p[k] += aa * a[k-1];

	    if (j != nm1) {
		DOloop(k,1,j) q[k] = -a[k];
		q[1] += alpha[nfcns - 1 - j];
	    }
	}
	DOloop(j,1,nfcns) alpha[j] = p[j];
    }

    if (nfcns <= 3) {
	  alpha[nfcns+1] = alpha[nfcns+2] = 0.0;
    }
    return 0;
}


/*
 *-----------------------------------------------------------------------
 * FUNCTION: eff
 *  FUNCTION TO CALCULATE THE DESIRED MAGNITUDE RESPONSE
 *  AS A FUNCTION OF FREQUENCY.
 *  AN ARBITRARY FUNCTION OF FREQUENCY CAN BE
 *  APPROXIMATED IF THE USER REPLACES THIS FUNCTION
 *  WITH THE APPROPRIATE CODE TO EVALUATE THE IDEAL
 *  MAGNITUDE.  NOTE THAT THE PARAMETER FREQ IS THE
 *  VALUE OF NORMALIZED FREQUENCY NEEDED FOR EVALUATION.
 *-----------------------------------------------------------------------
 */
static double eff(double freq, double *fx, int lband, int jtype)
{
      if (jtype != 2) return fx[lband];
      else            return fx[lband] * freq;
}

/*
 *-----------------------------------------------------------------------
 * FUNCTION: wate
 *  FUNCTION TO CALCULATE THE WEIGHT FUNCTION AS A FUNCTION
 *  OF FREQUENCY.  SIMILAR TO THE FUNCTION eff, THIS FUNCTION CAN
 *  BE REPLACED BY A USER-WRITTEN ROUTINE TO CALCULATE ANY
 *  DESIRED WEIGHTING FUNCTION.
 *-----------------------------------------------------------------------
 */
static double wate(double freq, double *fx, double *wtx, int lband, int jtype)
{
      if (jtype != 2)          return wtx[lband];
      if (fx[lband] >= 0.0001) return wtx[lband] / freq;
      return                          wtx[lband];
}

/*********************************************************/

/*  This routine accepts basic input information and puts it in 
 *  the form expected by remez.

 *  Adpated from main() by Travis Oliphant
 */

int pre_remez(double *h2, int numtaps, int numbands, double *bands,
                     double *response, double *weight, int type, int maxiter,
                     int grid_density) {
  
  int jtype, nbands, nfilt, lgrid, nz;
  int neg, nodd, nm1;
  int j, k, l, lband, dimsize;
  double delf, change, fup, temp;
  double *tempstor, *edge, *h, *fx, *wtx;
  double *des, *grid, *wt, *alpha, *work;
  double dev;
  int ngrid;
  int *iext;
  int nfcns, wrksize, total_dsize, total_isize;

  lgrid = grid_density;
  dimsize = (int) ceil(numtaps/2.0 + 2);
  wrksize = grid_density * dimsize;
  nfilt = numtaps;
  jtype = type; nbands = numbands;
  printf("  jtype=%02d\n", jtype);

  /* Note:  code assumes these arrays start at 1 */
  edge = bands-1; 
  h = h2 - 1;
  fx = response - 1;
  wtx = weight - 1;

  total_dsize = (dimsize+1)*7 + 3*(wrksize+1);
  total_isize = (dimsize+1);
  /* Need space for:  (all arrays ignore the first element).

     des  (wrksize+1)
     grid (wrksize+1)
     wt   (wrksize+1)
     iext (dimsize+1)   (integer)
     alpha (dimsize+1)
     work  (dimsize+1)*6 

  */
  tempstor = malloc((total_dsize)*sizeof(double)+(total_isize)*sizeof(int));
  if (tempstor == NULL) return -2;

  des = tempstor; grid = des + wrksize+1;
  wt = grid + wrksize+1; alpha = wt + wrksize+1;
  work = alpha + dimsize+1; iext = (int *)(work + (dimsize+1)*6);

  /* Set up problem on dense_grid */

  neg = 1;
  if (jtype == 1) neg = 0;
  nodd = nfilt % 2;
  nfcns = nfilt / 2;
  if (nodd == 1 && neg == 0) nfcns = nfcns + 1;

    /*
     * SET UP THE DENSE GRID. THE NUMBER OF POINTS IN THE GRID
     * IS (FILTER LENGTH + 1)*GRID DENSITY/2
     */
    grid[1] = edge[1];
    delf = lgrid * nfcns;
    delf = 0.5 / delf;
    if (neg != 0) {
	if (edge[1] < delf) grid[1] = delf;
    }
    j = 1;
    l = 1;
    lband = 1;

    /*
     * CALCULATE THE DESIRED MAGNITUDE RESPONSE AND THE WEIGHT
     * FUNCTION ON THE GRID
     */
    for (;;) {
	fup = edge[l + 1];
	do {
	    temp = grid[j];
	    des[j] = eff(temp,fx,lband,jtype);
	    wt[j] = wate(temp,fx,wtx,lband,jtype);
	    if (++j > wrksize) {
                /* too many points, or too dense grid */
                free(tempstor);
                return -1;
            } 
	    grid[j] = temp + delf;
	} while (grid[j] <= fup);

	grid[j-1] = fup;
	des[j-1] = eff(fup,fx,lband,jtype);
	wt[j-1] = wate(fup,fx,wtx,lband,jtype);
	++lband;
	l += 2;
	if (lband > nbands) break;
	grid[j] = edge[l];
    }

    ngrid = j - 1;
    if (neg == nodd) {
	if (grid[ngrid] > (0.5-delf)) --ngrid;
    }

    // PRINT GRID
    // DOloop(j,1,ngrid) {
    //    printf("  j=%02d: grid[j]=%g\n", j, grid[j]);
    //}

    /*
     * SET UP A NEW APPROXIMATION PROBLEM WHICH IS EQUIVALENT
     * TO THE ORIGINAL PROBLEM
     */
    if (neg <= 0) {
	if (nodd != 1) {
	    DOloop(j,1,ngrid) {
		change = cos(PI*grid[j]);
		des[j] = des[j] / change;
		wt[j]  = wt[j] * change;
	    }
	}
    } else {
	if (nodd != 1) {
	    DOloop(j,1,ngrid) {
		change = sin(PI*grid[j]);
		des[j] = des[j] / change;
		wt[j]  = wt[j]  * change;
	    }
	} else {
	    DOloop(j,1,ngrid) {
		change = sin(TWOPI * grid[j]);
		des[j] = des[j] / change;
		wt[j]  = wt[j]  * change;
	    }
	}
    }

    /*XX*/
    temp = (double)(ngrid-1) / (double)nfcns;
    DOloop(j,1,nfcns) {
	iext[j] = (int)((j-1)*temp) + 1; /* round? !! */
    }
    iext[nfcns+1] = ngrid;
    nm1 = nfcns - 1;
    nz  = nfcns + 1;

  //printf("pre_remez FOO!\n");
  //printf("ngrid = %d\n", ngrid);
  //printf("numbands = %d\n", numbands);
    if (remez(&dev, des, grid, edge, wt, ngrid, numbands, iext, alpha, nfcns,
              maxiter, work, dimsize) < 0) {
        free(tempstor);
        return -1;
    }

    /*
     * CALCULATE THE IMPULSE RESPONSE.
     */
    if (neg <= 0) {

	if (nodd != 0) {
	    DOloop(j,1,nm1) {
		h[j] = 0.5 * alpha[nz-j];
	    }
	    h[nfcns] = alpha[1];
	} else {
	    h[1] = 0.25 * alpha[nfcns];
	    DOloop(j,2,nm1) {
		h[j] = 0.25 * (alpha[nz-j] + alpha[nfcns+2-j]);
	    }
	    h[nfcns] = 0.5*alpha[1] + 0.25*alpha[2];
	}
    } else {
	if (nodd != 0) {
	    h[1] = 0.25 * alpha[nfcns];
	    h[2] = 0.25 * alpha[nm1];
	    DOloop(j,3,nm1) {
		h[j] = 0.25 * (alpha[nz-j] - alpha[nfcns+3-j]);
	    }
	    h[nfcns] = 0.5 * alpha[1] - 0.25 * alpha[3];
	    h[nz] = 0.0;
	} else {
	    h[1] = 0.25 * alpha[nfcns];
	    DOloop(j,2,nm1) {
		h[j] = 0.25 * (alpha[nz-j] - alpha[nfcns+2-j]);
	    }
	    h[nfcns] = 0.5 * alpha[1] - 0.25 * alpha[2];
	}
    }

    DOloop(j,1,nfcns){
        k = nfilt + 1 - j;
        if (neg == 0)
           h[k] = h[j];
        else
           h[k] = -h[j];
    }
    if (neg == 1 && nodd == 1) h[nz] = 0.0;

  free(tempstor);
  return 0;

}

/**************************************************************
 * End of remez routines 
 **************************************************************/


/****************************************************/
/* End of python-independent routines               */
/****************************************************/

/*
gcc -fPIC -c tkfork.c
gcc -shared  tkfork.o -o tkfork.so
LD_LIBRARY_PATH=. julia
x=[1.,2.,3.,4.]
ccall((:simple_test, "tkfork"), Void, (Cint, Ptr{Cfloat}), length(x), x)
*/
void simple_test(int xlen, double *x)
{
    printf("FOO!\n");
    for (int i = 0; i < xlen; i++) {
        printf("x[%d]=%f\n", i, x[i]);
    }
}


