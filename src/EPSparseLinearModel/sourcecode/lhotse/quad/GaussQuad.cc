/*
 * Copyright (C) 1999-2006 Matthias Seeger
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
 */
/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Library source file
 * Module: quad
 * Desc.:  Definition of class GaussQuad
 * ------------------------------------------------------------------- */

#include "lhotse/quad/GaussQuad.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_machine.h>

//BEGINNS(quad)
  // External Fortran code

#ifdef HAVE_FORTRAN
  extern "C" void
  F77_FUNC(gausq2,GAUSQ2) (const int* n,double* d,double* e,double* z,
			   int* ierr,const double* machep);
#endif

  // Static methods

  void GaussQuad::mygaussq(int kind,int n,double alpha,double* absc,
			   double* lwgts,double* b)
  {
#ifndef HAVE_FORTRAN
    throw NotImplemException(EXCEPT_MSG("NO EXTERNAL FORTRAN CODE"));
#else
    double lmuzero,machep=GSL_DBL_EPSILON;
    int i,ierr;

    myclass(kind,n,alpha,b,absc,lmuzero);
    ArrayUtils<double>::fill(lwgts+1,0.0,n-1); lwgts[0]=1.0;
    F77_FUNC(gausq2,GAUSQ2) (&n,absc,b,lwgts,&ierr,&machep);
    if (ierr!=0) {
      cout << "ERROR(GAUSQ2): Could not det. " << ierr
	   << "-th eigenvalue after 30 iters.!" << endl;
      throw NumericalException(EXCEPT_MSG(""));
    }
    for (i=0; i<n; i++)
      lwgts[i]=2.0*log(fabs(lwgts[i]))+lmuzero;
#endif
  }

  void GaussQuad::gaussHermite(int n,double* wgths,double* abscs,
			       double* sbuff)
  {
    ArrayHandle<double> shand;

    if (n<=0) throw InvalidParameterException(EXCEPT_MSG(""));
    if (sbuff==0) {
      shand.changeRep(n); sbuff=shand.p();
    }
    mygaussq(4,n,0.0,abscs,wgths,sbuff);
    ArrayUtils<double>::applyFunc(wgths,wgths,n,ptr_fun(exp));
  }

  void GaussQuad::gaussLaguerre(int n,double alpha,double* lwgths,
				double* abscs,double* sbuff)
  {
    ArrayHandle<double> shand;

    if (n<=0 || alpha<=-1.0) throw InvalidParameterException(EXCEPT_MSG(""));
    if (sbuff==0) {
      shand.changeRep(n); sbuff=shand.p();
    }
    //cout << "alpha=" << alpha << endl; // DEBUG!
    mygaussq(6,n,alpha,abscs,lwgths,sbuff);
  }

  /*
   * We use good starting guesses for the Hermite poly roots, then
   * use Newton to refine these. Newton is run for at most 10 iters.,
   * stopped early if the step is smaller than 5e-14. If it runs for
   * the full 10 iters., the last iterate is used (no error or
   * warning message).
   * The weights must sum to 1, so the middle weight is determined
   * in that way.
   */
  void GaussQuad::gaussHermiteOld(int n,double* wgths,double* abscs,
				  bool complete)
  {
    int i,j,k,m=n/2;
    double hp0,hpm1,hpm2,dhp0,x,temp,sq,lastsq,sum;

    if (n<=0) throw InvalidParameterException(EXCEPT_MSG(""));
    sum=0.0; // sum of weights
    for (i=0; i<m; i++) {
      // Initial guesses for root
      switch (i) {
      case 0:
	temp=(double) (2*n+1);
	x=sqrt(temp)-1.85575*exp(-log(temp)/6.0);
	break;
      case 1:
	x-=1.14*exp(0.426*log((double) n))/x;
	break;
      case 2:
	x=1.86*x-0.86*abscs[0];
	break;
      case 3:
	x=1.91*x-0.91*abscs[1];
	break;
      default:
	x=2.0*x-abscs[i-2];
      }
      // Use Newton's method to obtain refined asbscissa from starting
      // point. The orthonormal Hermite polys are evaluated using their
      // recursive definition
      for (j=0; j<10; j++) {
	hp0=1.0; hpm1=0.0;
	sq=0.0;
	for (k=1; k<=n; k++) {
	  hpm2=hpm1; // shift
	  hpm1=hp0;
	  lastsq=sq; sq=sqrt((double) k);
	  hp0=(x*sqrt(2)*hpm1-lastsq*hpm2)/sq;
	}
	dhp0=hpm1*sqrt((double) (2*n));
	temp=-hp0/dhp0; // Newton step
	x+=temp;
	if (fabs(temp)<5e-14) break; // accurate enough
      }
      abscs[i]=x; // abscissa
      sum+=(wgths[i]=2.0/(dhp0*dhp0)); // weight
    }
    if (n%2==1) {
      wgths[m]=1.0-2.0*sum; // weights must sum to 1
      abscs[m]=0.0;
    }
    if (complete) {
      for (i=0; i<m; i++) {
      wgths[n-i-1]=wgths[i];
      abscs[n-i-1]=-abscs[i];
      }
    }
  }
//ENDNS
