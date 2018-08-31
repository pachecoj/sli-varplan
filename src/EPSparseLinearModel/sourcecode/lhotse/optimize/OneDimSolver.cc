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
 * Module: optimize
 * Desc.:  Definition of class OneDimSolver
 * ------------------------------------------------------------------- */

#include "lhotse/optimize/OneDimSolver.h"

//BEGINNS(optimize)

  const int OneDimSolver::brackRightRegular ;
  const int OneDimSolver::brackRightBound   ;
  const int OneDimSolver::brackRightInfinite;

  // Public static methods

  /*
   * Losely related to NR routine 'rtsafe', which is however not
   * correct. We correct the mistakes and also allow for a missing right
   * bracket end ('rightInf'==true).
   */
#define MAXIT 100 // for 'newton'
  double OneDimSolver::newton(FuncOneDim* func,double l,double r,double acc,
			      int brRight,double boundR)
  {
    int j;
    double fh,fl,df,f;
    double dx,temp,xh,xl,rts,olds,alpha;
    bool nextBisect,didNewton,numerErr;

    if (!func->hasDerivative())
      throw InvalidParameterException("'func' must return derivatives!");
    func->eval(l,&fl,&df);
    if (fl==0.0) return l;
    if (brRight==brackRightRegular) {
      func->eval(r,&fh,&df);
      if ((fl>0.0 && fh>0.0) || (fl<0.0 && fh<0.0)) {
	cout << "l=" << l << ": fl=" << fl << endl
	     << "r=" << r << ": fr=" << fh << endl;
	throw InvalidParameterException("Root must be bracketed in [l,r]");
      }
    } else {
      /* Use method described in header comment in order to find right bracket
	 end. */
      bool isBound=(brRight==brackRightBound);
      j=0;
      dx=acc; // start with step of size 'acc'
      if (isBound && l+dx>=boundR-acc)
	throw InvalidParameterException("Cannot do any step left of 'boundR'");
      for (;;) {
	if (j++>MAXIT)
	  throw NumericalException("Maximum number of iterations exceeded in OneDimSolver::newton");
	do {
	  rts=l+dx;
	  numerErr=false;
	  try {
	    func->eval(rts,&f,&df);
	  } catch (...) {
	    if (dx==acc)
	      throw NumericalException("Newton failed: Cannot find right bracket end!");
	    dx=acc; numerErr=true; // try again with step size 'acc'
	  }
	} while (numerErr);
	if (f*fl<0.0) break; // OK, found right bracket end
	alpha=(fl-f)/(l-rts)-df;
	if (alpha*f<0.0) {
	  // OK, use quadratic approx.
	  alpha/=(l-rts);
	  dx=(sqrt(df*df-4.0*alpha*f)-df)/(2.0*alpha);
	} else if (fabs(df)>1e-8) {
	  // Use linear approx. (Newton step)
	  dx=-f/df;
	} else
	  dx=acc;
	if (isBound && rts+dx>boundR-acc)
	  dx=boundR-acc-rts; // bold!
	if (dx<acc) {
	  dx=acc; // minimum step size
	  if (isBound && rts+dx>boundR-acc)
	    throw NumericalException("Newton failed: Cannot find right bracket end!");
	}
	l=rts; // new left bracket end
      }
      r=rts; // right bracket end detected
      fh=f;
    }
    if (fh==0.0) return r;
    if (fabs(l-r)<acc) return l;
    if (fl<0.0) {
      xl=l; xh=r;
    } else {
      xh=l; xl=r;
    }
    rts=0.5*(l+r); // start with bisection step
    func->eval(rts,&f,&df);
    if (f<0.0) xl=rts;
    else xh=rts;
    if (fabs(xl-xh)<acc) return rts;
    nextBisect=false; // do Newton step next
    for (j=1; j<=MAXIT; j++) {
      /* We choose a bisection step if either 'nextBisect' is true or the
	 Newton step falls out of the bracket.
	 If we took the Newton step and find out that by doing so, the
	 bracket shrunk by a fraction less than 0.1, we set 'nextBisect'
	 to true, which leads to bisection in the next turn. We can therefore
	 guarantee that the bracket will shrink within two iterations by
	 at least a 9/20 fraction. */
      temp=fabs(df*(xl-xh));
      if (nextBisect || (((xh-rts)*df+f)*((xl-rts)*df+f))>=0.0 ||
	  fabs(df)<1e-15) {
	// Bisection step
	dx=0.5*(xh-xl);
	rts=xl+dx;
	didNewton=false;
      } else {
	// Newton step
	dx=f/df;
	rts-=dx;
	didNewton=true;
      }
      func->eval(rts,&f,&df);
      olds=fabs(xl-xh); // old bracket size
      if (f<0.0) xl=rts;
      else xh=rts;
      temp=fabs(xl-xh); // new bracket size
      if (temp<acc) return rts;
      nextBisect=didNewton && (temp>0.9*olds); // Bisection step next?
    }
    throw NumericalException("Maximum number of iterations exceeded in OneDimSolver::newton");
  }
#undef MAXIT
//ENDNS
