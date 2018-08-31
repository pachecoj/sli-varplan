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
 * Desc.:  Header class GaussQuad
 * ------------------------------------------------------------------- */

#ifndef GAUSSQUAD_H
#define GAUSSQUAD_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/quad/default.h"
#include "lhotse/matrix/ArrayUtils.h"
#include "lhotse/specfun/Specfun.h"

//BEGINNS(quad)
  /**
   * Static methods to compute weights, abscissas for Gaussian quadrature
   * rules.
   * The quadrature rule is
   *   \int_I W(x) f(x) dx \approx \sum_{i=1}^n w_i f(x_i).
   * Here, x_i are abscissas, w_i weights. Note that w_i > 0.
   * The methods here return x_1,...,x_n and w_1,...,w_n.
   * <p>
   * NOTE: Suitable only for smooth, "polynomial-like" f(x).
   * We partly re-implement the NETLIB procedure GAUSSQ from the Fortran
   * code. The subroutine GAUSQ2 is just taken from there.
   * The original Fortran code cannot handle large alpha values in the
   * Gauss-Laguerra rule, because it does not use log weights.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class GaussQuad
  {
  public:
    // Public static methods

    /**
     * Weights, abscissas for Gauss-Hermite quadrature, weight function
     *   W(x) = (1/\sqrt{\pi}) e^{-x^2}, I = (-infty,+infty).
     * <p>
     * NOTE: In order to do N(0,1) expectations over g(t), use the function
     *   f(x) = g(\sqrt{2} x) / \sqrt{\pi}
     * Needs scratch buffer of size n+2, can be given via 'sbuff', is
     * otherwise alloc. here.
     *
     * @param n        Number of points
     * @param wgths    Weights w_1,...,w_n
     * @param abscs    Abscissas x_1,...,x_n
     * @param sbuff    Scratch buffer of size n. Def.: 0
     */

    static void gaussHermite(int n,double* wgths,double* abscs,
			     double* sbuff=0);

    /**
     * Log weights, abscissas for Gauss-Laguerre quadrature, weight function
     *   W(x) = x^alpha e^{-x}, I = (0,+infty), alpha > -1.
     * <p>
     * Needs scratch buffer of size n+2, can be given via 'sbuff', is
     * otherwise alloc. here.
     *
     * @param n        Number of points
     * @param alpha    S.a. Must be > -1
     * @param lwgths   Log weights log(w_1),...,log(w_n)
     * @param abscs    Abscissas x_1,...,x_n
     * @param sbuff    Scratch buffer of size n. Def.: 0
     */
    static void gaussLaguerre(int n,double alpha,double* wgths,double* abscs,
			      double* sbuff=0);

    /**
     * DOWNW. COMPAT.! DO NOT USE
     *
     * Weights, abscissas for Gauss-Hermite quadrature, weight function
     *   W(x) = (1/\sqrt{\pi}) e^{-x^2}.
     * If m=floor((n+1)/2), only the first m of w_j, x_j are returned (see
     * header comment). All weights are positive.
     * NOTE: In order to do N(0,1) expectations over g(t), use the function
     *   f(x) = g(\sqrt{2} x).
     * If 'complete'==true, we return the full set of n weights/abscissas
     * in 'wgths', 'abscs' even though they are redundant.
     *
     * @param n        Number of points
     * @param wgths    Weights w_1,...,w_m
     * @param abscs    Abscissas x_1,...,x_m
     * @param complete S.a. Def.: false
     */
    static void gaussHermiteOld(int n,double* wgths,double* abscs,
				bool complete=false);

  protected:
    // Internal static methods

    /**
     * Re-implementation of NETLIB GAUSSQ.
     * The original routine supports 6 different rules, selected by
     * 'kind'. In the moment, we just do two of them:
     * - 4: Hermite quadrature, w(x)=exp(-x*x) on (-inf,+inf)
     * - 6: Generalized Laguerre quadrature,
     *      w(x)=exp(-x)*x**alpha on (0,+inf), alpha>-1
     * We also do not bother with the kpts, endpts arguments.
     * NOTE: We return log weights, not weights. For large 'alpha',
     * weights overflow immediately in GAUSSQ!
     *
     * @param kind  Kind of rule
     * @param n     Number of nodes
     * @param alpha Used only for Laguerre ('kind'==6)
     * @param absc  Abscissas written here (arg. t of gaussq)
     * @param lwgts Log weights written here (arg. log(w) of gaussq)
     * @param b     Scratch array of size n
     */
    static void mygaussq(int kind,int n,double alpha,double* absc,
			 double* lwgts,double* b);

    /**
     * Differs from NETLIB, in that 'lmuzero' returns log(muzero).
     */
    static void myclass(int kind,int n,double alpha,double* b,double* a,
			double& lmuzero);
  };

  // Inline methods

  inline void GaussQuad::myclass(int kind,int n,double alpha,double* b,
				 double* a,double& lmuzero)
  {
    int i;
    double temp;

    switch (kind) {
    case 4:
      // Hermite
      lmuzero=0.5*M_LNPI;
      ArrayUtils<double>::fill(a,0.0,n);
      for (i=0; i<n-1; i++) b[i]=sqrt(0.5*((double) (i+1)));
      break;
    case 6:
      // Generalized Laguerre
      lmuzero=Specfun::logGamma(alpha+1.0);
      for (i=0; i<n-1; i++) {
	temp=(double) (i+1);
	a[i]=alpha+2.0*temp-1.0;
	b[i]=sqrt(temp*(temp+alpha));
      }
      a[n-1]=alpha+2.0*((double) n)-1.0;
      break;
    }
  }
//ENDNS

#endif
