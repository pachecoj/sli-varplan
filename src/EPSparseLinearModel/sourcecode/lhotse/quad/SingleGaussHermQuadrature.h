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
 * Desc.:  Header abstract class SingleGaussHermQuadrature
 * ------------------------------------------------------------------- */

#ifndef SINGLEGAUSSHERMQUADRATURE_H
#define SINGLEGAUSSHERMQUADRATURE_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/quad/default.h"
#include "lhotse/quad/SingleQuadrature.h"
#include "lhotse/quad/GaussQuad.h"

//BEGINNS(quad)
  /**
   * Abstract subclass of 'SingleQuadrature'. Here, 'compLogT' has to be
   * supplied to compute log t(u) for given u. This is used as subroutine
   * to compute the service methods using Gauss-Hermite quadrature (see
   * 'GaussQuad' in module 'quad').
   * <p>
   * ATTENTION: This is not really a black box. Gauss-Hermite quadrature
   * requires that t(u), log t(u) can be well approximated by a low-order
   * polynomial over the range where the Gaussian integrated against has
   * significant mass.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class SingleGaussHermQuadrature : public virtual SingleQuadrature
  {
  protected:
    // Constants, Types

    static const int numQuadPoints=15; // Number of abscissas for quadrature

    // Members

    /*
     * 'quadAbs' stores the abscissas for quadrature (prov. by 'GaussQuad',
     * mult. by sqrt{2}), while 'quadWts' stores log w_j. We do not use
     * the symmetry in abscissas and weights.
     */
    StVector quadAbs,quadWts,equadWts;
    mutable StVector buffVec;

  public:
    // Public methods

    /**
     * Constructor. Initializes the parameters for the quadrature routine
     */
    SingleGaussHermQuadrature() : SingleQuadrature() {
      quadAbs.zeros(numQuadPoints);
      quadWts.zeros(numQuadPoints);
      equadWts.zeros(numQuadPoints);
      buffVec.zeros(numQuadPoints);
      GaussQuad::gaussHermiteOld(numQuadPoints,equadWts.getFlatBuff(),
				 quadAbs.getFlatBuff(),true);
      quadWts.apply1(equadWts,ptr_fun(log));
      quadAbs.prod(sqrt(2.0));
    }

    double compLogPart(double mean,double var,double* alpha,double* nu) const;

    double compExpectLog(double mean,double var,double* lmean=0,double* lvar=0)
      const;

  protected:
    // Internal methods

    /**
     * Returns log t(u), u=='u'.
     * Has to be implemented by subclasses.
     *
     * @param u Arg.
     * @return  log t(u)
     */
    virtual double compLogT(double u) const = 0;
  };

  // Inline methods

  inline double
  SingleGaussHermQuadrature::compLogPart(double mean,double var,double* alpha,
					 double* nu) const
  {
    int j;
    double sqa=sqrt(var),logz;

    for (j=0; j<numQuadPoints; j++)
      buffVec[j]=compLogT(sqa*quadAbs[j]+mean)+quadWts[j];
    logz=buffVec.stableLogSum();
    if (alpha!=0) {
      if (nu==0) throw InvalidParameterException(EXCEPT_MSG("alpha,nu"));
      double asum=0.0,nsum=0.0,temp,arg;
      for (j=0; j<numQuadPoints; j++) {
	asum+=(temp=exp(buffVec[j]-logz)*(arg=quadAbs[j]));
	nsum+=temp*arg;
      }
      *alpha=temp=asum/sqa;
      *nu=temp*temp+(1.0-nsum)/var;
    }

    return logz;
  }

  inline double
  SingleGaussHermQuadrature::compExpectLog(double mean,double var,
					   double* lmean,double* lvar) const
  {
    int j;
    double sqa=sqrt(var),sum;

    if (lmean==0) {
      for (j=0,sum=0.0; j<numQuadPoints; j++)
	sum+=equadWts[j]*compLogT(sqa*quadAbs[j]+mean);
    } else {
      if (lvar==0)
	throw InvalidParameterException(EXCEPT_MSG("lmean,lvar"));
      double hsum,asum,arg,temp;
      for (j=0,sum=hsum=asum=0.0; j<numQuadPoints; j++) {
	sum+=(temp=equadWts[j]*compLogT(sqa*(arg=quadAbs[j])+mean));
	hsum+=(temp*=arg);
	asum+=(temp*arg);
      }
      *lmean=hsum; *lvar=asum;
    }

    return sum;
  }
//ENDNS

#endif
