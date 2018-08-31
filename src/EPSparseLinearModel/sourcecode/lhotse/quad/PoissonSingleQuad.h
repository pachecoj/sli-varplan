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
 * Desc.:  Header abstract class PoissonSingleQuad
 * ------------------------------------------------------------------- */

#ifndef POISSONSINGLEQUAD_H
#define POISSONSINGLEQUAD_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/quad/default.h"
#include "lhotse/quad/SingleGaussHermQuadrature.h"

//BEGINNS(quad)
  /**
   * Quadrature tools for Poisson (exponential) likelihood:
   *   t(u) = exp(-tau exp(u+off))            ['state'==false]
   *   t(u) = exp(u+off) exp(-tau exp(u+off)) ['state'==true]
   * Here, tau=='tau' is positive, off=='off'.
   * We cannot compute 'compLogPart' analytically, therefore inherit from
   * 'SingleGaussHermQuadrature'.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class PoissonSingleQuad : public SingleGaussHermQuadrature
  {
  protected:
    // Members

    bool state; // See header comment
    double tau; // "
    double off; // "

  public:
    // Public methods

    /**
     * Default constructor.
     */
    PoissonSingleQuad() : SingleGaussHermQuadrature(),state(false),tau(1.0),
    off(0.0) {}

    /**
     * @param ptau   Value for 'tau' (>0)
     * @param pstate Value for 'state'
     * @param poff   Value for 'off''
     */
    void setPars(double ptau,bool pstate,double poff) {
      if (ptau<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
      state=pstate; tau=ptau; off=poff;
    }

    double compExpectLog(double mean,double var,double* lmean=0,double* lvar=0)
      const;

  protected:
    // Internal methods

    double compLogT(double u) const;
  };

  // Inline methods

  inline double
  PoissonSingleQuad::compExpectLog(double mean,double var,double* lmean,
				   double* lvar) const
  {
    double htil=mean+off;
    double rho=state?1.0:0.0,temp=tau*exp(0.5*var+htil);

    if (lmean!=0) {
      *lmean=sqrt(var)*(rho-temp);
      *lvar=rho*htil-temp*(1.0+var);
    }

    return rho*htil-temp;
  }

  inline double PoissonSingleQuad::compLogT(double u) const
  {
    double temp1=u+off;
    double temp2=-tau*exp(temp1);

    return state?(temp1+temp2):temp2;
  }
//ENDNS

#endif
