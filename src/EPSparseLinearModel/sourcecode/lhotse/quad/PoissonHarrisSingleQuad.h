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
 * Desc.:  Header abstract class PoissonHarrisSingleQuad
 * ------------------------------------------------------------------- */

#ifndef POISSONHARRISSINGLEQUAD_H
#define POISSONHARRISSINGLEQUAD_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/quad/default.h"
#include "lhotse/quad/SingleGaussHermQuadrature.h"

//BEGINNS(quad)
  /**
   * Same as 'PoissonSingleQuad', but the rate is obtained as f(u+off),
   * where f(.) is a function suggested by Harris et. al., Nature 424
   * (2003). Namely, f(x) = exp(x) for x<0, f(x) = 1+x for x>=0.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class PoissonHarrisSingleQuad : public SingleGaussHermQuadrature
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
    PoissonHarrisSingleQuad() : SingleGaussHermQuadrature(),state(false),
				tau(1.0),off(0.0) {}

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
      const {
      throw NotImplemException(EXCEPT_MSG(""));
    }

  protected:
    // Internal methods

    double compLogT(double u) const;
  };

  // Inline methods

  inline double PoissonHarrisSingleQuad::compLogT(double u) const
  {
    double temp1=u+off,ret;

    if (temp1<0.0) {
      ret=-tau*exp(temp1);
      if (state) ret+=temp1;
    } else {
      ret=-tau*(temp1+1.0);
      if (state) ret+=log1p(temp1);
    }

    return ret;
  }
//ENDNS

#endif
