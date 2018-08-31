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
 * Module: rando
 * Desc.:  Header abstract class DistrGamma
 * ------------------------------------------------------------------- */

#ifndef DISTRGAMMA_H
#define DISTRGAMMA_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/rando/default.h"
#include "lhotse/rando/DistribSamp.h"
#include "lhotse/specfun/Specfun.h"

//BEGINNS(rando)
  /**
   * Univariate gamma distribution
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class DistrGamma : public DistribSamp
  {
  protected:
    double a;      // shape parameter
    double b;      // scale parameter
    double abterm; // temp. value dep. on a,b (see pdf)

  public:
    // Constructors

    DistrGamma(double pa) : DistribSamp(),a(pa),b(1.0) {
      if (pa<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
    }

    DistrGamma(double av,double bv) : DistribSamp(),a(av),b(bv) {
      if (a<=0.0 || b<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
    }

    // Public methods

    double getShapePar() const {
      return a;
    }

    double getScalePar() const {
      return b;
    }

    void setShapePar(double av) {
      if (av<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
      a=av;
    }

    void setScalePar(double bv) {
      if (bv<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
      b=bv;
    }

    /**
     * Returns value of pdf at given point
     */
    double pdf(double x) const {
      return (x>0.0) ? exp((a-1.0)*log(x) - abterm - x/b) : 0.0;
    }

    /**
     * Returns value of cdf at given point
     */
    double cdf(double x) const {
      return Specfun::cdfGamma(a,b,x);
    }

    double sample(Generator* gen=0) const {
      if (gen==0) gen=generator();
      return (double) Random::devGamma(a,b,gen);
    }

  protected:
    // Internal methods

    void calcTemp() {
      abterm=Specfun::logGamma(a)+a*log(b);
    }
  };
//ENDNS

#endif
