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
 * Desc.:  Header abstract class DistrPoisson
 * ------------------------------------------------------------------- */

#ifndef DISTRNORMAL_H
#define DISTRNORMAL_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/rando/default.h"
#include "lhotse/rando/DistribSamp.h"
#include "lhotse/specfun/Specfun.h"
#include "lhotse/rando/Random.h"

////BEGINNS(rando)
  /**
   * Univariate Poisson distribution
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class DistrPoisson : public DistribSamp
  {
  protected:
    double rate; // rate parameter (mean)

  public:
    // Constructors

    DistrPoisson() : DistribSamp(),rate(1.0) {}

    DistrPoisson(double r) : DistribSamp(),rate(r) {
      if (r<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
    }

    // Public methods

    double getRate() const {
      return rate;
    }

    void setRate(double r) {
      if (r<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
      rate=r;
    }

    /**
     * Returns value of pdf at given point
     */
    double pdf(double x) const {
      double temp=floor(x);
      return (temp==x)?Specfun::pdfPoisson(rate,(int) x):0.0;
    }

    /**
     * Returns value of cdf at given point
     */
    double cdf(double x) const {
      return Specfun::cdfPoisson(rate,(int) floor(x));
    }

    double sample(Generator* gen=generator()) const {
      return (double) Random::devPoisson(rate,gen);
    }
  };
////ENDNS

#endif
