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
 * Desc.:  Header abstract class DistrBeta
 * ------------------------------------------------------------------- */

#ifndef DISTRBETA_H
#define DISTRBETA_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/rando/default.h"
#include "lhotse/rando/DistribSamp.h"
#include "lhotse/specfun/Specfun.h"
#include "lhotse/rando/Random.h"

//BEGINNS(rando)
  /**
   * Univariate Beta distribution
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class DistrBeta : public DistribSamp
  {
  protected:
    double a,b; // shape parameters

  public:
    // Constructors

    DistrBeta(double pa,double pb) : DistribSamp(),a(pa),b(pb) {
      if (a<=0.0 || b<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
    }

    // Public methods

    double getA() const {
      return a;
    }

    double getB() const {
      return b;
    }

    void setA(double pa) {
      if (pa<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
      a=pa;
    }

    void setB(double pb) {
      if (pb<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
      b=pb;
    }

    /**
     * Returns value of pdf at given point
     */
    double pdf(double x) const {
      return Specfun::pdfBeta(a,b,x);
    }

    /**
     * Returns value of cdf at given point
     */
    double cdf(double x) const {
      return Specfun::cdfBeta(a,b,x);
    }

    double sample(Generator* gen=0) const {
      if (gen==0) gen=generator();
      return (double) Random::devBeta(a,b,gen);
    }
  };
//ENDNS

#endif
