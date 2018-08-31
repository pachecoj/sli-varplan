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
 * Desc.:  Header abstract class DistrExp
 * ------------------------------------------------------------------- */

#ifndef DISTREXP_H
#define DISTREXP_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/rando/default.h"
#include "lhotse/rando/DistribSamp.h"
#include "lhotse/specfun/Specfun.h"

//BEGINNS(rando)
  /**
   * Univariate exponential distribution
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class DistrExp : public DistribSamp
  {
  protected:
    double mean;
    bool defPar; // Default parameter values (mean 1) active?

  public:
    // Constructors

    DistrExp() : DistribSamp(),mean(1.0),defPar(true) {}

    DistrExp(double m) : DistribSamp(),mean(m) {
      if (m<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
      defPar=(mean==1.0);
    }

    // Public methods

    double getMean() const {
      return mean;
    }

    void setMean(double m) {
      if (m<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
      mean=m;
      defPar=(mean==1.0);
    }

    /**
     * Returns value of pdf at given point
     */
    double pdf(double x) const {
      return defPar ? Specfun::pdfExp(x) :
	Specfun::pdfExp(mean,x);
    }

    /**
     * Returns value of cdf at given point
     */
    double cdf(double x) const {
      return defPar ? Specfun::cdfExp(x) :
	Specfun::cdfExp(mean,x);
    }

    double sample(Generator* gen=0) const {
      if (gen==0) gen=generator();
      return defPar ? (double) Random::devExp(gen) :
	(double) Random::devExp(mean,gen);
    }
  };
//ENDNS

#endif
