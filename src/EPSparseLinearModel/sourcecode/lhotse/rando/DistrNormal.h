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
 * Desc.:  Header abstract class DistrNormal
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

//BEGINNS(rando)
  /**
   * Univariate normal distribution
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class DistrNormal : public DistribSamp
  {
  protected:
    double mean;
    double stddev;
    bool defPar; // Default parameter values (mean 0,stddev 1) active?

  public:
    // Constructors

    DistrNormal() : DistribSamp(),mean(0.0),stddev(1.0),defPar(true) {}

    DistrNormal(double m,double sd) : DistribSamp(),mean(m),stddev(sd) {
      if (sd<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
      defPar=(mean==0.0 && stddev==1.0);
    }

    // Public methods

    double getMean() const {
      return mean;
    }

    double getStddev() const {
      return stddev;
    }

    void setMean(double m) {
      mean=m;
      defPar=(mean==0.0 && stddev==1.0);
    }

    void setStddev(double sd) {
      if (sd<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
      stddev=sd;
      defPar=(stddev==1.0 && mean==0.0);
    }

    /**
     * Returns value of pdf at given point
     */
    double pdf(double x) const {
      return defPar?Specfun::pdfNormal(x):
	Specfun::pdfNormal(mean,stddev,x);
    }

    /**
     * Returns value of cdf at given point
     */
    double cdf(double x) const {
      return defPar?Specfun::cdfNormal(x):
	Specfun::cdfNormal(mean,stddev,x);
    }

    double sample(Generator* gen=0) const {
      if (gen==0) gen=generator();
      return defPar?Random::devNormalRatUnif(gen):
	Random::devNormalRatUnif(mean,stddev,gen);
    }
  };
//ENDNS

#endif
