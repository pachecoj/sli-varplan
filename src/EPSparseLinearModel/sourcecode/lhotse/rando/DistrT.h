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
 * Desc.:  Header abstract class DistrT
 * ------------------------------------------------------------------- */

#ifndef DISTRT_H
#define DISTRT_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/rando/default.h"
#include "lhotse/rando/DistribSamp.h"
#include "lhotse/specfun/Specfun.h"
#include "lhotse/rando/Random.h"

//BEGINNS(rando)
  /**
   * Univariate t distribution
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class DistrT : public DistribSamp
  {
  protected:
    double degf; // degrees of freedom (>0)
    double mean;
    double stddev;
    bool defPar; // Default parameter values (mean 0,stddev 1) active?

  public:
    // Constructors

    DistrT() : DistribSamp(),degf(1.0),mean(0.0),stddev(1.0),defPar(true) {}

    DistrT(double df) : DistribSamp(),degf(df),mean(0.0),stddev(1.0),
      defPar(true) {
      if (df<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
    }

    DistrT(double df,double m,double sd) : DistribSamp(),degf(df),mean(m),
      stddev(sd) {
      if (df<=0.0 || sd<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
      defPar=(mean==0.0 && stddev==1.0);
    }

    // Public methods

    double getDegFreedom() const {
      return degf;
    }

    double getMean() const {
      return mean;
    }

    double getStddev() const {
      return stddev;
    }

    void setDegFreedom(double df) {
      if (df<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
      degf=df;
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
      return defPar?Specfun::pdfT(degf,x):
	Specfun::pdfT(degf,(x-mean)/stddev)/stddev;
    }

    /**
     * Returns value of cdf at given point
     */
    double cdf(double x) const {
      return defPar?Specfun::cdfT(degf,x):
	Specfun::cdfT(degf,(x-mean)/stddev);
    }

    double sample(Generator* gen=0) const {
      if (gen==0) gen=generator();
      return defPar?Random::devT(degf,gen):
	Random::devT(degf,mean,stddev,gen);
    }
  };
//ENDNS

#endif
