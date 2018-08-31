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
 * Desc.:  Header abstract class DistrUnif
 * ------------------------------------------------------------------- */

#ifndef DISTRUNIF_H
#define DISTRUNIF_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/rando/default.h"
#include "lhotse/rando/DistribSamp.h"

//BEGINNS(rando)
  /**
   * Univariate uniform distribution on (a,b)
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class DistrUnif : public DistribSamp
  {
  protected:
    double a,b; // interval limits (a,b)
    bool defPar; // Default parameter values (a=0,b=1) active?

  public:
    // Constructors

    DistrUnif() : DistribSamp(),a(0.0),b(1.0),defPar(true) {}

    DistrUnif(double pa,double pb) : DistribSamp(),a(pa),b(pb) {
      if (a>=b) throw InvalidParameterException(EXCEPT_MSG(""));
      defPar=(a==0.0 && b==1.0);
    }

    // Public methods

    double getLeftLimit() const {
      return a;
    }

    double getRightLimit() const {
      return b;
    }

    void setLeftLimit(double pa) {
      if (pa>=b) throw InvalidParameterException(EXCEPT_MSG(""));
      a=pa;
      defPar=(a==0.0 && b==1.0);
    }

    void setRightLimit(double pb) {
      if (a>=pb) throw InvalidParameterException(EXCEPT_MSG(""));
      b=pb;
      defPar=(b==1.0 && a==0.0);
    }

    /**
     * Returns value of pdf at given point
     */
    double pdf(double x) const {
      return (x>a && x<b) ? 1.0/(b-a) : 0.0;
    }

    /**
     * Returns value of cdf at given point
     */
    double cdf(double x) const {
      if (x<=a) return 0.0;
      else if (x<b) return (x-a)/(b-a);
      else return 1.0;
    }

    double sample(Generator* gen=0) const {
      if (gen==0) throw NoGeneratorException();
      double u=gen->get();
      return defPar ? u : u*(b-a)+a;
    }
  };
//ENDNS

#endif
