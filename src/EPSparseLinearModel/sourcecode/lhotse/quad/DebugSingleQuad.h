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
 * Desc.:  Header abstract class DebugSingleQuad
 * ------------------------------------------------------------------- */

#ifndef DEBUGSINGLEQUAD_H
#define DEBUGSINGLEQUAD_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/quad/default.h"
#include "lhotse/quad/SingleQuadrature.h"

//BEGINNS(quad)
  /**
   * FOR DEBUGGING!
   *
   * t(u) is the linear part of a normal noise model:
   *   t(u) = sigma^{-2} y u
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class DebugSingleQuad : public SingleQuadrature
  {
  protected:
    // Members

    double y,sigsq; // Parameters y (mean), sigma^2 (variance)

  public:
    // Public methods

    /**
     * Default constructor.
     */
    DebugSingleQuad() : SingleQuadrature(),y(0.0),sigsq(1.0) {}

    /**
     * @return Value of y
     */
    double getMean() const {
      return y;
    }

    /**
     * @param ny New value for y
     */
    void setMean(double ny) {
      y=ny;
    }

    /**
     * @return Value of sigma^2
     */
    double getVar() const {
      return sigsq;
    }

    /**
     * @param ns New value for sigma^2
     */
    void setVar(double ns) {
      if (ns<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
      sigsq=ns;
    }

    double compLogPart(double mean,double var,double* alpha=0,double* nu=0)
      const;
    
    double compExpectLog(double mean,double var,double* lmean=0,
			 double* lvar=0) const;
  };

  // Inline methods

  inline double
  DebugSingleQuad::compLogPart(double mean,double var,double* alpha,
				double* nu) const
  {
    double temp=y/sigsq,logpart;

    logpart=0.5*temp*(2.0*mean+var*temp);
    if (alpha!=0) {
      if (nu==0) throw InvalidParameterException(EXCEPT_MSG("alpha,nu"));
      *alpha=temp; *nu=0.0;
    }

    return logpart;
  }

  inline double
  DebugSingleQuad::compExpectLog(double mean,double var,double* lmean,
				 double* lvar) const
  {
    double temp=y/sigsq,elog;

    elog=temp*mean;
    if (lmean!=0) {
      if (lvar==0) throw InvalidParameterException(EXCEPT_MSG(""));
      *lmean=temp*sqrt(var);
      *lvar=elog;
    }

    return elog;
  }
//ENDNS

#endif
