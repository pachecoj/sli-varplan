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
 * Desc.:  Header abstract class SingleQuadrature
 * ------------------------------------------------------------------- */

#ifndef SINGLEQUADRATURE_H
#define SINGLEQUADRATURE_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/quad/default.h"
#include "lhotse/matrix/StVector.h"

//BEGINNS(quad)
  /**
   * Abstract base class for implementations of tools to compute Gaussian
   * expectations over u^k t(u) and u^k log t(u), k=0,1,2. Here, u\in \R
   * (univariate).
   * t(u) in typical applications is a likelihood factor. Likelihood factors
   * are maintained in 'SingleLikehoodFactor' classes, which uses classes
   * of the 'SingleQuadrature' hierarchy for their implementation.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class SingleQuadrature
  {
  public:
    // Public methods

    virtual ~SingleQuadrature() {}

    /**
     * Define the partition function
     *   Z = E[t(u)], E[.] over N(h,a), h=='mean',a=='var'.
     * The method returns log Z.
     * <p>
     * ADF terms (tilted moments):
     * If 'alpha', 'nu' are given (either both or none), the method also
     * returns
     *   alpha = (d/dh) log Z
     *   nu = -(d^2/dh^2) log Z
     * The moments hhat, ahat of the tilted distribution prop. to
     * t(u) N(u | h,a) are then computed as
     *   hhat = h + a alpha
     *   ahat = a (1 - nu a)
     *
     * @param mean  S.a.
     * @param var   S.a.
     * @param alpha S.a. Optional
     * @param nu    S.a. Optional
     * @return      log Z
     */
    virtual double compLogPart(double mean,double var,double* alpha=0,
			       double* nu=0) const = 0;

    /**
     * Computes
     *   E[log t(u)], E[.] over N(h,a), h=='mean', a=='var'.
     * If 'lmean', 'lvar' are given (either both or none), we also compute
     * the standardized moments
     *   lmean = E[log t(u) v], lvar = E[log t(u) v^2],
     *   v = (u-h)/sqrt(a)
     * NOTE: An implementation is not required by some users of this class.
     * This is why the def. implementation is given here, and the method is
     * not abstract. The def. implementation throws an exception.
     *
     * @param mean  S.a.
     * @param var   S.a.
     * @param lmean S.a. Optional
     * @param lvar  "
     * @return      S.a.
     */
    virtual double compExpectLog(double mean,double var,double* lmean=0,
				 double* lvar=0) const {
      throw NotImplemException(EXCEPT_MSG("Not implemented"));
    }
  };
//ENDNS

#endif
