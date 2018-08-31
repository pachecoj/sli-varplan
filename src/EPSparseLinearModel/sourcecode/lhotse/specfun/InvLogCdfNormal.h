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
 * Desc.:  Header class InvLogCdfNormal
 * ------------------------------------------------------------------- */

#ifndef INVLOGCDFNORMAL_H
#define INVLOGCDFNORMAL_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/specfun/default.h"
#include "lhotse/specfun/Specfun.h"
#include "lhotse/optimize/FuncOneDim.h"

//BEGINNS(specfun)
  /**
   * Helper class for 'InvLogCdfNormal::solve'. Implements the function
   * log Phi(x) - z.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class InvLogCdfNormal_Func : public FuncOneDim
  {
  public:
    double z;
    bool flipArg; // Repres. log Phi(-x) - z instead?

    InvLogCdfNormal_Func() : z(0.0),flipArg(false) {}

    bool hasDerivative() const {
      return true;
    }

    void eval(double x,double* f,double* df) {
      double xv=flipArg?(-x):x;

      *f=Specfun::logCdfNormal(xv)-z;
      *df=Specfun::derivLogCdfNormal(xv);
    }
  };

  /**
   * Provides static method for computing the inverse of log Phi(x), Phi
   * the c.d.f. of a standard Gaussian N(0,1). See 'Specfun::logCdfNormal'.
   * We use the Newton method in order to find the zero of
   *   log Phi(x) - z
   * pdf's, cdf's of standard distributions.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class InvLogCdfNormal
  {
  protected:
    // Static members

    static InvLogCdfNormal_Func func;

  public:
    // Public static methods

    /**
     * Finds zero of
     *   log Phi(x) - z.
     * An initial bracket [l,r] must be given. If l==-infty or r==+infty,
     * use 'brinf'. One bracket end must be finite.
     * <p>
     * We use 'Specfun::logCdfNormal' and 'Specfun::derivLogCdfNormal'
     * from class 'Specfun', which requires GSL. We use 'OneDimSolver'
     * from the 'optimize' module. A 'NumericalException' is thrown if
     * x cannot be determined.
     *
     * @param z     Value of z
     * @param l     Left bracket end
     * @param r     Right bracket end
     * @param acc   Absolute accuracy required
     * @param brinf 0: l,r finite. 1: r==+infty ('r' ignored).
     *              2: l==-infty ('l' ignored). Def.: 0
     * @return      Solution x
     */
    static double solve(double z,double l,double r,double acc,int brinf=0);
  };
//ENDNS

#endif
