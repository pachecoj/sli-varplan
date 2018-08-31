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
 * Module: optimize
 * Desc.:  Header abstract class FuncOneDim
 * ------------------------------------------------------------------- */

#ifndef FUNCONEDIM_H
#define FUNCONEDIM_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/optimize/default.h"

//BEGINNS(optimize)
  /**
   * Simple abstract class to represent a function from a subset of \R into
   * \R, together with its first derivative. This is used by the 1-D solvers
   * in class 'OneDimSolver'.
   * All subclasses must implement the virtual methods def. here. Most of them
   * will have additional methods, e.g. to chance other parameters than the
   * eval. point the function depends upon.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class FuncOneDim
  {
  public:
    // Destructor
    virtual ~FuncOneDim() {}

    // Public methods

    /**
     * If this method returns true, then 'eval' will return valid deriv.
     * values via 'df', otherwise the return of 'df' is undefined.
     *
     * @return See above
     */
    virtual bool hasDerivative() const = 0;

    /**
     * Returns function value and derivative at a point 'x' in the domain. If
     * 'x' is not in the domain, a 'InvalidParameterException' should be
     * thrown.
     * NOTE: The deriv. is only returned if 'hasDerivative' returns true,
     * otherwise 'df' is not touched.
     *
     * @param x  Eval. point
     * @param f  Function value ret. here
     * @param df Derivative ret. here (but see above!)
     */
    virtual void eval(double x,double* f,double* df)=0;
  };
//ENDNS

#endif
