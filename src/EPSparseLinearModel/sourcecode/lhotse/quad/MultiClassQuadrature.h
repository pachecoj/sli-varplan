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
 * Desc.:  Header abstract class MultiClassQuadrature
 * ------------------------------------------------------------------- */

#ifndef MULTICLASSQUADRATURE_H
#define MULTICLASSQUADRATURE_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/quad/default.h"
#include "lhotse/quad/MultiQuadrature.h"

//BEGINNS(quad)
  /**
   * Abstract subclass of 'MultiQuadrature' for multi-class likelihoods.
   * Has member 'cls' for current target value and vector 'bias' of
   * intercept (bias) parameters.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class MultiClassQuadrature : public virtual MultiQuadrature
  {
  public:
    // Public methods

    /**
     * @param cl New current class
     */
    virtual void setCls(int cl) = 0;

    /**
     * @return Current class
     */
    virtual int getCls() const = 0;

    /**
     * @param bvec New bias vector
     */
    virtual void setBias(const StVector& bvec) = 0;

    /**
     * @return Bias vector
     */
    virtual const StVector& getBias() const = 0;

    /**
     * Reset bias vector to 0
     */
    virtual void resetBias() = 0;
  };
//ENDNS

#endif
