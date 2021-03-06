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
 * Desc.:  Header class SoftmaxGaussQuadrature
 * ------------------------------------------------------------------- */

#ifndef SOFTMAXGAUSSQUADRATUREDEBUG_H
#define SOFTMAXGAUSSQUADRATUREDEBUG_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/quad/default.h"
#include "lhotse/quad/GaussProdQuadrature.h"
#include "lhotse/quad/SoftmaxADFMoments.h"
#include "lhotse/matrix/FastUtils.h"

//BEGINNS(quad)
  /**
   * Combines 'GaussProdQuadrature' and 'SoftmaxADFMoments'.
   * <p>
   * DEBUG variant, without special implementations.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class SoftmaxGaussQuadratureDebug : public GaussProdQuadrature,
				      public SoftmaxADFMoments
  {
  public:
    // Public methods

    /**
     * Constructor.
     *
     * @param dim Dimension of u
     */
    SoftmaxGaussQuadratureDebug(int dim) : GaussProdQuadrature(dim),
					   SoftmaxADFMoments(dim) {}
  };
//ENDNS

#endif
