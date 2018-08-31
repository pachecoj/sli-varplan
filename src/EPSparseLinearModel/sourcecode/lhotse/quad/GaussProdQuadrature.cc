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
 * Desc.:  Definition of class GaussProdQuadrature
 * ------------------------------------------------------------------- */

#include "lhotse/quad/GaussProdQuadrature.h"

//BEGINNS(quad)
  // Static members

  // off=sqrt(3), w0=2/3, w1=1/6
  const double GaussProdQuadrature::off(1.73205080756887729352);
  const double GaussProdQuadrature::logw0(-0.40546510810816438198);
  const double GaussProdQuadrature::logw1(-1.79175946922805500085);

//ENDNS
