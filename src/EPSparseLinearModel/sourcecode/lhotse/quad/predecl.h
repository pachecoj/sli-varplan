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
 * Desc.:  Module predeclarations
 * ------------------------------------------------------------------- */

#ifndef DE_QUAD_H
#define DE_QUAD_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

//BEGINNS(quad)
  class GaussQuad;
  class SingleQuadrature;
//  class DebugSingleQuad;
  class SingleGaussHermiteQuadrature;
  class PoissonSingleQuad;
  class PoissonHarrisSingleQuad;
  class SingleLikehoodFactor;
  class BinClassSingleLikehood;
  class RegressSingleLikehood;
  class ProbitSingleLikehood;
  class LogisticSingleLikehood;
  class NormalSingleLikehood;
  class LaplaceSingleLikehood;
  class MultiQuadrature;
  class MultiQuadParams;
  class MultiClassQuadrature;
  class MultiProbitQuadrature;
  class MultiProbitQuadParams;
  class GaussProdQuadrature;
  class ExactMonomQuadrature;
  class SoftmaxADFMoments;
  class SoftmaxExMonQuadrature;
  class SoftmaxGaussQuadrature;
  class SoftmaxGaussQuadratureDebug;
//ENDNS

#endif
