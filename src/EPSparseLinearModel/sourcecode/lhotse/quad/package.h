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
 * Desc.:  Package interface
 * ------------------------------------------------------------------- */

#ifndef PI_QUAD_H
#define PI_QUAD_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

// Module includes

#include "lhotse/quad/default.h"

#include "lhotse/quad/GaussQuad.h"
#include "lhotse/quad/SoftmaxExMonQuadrature.h"
#include "lhotse/quad/DebugSingleQuad.h"
#include "lhotse/quad/PoissonSingleQuad.h"
#include "lhotse/quad/SingleGaussHermQuadrature.h"
#include "lhotse/quad/SingleQuadrature.h"
#include "lhotse/quad/BinClassSingleLikehood.h"
#include "lhotse/quad/LaplaceSingleLikehood.h"
#include "lhotse/quad/LogisticSingleLikehood.h"
#include "lhotse/quad/NormalSingleLikehood.h"
#include "lhotse/quad/ProbitSingleLikehood.h"
#include "lhotse/quad/RegressSingleLikehood.h"
#include "lhotse/quad/SingleLikehoodFactor.h"
#include "lhotse/quad/MultiQuadrature.h"
#include "lhotse/quad/MultiClassQuadrature.h"
#include "lhotse/quad/MultiProbitQuadrature.h"
#include "lhotse/quad/GaussProdQuadrature.h"
#include "lhotse/quad/ExactMonomQuadrature.h"
#include "lhotse/quad/SoftmaxADFMoments.h"
#include "lhotse/quad/SoftmaxGaussQuadratureDebug.h"
#include "lhotse/quad/SoftmaxGaussQuadrature.h"

#endif
