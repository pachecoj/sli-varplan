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
 * Module: matrix
 * Desc.:  Package interface
 * ------------------------------------------------------------------- */

#ifndef PI_MATRIX_H
#define PI_MATRIX_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

// Module includes

#include "lhotse/matrix/default.h"

#include "lhotse/matrix/ArrayUtils.h"
#include "lhotse/matrix/BaseLinMat.h"
#include "lhotse/matrix/BaseLinVec.h"
#include "lhotse/matrix/BaseMatrix.h"
#include "lhotse/matrix/BaseVector.h"
#include "lhotse/matrix/FastUtils.h"
#include "lhotse/matrix/FinSet.h"
#include "lhotse/matrix/IncompleteCholesky.h"
#include "lhotse/matrix/IndexVector.h"
#include "lhotse/matrix/LogVector.h"
#include "lhotse/matrix/MatDefMembers.h"
#include "lhotse/matrix/MatTimeStamp.h"
#include "lhotse/matrix/Matrix.h"
#include "lhotse/matrix/StMatrix.h"
#include "lhotse/matrix/StVector.h"
#include "lhotse/matrix/TempBaseMatrix.h"
#include "lhotse/matrix/TempBaseVector.h"
#include "lhotse/matrix/TempIndexVector.h"
#include "lhotse/matrix/TempMatMethods.h"
#include "lhotse/matrix/TempStMatrix.h"
#include "lhotse/matrix/TempStVector.h"
#include "lhotse/matrix/TransInPlace.h"
#include "lhotse/matrix/Vector.h"
#include "lhotse/matrix/WriteBackMat.h"
#include "lhotse/matrix/WriteBackVec.h"
#include "lhotse/matrix/cblas_for_cpp.h"
#include "lhotse/matrix/SimpSparseMatrix.h"
#include "lhotse/matrix/ICholKernelMatrix.h"
#include "lhotse/matrix/IncompCholMatrix.h"

#endif
