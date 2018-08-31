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
 * Desc.:  Module predeclarations
 * ------------------------------------------------------------------- */

#ifndef DE_MATRIX_H
#define DE_MATRIX_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

//BEGINNS(matrix)
  class WrongDimensionException;
  class DimMismatchException;
  class EqualArgsException;
  class MatrixCreationException;
  class MaskObjectException;
  class MaskNotImplException;

  class Vector;
  template<class T> class ArrayUtilsBasic;
  template<class T> class ArrayUtils;
  template<class T> class BaseLinVec;
  template<class T> class BaseLinMat;
  class FastUtils;
  class TempMatMethods;
  template<class T> class TempMatMethods_AssignElement;
  template<class T> class MatDefMembers;
  class MatTimeStamp;
  template<class T> class WriteBackVec;
  template<class T> class TempBaseVector;
  template<class T> class BaseVector;
  class TempStVector;
  class StVector;
  class TempIndexVector;
  class IndexVector;
  template<class T> class FinSet;
  class LogVector;
  class Matrix;
  template<class T> class TransInPlace;
  class MatStrct;
  template<class T> class WriteBackMat;
  template<class T> class TempBaseMatrix;
  template<class T> class BaseMatrix;
  class TempStMatrix;
  class StMatrix;
  class IncompCholMatrix;
  class IncompleteCholesky;
  class ICholKernelMatrix;
  class SimpSparseMatrix;
//ENDNS

#endif
