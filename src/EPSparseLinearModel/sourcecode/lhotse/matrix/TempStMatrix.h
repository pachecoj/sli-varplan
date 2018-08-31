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
 * Desc.:  Header class TempStMatrix
 * ------------------------------------------------------------------- */

#ifndef TEMPSTMATRIX_H
#define TEMPSTMATRIX_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/matrix/default.h"
#include "lhotse/matrix/TempMatMethods.h"
#include "lhotse/matrix/TempMatMethods_AssignElement.h"

//BEGINNS(matrix)
  /**
   * Helper class for 'StMatrix'.
   * Much the same as 'TempBaseVector' for 'BaseVector', see comments there.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class TempStMatrix
  {
    friend class StMatrix;

  protected:
    // Members

    StMatrix* repCorr; // pointer to representation
    mutable TempMatMethods* rep;

    // Methods

  protected:
    /**
     * Must only be invoked by methods of 'StMatrix'!
     *
     * @param pp Pointer to newly dyn. created 'BaseMatrix' object
     */
    TempStMatrix(StMatrix* ppCorr,TempMatMethods* pp) : repCorr(ppCorr),
    rep(pp) {}

  public:
    TempStMatrix(const TempStMatrix& src) : rep(src.rep),repCorr(src.repCorr) {
      src.rep=0; // 'src' is not usable from now on!
    }

    ~TempStMatrix() {
      if (rep!=0) delete rep;
    }
    
    /**
     * Allows to call 'StMatrix' methods via dereferencing.
     *
     * @return Pointer to underlying reference
     */
    StMatrix* operator->() const {
      return repCorr;
    }

    /**
     * Drives automatic conversion to 'StMatrix&' (non-const), so objects
     * can be used as l- or r-value arguments.
     *
     * @return Ref. to repres.
     */
    operator StMatrix&() const {
      return *repCorr;
    }

    /**
     * See 'TempBaseMatrix'.
     *
     * @param arg Source argument
     */
    StMatrix& operator=(const BaseMatrix<double>& arg);

    StMatrix& operator=(const BaseVector<double>& arg);

    /**
     * See 'TempBaseMatrix'.
     *
     * @param arg Source argument
     */
    StMatrix& operator=(const TempStMatrix& arg);

    /**
     * See 'TempBaseMatrix'.
     *
     * @param arg Source argument
     */
    StMatrix& operator=(double arg);
  };
//ENDNS

#endif
