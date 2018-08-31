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
 * Desc.:  Header class TempBaseMatrix
 * ------------------------------------------------------------------- */

#ifndef TEMPBASEMATRIX_H
#define TEMPBASEMATRIX_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/matrix/default.h"
#include "lhotse/matrix/TempMatMethods.h"
#include "lhotse/matrix/TempMatMethods_AssignElement.h"

//BEGINNS(matrix)
  /**
   * Helper class for 'BaseMatrix'.
   * Much the same as 'TempBaseVector' for 'BaseVector', see comments there.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  template <class T> class TempBaseMatrix
  {
    friend class BaseMatrix<T>;

  protected:
    // Members

    BaseMatrix<T>* repCorr; // pointer to representation
    mutable TempMatMethods* rep;

    // Methods
  protected:
    /**
     * Must only be invoked by methods of 'BaseMatrix'!
     *
     * @param pp Pointer to newly dyn. created 'BaseMatrix' object
     */
    TempBaseMatrix(BaseMatrix<T>* ppCorr,TempMatMethods* pp) : repCorr(ppCorr),
    rep(pp) {}

  public:
    /**
     * Copy constructor.
     * NOTE: Not a CC in the usual sense! The repr. pointer of the 'src'
     * object is set to 0, so after the copying 'src' can only be destroyed.
     * ==> just to allow for creation and copying of temp. versions
     */
    TempBaseMatrix(const TempBaseMatrix<T>& src) : rep(src.rep),
    repCorr(src.repCorr) {
      src.rep=0; // 'src' is not usable from now on!
    }

    ~TempBaseMatrix() {
      if (rep!=0) delete rep;
    }
    
    /**
     * Allows to call 'BaseMatrix' methods via dereferencing.
     *
     * @return Pointer to underlying reference
     */
    BaseMatrix<T>* operator->() const {
      return repCorr;
    }

    /**
     * Drives automatic conversion to 'BaseMatrix&' (non-const), so objects
     * can be used as l- or r-value arguments.
     *
     * @return Ref. to repres.
     */
    operator BaseMatrix<T>&() const {
      return *repCorr;
    }

    /**
     * Implements '=' for user convenience. It allows you to write
     *   mat(rng1,rng2)=mat2;
     * instead of
     *   ((BaseMatrix<T>&) mat(rng1,rng2))=mat2;
     * NOTE: The latter is still valid, and does the same thing.
     *
     * @param arg Source argument (matrix or vector)
     */
    BaseMatrix<T>& operator=(const BaseMatrix<T>& arg) {
      //cout << "TempBaseMatrix::operator=" << endl;
      rep->assignVirtual(&arg);

      return *repCorr;
    }

    BaseMatrix<T>& operator=(const BaseVector<T>& arg) {
      //cout << "TempBaseMatrix::operator=" << endl;
      rep->assignVirtual(&arg);

      return *repCorr;
    }

    /**
     * We need this one as well, otherwise
     *   mat1(rng1,rng3)=mat2(rng2,rng4);
     * breaks! Why? In this case, the automatic conversion to
     * BaseMatrix<T>& is somehow not used, instead the compiler invents a
     * '=' operator which just copies members (which are pointers) across!
     *
     * @param arg Source argument (matrix or vector)
     */
    BaseMatrix<T>& operator=(const TempBaseMatrix<T>& arg) {
      //cout << "TempBaseMatrix-=TV: rep=" << rep << ",repc=" << repCorr << ",arg.repc=" << arg.repCorr << endl;
      rep->assignVirtual(arg.repCorr);

      return *repCorr;
    }

    /**
     * Implements '=' with scalar r.h.s. for user convenience. It allows
     * you to write
     *   mat(rng1,rng2)=scalar;
     * instead of
     *   ((BaseMatrix<T>&) mat(rng1,rng2))=scalar;
     * NOTE: The latter is still valid, and does the same thing.
     *
     * @param arg Source argument
     */
    BaseMatrix<T>& operator=(T arg) {
      //cout << "TempBaseMatrix::operator=" << endl;
      TempMatMethods_AssignElement<T> a(arg);
      rep->assignVirtual(&a);

      return *repCorr;
    }
  };
//ENDNS

#endif
