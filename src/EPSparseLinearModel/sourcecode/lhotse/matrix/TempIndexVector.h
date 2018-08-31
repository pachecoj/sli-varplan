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
 * Desc.:  Header class TempIndexVector
 * ------------------------------------------------------------------- */

#ifndef TEMPINDEXVECTOR_H
#define TEMPINDEXVECTOR_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/matrix/default.h"
#include "lhotse/matrix/TempMatMethods.h"

//BEGINNS(matrix)
  /**
   * Helper class for 'IndexVector'.
   * See header comments of 'TempBaseVector', just substitute 'BaseVector'
   * for 'IndexVector'.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class TempIndexVector
  {
    friend class IndexVector;

  protected:
    // Members

    IndexVector* repCorr; // pointer to representation
    mutable TempMatMethods* rep;

    // Methods
  protected:
    /**
     * Must only be invoked by methods of 'IndexVector'!
     *
     * @param ppCorr Pointer to newly dyn. created object
     * @param pp     "
     */
    TempIndexVector(IndexVector* ppCorr,TempMatMethods* pp) :
      repCorr(ppCorr),rep(pp) {}

  public:
    /**
     * Copy constructor.
     * NOTE: Not a CC in the usual sense! The repr. pointer of the 'src'
     * object is set to 0, to after the copying 'src' can only be destroyed.
     * ==> just to allow for creation and copying of temp. versions
     */
    TempIndexVector(const TempIndexVector& src) : rep(src.rep),
    repCorr(src.repCorr) {
      src.rep=0; // 'src' is not usable from now on!
    }

    ~TempIndexVector() {
      if (rep!=0) delete rep;
    }
    
    /**
     * Allows to call 'IndexVector' methods via dereferencing.
     *
     * @return Pointer to underlying reference
     */
    IndexVector* operator->() const {
      return repCorr;
    }

    /**
     * Drives automatic conversion to 'IndexVector&' (non-const), so objects
     * can be used as l- or r-value arguments.
     *
     * @return Ref. to repres.
     */
    operator IndexVector&() const {
      return *repCorr;
    }

    /**
     * See 'TempBaseVector'.
     *
     * @param arg Source argument
     */
    IndexVector& operator=(const BaseVector<int>& arg);

    /**
     * See 'TempBaseVector'.
     *
     * @param arg Source argument
     */
    IndexVector& operator=(const TempIndexVector& arg);

    /**
     * See 'TempBaseVector'.
     *
     * @param arg Source argument
     */
    IndexVector& operator=(int arg);
  };
//ENDNS

#endif
