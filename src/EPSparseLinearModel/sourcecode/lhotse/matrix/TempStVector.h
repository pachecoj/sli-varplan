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
 * Desc.:  Header class TempStVector
 * ------------------------------------------------------------------- */

#ifndef TEMPSTVECTOR_H
#define TEMPSTVECTOR_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/matrix/default.h"
#include "lhotse/matrix/TempMatMethods.h"

//BEGINNS(matrix)
  /**
   * Helper class for 'StVector'.
   * See header comments of 'TempBaseVector', just substitute 'BaseVector'
   * for 'StVector'.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class TempStVector
  {
    friend class StVector;
    friend class StMatrix;

  protected:
    // Members

    StVector* repCorr; // pointer to representation
    mutable TempMatMethods* rep;

    // Methods
  protected:
    /**
     * Must only be invoked by methods of 'StVector'!
     *
     * @param ppCorr Pointer to newly dyn. created 'StVector' object
     * @param pp     "
     */
    TempStVector(StVector* ppCorr,TempMatMethods* pp) : repCorr(ppCorr),
    rep(pp) {
      //cout << "this=" << this << endl;
      //cout << "TempStVector-C: rep=" << rep << ",repc=" << repCorr << endl;
    }

  public:
    /**
     * Copy constructor.
     * NOTE: Not a CC in the usual sense! The repr. pointer of the 'src'
     * object is set to 0, to after the copying 'src' can only be destroyed.
     * ==> just to allow for creation and copying of temp. versions
     */
    TempStVector(const TempStVector& src) : rep(src.rep),repCorr(src.repCorr) {
      //cout << "this=" << this << endl;
      //cout << "TempStVector-CC: rep=" << rep << ",repc=" << repCorr << endl;
      src.rep=0; // 'src' is not usable from now on!
    }

    ~TempStVector() {
      //cout << "this=" << this << endl;
      //cout << "TempStVector-D: rep=" << rep << ",repc=" << repCorr << endl;
      if (rep!=0) { delete rep; rep=0; repCorr=0; }
    }

    /**
     * Allows to call 'StVector' methods via dereferencing.
     *
     * @return Pointer to underlying reference
     */
    StVector* operator->() const {
      //cout << "this=" << this << endl;
      //cout << "TempStVector-->: rep=" << rep << ",repc=" << repCorr << endl;
      return repCorr;
    }

    /**
     * Drives automatic conversion to 'StVector&' (non-const), so objects
     * can be used as l- or r-value arguments.
     *
     * @return Ref. to repres.
     */
    operator StVector&() const {
      //cout << "this=" << this << endl;
      //cout << "TempStVector-&: rep=" << rep << ",repc=" << repCorr << endl;
      return *repCorr;
    }

    /**
     * See 'TempBaseVector'.
     *
     * @param arg Source argument
     */
    StVector& operator=(const BaseVector<double>& arg);

    /**
     * See 'TempBaseVector'.
     *
     * @param arg Source argument
     */
    StVector& operator=(const TempStVector& arg);

    /**
     * See 'TempBaseVector'.
     *
     * @param arg Source argument
     */
    StVector& operator=(double arg);
  };
//ENDNS

#endif
