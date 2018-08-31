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
 * Project source file
 * Module: eplin
 * Desc.:  Header class MaskEPPrior
 * ------------------------------------------------------------------- */

#ifndef EPLIN_MASKEPPRIOR_H
#define EPLIN_MASKEPPRIOR_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eplin/default.h"
#include "src/eplin/EPPrior.h"
#include "lhotse/matrix/StVector.h"
#include "lhotse/matrix/StMatrix.h"

//BEGINNS(eplin)
  /**
   * Standard implementation of 'EPPrior' for dense matrix X. Configured by
   * passing refs. for X, v (at construction), of which masks are maintained
   * (no copies are drawn!). X^T can be passed instead of X, in which case
   * 'trans'==true.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class MaskEPPrior : public EPPrior
  {
  protected:
    StMatrix xmsk;         // Mask to X (or X^T if 'trans'==true)
    bool trans;
    StVector vmsk;         // Mask to v
    double sigsq;          // Noise variance sigma^2 (-1 if not maintained)
    StVector b0;           // Vector b^(0)

  public:
    // Public methods

    /**
     * Constructor
     *
     * @param xm  Matrix X (or X^T), masked here
     * @param vv  Vector v, masked here
     * @param sig Value for sigma^2. Def: -1 (not maintained here)
     * @param trs Is 'xm'==X^T? Otherwise: X. Def.: false
     */
    MaskEPPrior(const StMatrix& xm,const StVector& vv,double sig=-1.0,
		bool trs=false);

    /**
     * Since a mask of X (or X^T) and v is kept here only, the underlying
     * objects can be changed (their size must not change!). If this is
     * the case, this object must be notified by calling this method.
     */
    void update();

    int rows() const {
      return trans?xmsk.cols():xmsk.rows();
    }

    int cols() const {
      return trans?xmsk.rows():xmsk.cols();
    }

    double getSigSq() const {
      return sigsq;
    }

    const StVector& getBVec() const {
      return b0;
    }

    void getVVec(StVector& vvec,bool accum=false) const {
      if (!accum) vvec=vmsk;
      else vvec.addprod(1.0,vmsk);
    }

    void getRow(int i,StVector& row) const;

    void getCol(int i,StVector& col) const;

    void getXMat(StMatrix& xmat) const;

    void multT(const StVector& v,StVector& res) const;

    void multT(const StMatrix& v,StMatrix& res) const;

    void mult(const StVector& v,StVector& res) const;

    void mult(const StMatrix& v,StMatrix& res) const;

    void outerProd(const StVector& dg,StMatrix& res) const;

    void outerTProd(StMatrix& res) const;
  };
//ENDNS

#endif
