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
 * Desc.:  Header abstract class EPPrior
 * ------------------------------------------------------------------- */

#ifndef EPLIN_EPPRIOR_H
#define EPLIN_EPPRIOR_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eplin/default.h"
#include "lhotse/matrix/predecl.h"

//BEGINNS(eplin)
  /**
   * Class representing prior towards 'ExpectPropLinear'. The prior has
   * density of a Gaussian with natural parameters:
   *   Pi^(0) = X^T X, X m-by-n
   *   b^(0)  = X^T v
   * for some v. This is not a normalizable distribution for m<n. Services
   * based on X have to be implemented, and b^(0) must be returned by
   * 'getBVec', v by 'getVVec'.
   * <p>
   * Noise variance sigma^2:
   * This is maintained here (as fixed value) iff 'getSigSq' ret. a positive
   * value. If sigma^2 is not maintained here, 'getSigSq' ret. -1.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class EPPrior
  {
  public:
    // Public methods

    /**
     * @return Number m of rows of X (size of v)
     */
    virtual int rows() const = 0;

    /**
     * @return Number n of cols. of X (size of b^(0))
     */
    virtual int cols() const = 0;

    /**
     * DOWNW. COMPAT.
     *
     * @return Noise stddev. sigma
     */
    virtual double getSigma() const {
      double sigsq=getSigSq();

      return (sigsq<0.0)?(-1.0):sqrt(sigsq);
    }

    /**
     * @return Noise variance sigma^2; or -1 if noise var. is not
     *         maintained here
     */
    virtual double getSigSq() const = 0;

    /**
     * @return Natural vector b^(0)
     */
    virtual const StVector& getBVec() const = 0;

    /**
     * @param vvec  Vector v ret. here, s.t. b^(0) = X v
     * @param accum Shall v be added to 'vvec'? Def: false
     */
    virtual void getVVec(StVector& vvec,bool accum=false) const = 0;

    /**
     * @param i   Row number
     * @param row Row i of X ret. here
     */
    virtual void getRow(int i,StVector& row) const = 0;

    /**
     * @param i   Col. number
     * @param col Col i of X ret. here
     */
    virtual void getCol(int i,StVector& col) const = 0;

    /**
     * NOTE: This is just a fix for the moment. For large sparse matrices X,
     * this will not be feasible. Have to replace methods relying on
     * 'getXMat' by s.th. smarter then!
     *
     * @param xmat Matrix X (dense) ret. here
     */
    virtual void getXMat(StMatrix& xmat) const = 0;

    /**
     * @param a   Input vector a
     * @param res Result vector res = X^T a ret. here
     */
    virtual void multT(const StVector& a,StVector& res) const = 0;

    /**
     * @param a   Input matrix A
     * @param res Result matrix Res = X^T A ret. here
     */
    virtual void multT(const StMatrix& a,StMatrix& res) const = 0;

    /**
     * @param a   Input vector a
     * @param res Result vector res = X a ret. here
     */
    virtual void mult(const StVector& a,StVector& res) const = 0;

    /**
     * @param a   Input matrix A
     * @param res Result matrix Res = X A ret. here
     */
    virtual void mult(const StMatrix& a,StMatrix& res) const = 0;

    /**
     * NOTE: 'dg' may contain negative entries as well.
     *
     * @param dg  Diagonal matrix D (as vector)
     * @param res Matrix X D X^T ret. here
     */
    virtual void outerProd(const StVector& dg,StMatrix& res) const = 0;

    /**
     * @param res Matrix X^T X ret. here
     */
    virtual void outerTProd(StMatrix& res) const = 0;
  };
//ENDNS

#endif
