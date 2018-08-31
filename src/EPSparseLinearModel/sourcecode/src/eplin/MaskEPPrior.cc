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
 * Desc.:  Definition of class MaskEPPrior
 * ------------------------------------------------------------------- */

#include "src/eplin/MaskEPPrior.h"

//BEGINNS(ep)

  MaskEPPrior::MaskEPPrior(const StMatrix& xm,const StVector& vv,
			   double sig,bool trs) : trans(trs)
  {
    sigsq=(sig<=0.0)?(-1.0):sig;
    if (vv.size()!=(trs?xm.cols():xm.rows()))
      throw WrongDimensionException(EXCEPT_MSG(""));
    xmsk.reassign(xm); vmsk.reassign(vv);
    update(); // Compute inner repres.
  }

  void MaskEPPrior::update()
  {
    xmsk.mulVec(b0,vmsk,!trans);
  }

  void MaskEPPrior::getRow(int i,StVector& row) const
  {
    if (i<0 || i>=rows()) throw OutOfRangeException(EXCEPT_MSG(""));
    if (!trans) row=xmsk(i);
    else row=xmsk(full(),i);
  }

  void MaskEPPrior::getCol(int i,StVector& col) const
  {
    if (i<0 || i>=cols()) throw OutOfRangeException(EXCEPT_MSG(""));
    if (trans) col=xmsk(i);
    else col=xmsk(full(),i);
  }

  void MaskEPPrior::getXMat(StMatrix& xmat) const
  {
    if (!trans) xmat=xmsk;
    else xmat.trans(xmsk);
  }

  void MaskEPPrior::multT(const StVector& v,StVector& res) const
  {
    if (v.size()!=rows()) throw WrongDimensionException(EXCEPT_MSG(""));
    xmsk.mulVec(res,v,!trans);
  }

  void MaskEPPrior::multT(const StMatrix& v,StMatrix& res) const
  {
    if (v.rows()!=rows()) throw WrongDimensionException(EXCEPT_MSG(""));
    res.mul(xmsk,!trans,v,false);
  }

  void MaskEPPrior::mult(const StVector& v,StVector& res) const
  {
    if (v.size()!=cols()) throw WrongDimensionException(EXCEPT_MSG(""));
    xmsk.mulVec(res,v,trans);
  }

  void MaskEPPrior::mult(const StMatrix& v,StMatrix& res) const
  {
    if (v.rows()!=cols()) throw WrongDimensionException(EXCEPT_MSG(""));
    res.mul(xmsk,trans,v,false);
  }

  void MaskEPPrior::outerProd(const StVector& dg,StMatrix& res) const
  {
    if (dg.size()!=cols()) throw WrongDimensionException(EXCEPT_MSG(""));
    res.symMulDiag(xmsk,trans,dg); // X D X^T
    res.makeSymmAuto();
  }

  void MaskEPPrior::outerTProd(StMatrix& res) const
  {
    res.symMul(xmsk,!trans); // X^T X
    res.makeSymmAuto();
  }
//ENDNS
