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
 * Module: specfun
 * Desc.:  Definition of class InvLogCdfNormal
 * ------------------------------------------------------------------- */

#include "lhotse/specfun/InvLogCdfNormal.h"
#include "lhotse/optimize/OneDimSolver.h"

//BEGINNS(specfun)
  // Static members

  InvLogCdfNormal_Func InvLogCdfNormal::func;

  // Static methods

  double InvLogCdfNormal::solve(double z,double l,double r,double acc,
				int brinf)
  {
    double x;

    func.flipArg=(brinf==2);
    if (brinf==0 && l>=r) throw InvalidParameterException(EXCEPT_MSG(""));
    else if (func.flipArg) l=-r;
    func.z=z;
    x=OneDimSolver::newton(&func,l,r,acc,
			   (brinf==0)?OneDimSolver::brackRightRegular:
			   OneDimSolver::brackRightInfinite);
    if (func.flipArg) x=-x;

    return x;
  }
//ENDNS
