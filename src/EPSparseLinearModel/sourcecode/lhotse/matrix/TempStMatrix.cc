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
 * Desc.:  Definition of class TempStMatrix
 * ------------------------------------------------------------------- */

#include "lhotse/matrix/TempStMatrix.h"
#include "lhotse/matrix/TempMatMethods_AssignElement.h"
#include "lhotse/matrix/StMatrix.h"
#include "lhotse/matrix/BaseVector.h"

//BEGINNS(matrix)
  // Public methods

  StMatrix& TempStMatrix::operator=(const BaseMatrix<double>& arg)
  {
    //cout << "this=" << this << endl;
    //cout << "TempStMatrix-=BM: rep=" << rep << ",repc=" << repCorr << ",arg=" << &arg << endl;
    rep->assignVirtual(&arg);

    return *repCorr;
  }

  StMatrix& TempStMatrix::operator=(const BaseVector<double>& arg)
  {
    //cout << "this=" << this << endl;
    //cout << "TempStMatrix-=BV: rep=" << rep << ",repc=" << repCorr << ",arg=" << &arg << endl;
    rep->assignVirtual(&arg);

    return *repCorr;
  }

  StMatrix& TempStMatrix::operator=(const TempStMatrix& arg)
  {
    //cout << "this=" << this << endl;
    //cout << "TempStMatrix-=TV: rep=" << rep << ",repc=" << repCorr << ",arg.repc=" << arg.repCorr << endl;
    rep->assignVirtual(arg.repCorr);

    return *repCorr;
  }

  StMatrix& TempStMatrix::operator=(double arg)
  {
    //cout << "this=" << this << endl;
    //cout << "TempStMatrix-=S: rep=" << rep << ",repc=" << repCorr << ",arg=" << arg << endl;
    TempMatMethods_AssignElement<double> a(arg);
    rep->assignVirtual(&a);

    return *repCorr;
  }
//ENDNS
