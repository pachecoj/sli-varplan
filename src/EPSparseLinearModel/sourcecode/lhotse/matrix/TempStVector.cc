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
 * Desc.:  Definition of class TempStVector
 * ------------------------------------------------------------------- */

#include "lhotse/matrix/TempStVector.h"
#include "lhotse/matrix/TempMatMethods_AssignElement.h"
#include "lhotse/matrix/StVector.h"

//BEGINNS(matrix)
  // Public methods

  StVector& TempStVector::operator=(const BaseVector<double>& arg)
  {
    //cout << "this=" << this << endl;
    //cout << "TempStVector-=BV: rep=" << rep << ",repc=" << repCorr << ",arg=" << &arg << endl;
    rep->assignVirtual(&arg);

    return *repCorr;
  }

  StVector& TempStVector::operator=(const TempStVector& arg)
  {
    //cout << "this=" << this << endl;
    //cout << "TempStVector-=TV: rep=" << rep << ",repc=" << repCorr << ",arg.repc=" << arg.repCorr << endl;
    rep->assignVirtual(arg.repCorr);

    return *repCorr;
  }

  StVector& TempStVector::operator=(double arg)
  {
    //cout << "this=" << this << endl;
    //cout << "TempStVector-=S: rep=" << rep << ",repc=" << repCorr << ",arg=" << arg << endl;
    TempMatMethods_AssignElement<double> a(arg);
    rep->assignVirtual(&a);

    return *repCorr;
  }
//ENDNS
