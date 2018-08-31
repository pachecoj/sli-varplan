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
 * Desc.:  Definition of class TempIndexVector
 * ------------------------------------------------------------------- */

#include "lhotse/matrix/TempIndexVector.h"
#include "lhotse/matrix/TempMatMethods_AssignElement.h"
#include "lhotse/matrix/IndexVector.h"

//BEGINNS(matrix)
  // Public methods

  IndexVector& TempIndexVector::operator=(const BaseVector<int>& arg)
  {
    //cout << "this=" << this << endl;
    //cout << "TempIndexVector-=BV: rep=" << rep << ",repc=" << repCorr << ",arg=" << &arg << endl;
    rep->assignVirtual(&arg);

    return *repCorr;
  }

  IndexVector& TempIndexVector::operator=(const TempIndexVector& arg)
  {
    //cout << "this=" << this << endl;
    //cout << "TempIndexVector-=TV: rep=" << rep << ",repc=" << repCorr << ",arg.repc=" << arg.repCorr << endl;
    rep->assignVirtual(arg.repCorr);

    return *repCorr;
  }

  IndexVector& TempIndexVector::operator=(int arg)
  {
    //cout << "this=" << this << endl;
    //cout << "TempIndexVector-=S: rep=" << rep << ",repc=" << repCorr << ",arg=" << arg << endl;
    TempMatMethods_AssignElement<int> a(arg);
    rep->assignVirtual(&a);

    return *repCorr;
  }
//ENDNS
