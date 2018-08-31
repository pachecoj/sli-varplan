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
 * Desc.:  Definition of class LogVector
 * ------------------------------------------------------------------- */

#include "lhotse/matrix/LogVector.h"
#include "lhotse/matrix/StVector.h"
#include "lhotse/FileUtils.h"
#include "lhotse/NumberFormats.h"

//BEGINNS(matrix)
  const int LogVector::maxEmptySize;

  LogVector& LogVector::operator=(const LogVector& src) {
    if (&src!=this) {
      ensureCapacity(src.n,false);
      memmove(buff.p(),src.buff.p(),n*sizeof(double));
      window=src.window;
      active=src.active;
    }
    return *this;
  }

  LogVector& LogVector::operator=(const StVector& src) {
    ensureCapacity(src.size(),false);
    ((StVector&) StVector::mask(buff))=src;

    return *this;
  }

  ofstream& LogVector::save(ofstream& os) const
  {
    StVector::mask(buff)->save(os);
    return os;
  }

  ifstream& LogVector::load(ifstream& is)
  {
    StVector tempV;

    tempV.load(is);
    int sz=tempV.size();
    if (window==-1 || window>=sz) {
      ensureCapacity(sz,false);
      ((StVector&) StVector::mask(buff))=tempV;
    } else {
      ensureCapacity(window,false);
      ((StVector&) StVector::mask(buff))=
	(StVector&) tempV(Range(sz-window,sz-1));
    }

    return is;
  }
//ENDNS
