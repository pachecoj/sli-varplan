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
 * Module: matif
 * Desc.:  Header class MatlabString
 * ------------------------------------------------------------------- */

#ifndef MATLABSTRING_H
#define MATLABSTRING_H

#include "lhotse/matif/default.h"

//BEGINNS(matif)
  /**
   * Reads string from Matlab char row vector and stores non-persistent
   * copy in 'buff'.
   * NOTE: Present version does NOT allow for persistent objects.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class MatlabString
  {
  protected:
    // Members

    char* buff;  // String buffer (zero-terminated)

  public:
    // Constructors

    /**
     * Input is 'mxArray' which must be a char row vector
     *
     * @param arr Matlab char row vector
     */
    MatlabString(const mxArray* arr) {
      if (mxIsChar(arr)!=1)
	throw MatIFException("String argument (row vector) expected");
      if (mxGetM(arr)!=1)
	throw MatIFException("String argument (row vector) expected");
      int len=mxGetN(arr)+1;
      buff=(char*) mxMalloc(len); // non-persistent
      mxGetString(arr,buff,len); // convert
    }

    virtual ~MatlabString() {
      mxFree(buff);
    }

    // Overload new, delete

    void* operator new(size_t sz,bool pflag=false) {
      void* ptr=mxMalloc(sz); // non-persistent
      return ptr;
    }

    void operator delete(void* ptr) {
      mxFree(ptr);
    }

    // Public methods

    /**
     * @return String buffer (0-term.)
     */
    const char* getStr() const {
      return buff;
    }
  };
//ENDNS

#endif
