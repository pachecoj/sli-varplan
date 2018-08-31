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
 * Desc.:  Header class BaseLinMat
 * ------------------------------------------------------------------- */

#ifndef BASELINMAT_H
#define BASELINMAT_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/matrix/default.h"
#include "lhotse/matrix/WriteBackMat.h"

//BEGINNS(matrix)
  /**
   * Helper class used to interface code which requires matrices which are:
   * - contiguous in columns
   * - strided in rows, with a positive striding constant >= number of
   *   rows
   * Call this a "linear" matrix (or "flat").
   * <p>
   * The method 'getLinMat' of 'BaseMatrix' creates a linear matrix based on
   * its current object. If the latter is not an indexed mask, the linear
   * matrix simply refers into the "parent" matrix. Otherwise, a copy is
   * drawn. The buffer watcher stored here will have ref. count == 1 only in
   * this case. If the matrix is to be written back (flag for 'getLinMat'), we
   * store a pointer to the parent object in this case.
   * <p>
   * Structure patterns:
   * 'strpatt' contains the structure pattern for the matrix (see
   * 'WriteBackMat<T>::strctXXX' constants). If 'strpatt'!='strctNormal', only
   * a subpattern in the rect. matrix frame is to be used for the matrix
   * content (for both read/write-back).
   * NOTE: The structure pattern need not be identical with the concrete matrix
   * form, f.ex. a lower-triang. matrix can be stored with an upper-triang.
   * structure pattern. The concrete matrix form is NOT stored in this class.
   * <p>
   * Write back:
   * Happens in the destructor iff 'parMat'!=0. In this case, we write back
   * the flat copy into the matrix 'parMat'. We use 'parMat's 'writeBackBuff'
   * for this.
   * NOTE: 'writeBackBuff' has to conform to the structure pattern.
   * <p>
   * See 'BaseMatrix::getLinMat' and 'BaseMatrix::writeBackBuff' for more
   * details.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  template<class T> class BaseLinMat
  {
  public:
    int m,n;                 // size
    ArrayHandle<T> buff;     // buffer
    int stride;              // striding constant
    BaseMatrix<T>* parMat;   // see header comment
    unsigned char strpatt;   // "

    BaseLinMat() : m(0),n(0),buff(ArrayHandleZero<T>::get()),stride(0),
		   parMat(0),strpatt(WriteBackMat<T>::strctNormal) {}

    void init(int pm,int pn,const ArrayHandle<T>& pbuff,int pstride,
	      unsigned char pstrpatt=WriteBackMat<T>::strctNormal,
	      BaseMatrix<T>* ppmat=0) {
      m=pm; n=pn; buff=pbuff; stride=pstride; parMat=ppmat;
      strpatt=pstrpatt;
    }

    /**
     * Destructor. If 'parMat'!=0, we do a copy-back into this matrix.
     */
    ~BaseLinMat() {
      if (parMat!=0) {
	// 'buff' repres. contiguous matrix of same size as 'parMat'
	parMat->writeBackBuff(buff,strpatt);
      }
    }
  };
//ENDNS

#endif
