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
 * Desc.:  Header class BaseLinVec
 * ------------------------------------------------------------------- */

#ifndef BASELINVEC_H
#define BASELINVEC_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/matrix/default.h"
#include "lhotse/matrix/WriteBackVec.h"
#include "lhotse/matrix/ArrayUtils.h"

//BEGINNS(matrix)
  /**
   * Helper class used to interface code which does not allow for indexing
   * (like in indexed masks), i.e. vectors either have to be flat or linear
   * (arb. but fixed step size).
   * <p>
   * The method 'getLinVec' of 'BaseVector' creates a linear vector based on
   * its current object. If the latter is not an indexed mask, the linear
   * vector simply refers into the "parent" vector. Otherwise, a copy is
   * drawn (which is always flat; step size = 1). The buffer watcher stored
   * here will have ref. count == 1 only in this case. If the vector is to
   * be written back (flag for 'getLinVec'), we store a pointer to the
   * parent object in this case.
   * <p>
   * Write back:
   * Happens in the destructor iff 'parVec'!=0. In this case, we write back
   * the flat copy into the vector 'parVec'. We use 'parVec's 'writeBackBuff'
   * for this.
   * <p>
   * Negative step size:
   * There are different conventions for 'step'<0:
   * - Fortran (BLAS): Buffer pointer ref. to last element of the vector.
   *   This means absolute indexes into the buffer are non-neg. Buffer ptr.
   *   always ref. to leftmost used element in memory
   *     x_i --> buff[(i-n+1)*step], i=0,...,n-1
   * - Ours: Buffer pointer ref. to first element of vector, no matter what
   *   'step' is. Absolute indexes are negative if 'step'<0
   *     x_i -> buff[i*step], i=0,...,n-1
   * In this class, 'buff' follows the Fortran convention. The user should
   * use 'getBuff' to obtain the buffer pointer which automatically choses
   * the correct conv. (dep. on 'negFort': true -> Fortran; false -> ours).
   * <p>
   * See 'BaseVector::getLinVec' for more details.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  template<class T> class BaseLinVec
  {
  public:
    int n;                   // size
    ArrayHandle<T> buff;     // vector buffer
    int step;                // step size
    bool negFort;            // see header comm.
    BaseVector<T>* parVec;   // see header comment (for write-back)

    BaseLinVec() : n(0),buff(ArrayHandleZero<T>::get()),step(1),parVec(0),
		   negFort(true) {}

    /**
     * NOTE: Indep. of 'negFrt', 'pbuff' must (for neg. step size) conform
     * to the Fortran convention! 'negFort' is just used in 'getBuff'.
     */
    void init(int pn,const ArrayHandle<T>& pbuff,int pstep=1,bool negFrt=true,
	      BaseVector<T>* ppvec=0) {
      n=pn; buff=pbuff; step=pstep; parVec=ppvec;
      negFort=negFrt;
    }

    /**
     * Destructor. If 'parVec'!=0, we do a write-back into this vector.
     */
    ~BaseLinVec() {
      if (parVec!=0) {
	// NOTE: Vector must be flat here and has the same size as 'parVec'
	parVec->writeBackBuff(buff);
      }
    }

    /**
     * Returns buffer pointer, which for 'step'<0 depends on the convention
     * used (see header comm., 'negFort').
     * If 'convFort' is given, it overrides 'negFort' (true -> Fortran conv.;
     * false -> ours).
     *
     * @param convFort S.a. Optional
     * @return         Buffer pointer
     */
    T* getBuff(bool convFort) const {
      if (step>0 || convFort) return buff;
      else return buff.p()+((1-n)*step);
    }

    T* getBuff() const {
      return getBuff(negFort);
    }

    const T& operator[](int pos) const {
      if (pos<0 || pos>=n) throw OutOfRangeException(EXCEPT_MSG(""));
      // Fortran convention for 'step'<0:
      return (step>0)?buff[pos*step]:buff[(pos-n+1)*step];
    }

    T& operator[](int pos) {
      if (pos<0 || pos>=n) throw OutOfRangeException(EXCEPT_MSG(""));
      return (step>0)?buff[pos*step]:buff[(pos-n+1)*step];
    }

    /**
     * Copies vector into given flat buffer 'arr' (must have place for
     * 'n' entries).
     *
     * @param arr Target buffer
     */
    void copy(T* arr) const {
      // ATTENTION: 'ArrayUtils' is using our convention, not Fortran's
      ArrayUtils<T>::copy(arr,getBuff(false),n,1,step);
    }
  };
//ENDNS

#endif
