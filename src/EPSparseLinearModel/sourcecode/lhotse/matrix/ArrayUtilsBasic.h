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
 * Desc.:  Header class ArrayUtilsBasic
 * ------------------------------------------------------------------- */

#ifndef ARRAYUTILSBASIC_H
#define ARRAYUTILSBASIC_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/matrix/default.h"

//BEGINNS(matrix)
  /**
   * This class contains some elementrary inline methods of 'ArrayUtils',
   * which are used by low-level classes such as 'ArrayHandle'. We need
   * this class to break a mutual dependency between 'ArrayUtils' and
   * low-level classes in the global module.
   * <p>
   * This class does not rely on any of the low-level global classes,
   * except for 'StandardException' and the global exceptions. These must
   * be included in 'global.h' (global module) BEFORE this class.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  template<class T> class ArrayUtilsBasic
  {
  public:
    // Public static methods

    /**
     * a = f(b), b vector, f unary function applied to each element, start
     * to end. f : T1 -> T.
     *
     * @param a     Target vector a (element type T)
     * @param b     Vector b
     * @param n     Length
     * @param func  Function f (STL object)
     * @param aStep Step size for a (ignored if 'aInd' given)
     * @param bStep Step size for b (ignored if 'bInd' given)
     * @param aInd  Position index for a. 0 -> none
     * @param bInd  Position index for b. 0 -> none
     */
    template<class UnOp,class T2> static void
    applyFunc(T* a,const T2* b,int n,const UnOp& func,int aStep=1,int bStep=1,
	      const int* aInd=0,const int* bInd=0) {
      typedef typename UnOp::argument_type Arg;
      int i;
      if (aInd==0) {
	if (aStep==0) throw InvalidParameterException(EXCEPT_MSG("aStep"));
	if (bInd==0) {
	  if (bStep==0) throw InvalidParameterException(EXCEPT_MSG("bStep"));
	  for (i=0; i<n; i++,a+=aStep,b+=bStep) *a=(T) func((Arg) *b);
	} else {
	  for (i=0; i<n; i++,a+=aStep) *a=(T) func((Arg) b[*(bInd++)]);
	}
      } else {
	if (bInd==0) {
	  if (bStep==0) throw InvalidParameterException(EXCEPT_MSG("bStep"));
	  for (i=0; i<n; i++,b+=bStep) a[*(aInd++)]=(T) func((Arg) *b);
	} else {
	  for (i=0; i<n; i++) a[*(aInd++)]=(T) func((Arg) b[*(bInd++)]);
	}
      }
    }

    /**
     * a = f(b,c), b,c vectors, f binary function applied to each element
     * pair (b,c must have same length 'n'), start to end. f (T1,T2) -> T.
     *
     * @param a     Target vector a (element type T)
     * @param b     Vector b
     * @param c     Vector c
     * @param n     Length
     * @param func  Function object for f
     * @param aStep Step size for a (ignored if 'aInd' given)
     * @param bStep Step size for b (ignored if 'bInd' given)
     * @param cStep Step size for c (ignored if 'cInd' given)
     * @param aInd  Position index for a. 0 -> none
     * @param bInd  Position index for b. 0 -> none
     * @param cInd  Position index for c. 0 -> none
     */
    template<class BinOp,class T1,class T2> static void
    applyBinFunc(T* a,const T1* b,const T2* c,int n,const BinOp& func,
		 int aStep=1,int bStep=1,int cStep=1,const int* aInd=0,
		 const int* bInd=0,const int* cInd=0) {
      typedef typename BinOp::first_argument_type Arg1;
      typedef typename BinOp::second_argument_type Arg2;
      int i;
      if (aInd==0) {
	if (aStep==0) throw InvalidParameterException(EXCEPT_MSG("aStep"));
	if (bInd==0) {
	  if (bStep==0) throw InvalidParameterException(EXCEPT_MSG("bStep"));
	  if (cInd==0) {
	    if (cStep==0) throw InvalidParameterException(EXCEPT_MSG("cStep"));
	    for (i=0; i<n; i++,a+=aStep,b+=bStep,c+=cStep)
	      *a=(T) func((Arg1) *b,(Arg2) *c);
	  } else {
	    for (i=0; i<n; i++,a+=aStep,b+=bStep)
	      *a=(T) func((Arg1) *b,(Arg2) c[*(cInd++)]);
	  }
	} else {
	  if (cInd==0) {
	    if (cStep==0) throw InvalidParameterException(EXCEPT_MSG("cStep"));
	    for (i=0; i<n; i++,a+=aStep,c+=cStep)
	      *a=(T) func((Arg1) b[*(bInd++)],(Arg2) *c);
	  } else {
	    for (i=0; i<n; i++,a+=aStep)
	      *a=(T) func((Arg1) b[*(bInd++)],(Arg2) c[*(cInd++)]);
	  }
	}
      } else {
	if (bInd==0) {
	  if (bStep==0) throw InvalidParameterException(EXCEPT_MSG("bStep"));
	  if (cInd==0) {
	    if (cStep==0) throw InvalidParameterException(EXCEPT_MSG("cStep"));
	    for (i=0; i<n; i++,b+=bStep,c+=cStep)
	      a[*(aInd++)]=(T) func((Arg1) *b,(Arg2) *c);
	  } else {
	    for (i=0; i<n; i++,b+=bStep)
	      a[*(aInd++)]=(T) func((Arg1) *b,(Arg2) c[*(cInd++)]);
	  }
	} else {
	  if (cInd==0) {
	    if (cStep==0) throw InvalidParameterException(EXCEPT_MSG("cStep"));
	    for (i=0; i<n; i++,c+=cStep)
	      a[*(aInd++)]=(T) func((Arg1) b[*(bInd++)],(Arg2) *c);
	  } else {
	    for (i=0; i<n; i++)
	      a[*(aInd++)]=(T) func((Arg1) b[*(bInd++)],(Arg2) c[*(cInd++)]);
	  }
	}
      }
    }
  };

//ENDNS

#endif
