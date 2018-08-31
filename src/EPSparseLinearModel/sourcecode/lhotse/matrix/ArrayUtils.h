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
 * Desc.:  Header class ArrayUtils
 * ------------------------------------------------------------------- */

#ifndef ARRAYUTILS_H
#define ARRAYUTILS_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/matrix/default.h"
#include "lhotse/matrix/ArrayUtilsBasic.h"
#include <algorithm> // for sorting
#include <functional> // for sorting

//BEGINNS(matrix)
/* OLD CODE:
  /**
   * Helper function object, used by 'sort'. Implements '<' on 'pair<T,T2>',
   * by just looking at the first element.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   *
  template<class T,class T2> class ArrayUtils_CompareLess
  {
  public:
    bool operator() (const pair<T,T2>& a,const pair<T,T2>& b) {
      return (a.first<b.first);
    }
  };

  /**
   * Helper function object, used by 'sort'. Implements '>' on 'pair<T,T2>',
   * by just looking at the first element.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   *
  template<class T,class T2> class ArrayUtils_CompareGreater
  {
  public:
    bool operator() (const pair<T,T2>& a,const pair<T,T2>& b) {
      return (a.first>b.first);
    }
  };
*/

  /**
   * Some static inline methods doing basic non-arithmetic operations on
   * flat or indexed arrays.
   * Define freq. used methods which do not require arithmetical operations
   * here (i.e. hold for any T) to de-clutter the vector/matrix classes.
   * <p>
   * Negative step size convention:
   * ATTENTION: If step size parameter (say: step) is <0, we use our
   * internal convention:
   *   x_i -> x[i*step],
   * meaning that negative indexes are applied to x. This is DIFFERENT
   * from the Fortran convention typically used in 'BaseLinVec', and also
   * in 'FastUtils'. If xlin is a 'BaseLinVec' version of x, use
   * 'xlin.getBuff(false)' to obtain the buffer pointer corr. to our
   * convention used here.
   * NOTE: 'FastUtils' behaves differently, uses the Fortran convention!
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  template<class T> class ArrayUtils : public ArrayUtilsBasic<T>
  {
  public:
    // Public static methods

    /**
     * a = f(), f nullary function applied once for each element, from
     * start to end.
     *
     * @param a     Target vector a
     * @param n     Length
     * @param func  Function f
     * @param aStep Step size for a (ignored if 'aInd' given)
     * @param aInd  Position index for a. 0 -> none
     */
    template<class T2> static void applyNulFunc(T* a,int n,
						const NullaryFunc<T2>& func,
						int aStep=1,
						const int* aInd=0) {
      int i;
      if (aInd==0) {
	if (aStep==0) throw InvalidParameterException(EXCEPT_MSG("aStep"));
	for (i=0; i<n; i++,a+=aStep) *a=(T) func();
      } else {
	for (i=0; i<n; i++) a[*(aInd++)]=(T) func();
      }
    }

    /**
     * Applies accum. object 'acc' to each element of 'a', start to end.
     *
     * @param a     Source vector a
     * @param n     Length
     * @param acc   Accumulator object
     * @param aStep Step size for a (ignored if 'aInd' given)
     * @param aInd  Position index for a. 0 -> none
     */
    template<class T1,class T2> static
    void applyAcc(T* a,int n,const AccumulFunc<T1,T2>& acc,int aStep=1,
		  const int* aInd=0) {
      int i;
      if (aInd==0) {
	if (aStep==0) throw InvalidParameterException(EXCEPT_MSG("aStep"));
	for (i=0; i<n; i++,a+=aStep) acc((T1) *a);
      } else {
	for (i=0; i<n; i++) acc((T1) a[*(aInd++)]);
      }
    }

    /**
     * Applies accum. object 'acc' to each element of f(a), start to end,
     * where f=='func' (f: T -> T1, accum.: T1 -> T2).
     * <p>
     * NOTE: Could also be done by config. 'acc' with
     *   compose22(g,func)
     * instead of g, then calling 'applyAcc'.
     *
     * @param a     Source vector a
     * @param n     Length
     * @param func  Function f
     * @param acc   Accumulator object
     * @param aStep Step size for a (ignored if 'aInd' given)
     * @param aInd  Position index for a. 0 -> none
     */
    template<class UnOp,class T2,class T3> static void
    applyAccFunc(const T* a,int n,const UnOp& func,
		 const AccumulFunc<T3,T2>& acc,
		 int aStep=1,const int* aInd=0) {
      typedef typename UnOp::argument_type Arg;
      int i;
      if (aInd==0) {
	if (aStep==0) throw InvalidParameterException(EXCEPT_MSG("aStep"));
	for (i=0; i<n; i++,a+=aStep) acc((T3) func((Arg) *a));
      } else {
	for (i=0; i<n; i++) acc((T3) func((Arg) a[*(aInd++)]));
      }
    }

    /**
     * Applies accum. object 'acc' to each element of f(a,b), start to end,
     * where f=='func' (f: (T,T1) -> T2, accum.: T2 -> T3).
     *
     * @param a     Source vector a (element type T)
     * @param b     Source vector b
     * @param n     Length
     * @param func  Function f
     * @param acc   Accumulator object
     * @param aStep Step size for a (ignored if 'aInd' given)
     * @param bStep Step size for b (ignored if 'bInd' given)
     * @param aInd  Position index for a. 0 -> none
     * @param bInd  Position index for b. 0 -> none
     */
    template<class BinOp,class T2,class T3,class T4> static void
    applyAccBinFunc(const T* a,const T4* b,int n,const BinOp& func,
		    const AccumulFunc<T2,T3>& acc,int aStep=1,int bStep=1,
		    const int* aInd=0,const int* bInd=0) {
      typedef typename BinOp::first_argument_type Arg1;
      typedef typename BinOp::second_argument_type Arg2;
      int i;
      if (aInd==0) {
	if (aStep==0) throw InvalidParameterException(EXCEPT_MSG("aStep"));
	if (bInd==0) {
	  if (bStep==0) throw InvalidParameterException(EXCEPT_MSG("bStep"));
	  for (i=0; i<n; i++,a+=aStep,b+=bStep)
	    acc((T2) func((Arg1) *a,(Arg2) *b));
	} else {
	  for (i=0; i<n; i++,a+=aStep)
	    acc((T2) func((Arg1) *a,(Arg2) b[*(bInd++)]));
	}
      } else {
	if (bInd==0) {
	  if (bStep==0) throw InvalidParameterException(EXCEPT_MSG("bStep"));
	  for (i=0; i<n; i++,b+=bStep)
	    acc((T2) func((Arg1) a[*(aInd++)],(Arg2) *b));
	} else {
	  for (i=0; i<n; i++)
	    acc((T2) func((Arg1) a[*(aInd++)],(Arg2) b[*(bInd++)]));
	}
      }
    }

    /**
     * Fills vector a with element s.
     *
     * @param a     Vector a
     * @param s     Fill value
     * @param n     Length
     * @param aStep Step size for a (ignored if 'aInd' given)
     * @param aInd  Position index for a. 0 -> none
     */
    static void fill(T* a,T s,int n,int aStep=1,const int* aInd=0);

    /**
     * Copies vector b to vector a (must have enough space!). If both a,
     * b are flat (step size 1), 'copy' can be used with overlapping
     * buffers. Otherwise, copying is done from start to end for both.
     * NOTE: Use 'convert' if source/target type are different.
     *
     * @param a     Target vector a
     * @param b     Source vector b
     * @param n     Length
     * @param aStep Step size for a (ignored if 'aInd' given)
     * @param bStep Step size for b (ignored if 'bInd' given)
     * @param aInd  Position index for a. 0 -> none
     * @param bInd  Position index for b. 0 -> none
     */
    static void copy(T* a,const T* b,int n,int aStep=1,int bStep=1,
		     const int* aInd=0,const int* bInd=0);

    /**
     * Swaps contents of vectors a,b (must have same length), both proc.
     * start to end.
     *
     * @param a     Vector a
     * @param b     Vector b
     * @param n     Length
     * @param aStep Step size for a (ignored if 'aInd' given)
     * @param bStep Step size for b (ignored if 'bInd' given)
     * @param aInd  Position index for a. 0 -> none
     * @param bInd  Position index for b. 0 -> none
     */
    static void swap(T* a,T* b,int n,int aStep=1,int bStep=1,
		     const int* aInd=0,const int* bInd=0);

    /**
     * Converts vector b (elem. type T2) to elem. type T and writes it into
     * a. Uses static cast.
     *
     * @param a     Target vector a
     * @param b     Source vector b (elem. type T2)
     * @param n     Length
     * @param aStep Step size for a (ignored if 'aInd' given)
     * @param bStep Step size for b (ignored if 'bInd' given)
     * @param aInd  Position index for a. 0 -> none
     * @param bInd  Position index for b. 0 -> none
     */
    template<class T2> static void convert(T* a,const T2* b,int n,int aStep=1,
					   int bStep=1,const int* aInd=0,
					   const int* bInd=0) {
      int i;
      if (aInd==0) {
	if (aStep==0) throw InvalidParameterException(EXCEPT_MSG("aStep"));
	if (bInd==0) {
	  if (bStep==0) throw InvalidParameterException(EXCEPT_MSG("bStep"));
	  for (i=0; i<n; i++,a+=aStep,b+=bStep) *a=static_cast<T>(*b);
	} else {
	  for (i=0; i<n; i++,a+=aStep) *a=static_cast<T>(b[*(bInd++)]);
	}
      } else {
	if (bInd==0) {
	  if (bStep==0) throw InvalidParameterException(EXCEPT_MSG("bStep"));
	  for (i=0; i<n; i++,b+=bStep) a[*(aInd++)]=static_cast<T>(*b);
	} else {
	  for (i=0; i<n; i++) a[*(aInd++)]=static_cast<T>(b[*(bInd++)]);
	}
      }
    }

    /**
     * Sorts vector a (of length n), in ascending order if 'asc'==true,
     * in descending order otherwise.
     * NOTE: No step size !=1 or indexing allowed here.
     *
     * @param a     Vector a
     * @param n     Length
     * @param asc   In ascending order? Def.: true
     */
    static void sort(T* a,int n,bool asc=true) {
      if (asc)
	std::sort(a,a+n,std::less<T>());
      else
	std::sort(a,a+n,std::greater<T>());
    }

    /**
     * Sorts vector a (of length n), in ascending order if asc==true,
     * in descending order otherwise. In addition, a T2 index has to be
     * passed in ind, of the same length n. We apply the same permutation
     * to a and ind.
     * <p>
     * NOTE: A flat copy is created in any case, pairing elements of a and
     * ind, then STL's 'sort' is used.
     *
     * @param a     Vector a
     * @param ind   Index ind
     * @param n     Length
     * @param aStep Step size for a (ignored if 'aInd' given)
     * @param iStep Step size for ind (ignored if 'iInd' given)
     * @param aInd  Position index for a. 0 -> none
     * @param iInd  Position index for ind. 0 -> none
     * @param asc   In ascending order? Def.: true
     */
    template<class T2> static void sortInd(T* a,T2* ind,int n,int aStep=1,
					   int iStep=1,const int* aInd=0,
					   const int* iInd=0,bool asc=true) {
      // Draw flat copy (paired elements)
      ArrayHandle<pair<T,T2> > flatCp(n);
      int i;
      T* aB;
      T2* iB;
      if (aInd==0) {
	if (aStep==0) throw InvalidParameterException(EXCEPT_MSG("aStep"));
	for (i=0,aB=a; i<n; i++,aB+=aStep) flatCp[i].first=*aB;
      } else {
	for (i=0; i<n; i++) flatCp[i].first=a[*(aInd++)];
	aInd-=n;
      }
      if (iInd==0) {
	if (iStep==0) throw InvalidParameterException(EXCEPT_MSG("iStep"));
	for (i=0,iB=ind; i<n; i++,iB+=iStep) flatCp[i].second=*iB;
      } else {
	for (i=0; i<n; i++) flatCp[i].second=ind[*(iInd++)];
	iInd-=n;
      }
      // Sort copy using STL
      if (asc)
        std::sort(flatCp.p(),flatCp.p()+n,pair1st<T2>(less<T>()));
	//std::sort(flatCp.p(),flatCp.p()+n,ArrayUtils_CompareLess<T,T2>());
      else
        std::sort(flatCp.p(),flatCp.p()+n,pair1st<T2>(greater<T>()));
	//std::sort(flatCp.p(),flatCp.p()+n,ArrayUtils_CompareGreater<T,T2>());
      // Copy back
      if (aInd==0) {
	for (i=0,aB=a; i<n; i++,aB+=aStep) *aB=flatCp[i].first;
      } else {
	for (i=0; i<n; i++) a[*(aInd++)]=flatCp[i].first;
      }
      if (iInd==0) {
	for (i=0,iB=ind; i<n; i++,iB+=iStep) *iB=flatCp[i].second;
      } else {
	for (i=0; i<n; i++) ind[*(iInd++)]=flatCp[i].second;
      }
    }

    /**
     * Selection and accumulation.
     * 'imap' has 2*k-1 int entries, ie. k consec. tupels i(s),j(s), where
     * j(k-1) is not used. Works by running a counter pos(s), with pos(0)=0.
     * For s=0,...,k-1:
     *   a[pos(s)] += b[i(s)],
     *   pos(s+1) = pos(s) + j(s),
     * 'a' must have suff. size, and be init. properly. i(s) must lie within
     * 'b'.
     *
     * @param a     Target vector a
     * @param b     Vector b
     * @param imap  Int vector (must be flat)
     * @param k     Number of tupels in 'imap'
     * @param aStep Step size for a
     * @param bStep Step size for b
     * @param aInd  Position index for a. 0 -> none
     * @param bInd  Position index for b. 0 -> none
     */
    static void accuMap(T* a,const T* b,const int* imap,int k,int aStep=1,
			int bStep=1,const int* aInd=0,const int* bInd=0);
  };

  template<class T> inline
  void ArrayUtils<T>::fill(T* a,T s,int n,int aStep,const int* aInd)
  {
    int i;
    if (aInd==0) {
      if (aStep==0) throw InvalidParameterException(EXCEPT_MSG("aStep"));
      for (i=0; i<n; i++,a+=aStep) *a=s;
    } else {
      for (i=0; i<n; i++) a[*(aInd++)]=s;
    }
  }

  template<class T> inline
  void ArrayUtils<T>::copy(T* a,const T* b,int n,int aStep,int bStep,
			   const int* aInd,const int* bInd)
  {
    int i;
    if (aInd==0) {
      if (aStep==0) throw InvalidParameterException(EXCEPT_MSG("aStep"));
      if (bInd==0) {
	if (bStep==0) throw InvalidParameterException(EXCEPT_MSG("bStep"));
	if (aStep==1 && bStep==1 && a>b) {
	  // From right to left (reverse)
	  a+=(n-1); b+=(n-1); aStep=bStep=-1;
	}
	for (i=0; i<n; i++,a+=aStep,b+=bStep) *a=(*b);
      } else {
	for (i=0; i<n; i++,a+=aStep) *a=b[*(bInd++)];
      }
    } else {
      if (bInd==0) {
	if (bStep==0) throw InvalidParameterException(EXCEPT_MSG("bStep"));
	for (i=0; i<n; i++,b+=bStep) a[*(aInd++)]=*b;
      } else {
	for (i=0; i<n; i++) a[*(aInd++)]=b[*(bInd++)];
      }
    }
  }

  template<class T> inline
  void ArrayUtils<T>::swap(T* a,T* b,int n,int aStep, int bStep,
			   const int* aInd,const int* bInd)
  {
    int i;
    T temp;

    if (aInd==0) {
      if (aStep==0) throw InvalidParameterException(EXCEPT_MSG("aStep"));
      if (bInd==0) {
	if (bStep==0) throw InvalidParameterException(EXCEPT_MSG("bStep"));
	for (i=0; i<n; i++,a+=aStep,b+=bStep) {
	  temp=*a; *a=*b; *b=temp;
	}
      } else {
	T* bP;
	for (i=0; i<n; i++,a+=aStep) {
	  bP=b+*(bInd++);
	  temp=*a; *a=*bP; *bP=temp;
	}
      }
    } else {
      if (bInd==0) {
	if (bStep==0) throw InvalidParameterException(EXCEPT_MSG("bStep"));
	T* aP;
	for (i=0; i<n; i++,b+=bStep) {
	  aP=a+*(aInd++);
	  temp=*aP; *aP=*b; *b=temp;
	}
      } else {
	T* aP,*bP;
	for (i=0; i<n; i++) {
	  aP=a+*(aInd++); bP=b+(*bInd++);
	  temp=*aP; *aP=*bP; *bP=temp;
	}
      }
    }
  }

  template<class T> inline
  void ArrayUtils<T>::accuMap(T* a,const T* b,const int* imap,int k,int aStep,
			      int bStep,const int* aInd,const int* bInd)
  {
    int s;

    if (aInd==0) {
      if (aStep==0) throw InvalidParameterException(EXCEPT_MSG("aStep"));
      if (bInd==0) {
	if (bStep==0) throw InvalidParameterException(EXCEPT_MSG("bStep"));
	for (s=0; s<k-1; s++) {
	  (*a)+=*(b+(*(imap++))*bStep);
	  a+=(aStep*(*(imap++)));
	}
	(*a)+=*(b+(*imap)*bStep);
      } else {
	for (s=0; s<k-1; s++) {
	  (*a)+=b[bInd[*(imap++)]];
	  a+=(aStep*(*(imap++)));
	}
	(*a)+=b[bInd[*imap]];
      }
    } else {
      if (bInd==0) {
	if (bStep==0) throw InvalidParameterException(EXCEPT_MSG("bStep"));
	for (s=0; s<k-1; s++) {
	  a[*aInd]+=*(b+(*(imap++))*bStep);
	  aInd+=(*(imap++));
	}
	a[*aInd]+=*(b+(*imap)*bStep);
      } else {
	for (s=0; s<k-1; s++) {
	  a[*aInd]+=b[bInd[*(imap++)]];
	  aInd+=(*(imap++));
	}
	a[*aInd]+=b[bInd[*imap]];
      }
    }
  }

//ENDNS

#endif
