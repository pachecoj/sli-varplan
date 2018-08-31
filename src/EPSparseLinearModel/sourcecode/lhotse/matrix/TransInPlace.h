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
 * Desc.:  Header class TransInPlace
 * ------------------------------------------------------------------- */

#ifndef TRANSINPLACE_H
#define TRANSINPLACE_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/matrix/default.h"
#include "lhotse/matrix/ArrayUtils.h"

//BEGINNS(matrix)
  /**
   * Static method for transition of a rectangular T matrix in situ. The
   * matrix is stored in a flat T* buffer, where we assume column-major
   * ordering. The method can be use for row-major ordering just as well,
   * by switching row and column dimension.
   * <p>
   * This is just a translation of Algorithm 467 (CACM):
   *   Brenner: Matrix Transposition in Place, CACM 16(11), 1973
   *   Available: portal.acm.org
   * <p>
   * Translation from Fortran is made hard by several things:
   * - array indexing starts with 1, not 0!
   * - lots of goto instead of proper loops
   * - if <expr> <p1>,<p2>,<p3>
   *   ==> <expr> <  0: goto <p1>
   *       <expr> == 0: goto <p2>
   *       <expr> >  0: goto <p3>
   * - do <p1> <var>=<start>,<end>,<step>
   *   ...
   *   <p1>: continue
   *   ==> like for loop from <start> to <end> in steps of <step>,
   *       executing stuff in ...
   * <p>
   * The rough idea is to perform the permutation onto the flat buffer
   * one subcycle after the other (each permutation decomposes into a
   * unique set of subcycles). The game is how to ensure that each
   * subcycle is performed once and only once. A simple but inefficient
   * way is to count up from small numbers and then search for the
   * smallest position in each subcycle, always starting from such
   * minimal one. A number can be skipped as starting position if a
   * smaller pos. in its subcycle is detected. The algorithm here is
   * faster and uses some clever number theory.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  template<class T> class TransInPlace
  {
  public:
    // Public static methods

    /**
     * Transposes matrix (or tensor) A in place. A is given via the flat array
     * 'a', in column major ordering and has size 'nn'-by-'n1'-by-'n2'. The
     * transposition flips the last two dimensions. In order to transpose
     * A in row-major ordering, just flip n1, n2 before calling.
     * NOTE: A must fill 'A' without gaps!
     * <p>
     * The algorithm uses a temp. boolean vector of size (n1+n2)/2, but no
     * other temp. storage. If 'nn'>1, we also req. O(nn) temp. storage.
     *
     * @param a  A in column-major format
     * @param n1 Number rows of A (or cols if row-major)
     * @param n2 Number cols of A (or rows if row-major)
     * @param nn S.a. Def.: 1
     */
    static void trans(T* a,int n1,int n2,int nn=1);

  protected:
    // Internal static methods

    /**
     * Factors 'n' into prime factors: 'ifact' contains the factors, 'ipower'
     * the factors up their multiplicities, 'nexp' their multiplicities.
     * NOTE: Works only for up to 8 distinct factors! Arrays must be of size
     * >= 8.
     *
     * @param n      S.a.
     * @param ifact  S.a.
     * @param ipower S.a.
     * @param nexp   S.a.
     * @return       Number of dist. prime factors
     */
    static int factor(int n,int* ifact,int* ipower,int* nexp);
  };

  /*
   * Protection against overflows:
   * The permutation is given by k -> k' = k*n (mod m), but k*n can be as
   * large as (m-1)*n and might overflow. We can also write
   *   k = i n2 + j ==> k' = j n + i,
   * so k' can be computed safely with 1 division, 2 multiplications.
   */
  template<class T> inline
  void TransInPlace<T>::trans(T* a,int n1,int n2,int nn)
  {
    int i,j,ip;
    int n=n1,n12=n1*n2;
    int m=n12-1;
    T temp;

    if (n1<2 || n2<2) return; // nothing to do
    if (n1==n2) {
      // Square matrix: simple special case
      T* a1P=a,*a2P;
      int num,off;
      for (i=1,off=n; i<=n; i++,off+=(n+1)) {
	a1P+=nn*i; a2P=a+nn*off;
	num=n-i;
	if (nn==1) {
	  for (j=0; j<num; j++) {
	    temp=*a1P;
	    *(a1P++)=*a2P;
	    *a2P=temp; a2P+=n;
	  }
	} else {
	  for (j=0; j<num; j++) {
	    ArrayUtils<T>::swap(a1P,a2P,nn);
	    a1P+=nn; a2P+=nn*n;
	  }
	}
      }
    } else {
      int nwork=(n1+n2)/2;
      ArrayHandle<bool> moved(nwork); // temp. array
      int ifact[8],ipower[8],nexp[8],iexp[8]; // for prime fact.
      // Modulus m factored into prime powers. Eight factors sufficient
      // 30
      int npower=factor(m,ifact,ipower,nexp);
      for (ip=0; ip<npower; ip++) iexp[ip]=0;
      // Generate every divisor of m less than m/2
      int idiv=1;
      int ncount,istart,mmist,isoid,itest,itemp;
      int ia1,ia2,mmia1,mmia2;
      bool cont;
      T atemp,btemp;
      ArrayHandle<T> atarr,btarr;
      if (nn>1) {
	atarr.changeRep(nn); btarr.changeRep(nn);
      }
      while (idiv<m/2) { // 50
	// The number of elements whose index is div. by idiv and by no other
	// divisor of m is the Euler totient function, phi(m/idiv)
	ncount=m/idiv;
	for (ip=0; ip<npower; ip++)
	  if (iexp[ip]!=nexp[ip])
	    ncount=(ncount/ifact[ip])*(ifact[ip]-1);
	for (i=0; i<nwork; i++) moved[i]=false;
	// The starting point of a subcycle is div. only by idiv and must
	// not appear in any other subcycle
	for (istart=idiv; ncount>0; istart+=idiv) { // 160
	  // 80
	  mmist=m-istart;
	  if (istart!=idiv) {
	    if (istart<=nwork && moved[istart-1]) continue;
	    // 90
	    isoid=istart/idiv;
	    for (ip=0,cont=false; ip<npower; ip++)
	      if (iexp[ip]!=nexp[ip] && (isoid%ifact[ip])==0) {
		cont=true; break;
	      }
	    if (cont) continue;
	    if (istart>nwork) {
	      itest=istart; cont=false;
	      do {
		// 110
		itemp=itest/n2;
		itest=itemp+(itest-itemp*n2)*n; // succ. in permut.
		if (itest<istart || itest>mmist) {
		  cont=true; break;
		}
	      } while (itest>istart && itest<mmist);
	      if (cont) continue;
	    }
	  }
	  // 120
	  if (nn==1) {
	    atemp=a[istart]; btemp=a[mmist];
	  } else {
	    ArrayUtils<T>::copy(atarr,a+nn*istart,nn);
	    ArrayUtils<T>::copy(btarr,a+nn*mmist,nn);
	  }
	  for (ia1=istart; ; ) {
	    // 130
	    itemp=ia1/n2;
	    ia2=itemp+(ia1-itemp*n2)*n; // succ. in permut.
	    mmia1=m-ia1; mmia2=m-ia2;
	    if (ia1<=nwork) moved[ia1-1]=true;
	    if (mmia1<=nwork) moved[mmia1-1]=true;
	    ncount-=2;
	    // Move two elements, the second from the negative subcycle.
	    // Check first for subcycle closure
	    if (nn==1) {
	      if (ia2==istart) {
		// 140
		a[ia1]=atemp; a[mmia1]=btemp;
		break;
	      } else if (mmia2==istart) {
		// 150
		a[ia1]=btemp; a[mmia1]=atemp;
		break;
	      } else {
		a[ia1]=a[ia2];
		a[mmia1]=a[mmia2];
		ia1=ia2;
	      }
	    } else {
	      if (ia2==istart) {
		// 140
		ArrayUtils<T>::copy(a+nn*ia1,atarr,nn);
		ArrayUtils<T>::copy(a+nn*mmia1,btarr,nn);
		break;
	      } else if (mmia2==istart) {
		// 150
		ArrayUtils<T>::copy(a+nn*ia1,btarr,nn);
		ArrayUtils<T>::copy(a+nn*mmia1,atarr,nn);
		break;
	      } else {
		ArrayUtils<T>::copy(a+nn*ia1,a+nn*ia2,nn);
		ArrayUtils<T>::copy(a+nn*mmia1,a+nn*mmia2,nn);
		ia1=ia2;
	      }
	    }
	  }
	}
	for (ip=0,cont=false; ip<npower; ip++) {
	  if (iexp[ip]!=nexp[ip]) {
	    iexp[ip]++;
	    idiv*=ifact[ip];
	    cont=true; break;
	  } else {
	    // 170
	    iexp[ip]=0;
	    idiv/=ipower[ip];
	  }
	}
	if (!cont) break; // leave loop
      }
    }
  }

  /*
   * Very simple (does not even use a sieve), but sufficient for our
   * small numbers (will not break RSA :-))
   */
  template<class T> inline
  int TransInPlace<T>::factor(int n,int* ifact,int* ipower,int* nexp)
  {
    int ip=-1,ifcur=0,npart=n,idiv=2;
    int iquot;

    for (;;) {
      // 10
      iquot=npart/idiv;
      if ((npart-iquot*idiv)==0) {
	// 20
	if (idiv>ifcur) {
	  ip++;
	  ifact[ip]=ipower[ip]=idiv; nexp[ip]=1;
	  ifcur=idiv;
	} else {
	  ipower[ip]*=idiv; nexp[ip]++;
	}
	npart=iquot;
      } else {
	// 60
	if (iquot<=idiv) break; // leave loop -> 100
	if (idiv==2)
	  idiv=3;
	else
	  idiv+=2;
      }
    }
    // 100
    if (npart>1) {
      if (npart>ifcur) {
	// 120
	ip++;
	ifact[ip]=ipower[ip]=npart; nexp[ip]=1;
      } else {
	// 130
	ipower[ip]*=npart; nexp[ip]++;
      }
    }
    // 140
    return (ip+1);
  }
//ENDNS

#endif
