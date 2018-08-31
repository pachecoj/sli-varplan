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
 * Desc.:  Definition of class IndexVector
 * ------------------------------------------------------------------- */

#include "lhotse/matrix/IndexVector.h"
#ifdef MATLAB_MEX
#include "lhotse/matif/MatlabMatrix.h"
#endif

// Local macro for args to 'checkTS'
#define TSARG __FILE__,__LINE__

//BEGINNS(matrix)
  // Public methods

#ifdef MATLAB_MEX
  IndexVector::IndexVector(const MatlabMatrix& mat,bool decr) :
    BaseVector<int>(mat)
  {
    setDefValues();
    if (decr)
      apply1(bind2nd(std::minus<int>(),1));
  }
#endif

  void IndexVector::complement(const BaseVector<int>& vec,int sz,
			       bool strict)
  {
    BaseVector<bool> inComp(sz,true);
    int i,j,lim=vec.size();
    if (sz<=0) throw InvalidParameterException("sz");

    checkTS(TSARG); vec.checkTS(TSARG);
    if (strict) {
      if (!vec.checkBounds(Interval<int>(0,sz,IntVal::ivClosed,
					 IntVal::ivOpen)))
	throw InvalidParameterException(EXCEPT_MSG("'vec' has entries out of range"));
    } else {
      if (!isValidVec(vec))
	throw InvalidParameterException(EXCEPT_MSG("'vec' has invalid entries"));
    }
    if (strict) {
      for (i=0; i<lim; i++) {
	j=vec[i];
	if (!inComp[j])
	  throw InvalidParameterException(EXCEPT_MSG("'vec' contains duplicate entries"));
	inComp[j]=false;
      }
    } else {
      for (i=0; i<lim; i++) inComp[vec[i]]=false;
    }
    ensureCapacity(sz-lim);
    for (i=j=0; i<sz; i++)
      if (inComp[i]) set(j++,i);
  }

  void IndexVector::invIndInsert(IndexVector& indJ,int elem,int pos)
  {
    if (elem<0 || elem>=indJ.size()) throw InvalidParameterException("elem");
    if (pos==-1) pos=n;
    else if (pos<0 || pos>n) throw InvalidParameterException("pos");

    checkTS(TSARG); indJ.checkTS(TSARG);
    insert(elem,1,pos); // insert 'elem' into I
    // Correct J and insert 'elem'
    for (int i=pos; i<n; i++) indJ.set(operator[](i),i);
  }

  void IndexVector::invIndRemove(IndexVector& indJ,int pos,bool presOrd)
  {
    if (pos<0 || pos>=n) throw InvalidParameterException("pos");
    checkTS(TSARG); indJ.checkTS(TSARG);
    if (!presOrd) {
      invIndExchange(indJ,pos,n-1);
      invIndRemove(indJ,n-1,true); // cheap
    } else {
      indJ.set(operator[](pos),-1);
      remove(pos);
      // Correct J
      for (int i=pos; i<n; i++)
	indJ.set(operator[](i),i);
    }
  }

  void IndexVector::histogram(BaseVector<int>& hist) const
  {
    int i;

    hist.fill(max()+1,0);
    for (i=0; i<n; i++) hist[operator[](i)]++;
  }
//ENDNS

#undef TSARG
