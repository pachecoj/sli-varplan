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
 * Module: GLOBAL
 * Desc.:  Definition of class Range
 * ------------------------------------------------------------------- */

#include <algorithm>

#include "lhotse/global.h"
#include "lhotse/Range.h"
#include "lhotse/matrix/BaseVector.h"

// Static members

const int Range::statFlat;
const int Range::statLinear;
const int Range::statIndex;
Range RangeFull::defR;

// Public methods

Range::Range(int pstart,int pend,int pstep) : start(pstart),step(pstep)
{
  if (pstart<0) throw InvalidParameterException(EXCEPT_MSG(""));
  if (pstep==1) {
    status=statFlat;
    if (pend!=-1) {
      if (pend<pstart)
	throw InvalidParameterException(EXCEPT_MSG("pend"));
      sz=pend-pstart+1;
    } else sz=-1; // open range
  } else if (pstep==0)
    throw InvalidParameterException(EXCEPT_MSG("pstep"));
  else if (pend<0 || (pend-pstart)%pstep!=0)
    throw InvalidParameterException(EXCEPT_MSG("Invalid linear range"));
  else {
    sz=(pend-pstart)/pstep+1;
    if (sz<=0)
      throw InvalidParameterException(EXCEPT_MSG("Invalid linear range"));
    status=statLinear;
  }
}

Range::Range(const ArrayHandle<int>& pindex) : status(statIndex),
					       sz(pindex.size())
{
  if (sz==0)
    throw InvalidParameterException(EXCEPT_MSG("'pindex' must not be empty"));
  if (DefIVal<int>::nonneg().check(pindex)!=0)
    throw InvalidParameterException(EXCEPT_MSG("'pindex' must be non-neg."));
  index=pindex; // copy handle
}

Range::Range(const BaseVector<int>& pindex) : status(statIndex),
					      sz(pindex.size())
{
  if (sz==0)
    throw InvalidParameterException(EXCEPT_MSG("'pindex' must not be empty"));
  if (!pindex.checkBounds(DefIVal<int>::nonneg()))
    throw InvalidParameterException(EXCEPT_MSG("'pindex' must be non-neg."));
  index=pindex.getFlatBuff();
}

int Range::getMaxPos(int n) const {
  if (status!=statIndex)
    return (step<0)?start:getEnd(n);
  else
    return BaseVector<int>::mask(index)->max();
}

bool Range::checkRange(int n) const {
  if (isFullRange(n)) return false;
  else if (status!=statIndex)
    return (start>=n || (sz!=-1 && start+(sz-1)*step>=n));
  else
    // already know that 'index' is non-neg.
    return (Interval<int>(0,n,IntVal::ivInf,IntVal::ivOpen).check(index)!=0);
}

bool Range::isUniqueMap() const
{
  if (status!=statIndex)
    return true;
  else
    return BaseVector<int>::mask(index)->areElemUnique();
}

void Range::mapIndex(const int* src,int n,int* trg,bool doCheck) const
{
  int i;

  if (doCheck && checkRange(n))
    throw OutOfRangeException(EXCEPT_MSG("Range violation"));
  if (status==statFlat)
    memmove(trg,src+start,size(n)*sizeof(int));
  else if (status==statLinear) {
    const int* sB=src+start;
    for (i=0; i<sz; i++,sB+=step) *(trg++)=*sB;
  } else {
    const int* iB=index;
    for (i=0; i<sz; i++)
      *(trg++)=src[*(iB++)];
  }
}

void Range::mapIndex(const ArrayHandle<int>& src,ArrayHandle<int>& trg) const
{
  bool doExt=(size(src.size())>trg.size());
  ArrayHandle<int> newTrg;

  if (!doExt) newTrg=trg;
  else
    newTrg.changeRep(size(src.size()));
  mapIndex(src.p(),src.size(),(int*) newTrg.p());
  trg=newTrg;
}

Range Range::subrange(const Range& rng) const
{
  int k;

  if (!isOpen()) {
    if (rng.checkRange(sz))
      throw InvalidParameterException(EXCEPT_MSG("'rng' violates range"));
    if (rng.isFullRange(sz))
      return Range(*this); // just a copy of this range
  } else if (rng.isOpen()) {
    // new one open too
    return Range(start+rng.getStart());
  }
  if (isFlatRange() && start==0)
    return Range(rng); // just a copy of 'rng'
  if (isLiteralRange()) {
    // this is literal
    if (rng.isLiteralRange()) {
      // new one literal as well
      // NOTE: Not both can be open, so rng.getEnd(sz-1) works
      return Range(start+step*rng.getStart(),
		   start+step*rng.getEnd(sz-1),step*rng.getStep());
    } else {
      // 'rng' is indexed. Need new index
      int rngSize=rng.size(0); // 'rng' not open
      ArrayHandle<int> newInd(rngSize);
      const int* rngInd=rng.getIndex();
      for (k=0; k<rngSize; k++)
	newInd[k]=start+step*rngInd[k];
      return Range(newInd);
    }
  } else {
    // this is indexed, but 'rng' is not the full range
    int rngSize=rng.size(sz);
    ArrayHandle<int> newInd(rngSize);
    rng.mapIndex(index.p(),sz,newInd.p(),false);
    return Range(newInd);
  }
}

Range Range::translate(int off) const
{
  if (off==0) return *this;
  switch (status) {
  case statFlat:
    if (start+off<0 || (sz!=-1 && start+sz-1+off<0))
      throw InvalidParameterException(EXCEPT_MSG(""));
    return Range(start+off,(sz==-1)?(-1):(start+off+sz-1));
  case statLinear:
    if (start+off<0 || start+off+(sz-1)*step<0)
      throw InvalidParameterException(EXCEPT_MSG(""));
    return Range(start+off,start+off+(sz-1)*step,step);
  default:
    BaseVector<int> msk;
    msk.reassign(index);
    if (!msk.checkBounds(Interval<int>(-off,0,IntVal::ivClosed)))
      throw InvalidParameterException(EXCEPT_MSG(""));
    ArrayHandle<int> tind(index.size());
    BaseVector<int>::mask(tind)->apply1(msk,bind2nd(std::plus<int>(),off));
    return Range(tind);
  }
}

// Global functions

const Range& full()
{
  return RangeFull::get();
}
