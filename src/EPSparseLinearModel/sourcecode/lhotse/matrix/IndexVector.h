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
 * Desc.:  Header class IndexVector
 * ------------------------------------------------------------------- */

#ifndef INDEXVECTOR_H
#define INDEXVECTOR_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/matrix/default.h"
#include "lhotse/matrix/BaseVector.h"
#include "lhotse/matrix/TempIndexVector.h"
#ifdef MATLAB_MEX
#include "lhotse/matif/predecl.h"
#endif

// Local macro for args to 'checkTS'
#define TSARG __FILE__,__LINE__

//BEGINNS(matrix)
  /**
   * Non-template specialization of 'BaseVector<int>'. Provides services
   * which are specific to element type 'int', most of them rel. to the
   * use of this object as index vector.
   * <p>
   * NOTE: This class does not have an own file format, but uses the one of
   * 'BaseVector<int>'.
   * <p>
   * NOTE: This class introduces element restrictions. All elements have to
   * be non-negative!
   * <p>
   * Index/inverse index services:
   * In many situations, an index I mapping {0,...,k-1} into {0,...,n-1}
   * (injective) is managed together with its inverse J mapping
   * {0,...,n-1} into {-1,0,...,k-1}:
   *   J(I(i)) = i  for 0<=i<k
   *   J(j)    = -1 if j\not\in I
   * This class offers corr. service methods 'invIndXXX'.
   * NOTE: In general, 'IndexVector' objects can contain duplicate entries.
   * <p>
   * NOTE: One of the main applications of 'IndexVector' is to build up and
   * maintain a selection index which can be used to construct a 'Range'
   * object. However, arithmetic progressions should not be constr. as 'Range'
   * from 'IndexVector', but directly as linear range.
   * <p>
   * Persistance:
   * The file format is the same as of the superclass BaseVector<int>.
   * ==> If the superclass is changed to BaseVector<int>, need to make sure
   *     'load' can load old BaseVector<int> formats!
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class IndexVector : public BaseVector<int>
  {
  public:
    // Constructors

    /**
     * Constructs empty vector
     */
    IndexVector() : BaseVector<int>() {
      setDefValues();
    }

    /**
     * Constructor. The vector is initialized with a.
     *
     * @param num Number of entries
     * @param a   Element to fill. Def.: def. fill value
     */
    IndexVector(int num,int a) : BaseVector<int>() {
      setDefValues();
      init(num,a);
    }

    explicit IndexVector(int num) : BaseVector<int>(num) {
      setDefValues();
    }

    /**
     * Copy constructor.
     *
     * @param vec Vector to copy
     */
    IndexVector(const BaseVector<int>& vec) : BaseVector<int>() {
      copyDefValues(vec); // def. value members
      assignInt(vec);
    }

    IndexVector(const BaseVecWrapper<int>& arg) : BaseVector<int>() {
      // use 'arg' fields to init. ours, then dispose of 'arg'
      convertWrapped((BaseVecWrapper<int>&) arg);
    }

#ifdef MATLAB_MEX
    /**
     * Special copy constructor from 'MatlabMatrix' object. This uses
     * static cast to convert from double to int. Matrix read in column
     * major order.
     * If 'decr'==true (def.: false), each element is decrem. by 1
     * to convert between Matlab's base 1 and C's base 0 for indexes.
     *
     * @param mat Matlab matrix
     */
    IndexVector(const MatlabMatrix& mat,bool decr=false);
#endif

    static TempIndexVector mask(int* pbuff,int pn,
				const Range& rng=RangeFull::get(),
				MemWatchBase* watch=0) {
      IndexVector* vec=new IndexVector();
      vec->reassign(pbuff,pn,rng,watch);
      return TempIndexVector(vec,vec);
    }

    static TempIndexVector mask(const ArrayHandle<int>& parr,
				const Range& rng=RangeFull::get()) {
      IndexVector* vec=new IndexVector();
      vec->reassign(parr.p(),parr.size(),rng,parr.getMemWatch());
      return TempIndexVector(vec,vec);
    }

    static TempIndexVector mask(const BaseVector<int>& vec,
				const Range& rng=RangeFull::get()) {
      IndexVector* mvec=new IndexVector();
      mvec->reassign(vec,rng);
      return TempIndexVector(mvec,mvec);
    }

    const TempIndexVector operator()(const Range& rng=RangeFull::get()) const {
      IndexVector* vec=DYNCAST(IndexVector,subrvalInt(rng));
      return TempIndexVector(vec,vec);
    }

    const int& operator()(int rng) const {
      return operator[](rng);
    }

    TempIndexVector operator()(const Range& rng=RangeFull::get()) {
      IndexVector* vec=DYNCAST(IndexVector,sublvalInt(rng));
      return TempIndexVector(vec,vec);
    }

    int operator()(int rng) {
      if (doCheckValid())
	throw TypeNotSuppException(EXCEPT_MSG("Not supported for this type. Use 'set'"));
      if (rng>=n)
	expand(rng+1);
      return operator[](rng);
    }

    /*
     * ATTENTION:
     * 'assignVirtual' must do the same as 'operator='. If you change
     * 'operator=' below, change them as well!
     */

    IndexVector& operator=(const BaseVector<int>& vec) {
      assignInt(vec);
      return *this;
    }

    IndexVector& operator=(int elem) {
      assignInt(elem);
      return *this;
    }

    void setDefValues() {
      BaseVector<int>::setDefValues();
      setDefFillValue(0); // Make sure 'defFill' is correct
    }

    /**
     * Constructs arithmetic progression a + k*s, k=0,1,...,K, where K is the
     * smallest k>=0 s.t. s*(a + k*s - b) >= 0. Here, a=='start', b=='end' and
     * s=='step'. s must not be 0 and is 1 by default.
     * <p>
     * NOTE: Do not use arith. progr. index vectors to build a 'Range'
     * (inefficient). ==> Method somewhat outdated!
     *
     * @param start See above
     * @param end   "
     * @param step  ". Optional. Def.: 1
     */
    void fillRange(int start,int end,int step=1) {
      register int k,limit,act;

      checkTS(TSARG);
      if (step==0) throw InvalidParameterException("step");
      limit=std::max(0,(int) ceil((end-start)/step));
      if (start<0 || start+limit*step<0)
	throw InvalidParameterException(EXCEPT_MSG(""));
      ensureCapacity(limit+1);
      for (k=0,act=start; k<=limit; k++) {
	set(k,act); act+=step;
      }
    }

    /**
     * Given an index vector 'vec' whose elements are in 0,...,'sz'-1, we
     * determine the complement of 'vec', i.e. the numbers between 0 and
     * 'sz'-1 which are not in 'vec' and write it into this vector. The
     * entries are sorted in ascending order, even if the 'vec' entries are
     * not. This is done by using a temp. boolean vector and ticking off the
     * entries of 'vec'.
     * NOTE: In general, 'vec' may contain duplicate entries and even entries
     * outside of 0,...,'sz'-1. However, if 'strict' is true, it is checked
     * whether there are duplicates or entries out of bounds in 'vec'. If so,
     * an exception is thrown.
     *
     * @param vec    Index vector
     * @param sz     See above
     * @param strict See above. Def.: false
     */
    void complement(const BaseVector<int>& vec,int sz,bool strict=false);

    /**
     * Index/inverse index method. This object is I, the inverse index J has
     * to be passed in 'indJ'.
     * Exchanges elements of I at i2=='pos1' and i2=='pos2'.
     *
     * @param indJ See above
     * @param pos1 "
     * @param pos2 "
     */
    void invIndExchange(IndexVector& indJ,int pos1,int pos2) {
      if (pos1<0 || pos1>=n || pos2<0 || pos2>=n)
	throw InvalidParameterException(EXCEPT_MSG(""));
      checkTS(TSARG); indJ.checkTS(TSARG);
      int j1=operator[](pos1),j2=operator[](pos2);
      set(pos2,j1); set(pos1,j2);
      indJ.set(j1,pos2); indJ.set(j2,pos1);
    }

    /**
     * Index/inverse index method. This object is I, the inverse index J has
     * to be passed in 'indJ'.
     * Inserts element 'elem' at position 'pos' into I. If 'pos'==-1 or
     * 'pos' is equal to the length of I, the element is appended at the end
     * (default). The ordering of I is preserved, i.e. all elements below 'pos'
     * are moved down by one position.
     *
     * @param indJ See above
     * @param elem Element to be inserted
     * @param pos  Insert position
     */
    void invIndInsert(IndexVector& indJ,int elem,int pos=-1);

    /**
     * Index/inverse index method. This object is I, the inverse index J has
     * to be passed in 'indJ'.
     * Removes element at position 'pos' within I. If 'presOrd'==true, the
     * ordering of I is preserved, i.e. all elements below 'pos' are moved up
     * by one position (default). Otherwise, the last element in I is moved to
     * position 'pos'. The latter is OK if the index represents an unordered
     * set.
     *
     * @param indJ    See above
     * @param pos     Remove position
     * @param presOrd Optional. See above. Def.: true
     */
    void invIndRemove(IndexVector& indJ,int pos,bool presOrd=true);

    /**
     * Same as 'invIndRemove', but instead of a position in I, an element in
     * I has to be passed via 'elem'. If 'elem' is not in I, the method does
     * nothing.
     */
    void invIndRemoveElem(IndexVector& indJ,int elem,bool presOrd=true) {
      if (elem>=0 && elem<indJ.size() && indJ[elem]!=-1)
	invIndRemove(indJ,indJ[elem],presOrd);
    }

    /**
     * Returns histrogram of this vector in 'hist', so that 'hist[i]'
     * contains the number of entries ==i in this vector.
     *
     * @param hist Historgram ret. here
     */
    void histogram(BaseVector<int>& hist) const;

  protected:
    /*
     * Internal methods
     * There is an element constraint: elements must not be negative
     */

    bool doCheckValid() const {
      return true;
    }

    /**
     * Elements must be nonnegative.
     */
    bool isValidElement(const int& elem) const {
      return (elem>=0);
    }

    /*
     * Only 'IndexVector' or subclasses enforce non-neg. constraint.
     * ATTENTION: dynamic_cast does not work in MEX scripts under Linux.
     * If MATLAB_MEX is set, this method always returns false.
     */
    bool isValidVecType(const BaseVector<int>* vecP) const {
#ifdef MATLAB_MEX
      return false; // cannot use dynamic_cast!!
#else
      return (DYNCAST(const IndexVector,vecP)!=0);
#endif
    }

    BaseVector<int>* newEmpty() const {
      return new IndexVector();
    }

    BaseVector<int>* copy() const {
      return new IndexVector(*this);
    }
  };
//ENDNS

#undef TSARG

#endif
