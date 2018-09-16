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
 * Desc.:  Header class StVector
 * ------------------------------------------------------------------- */

#ifndef STVECTOR_H
#define STVECTOR_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/matrix/default.h"
#include "lhotse/matrix/BaseVector.h"
#include "lhotse/matrix/FastUtils.h"
#include "lhotse/matrix/ArrayUtils.h"
#include "lhotse/matrix/TempStVector.h"
#include "lhotse/optimize/predecl.h"
#ifdef MATLAB_MEX
#include "lhotse/matif/predecl.h"
#endif

//USING(rando)

// Local macro for args to 'checkTS'
#define TSARG __FILE__,__LINE__

//BEGINNS(matrix)
  /**
   * Implements standard vector of double entries with common mathematical
   * operations and methods to store/retrieve.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class StVector : public BaseVector<double>
  {
    friend class StMatrix;
    friend class LineSearch; // still required??

  public:
    // Constructors

    /**
     * Constructs empty vector.
     */
    StVector() : BaseVector<double>() {
      setDefValues();
    }

    /**
     * Constructor. The vector is initialized with 'a'
     *
     * @param num Number of entries
     * @param a   Fill value. Def.: def. fill value
     */
    explicit StVector(int num) : BaseVector<double>(num) {
      setDefValues();
    }

    StVector(int num,double a) : BaseVector<double>() {
      setDefValues();
      init(num,a);
    }

    /**
     * Copy constructor.
     *
     * @param vec    Vector to copy/clone
     */
    StVector(const BaseVector<double>& vec) : BaseVector<double>() {
      copyDefValues(vec); // def. value members
      assignInt(vec); // copy
    }

    StVector(const BaseVecWrapper<double>& arg) : BaseVector<double>() {
      // use 'arg' fields to init. ours, then dispose of 'arg'
      convertWrapped((BaseVecWrapper<double>&) arg);
    }

    static TempStVector mask(double* pbuff,int pn,
			     const Range& rng=RangeFull::get(),
			     MemWatchBase* watch=0) {
      StVector* vec=new StVector();
      vec->reassign(pbuff,pn,rng,watch);
      return TempStVector(vec,vec);
    }

    static TempStVector mask(const ArrayHandle<double>& parr,
			     const Range& rng=RangeFull::get()) {
      StVector* vec=new StVector();
      vec->reassign(parr.p(),parr.size(),rng,parr.getMemWatch());
      return TempStVector(vec,vec);
    }

    static TempStVector mask(const BaseVector<double>& vec,
			     const Range& rng=RangeFull::get()) {
      StVector* mvec=new StVector();
      mvec->reassign(vec,rng);
      return TempStVector(mvec,mvec);
    }

    const TempStVector operator()(const Range& rng=RangeFull::get()) const {
      //cout << "StVector-()R: this=" << this << endl;
      StVector* vec=DYNCAST(StVector,subrvalInt(rng));
      return TempStVector(vec,vec);
    }

    TempStVector operator()(const Range& rng=RangeFull::get()) {
      //cout << "StVector-()L: this=" << this << endl;
      StVector* vec=DYNCAST(StVector,sublvalInt(rng));
      return TempStVector(vec,vec);
    }

    /*
     * ATTENTION:
     * 'assignVirtual' must do the same as 'operator='. If you change
     * 'operator=' below, change them as well!
     */

    StVector& operator=(const BaseVector<double>& vec) {
      //cout << "StVector-=: this=" << this << ",arg=" << &vec << endl;
      assignInt(vec);
      return *this;
    }

    StVector& operator=(double elem) {
      //cout << "StVector-=: this=" << this << ",arg=" << elem << endl;
      assignInt(elem);
      return *this;
    }

    void setDefValues() {
      BaseVector<double>::setDefValues();
      setDefFillValue(0.0); // Make sure 'defFill' is correct
    }

#ifdef MATLAB_MEX
    /**
     * Variant of 'reassign' to mask the vector stored in a
     * 'MatlabMatrix' object. The dimensions of 'mat' are ignored, the
     * whole matrix buffer (column-major) is masked as single vector.
     * A range 'rng' can be applied to this flat buffer (see
     * 'BaseVector::reassign' from flat buffer).
     * <p>
     * NOTE: There is no buffer watcher and memory control. Create and
     * use mask locally only, do not apply any Matlab functions on the
     * matrix behind 'mat' during the lifetime of the mask!
     * NOTE: Cannot call this method 'reassign', otherwise all superclass
     * 'reassign' are overwritten.
     *
     * @param mat Matlab matrix
     * @param rng Range. Def: full
     */
    virtual void maskMatlab(const MatlabMatrix& mat,
			    const Range& rng=RangeFull::get());
#endif

    // Basic arithmetic operations (public)

    /**
     * this = a    + b*c
     * this = this + b*c
     *
     * Here, a is a vector (this if not given), either b or c can be
     * a scalar (but not both). For a scalar, the cases +1/-1 are
     * treated seperately by our implementation.
     *
     * @param a Optional. Def.: this
     * @param b S.a.
     * @param c S.a.
     */
    virtual void addprod(const BaseVector<double>& a,
			 const BaseVector<double>& b,
			 const BaseVector<double>& c);
    virtual void addprod(const BaseVector<double>& a,double b,
			 const BaseVector<double>& c);
    virtual void addprod(const BaseVector<double>& a,
			 const BaseVector<double>& b,double c) {
      addprod(a,c,b);
    }

    virtual void addprod(const BaseVector<double>& b,
			 const BaseVector<double>& c);
    virtual void addprod(double b,const BaseVector<double>& c);
    virtual void addprod(const BaseVector<double>& b,double c) {
      addprod(c,b);
    }

    /**
     * this = a    + b/c
     * this = this + b/c
     *
     * Here, a is a vector (this if not given), either b or c can be
     * a scalar (but not both). For a scalar, the cases +1/-1 are
     * treated seperately by our implementation.
     *
     * @param a Optional. Def.: this
     * @param b S.a.
     * @param c S.a.
     */
    virtual void adddiv(const BaseVector<double>& a,
			const BaseVector<double>& b,
			const BaseVector<double>& c);
    virtual void adddiv(const BaseVector<double>& a,double b,
			const BaseVector<double>& c);
    virtual void adddiv(const BaseVector<double>& a,
			const BaseVector<double>& b,double c) {
      addprod(a,1.0/c,b);
    }

    virtual void adddiv(const BaseVector<double>& b,
			const BaseVector<double>& c);
    virtual void adddiv(double b,const BaseVector<double>& c);
    virtual void adddiv(const BaseVector<double>& b,double c) {
      addprod(1.0/c,b);
    }

    /**
     * this = b*c
     *
     * Either b or c can be a scalar (but not both).
     *
     * @param b S.a. Optional. Def.: this
     * @param c S.a.
     */
    virtual void prod(const BaseVector<double>& b,
		      const BaseVector<double>& c);
    virtual void prod(const BaseVector<double>& b,double c);
    virtual void prod(double b,const BaseVector<double>& c) {
      prod(c,b);
    }

    virtual void prod(const BaseVector<double>& c);
    virtual void prod(double c);

    /**
     * this = b/c
     *
     * Either b or c can be a scalar (but not both).
     * NOTE: There is no special method for b/this, just pass this for c.
     *
     * @param b S.a. Optional. Def.: this
     * @param c S.a.
     */
    virtual void div(const BaseVector<double>& b,const BaseVector<double>& c);
    virtual void div(const BaseVector<double>& b,double c) {
      prod(b,1.0/c);
    }

    virtual void div(double b,const BaseVector<double>& c);
    virtual void div(const BaseVector<double>& c);
    virtual void div(double c) {
      prod(1.0/c);
    }

    /**
     * this = a + s*(b.*c)
     *
     * @param a S.a.
     * @param b "
     * @param c "
     * @param s "
     */
    virtual void addsprod(const BaseVector<double>& a,
			  const BaseVector<double>& b,
			  const BaseVector<double>& c,double s);

    /**
     * this = a + s*(b./c)
     *
     * @param a S.a.
     * @param b "
     * @param c "
     * @param s "
     */
    virtual void addsdiv(const BaseVector<double>& a,
			 const BaseVector<double>& b,
			 const BaseVector<double>& c,double s);

    /**
     * this = a    + s*1
     * this = this + s*1
     *
     * @param a Def.: this
     * @param s Scalar
     */
    virtual void addscal(const BaseVector<double>& a,double s);
    virtual void addscal(double s);

    /**
     * ret = this^T*a
     * Optional: Weighted IP: ret=this^t*D*a, D diag. matrix
     *
     * @param a     A
     * @param d     Optional. Diag. matrix D
     * @return      Inner product
     */
    virtual double inner(const BaseVector<double>& a,
			 const BaseVector<double>* d=0) const;

    /**
     * this = a + b
     * Short for: addprod(a,1.0,b).
     *
     * @param a S.a. Def.: this
     * @param b "
     */
    virtual void add(const BaseVector<double>& a,const BaseVector<double>& b) {
      addprod(a,1.0,b);
    }

    virtual void add(const BaseVector<double>& b) {
      addprod(1.0,b);
    }

    /**
     * this = a - b
     * Short for: addprod(a,-1.0,b).
     *
     * @param a S.a.: Def.: this
     * @param b "
     */
    virtual void sub(const BaseVector<double>& a,const BaseVector<double>& b) {
      addprod(a,-1.0,b);
    }

    virtual void sub(const BaseVector<double>& b) {
      addprod(-1.0,b);
    }

    /**
     * L1 norm of this vector. 0 for empty vector
     *
     * @return Norm
     */
    virtual double norm1() const {
      return accumulate(accum_fun(std::plus<double>(),0.0),ptr_fun<double, double>(fabs));
    }

    /**
     * L2 norm of this vector
     *
     * @return Norm
     */
    virtual double norm2() const {
      return sqrt(inner(*this));
    }

    /**
     * Maximum norm of this vector (max. abs. value). 0 for empty vector
     *
     * @return Norm
     */
    virtual double normInf() const {
      return accumulate(accum_fun(BinFuncMax<double>(),0.0),ptr_fun<double, double>(fabs));
    }

    /**
     * @return Sum of elements
     */
    virtual double sum() const;

    /**
     * If the vector represents a diagonal matrix with positive entries, this
     * method computes the log determinant.
     *
     * @return Log determinant of pos. diagonal matrix
     */
    virtual double logDet() const {
      return accumulate(accum_fun(std::plus<double>(),0.0),ptr_fun<double, double>(log));
    }

    // Initialization methods

    /**
     * Initializes vector with zeros. With no arguments, the original
     * dimensions are kept
     *
     * @param num Optional. Number of entries
     */
    virtual void zeros() {
      fill(n,0.0);
    }

    virtual void zeros(int num) {
      fill(num,0.0);
    }

    /**
     * Initializes vector with standard unit vector k (1 at position k, 0
     * elsewhere). Note: k starts from 0. If no number of entries is specified,
     * the vector size is not changed
     *
     * @param k   See above
     * @param num Optional. Number of entries
     */
    virtual void eye(int k);
    virtual void eye(int k,int num);

    /**
     * Initializes vector with 1s.
     *
     * @param num Optional. Number of entries
     */
    virtual void ones() {
      fill(n,1.0);
    }

    virtual void ones(int num) {
      fill(num,1.0);
    }

    // Other methods

    /**
     * Method useful in the context of computing posterior or marginal
     * probabilities. Given vector (a_i)_i in this object, we compute
     *   s = \log \sum_i \exp(a_i)
     * in a stable way. This is done by first computing M = \max_i a_i, then
     *   s = M + \log \sum_i \exp(a_i-M).
     * NOTE: The function is called LOGSUMEXP.
     *
     * @return Sum s
     */
    virtual double stableLogSum() const;

    /**
     * Alias for 'stableLogSum'
     */
    virtual double logsumexp() const {
      return stableLogSum();
    }

    /**
     * Computes cumulative partial sums of the elements of 'vec' and stores
     * them into this vector. Namely, if the elements of 'vec' are x_0,x_1,
     * ..., then the elements of this vector will be s_0,s_1,..., with
     *   s_i = x_0 + ... + x_i.
     * NOTE: This vector and 'vec' may be the same object.
     *
     * @param vec Source vector (def.: this vector)
     */
    virtual void cumulSum(const BaseVector<double>& vec);

    virtual void cumulSum() {
      cumulSum(*this);
    }

    /**
     * Computes maximum symmetric relative difference between A (this) and
     * B ('b'), for each component
     *   max( |(b-a)/a|, |(a-b)/b| ),
     * where the expressions are 0 for denom. 0. Returns maximum of all these.
     *
     * @param b Source vector
     * @return  S.a.
     */
    virtual double maxRelDiff(const BaseVector<double>& b) const;

    // IO methods (public)

    /**
     * Loads vector from (binary) input file stream.
     * The current format is the same as for BaseVector<double>, but this
     * method can load older StVector formats as well.
     *
     * @param is    Binary input file stream
     * @param noTag If true, we do not load the tag and assume the current
     *              FF
     */
    ifstream& load(ifstream& is,bool noTag=false);

  protected:
    // Internal methods

    BaseVector<double>* copy() const {
      return new StVector(*this);
    }

    BaseVector<double>* newEmpty() const {
      return new StVector();
    }

    /**
     * Helper for arithmetic methods.
     * If b,c not given: 'checkTS' called for this vector and a, then
     * this vector is brought to length of a.
     * If b or b,c given: 'checkTS' called for this vector,a,b. Then, if a,b,
     * (c) have diff. sizes, DimMismatchException is thrown. Otherwise, this
     * vector is brought to length of a.
     * If 'doThis'==true, instead of changing the size of this, its length
     * is compared to the others just as well.
     *
     * @param a      S.a.
     * @param b      S.a. Optional
     * @param c      S.a. Optional
     * @param doThis S.a. Def.: false
     */
    void allSameSize(const BaseVector<double>& a,bool doThis=false) {
      checkTS(TSARG); a.checkTS(TSARG);
      if (!doThis)
	ensureCapacity(a.size());
      else if (n!=a.size()) throw DimMismatchException(EXCEPT_MSG(""));
    }

    void allSameSize(const BaseVector<double>& a,const BaseVector<double>& b,
		     bool doThis=false) {
      checkTS(TSARG); a.checkTS(TSARG); b.checkTS(TSARG);
      int pn=a.size();
      if (pn!=b.size()) throw DimMismatchException(EXCEPT_MSG(""));
      if (!doThis)
	ensureCapacity(pn);
      else if (n!=pn) throw DimMismatchException(EXCEPT_MSG(""));
    }

    void allSameSize(const BaseVector<double>& a,const BaseVector<double>& b,
		     const BaseVector<double>& c,bool doThis=false) {
      checkTS(TSARG); a.checkTS(TSARG); b.checkTS(TSARG); c.checkTS(TSARG);
      int pn=a.size();
      if (pn!=b.size() || pn!=c.size())
	throw DimMismatchException(EXCEPT_MSG(""));
      if (!doThis)
	ensureCapacity(pn);
      else if (n!=pn) throw DimMismatchException(EXCEPT_MSG(""));
    }
  };

  // Public methods

  inline void StVector::addprod(const BaseVector<double>& a,
				const BaseVector<double>& b,
				const BaseVector<double>& c)
  {
    allSameSize(a,b,c);
    operator=(a);
    FastUtils::addprod(buff,b.getBuffPtr_INT(),c.getBuffPtr_INT(),
		       n,1.0,step,b.getStep_INT(),c.getStep_INT(),
		       getMindex_INT(),b.getMindex_INT(),b.getMindex_INT());
  }

  inline void StVector::addprod(const BaseVector<double>& a,double b,
				const BaseVector<double>& c)
  {
    allSameSize(a,c);
    FastUtils::addsmul(buff,a.getBuffPtr_INT(),c.getBuffPtr_INT(),
		       b,n,step,a.getStep_INT(),c.getStep_INT(),
		       getMindex_INT(),a.getMindex_INT(),c.getMindex_INT());
  }

  inline void StVector::addprod(const BaseVector<double>& b,
				const BaseVector<double>& c)
  {
    allSameSize(b,c,true);
    FastUtils::addprod(buff,b.getBuffPtr_INT(),c.getBuffPtr_INT(),
		       n,1.0,step,b.getStep_INT(),c.getStep_INT(),
		       getMindex_INT(),b.getMindex_INT(),c.getMindex_INT());
  }

  inline void StVector::addprod(double b,const BaseVector<double>& c)
  {
    allSameSize(c,true);
    FastUtils::addsmul(buff,c.getBuffPtr_INT(),b,n,step,c.getStep_INT(),
		       getMindex_INT(),c.getMindex_INT());
  }

  inline void StVector::adddiv(const BaseVector<double>& a,
			       const BaseVector<double>& b,
			       const BaseVector<double>& c)
  {
    allSameSize(a,b,c);
    operator=(a);
    FastUtils::adddiv(buff,b.getBuffPtr_INT(),c.getBuffPtr_INT(),
		      n,1.0,step,b.getStep_INT(),c.getStep_INT(),
		      getMindex_INT(),b.getMindex_INT(),b.getMindex_INT());
  }

  inline void StVector::adddiv(const BaseVector<double>& a,double b,
			       const BaseVector<double>& c)
  {
    allSameSize(a,c);
    FastUtils::addinv(buff,a.getBuffPtr_INT(),c.getBuffPtr_INT(),
		      b,n,step,a.getStep_INT(),c.getStep_INT(),getMindex_INT(),
		      a.getMindex_INT(),c.getMindex_INT());
  }

  inline void StVector::adddiv(const BaseVector<double>& b,
			       const BaseVector<double>& c)
  {
    allSameSize(b,c,true);
    FastUtils::adddiv(buff,b.getBuffPtr_INT(),c.getBuffPtr_INT(),
		      n,1.0,step,b.getStep_INT(),c.getStep_INT(),
		      getMindex_INT(),b.getMindex_INT(),c.getMindex_INT());
  }

  inline void StVector::adddiv(double b,const BaseVector<double>& c)
  {
    allSameSize(c,true);
    FastUtils::addinv(buff,c.getBuffPtr_INT(),b,n,step,c.getStep_INT(),
		      getMindex_INT(),c.getMindex_INT());
  }

  inline void StVector::prod(const BaseVector<double>& b,
			     const BaseVector<double>& c)
  {
    allSameSize(b,c);
    FastUtils::prod(buff,b.getBuffPtr_INT(),c.getBuffPtr_INT(),
		    n,step,b.getStep_INT(),c.getStep_INT(),getMindex_INT(),
		    b.getMindex_INT(),c.getMindex_INT());
  }

  inline void StVector::prod(const BaseVector<double>& b,double c)
  {
    allSameSize(b);
    FastUtils::smul(buff,b.getBuffPtr_INT(),c,n,step,b.getStep_INT(),
		    getMindex_INT(),b.getMindex_INT());
  }

  inline void StVector::prod(const BaseVector<double>& c)
  {
    allSameSize(c,true);
    FastUtils::prod(buff,c.getBuffPtr_INT(),n,step,c.getStep_INT(),
		    getMindex_INT(),c.getMindex_INT());
  }

  inline void StVector::prod(double c)
  {
    FastUtils::smul(buff,c,n,step,getMindex_INT());
  }

  inline void StVector::div(const BaseVector<double>& b,
			    const BaseVector<double>& c)
  {
    allSameSize(b,c);
    FastUtils::div(buff,b.getBuffPtr_INT(),c.getBuffPtr_INT(),
		   n,step,b.getStep_INT(),c.getStep_INT(),getMindex_INT(),
		   b.getMindex_INT(),c.getMindex_INT());
  }

  inline void StVector::div(double b,const BaseVector<double>& c)
  {
    allSameSize(c);
    FastUtils::inv(buff,c.getBuffPtr_INT(),b,n,step,c.getStep_INT(),
		   getMindex_INT(),c.getMindex_INT());
  }

  inline void StVector::div(const BaseVector<double>& c)
  {
    allSameSize(c,true);
    FastUtils::div(buff,c.getBuffPtr_INT(),n,step,c.getStep_INT(),
		   getMindex_INT(),c.getMindex_INT());
  }

  inline void StVector::addscal(const BaseVector<double>& a,double s)
  {
    allSameSize(a);
    FastUtils::addscal(buff,a.getBuffPtr_INT(),s,n,step,
		       a.getStep_INT(),getMindex_INT(),a.getMindex_INT());
  }

  inline void StVector::addscal(double s)
  {
    FastUtils::addscal(buff,s,n,step,getMindex_INT());
  }

  inline double StVector::inner(const BaseVector<double>& a,
				const BaseVector<double>* d) const
  {
    checkTS(TSARG); a.checkTS(TSARG);
    if (a.size()!=n) throw DimMismatchException(EXCEPT_MSG(""));
    if (d==0) {
      if (&a==this)
	return FastUtils::inner(buff,n,step,getMindex_INT());
      else
	return FastUtils::inner(buff,a.getBuffPtr_INT(),n,step,
				a.getStep_INT(),getMindex_INT(),
				a.getMindex_INT());
    } else {
      if (&a==this)
	return FastUtils::mahal(buff,d->getBuffPtr_INT(),n,step,
				d->getStep_INT(),getMindex_INT(),
				d->getMindex_INT());
      else
	return FastUtils::mahal(buff,a.getBuffPtr_INT(),d->getBuffPtr_INT(),n,
				step,a.getStep_INT(),d->getStep_INT(),
				getMindex_INT(),a.getMindex_INT(),
				d->getMindex_INT());
    }
  }

  inline void StVector::addsprod(const BaseVector<double>& a,
				 const BaseVector<double>& b,
				 const BaseVector<double>& c,double s)
  {
    allSameSize(a,b,c);
    if (&a==this)
      FastUtils::addprod(buff,b.getBuffPtr_INT(),
			 c.getBuffPtr_INT(),n,s,step,b.getStep_INT(),
			 c.getStep_INT(),getMindex_INT(),b.getMindex_INT(),
			 c.getMindex_INT());
    else {
      FastUtils::prod(buff,b.getBuffPtr_INT(),c.getBuffPtr_INT(),n,
		      step,b.getStep_INT(),c.getStep_INT(),getMindex_INT(),
		      b.getMindex_INT(),c.getMindex_INT());
      FastUtils::addsmul(buff,a.getBuffPtr_INT(),s,n,step,a.getStep_INT(),
			 getMindex_INT(),a.getMindex_INT());
    }
  }

  inline void StVector::addsdiv(const BaseVector<double>& a,
				const BaseVector<double>& b,
				const BaseVector<double>& c,double s)
  {
    allSameSize(a,b,c);
    if (&a==this)
      FastUtils::adddiv(buff,b.getBuffPtr_INT(),c.getBuffPtr_INT(),
			n,s,step,b.getStep_INT(),c.getStep_INT(),getMindex_INT(),
			b.getMindex_INT(),c.getMindex_INT());
    else {
      FastUtils::div(buff,b.getBuffPtr_INT(),c.getBuffPtr_INT(),n,
		     step,b.getStep_INT(),c.getStep_INT(),getMindex_INT(),
		     b.getMindex_INT(),c.getMindex_INT());
      FastUtils::addsmul(buff,a.getBuffPtr_INT(),s,n,step,a.getStep_INT(),
			 getMindex_INT(),a.getMindex_INT());
    }
  }

  inline void StVector::eye(int k)
  {
    checkTS(TSARG);
    if (k>=n || k<0) throw WrongDimensionException(EXCEPT_MSG(""));
    zeros();
    operator[](k)=1.0;
  }

  inline void StVector::eye(int k,int num)
  {
    checkTS(TSARG);
    if (k>=num || k<0) throw WrongDimensionException(EXCEPT_MSG(""));
    zeros(num);
    operator[](k)=1.0;
  }

  inline void StVector::cumulSum(const BaseVector<double>& vec)
  {
    checkTS(TSARG); vec.checkTS(TSARG);
    ensureCapacity(vec.size());
    FastUtils::cumulSum(buff,vec.getBuffPtr_INT(),n,step,vec.getStep_INT(),
			getMindex_INT(),vec.getMindex_INT());
  }

  inline double StVector::stableLogSum() const
  {
    checkTS(TSARG);

    return FastUtils::logsumexp(buff,n,step,getMindex_INT());
  }

  inline double StVector::sum() const
  {
    checkTS(TSARG);

    return FastUtils::sum(buff,n,step,getMindex_INT());
  }
//ENDNS

#undef TSARG

#endif
