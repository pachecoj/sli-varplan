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
 * Desc.:  Header class StMatrix
 * ------------------------------------------------------------------- */

#ifndef STMATRIX_H
#define STMATRIX_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

/*
 * TODO:
 * 0. Structure pattern support for more methods (e.g. 'apply1')
 * 1. IMPORTANT: Clean up old Cholesky mess (no more separ. diagonal).
 *    ==> Dangerous, some old methods require 'lowNDg' as pattern
 *        and modify 'strpatt' to that value!
 *    ==> Kick out old methods!
 * 4. LAPACK support:
 *    Check how LAPACK reports errors, catch those. Then, bind in as
 *    required.
 *    ==> NOTE: Some LAPACK routines already used, no error checking!
 *    In general, errors reported via the INFO argument, but if args
 *    are wrong (negative INFO), also a special function is called
 *    which typically leads to program termination.
 */

/*
 * About LAPACK
 *
 * There is no official C interface to LAPACK. clapack.h supplied by ATLAS
 * is incomplete and has errors, we do not use it. We call the Fortran
 * routines of LAPACK directly. At present, only a few LAPACK routines are
 * actually used here.
 *
 * About BLAS
 *
 * We use the CBLAS C interface. The most important BLAS functions are
 * wrapped in class 'FastUtils'.
 *
 * NOTE: Man pages for most BLAS, LAPACK functions can be found at:
 *   http://www.mathkeisan.com/man.html  [copyright NEC]
 */

#include <algorithm>
#include "lhotse/matrix/default.h"
#include "lhotse/matrix/BaseMatrix.h"
#include "lhotse/matrix/FastUtils.h"
#include "lhotse/matrix/StVector.h"
#include "lhotse/matrix/TempStMatrix.h"
#ifdef MATLAB_MEX
#include "lhotse/matif/predecl.h"
#endif
#ifdef MATLAB_DEBUG_OLD
#include "lhotse/MatlabDebug.h"
#endif

//USING(rando);

// Local macro for args to 'checkTS'
#define TSARG __FILE__,__LINE__

//BEGINNS(matrix)
  /**
   * Subclass of 'BaseMatrix<double>', specialized to double elements
   * (unconstrained). Numerical methods are implemented here.
   * <p>
   * Structure pattern:
   * Some arithmetic methods work on matrices with special structure,
   * which require only a part of the full rectangular matrix frame
   * (f.ex. just the upper triangle). For these methods, the member
   * 'strpatt' has to be set correctly before they are called
   * ('setStrctPatt').
   * The def. and initial pattern is the full matrix: 'normal'
   * (see 'MatStrct'). The other structure patterns require square
   * matrices.
   *
   * For such a method, a symmetric matrix can be repres. by a lower or
   * upper triangular matrix (either in most cases). Typically, it must
   * have a diagonal, in spec. cases the diag. can be passed as extra
   * vector.
   *
   * NOTE: If a matrix is re-sized to non-square with a square 'strpatt',
   * this mismatch is noticed only once a method is called which needs
   * 'strpatt'. It will reset 'strpatt' to 'normal' for a non-square
   * matrix (ATTENTION!).
   * NOTE: 'strpatt' is NOT stored in file. 'load' init. 'strpatt' to
   * 'normal'. There is no special storage format for structured
   * matrices (check subclasses).
   * NOTE: 'strpatt' is used only by methods marked (in comments) as
   * "Uses structure pattern", all others ignore 'strpatt' (work as if
   * 'normal')!
   *
   * ATTENTION: The standard methods which copy matrices (copy constructor,
   * operator=) do NOT use 'strpatt' in the sense that they always copy
   * all elements. If the source matrix is of type 'StMatrix', the target
   * matrix inherits its str. patt., otherwise the target matrix has str.
   * patt. 'normal'. The method 'assignStrct' honours the str. patt. of the
   * source matrix, which the target matrix inherits (and only elem. acc.
   * to the str. patt. are actually copied).
   * Also, none of the standard methods to access matrix parts use 'strpatt',
   * e.g. masking the row of a matrix will always use the full row, indep.
   * of 'strpatt'. The mask matrices generated have str. patt. 'normal'.
   *
   * NOTE: 'strpatt' is not conserved subindexing. If the matrix is resized
   * to become non-square, it is reset to 'normal'.
   * ==> Best to set 'strpatt' to the desired value just before the
   *     arithm. operation is done
   * The structure pattern is conserved in assignment (if the r-value is
   * a 'StMatrix' itself) and in copy construction.
   *
   * In general, the str. patt. should be regarded "volatile" and applic.
   * only to specific arithmetic methods (these are marked). Other methods,
   * esp. the ones inherited from 'BaseMatrix<double>', ignore it. To be
   * safe, reassign the str. patt. just before the arithmetic method which
   * requires it.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class StMatrix : public BaseMatrix<double>
  {
  protected:
    // Additional members

    mutable uchar strpatt; // See header comment

  public:
    // Constructors

    StMatrix() : BaseMatrix<double>(),strpatt(normal) {
      setDefValues();
    }

    StMatrix(int rows,int cols) : BaseMatrix<double>(rows,cols),
				  strpatt(normal) {
      setDefValues();
    }

    explicit StMatrix(int rows) : BaseMatrix<double>(rows,rows),
				  strpatt(normal) {
      setDefValues();
    }

    StMatrix(const StMatrix& mat) : BaseMatrix<double>(),strpatt(normal) {
      mat.checkTS(TSARG);
      copyDefValues(mat); // copy def. values
      assignInt(mat);
      mat.updateStrctPatt();
      strpatt=mat.strpatt; // does copy str. patt.
    }

    StMatrix(const BaseMatWrapper<double>& arg) : BaseMatrix<double>(),
						  strpatt(normal) {
      // use 'arg' fields to init. ours, then dispose of 'arg'
      convertWrapped((BaseMatWrapper<double>&) arg);
    }

    static TempStMatrix mask(double* pbuff,int pm,int pn,int strid,
			     MemWatchBase* watch=0) {
      StMatrix* mat=new StMatrix();
      mat->reassign(pbuff,pm,pn,strid,watch);
      return TempStMatrix(mat,mat);
    }

    static TempStMatrix mask(const ArrayHandle<double>& parr,int pm,int pn,
			     int strid) {
      StMatrix* mat=new StMatrix();
      mat->reassign(parr.p(),pm,pn,strid,parr.getMemWatch());
      return TempStMatrix(mat,mat);
    }

    static TempStMatrix mask(double* pbuff,int strid,const Range& rngR,
			     const Range& rngC,MemWatchBase* watch=0) {
      StMatrix* mat=new StMatrix();
      mat->reassign(pbuff,strid,rngR,rngC,watch);
      return TempStMatrix(mat,mat);
    }

    static TempStMatrix mask(const ArrayHandle<double>& pbuff,int strid,
			     const Range& rngR,const Range& rngC) {
      StMatrix* mat=new StMatrix();
      mat->reassign(pbuff.p(),strid,rngR,rngC,pbuff.getMemWatch());
      return TempStMatrix(mat,mat);
    }

    static TempStMatrix mask(const BaseMatrix<double>& mat,
			     const Range& rngR=RangeFull::get(),
			     const Range& rngC=RangeFull::get()) {
      StMatrix* mmat=new StMatrix();
      mmat->reassign(mat,rngR,rngC);
      return TempStMatrix(mmat,mmat);
    }

    static TempStMatrix mask(const BaseVector<double>& vec,bool col=true) {
      StMatrix* mmat=new StMatrix();
      mmat->reassign(vec,col);
      return TempStMatrix(mmat,mmat);
    }

    void setDefValues() {
      BaseMatrix<double>::setDefValues();
      setDefFillValue(0.0); // Make sure 'defFill' is correct
    }

#ifdef MATLAB_MEX
    /**
     * Assigns mask matrix to wrap Matlab matrix 'mat'.
     * ATTENTION: There is no memory watcher for 'mat', so this vector
     * won't use one either
     * ==> Only mask 'MatlabMatrix' objects temporarily!
     *
     * @param mat Matlab matrix
     * @param rngR Row range. Def: all
     * @param rngC Col range. Def: all
     */
    virtual void maskMatlab(MatlabMatrix& mat,
			    const Range& rngR=RangeFull::get(),
			    const Range& rngC=RangeFull::get());
#endif

    const TempStMatrix operator()(const Range& rngR=RangeFull::get(),
				  const Range& rngC=RangeFull::get()) const {
      StMatrix* mat=DYNCAST(StMatrix,subrvalInt(rngR,rngC));
      return TempStMatrix(mat,mat);
    }

    const TempStVector operator()(const Range& rngR,int rngC) const {
      StVector* vec=DYNCAST(StVector,subrvalInt(rngR,rngC));
      return TempStVector(vec,vec);
    }

    const TempStVector operator()(int rngR,const Range& rngC=
				  RangeFull::get()) const {
      StVector* vec=DYNCAST(StVector,subrvalInt(rngR,rngC));
      return TempStVector(vec,vec);
    }

    TempStMatrix operator()(const Range& rngR=RangeFull::get(),
			    const Range& rngC=RangeFull::get()) {
      StMatrix* mat=DYNCAST(StMatrix,sublvalInt(rngR,rngC));
      return TempStMatrix(mat,mat);
    }

    TempStVector operator()(const Range& rngR,int rngC) {
      StVector* vec=DYNCAST(StVector,sublvalInt(rngR,rngC));
      return TempStVector(vec,vec);
    }

    TempStVector operator()(int rngR,const Range& rngC=RangeFull::get()) {
      StVector* vec=DYNCAST(StVector,sublvalInt(rngR,rngC));
      return TempStVector(vec,vec);
    }

    TempStVector diag(int off=0,const Range& rng=RangeFull::get()) const {
      StVector* vec=DYNCAST(StVector,diagInt(off,rng));
      return TempStVector(vec,vec);
    }

    TempStVector operator[](int row) const {
      StVector* vec=DYNCAST(StVector,subrvalInt(row,RangeFull::get()));
      return TempStVector(vec,vec);
    }

    /**
     * The structure pattern is inherited from 'mat' if 'mat' is of type
     * 'StMatrix', otherwise it is set to 'normal'.
     * ATTENTION: This does not mean that the str. patt. is properly
     * supported here. This method copies all elements from 'mat'. To
     * assign only the elements specified by the str. patt., use
     * 'assignStrct'.
     *
     * @param mat Source argument
     */
    StMatrix& operator=(const BaseMatrix<double>& mat) {
      assignInt(mat);
      strpatt=mat.getStrctPatt();
      return *this;
    }

    StMatrix& operator=(const BaseVector<double>& vec) {
      assignInt(vec);
      strpatt=normal;
      return *this;
    }

    StMatrix& operator=(double elem) {
      assignInt(elem);
      strpatt=normal;
      return *this;
    }

    /*
     * ATTENTION:
     * 'assignVirtual' must do the same as 'operator='. If you change
     * 'operator=' below, change them as well!
     */

    void assignVirtual(const TempMatMethods* arg) {
      //cout << "StMatrix::assignVirtual" << endl;
      // Is 'arg' a matrix?
      const BaseMatrix<double>* aptr=DYNCAST(const BaseMatrix<double>,arg);
      if (aptr==0) {
	// Is 'arg' a vector?
	const BaseVector<double>* avptr=DYNCAST(const BaseVector<double>,arg);
	if (avptr==0) throw InternalException(EXCEPT_MSG(""));
	operator=(*avptr);
      } else
	operator=(*aptr);
    }

    void assignVirtual(const TempMatMethods_AssignElement_Generic* arg) {
      //cout << "StMatrix::assignVirtual" << endl;
      const TempMatMethods_AssignElement<double>* aptr=
	DYNCAST(const TempMatMethods_AssignElement<double>,arg);
      if (aptr==0) throw InternalException(EXCEPT_MSG(""));
      operator=(aptr->elem);
    }
    
    /**
     * Same as 'operator=', but structure pattern of argument 'mat' is
     * observed ('mat' must be square for a square pattern, otherwise its
     * pattern  is changed to 'normal' before the assignment), in that
     * only elements corr. to the str. patt. are in fact copied.
     * <p>
     * If this matrix is normal, it inherits the structure pattern of 'mat'.
     * If a mask (and space for 'mat'), its str. patt. is not changed.
     * If this matrix has a size different from 'mat', it is resized and
     * all elements are set to the def. fill value before the assignment.
     * Otherwise, elements not covered by the assignment pattern are not
     * accessed. This is true also if this matrix is a mask and large enough
     * to store the new pattern.
     * <p>
     * NOTE: 'mat' and this matrix can have diff. structure patt. !=
     * 'normal', if they corr. to parts of the same size. If this
     * matrix has str. patt. != 'normal', it is the target pattern
     * for 'assignInt'. An exc. is thrown if source and target pattern are
     * not compatible.
     *
     * @param mat Source matrix
     */
    virtual void assignStrct(const StMatrix& mat) {
      mat.updateStrctPatt();
      if (strpatt==normal) {
	assignInt(mat,mat.strpatt);
	if (~isMask()) strpatt=mat.strpatt; // inherit pattern
      } else
	assignInt(mat,mat.strpatt,strpatt);
    }

    StVector max(bool rowMaj=false) const {
      BaseVecWrapper<double> resWrap(newEmptyVector());
      BaseMatrix<double>::max(*(resWrap.getRep()),rowMaj);

      return resWrap;
    }

    StVector min(bool rowMaj=false) const {
      BaseVecWrapper<double> resWrap(newEmptyVector());
      BaseMatrix<double>::min(*(resWrap.getRep()),rowMaj);

      return resWrap;
    }

    /**
     * Sets structure pattern flag to value 'spatt'. See
     * 'MatStrct::XXX' for poss. values. See header comment.
     * NOTE: Values != 'normal' require a square matrix!
     *
     * @param spatt S.a.
     */
    virtual void setStrctPatt(uchar spatt) const {
      if (spatt>last)
	throw InvalidParameterException(EXCEPT_MSG("spatt"));
      if (spatt!=normal && m!=n)
	throw InvalidParameterException(EXCEPT_MSG("Need square matrix"));
      strpatt=spatt;
    }

    /**
     * Returns structure pattern flag.
     * NOTE: If this matrix is not square, 'strpatt' is reset to
     * 'normal'.
     *
     * @return S.a.
     */
    uchar getStrctPatt() const {
      updateStrctPatt();

      return strpatt;
    }

    /**
     * Same as 'fill', but the structure pattern 'strpatt' is observed.
     * If 'rows'!='cols', 'strpatt' is set to 'normal'.
     * NOTE: If the matrix size changes, 'fill' is used (which ignores
     * 'strpatt'). To avoid that, use 'expand' then 'fillStrct'.
     * NOTE: Uses structure pattern.
     *
     * @param rows Optional. Number of rows
     * @param cols Optional. Number of cols. Def.: same as rows
     * @param a    Optional. Fill element. Def.: def. fill value
     */
    virtual void fillStrct() {
      fillStrct(m,n,defFill);
    }

    virtual void fillStrct(int rows,int cols,double a);

    virtual void fillStrct(int rows,int cols=-1) {
      if (cols==-1)
	fillStrct(rows,rows,defFill);
      else
	fillStrct(rows,cols,defFill);
    }

    // Arithmetic operations (matrix-vector, BLAS-II)

    /**
     * a = alpha*op(M)*b + beta*a, M==this
     *
     * Form of M dep. on 'strpatt': normal for 'normal', symmetric
     * for 'upper' or 'lower'. No other values allowed.
     * For triangular M, use 'triMulVec'.
     * <p>
     * op(M) is M for 'trans'==false (def.), M^T otherwise.
     * NOTE: 'a' has to have correct size if 'beta'!=0.
     * NOTE: Uses structure pattern.
     * <p>
     * NOTE: Special implementation (only for str. patt. 'normal'): This
     * matrix can be indexed mask, but no copy is drawn here.
     *
     * @param avec  Input/output vector
     * @param bvec  Input vector
     * @param trans Transpose M? Def.: false
     * @param alpha Scalar. Def: 1
     * @param beta  Scalar. Def.: 0
     */
    virtual void mulVec(StVector& avec,const StVector& bvec,bool trans=false,
			double alpha=1.0,double beta=0.0) const;

    /**
     * a = op(M)*a, M==this (M triangular)
     *
     * M is triangular, form dep. on structure pattern ('normal' not
     * allowed).
     * op(M) is M for 'trans'==false (def.), M^T otherwise.
     * NOTE: Uses structure pattern.
     *
     * @param avec  Input/output vector
     * @param trans Transpose M? Def.: false
     */
    virtual void triMulVec(StVector& avec,bool trans=false) const;

    /**
     * a = op(M^{-1})*a, M==this (M triangular)
     *
     * M is triangular, form dep. on structure pattern ('normal' not
     * allowed).
     * op(X) is X for 'trans'==false (def.), X^T otherwise.
     * NOTE: Uses structure pattern.
     * NOTE: For multiple r.h.s., use BLAS-III 'backsubst' version.
     *
     * @param avec  Input/output vector
     * @param trans S.a. Def.: false
     */
    virtual void backsubst(StVector& avec,bool trans=false) const;

    /**
     * a = op(M^{-1}*D^{-1})*a, M==this (M triangular)
     *
     * M is triangular, form dep. on structure pattern ('normal' not
     * allowed). D diagonal, positive.
     * op(X) is X for 'trans'==false (def.), X^T otherwise.
     * NOTE: Uses structure pattern.
     * NOTE: For multiple r.h.s., use BLAS-III 'backsubst' version.
     *
     * @param avec  Input/output vector
     * @param dvec  Diag. matrix, positive
     * @param trans S.a. Def.: false
     */
    virtual void backsubst(StVector& avec,const StVector& dvec,
			   bool trans=false) const;

    /**
     * M = M + alpha*b*c^T, M==this (rank 1 update)
     * M = M + alpha*b*b^T (if c not given)
     *
     * If c is given, M must have structure 'normal'. Otherwise, it can
     * be symmetric ('upper' / 'lower').
     * NOTE: Uses structure pattern.
     * NOTE: Use 'symMul' for updates larger than rank 1.
     *
     * @param bvec  Input vector
     * @param alpha Scalar. Def.: 1
     * @param cvec  Optional. Input vector
     */
    virtual void rankOne(const StVector& bvec,double alpha=1.0);
    virtual void rankOne(const StVector& bvec,double alpha,
			 const StVector& cvec);

    /**
     * M = M + alpha*a*1^T, M==this ('col'==true)
     * M = M + alpha*1*a^T, M==this ('col'==false)
     *
     * Special case of 'rankOne'. M must have structure 'normal'.
     * NOTE: Uses structure pattern.
     *
     * @param avec  Input vector a
     * @param alpha Scalar. Def.: 1
     * @param col   S.a. Def.: true
     */
    virtual void addVec(const StVector& avec,double alpha=1.0,bool col=true);

    /**
     * M = M + alpha*(b*c^T + c*b^T), M==this (rank 2 update, M symmetric)
     *
     * M must be symmetric (str. pattern 'upper' / 'lower').
     * NOTE: Uses structure pattern.
     *
     * @param bvec  Input vector
     * @param cvec  Input vector
     * @param alpha Scalar. Def.: 1
     */
    virtual void symRankTwo(const StVector& bvec,const StVector& cvec,
			    double alpha=1.0);

    // Arithmetic methods (matrix-matrix, BLAS-III)

    /**
     * A = alpha*A*op(M^{-1}), M==this (M triangular, 'left'==false)
     * A = alpha*op(M^{-1})*A, M==this (M triangular, 'left'==true)
     *
     * M is triangular, form dep. on structure pattern ('normal' not
     * allowed). A must be normal ('normal').
     * op(X) is X for 'trans'==false (def.), X^T otherwise.
     * NOTE: Uses structure pattern.
     * NOTE: The def. Cholesky usage leads to upper triang. M, then
     * A = A*M^{-1}. This is more efficient than A = A*M^{-T}.
     *
     * @param amat  Input/output matrix
     * @param left  S.a. Def.: false
     * @param trans S.a. Def.: false
     * @param alpha S.a. Def.: 1
     */
    virtual void backsubst(StMatrix& amat,bool left=false,bool trans=false,
			   double alpha=1.0) const;

    /**
     * A = A*op(M^{-1}*D^{-1}), M==this (M triangular, 'left'==false)
     * A = op(M^{-1}*D^{-1})*A, M==this (M triangular, 'left'==true)
     *
     * M is triangular, form dep. on structure pattern ('normal' not
     * allowed). A must be normal ('normal'). D is diagonal, positive.
     * op(X) is X for 'trans'==false (def.), X^T otherwise.
     * NOTE: Uses structure pattern.
     * NOTE: The def. Cholesky usage leads to upper triang. M, then
     * A = A*M^{-1}*D^{-1}.
     *
     * @param amat  Input/output matrix
     * @param dvec  Diag. matrix, positive
     * @param left  S.a. Def.: false
     * @param trans S.a. Def.: false
     */
    virtual void backsubst(StMatrix& amat,const StVector& dvec,
			   bool left=false,bool trans=false) const;

    /**
     * M = alpha*op(A)*op(B) + beta*M, M==this
     *
     * M must be normal ('normal') and have the right size if
     * 'beta'!=0. If 'beta'==0, the str. patt. is set to 'normal'.
     * op(X)==X or X^T, dep. on 'xtrans'. A,B can be normal
     * ('normal') or symmetric ('upper' / 'lower').
     * For triangular A or B, use 'triMul'.
     * NOTE: This object must be different from 'amat', 'bmat'.
     * <p>
     * NOTE: Uses structure pattern.
     * NOTE: The following cases are not supported by BLAS-III and require
     * temp. copies:
     * - both A,B symmetric: one of them is copied
     * - A or B symmetric, 'xtrans' flag for other ==true: other is explic.
     *   transposed
     * ==> try to avoid these.
     *
     * @param amat   Input matrix
     * @param atrans S.a. Def: false
     * @param bmat   Input matrix
     * @param btrans S.a. Def: false
     * @param alpha  Scalar. Def.: 1
     * @param beta   Scalar. Def.: 0
     */
    virtual void mul(const StMatrix& amat,bool atrans,
		     const StMatrix& bmat,bool btrans=false,double alpha=1.0,
		     double beta=0.0);
    virtual void mul(const StMatrix& amat,const StMatrix& bmat,
		     double alpha=1.0,double beta=0.0) {
      mul(amat,false,bmat,false,alpha,beta);
    }

    /**
     * M = alpha*A*A^T + beta*M, M==this (M symmetric, 'atrans'==false)
     * M = alpha*A^T*A + beta*M, M==this (M symmetric, 'atrans'==true)
     *
     * M must be symmetric ('upper' / 'lower'). Must have right
     * size if beta!=0. A must be normal ('normal').
     * NOTE: If beta==0 and this matrix is not symmetric, 'strpatt' is set
     * to 'lower'. To force 'upper', set 'strpatt'='upper'
     * before call.
     * NOTE: Uses structure pattern.
     *
     * @param amat   Input matrix
     * @param atrans S.a.
     * @param alpha  Scalar. Def.: 1
     * @param beta   Scalar. Def.: 0
     */
    virtual void symMul(const StMatrix& amat,bool atrans,double alpha=1.0,
			double beta=0.0);

    /**
     * M = alpha*A*D*A^T + beta*M, M==this (M symmetric, 'atrans'==false)
     * M = alpha*A^T*D*A + beta*M, M==this (M symmetric, 'atrans'==true)
     *
     * M must be symmetric ('upper' / 'lower'). Must have right
     * size if beta!=0. A must be normal ('normal'). D is a diagonal
     * matrix given by 'dvec'.
     * NOTE: If beta==0 and this matrix is not symmetric, 'strpatt' is set
     * to 'lower'. To force 'upper', set 'strpatt'='upper'
     * before call.
     * NOTE: Uses structure pattern.
     * NOTE: Uses BLAS-I only, not BLAS-III ==> Slower than 'symMul'!
     * If D nonneg., it is better to compute A*D^(1/2), then use 'symMul'.
     *
     * @param amat   Input matrix
     * @param atrans S.a.
     * @param dvec   S.a.
     * @param alpha  Scalar. Def.: 1
     * @param beta   Scalar. Def.: 0
     */
    virtual void symMulDiag(const StMatrix& amat,bool atrans,
			    const StVector& dvec,double alpha=1.0,
			    double beta=0.0);

    /**
     * M = alpha*(A*B^T + B*A^T) + beta*M, M==this (M symmetric,
     *                                              'atrans'==false)
     * M = alpha*(A^T*B + B^T*A) + beta*M, M==this (M symmetric,
     *                                              'atrans'==true)
     *
     * M must be symmetric ('upper' / 'lower'). Must have right
     * size if beta!=0. A,B must be normal ('normal').
     * NOTE: If beta==0 and this matrix is not symmetric, 'strpatt' is set
     * to 'lower'. To force 'upper', set 'strpatt'='upper'
     * before call.
     * NOTE: Uses structure pattern.
     *
     * @param amat   Input matrix
     * @param atrans S.a.
     * @param bmat   Input matrix
     * @param alpha  Scalar. Def.: 1
     * @param beta   Scalar. Def.: 0
     */
    virtual void symMul(const StMatrix& amat,bool atrans,
			const StMatrix& bmat,double alpha=1.0,
			double beta=0.0);

    /**
     * A = alpha*op(M)*A, M==this (M triangular, 'left'==true)
     * A = alpha*A*op(M), M==this (M triangular, 'left'==false)
     *
     * M must triangular (!='normal'). A must be normal ('normal').
     * op(M)==M or M^T, dep. on 'mtrans'.
     * NOTE: Uses structure pattern.
     *
     * @param amat   Input/output matrix
     * @param mtrans S.a.
     * @param left   S.a.
     * @param alpha  Scalar. Def.: 1
     */
    virtual void triMul(StMatrix& amat,bool mtrans,bool left,
			double alpha=1.0) const;

    // Other arithmetic methods

    /**
     * M = diag(a)*B,      M==this ('inv'==false)
     * M = diag(a)^{-1}*B, M==this ('inv'==true)
     *
     * M,B can be the same matrix. B's str. pattern must not be
     * 'lowNDg', 'uppNDg'. M inherits B's str. patt.
     * NOTE: Uses structure pattern.
     * If 'incr'==true, result is added to M (must have same size and
     * str. patt.).
     * <p>
     * NOTE: Special implementation: no matrices are copied even if indexed
     * masks.
     *
     * @param avec Input vector
     * @param bmat Input matrix
     * @param inv  S.a. Def.: false
     * @param incr S.a. Def.: false
     */
    virtual void mulDiag(const StVector& avec,const StMatrix& bmat,
			 bool inv=false,bool incr=false);

    /**
     * M = B*diag(a),      M==this ('inv'==false)
     * M = B*diag(a)^{-1}, M==this ('inv'==true)
     *
     * See 'mulDiag' above.
     * NOTE: Special implementation: no matrices are copied even if indexed
     * masks.
     *
     * @param bmat Input matrix
     * @param avec Input vector
     * @param inv  S.a. Def.: false
     * @param incr S.a. Def.: false
     */
    virtual void mulDiag(const StMatrix& bmat,const StVector& avec,
			 bool inv=false,bool incr=false);

    /**
     * M = alpha*D*A*D + beta*M, D diagonal
     *
     * Both M, A must be symmetric ('lower' / 'upper'). If
     * 'beta'==0, M inherits str. patt. from A. If 'beta'!=0, A, M must be
     * different objects.
     * NOTE: Uses structure pattern.
     * NOTE: Uses BLAS-I only. Matrices are not copied.
     *
     * @param amat  Input matrix A
     * @param dvec  Diag. matrix D
     * @param alpha S.a. Def.: 1
     * @param beta  S.a. Def.: 0
     */
    virtual void mulDiagSym(const StMatrix& amat,const StVector& dvec,
			    double alpha=1.0,double beta=0.0);

    /**
     * M = alpha*(D*A + A*D) + beta*M, D diagonal
     *
     * Both M, A must be symmetric ('lower' / 'upper'). If
     * 'beta'==0, M inherits str. patt. from A. If 'beta'!=0, A, M must be
     * different objects.
     * NOTE: Uses structure pattern.
     * NOTE: Uses BLAS-I only. Matrices are not copied.
     *
     * @param amat  Input matrix A
     * @param dvec  Diag. matrix D
     * @param alpha S.a. Def.: 1
     * @param beta  S.a. Def.: 0
     */
    virtual void addMulDiagSym(const StMatrix& amat,const StVector& dvec,
			       double alpha=1.0,double beta=0.0);

    /**
     * M = A + alpha*B, M==this
     *
     * M,A can be the same matrix. If A not given, A==M. A,B must have same
     * structure pattern, this is inherited by M.
     * NOTE: Uses structure pattern.
     * NOTE: 'bmat' can be a scalar, in which case it repres.
     * the matrix s*1*1^T, s the scalar (all s).
     *
     * @param amat  Input matrix. Def.: this matrix
     * @param bmat  Input matrix. Can be scalar
     * @param alpha Scalar
     */
    virtual void addsmul(const StMatrix& bmat,double alpha);

    virtual void addsmul(const StMatrix& amat,const StMatrix& bmat,
			 double alpha) {
      if (amat.m!=bmat.m || amat.n!=bmat.n)
	throw DimMismatchException(EXCEPT_MSG(""));
      amat.updateStrctPatt(); bmat.updateStrctPatt();
      if (amat.strpatt!=bmat.strpatt)
	throw WrongStatusException(EXCEPT_MSG("amat/bmat.strpatt"));
      // Need to observe A's strpatt!
      assignInt(amat,amat.strpatt);
      strpatt=amat.strpatt;
      addsmul(bmat,alpha);
    }

    virtual void addsmul(double bmat,double alpha);

    virtual void addsmul(const StMatrix& amat,double bmat,double alpha) {
      amat.updateStrctPatt();
      // Need to observe A's strpatt!
      assignInt(amat,amat.strpatt);
      strpatt=amat.strpatt;
      addsmul(bmat,alpha);
    }

    /**
     * tr(M), M==this
     *
     * Sum of diagonal elements (also for non-square matrix).
     * NOTE: Uses structure pattern.
     *
     * @return Trace
     */
    virtual double trace() const {
      updateStrctPatt();
      if (strpatt==uppNDg || strpatt==lowNDg)
	return (double) n;
      else {
	Handle<StVector> msk(DYNCAST(StVector,diagInt(0,RangeFull::get())));
	return msk->sum();
      }
    }

    /**
     * tr(A^T*M), M==this
     *
     * Sum over comp.wise product. A,M must have same size. If any of A,M has
     * structure pattern 'upper' / 'lower', it is taken to be
     * symmetric. A,M must have the same structure pattern ('strpatt').
     * NOTE: Uses structure pattern.
     *
     * @param amat Input matrix A
     * @return     S.a.
     */
    virtual double traceProd(const StMatrix& a) const;

    /**
     * Sums columns, ret. in col. vector 'vec'. If 'rowMaj'==true, sums rows
     * and ret. in row. vector 'vec'.
     * If 'increm'==true and 'vec' has correct size, it is not reset to 0.
     *
     * @param vec    Output vector
     * @param rowMaj S.a. Def.: false
     * @param increm S.a. Def.: false
     */
    virtual void sum(StVector& vec,bool rowMaj=false,bool increm=false) const;

    /**
     * Same as above, but vector is returned.
     *
     * @param rowMaj S.a. Def.: false
     * @return       Sum vector
     */
    virtual StVector sum(bool rowMaj=false) const {
      BaseVecWrapper<double> resWrap(newEmptyVector());
      sum(*DYNCAST(StVector,resWrap.getRep()),rowMaj);

      return resWrap;
    }

    /**
     * M = A + diag(d), M==this
     *
     * If A not given, A==M. Size of d must be the smaller of rows, cols
     * of A. M inherits 'strpatt' of A (must not be unit diagonal structure).
     *
     * @param a A. Def.: this matrix M
     * @param d D
     */
    void addDiag(const StVector& d) {
      if (d.size()!=std::min(m,n))
	throw DimMismatchException(EXCEPT_MSG("d"));
      updateStrctPatt();
      if (strpatt==uppNDg || strpatt==lowNDg)
	throw InvalidParameterException(EXCEPT_MSG("strpatt"));
      Handle<StVector> msk(DYNCAST(StVector,diagInt(0,RangeFull::get())));
      msk->addprod(1.0,d);
    }

    void addDiag(const StMatrix& a,const StVector& d) {
      a.updateStrctPatt();
      if (a.strpatt==uppNDg || a.strpatt==lowNDg)
	throw InvalidParameterException(EXCEPT_MSG("a.strpatt"));
      if (d.size()!=std::min(a.rows(),a.cols()))
	throw DimMismatchException(EXCEPT_MSG("a,d"));
      // Need to observe A's strpatt!
      assignInt(a,a.strpatt);
      strpatt=a.strpatt;
      addDiag(d);
    }

    /**
     * M = alpha*A, M==this
     *
     * If A not given, A==M. For unit diagonal strpatt, the (unit) diagonal
     * is NOT multiplied. M inherits A's strpatt.
     * NOTE: Uses structure pattern.
     *
     * @param amat  A. Def: this matrix M
     * @param alpha alpha
     */
    void smul(const StMatrix& amat,double alpha);
    void smul(double alpha);

    /**
     * M = A + alpha*I, M==this
     *
     * If A not given, A==M. Structure pattern must not be unit diag. M
     * inherits A's strpatt.
     * NOTE: Uses structure pattern.
     *
     * @param a     A. Def: this matrix M
     * @param alpha alpha
     */
    void addseye(double alpha) {
      updateStrctPatt();
      if (strpatt==uppNDg || strpatt==lowNDg)
	throw WrongStatusException(EXCEPT_MSG("Wrong strpatt"));
      Handle<StVector> dg(DYNCAST(StVector,diagInt(0,RangeFull::get())));
      dg->addscal(alpha);
    }

    void addseye(const StMatrix& a,double alpha) {
      a.updateStrctPatt();
      if (a.strpatt==uppNDg || a.strpatt==lowNDg)
	throw WrongStatusException(EXCEPT_MSG("Wrong a.strpatt"));
      assignInt(a,a.strpatt);
      strpatt=a.strpatt;
      addseye(alpha);
    }

    /**
     * M = A \circ B (Schur prod., componentwise)
     *
     * If A not given, A==M. A,B must have same structure pattern. M
     * inherits this pattern.
     * NOTE: Uses structure pattern.
     *
     * @param a     A. Def.: this matrix M
     * @param b     B
     */
    void prod(const StMatrix& b);
    void prod(const StMatrix& a,const StMatrix& b);

    /**
     * Addition: M = A+B
     * Short for 'addsmul(a,b,1.0)'
     *
     * @param a A. Def.: this matrix
     * @param b B
     */
    void add(const StMatrix& a,const StMatrix& b) {
      addsmul(a,b,1.0);
    }

    void add(const StMatrix& b) {
      addsmul(b,1.0);
    }

    /**
     * Subtraction: M = A-B
     * Short for 'addsmul(a,b,-1.0)'
     *
     * @param a A. Def.: this matrix
     * @param b B
     */
    void sub(const StMatrix& a,const StMatrix& b) {
      addsmul(a,b,-1.0);
    }

    void sub(const StMatrix& b) {
      addsmul(b,-1.0);
    }

    /**
     * diag( M M^T )    ['rowMaj'==false]
     * diag( M^T M )    ['rowMaj'==true]
     *
     * Squares all elements, adds cols and returns result in 'vec' (col.
     * vector). If 'increm' is true and 'vec' has the correct size, it is not
     * reset to 0. If 'rowMaj'==true, the rows are added (after squaring).
     * If 'wvec' is given, col/row i is squared and multiplied by 'wvec[i]'
     * before being added up.
     * <p>
     * NOTE: Special implementation. This matrix can be indexed mask without
     * it being copied here.
     *
     * @param vec    Result
     * @param wvec   Weight vector. Optional
     * @param rowMaj S.a. Def.: false
     * @param increm Optional. Def.: false
     */
    void sumSquares(StVector& vec,bool rowMaj=false,bool increm=false) const;

    void sumSquares(StVector& vec,const StVector& wvec,bool rowMaj=false,
		    bool increm=false) const;

    /**
     * Comp. wise product this matrix and A, then sum of cols. If
     * 'rowMaj'==true, sum of rows. If 'increm'==true, the result is
     * added to 'vec' (not reset to 0).
     * If 'atrans'==true, the comp. wise product with A^T is taken.
     * <p>
     * NOTE: Special implementation. Both this matrix and A can be
     * indexed masks without the matrices being copied.
     *
     * @param a      A
     * @param vec    Result
     * @param atrans S.a. Def.: false
     * @param rowMaj S.a. Def.: false
     * @param increm Optional. Def.: false
     */
    void sumProd(const StMatrix& a,StVector& vec,bool atrans=false,
		 bool rowMaj=false,bool increm=false) const;

    // Special methods

    /**
     * Structure pattern is set to 'normal'.
     */
    void makeSymm(bool low=true) {
      BaseMatrix<double>::makeSymm(low);
      setStrctPatt(normal);
    }

    /**
     * Automatic version of 'makeSymm', uses structure pattern. This matrix
     * must be triangular (!='normal'). It is converted into full
     * symmetric ('normal'). If the structure has a unit diagonal and
     * 'dvec' is given, 'dvec' is used as diagonal.
     * NOTE: Uses structure pattern.
     *
     * @param dvec S.a. Optional
     */
    virtual void makeSymmAuto();

    virtual void makeSymmAuto(const StVector& dvec) {
      if (dvec.size()!=m) throw DimMismatchException(EXCEPT_MSG("dvec"));
      bool docp=(strpatt==uppNDg || strpatt==lowNDg);
      makeSymmAuto();
      if (docp) ((StVector&) diag())=dvec;
    }

    /**
     * Replaces this square matrix A by sym(A) = (1/2) (A+A^T).
     */
    void sym();

    /**
     * Computes maximum symmetric relative difference between A (this) and
     * B ('b'), for each component
     *   max( |(b-a)/a|, |(a-b)/b| ),
     * where the expressions are 0 for denom. 0. Returns maximum of all these.
     * If 'pos' is given, the pos. for the first maximum is returned.  Here,
     * (i,j) is coded as i+j*m.
     * Depending on the str. pattern, we loop over only the rel. part of
     * the matrix. The str. pattern of this matrix det. the part, the patt.
     * of 'b' is ignored.
     * ==> Cannot be used if this and 'b' have diff. patterns with the same
     *     size!
     * NOTE: Uses structure pattern.
     *
     * @param b   Source matrix
     * @param pos S.a. Def.: 0
     * @return    S.a.
     */
    virtual double maxRelDiff(const BaseMatrix<double>& b,int* pos=0) const {
      double mxel=0.0,aelem,belem,temp;
      int i,j,mxpos,fst,lst;

      if (b.rows()!=m || b.cols()!=n)
	throw WrongDimensionException(EXCEPT_MSG(""));
      updateStrctPatt();
      for (j=0; j<n; j++) {
	switch (strpatt) {
	case normal:
	  fst=0; lst=m-1;
	  break;
	case lower:
	  fst=j; lst=m-1;
	  break;
	case lowNDg:
	  fst=j+1; lst=m-1;
	  break;
	case upper:
	  fst=0; lst=j;
	  break;
	case uppNDg:
	  fst=0; lst=j-1;
	  break;
	}
	for (i=fst; i<=lst; i++) {
	  belem=b.get(i,j); aelem=get(i,j);
	  if (aelem!=0.0)
	    temp=fabs((belem-aelem)/aelem);
	  else
	    temp=0.0;
	  if (belem!=0.0)
	    temp=std::max(temp,fabs((aelem-belem)/belem));
	  if (temp>mxel) {
	    mxel=temp; mxpos=i+j*m;
	  }
	}
      }

      if (pos!=0) *pos=mxpos;
      return mxel;
    }

    // Initialization methods (public)

    /**
     * Initializes matrix with zeros. Called with one argument, the method
     * will produce a square matrix. With no arguments, the original
     * dimensions are kept. If one of the arguments is 0, the matrix is set
     * to the empty matrix.
     *
     * @param rows Optional. Number of rows
     * @param cols Optional. Number of columns
     */
    void zeros() {
      fill(m,n,0.0);
    }

    void zeros(int rows) {
      fill(rows,rows,0.0);
    }

    void zeros(int rows,int cols) {
      fill(rows,cols,0.0);
    }

    /**
     * Initializes matrix with ones. Called with one argument, the method
     * will produce a square matrix. With no arguments, the original
     * dimensions are kept. If one of the arguments is 0, the matrix is set
     * to the empty matrix.
     *
     * @param rows Optional. Number of rows
     * @param cols Optional. Number of columns
     */
    void ones() {
      fill(m,n,1.0);
    }

    void ones(int rows) {
      fill(rows,rows,1.0);
    }

    void ones(int rows,int cols) {
      fill(rows,cols,1.0);
    }

    /**
     * Initializes matrix with identity matrix. With no arguments, the
     * original dimensions are kept.
     * If 'rows'!='cols', I occupies upper left, rest is zeros.
     *
     * @param rows Optional. Number of rows
     * @param cols Optional. Number of columns
     */
    void eye() {
      eye(m,n);
    }

    void eye(int rows) {
      eye(rows,rows);
    }

    void eye(int rows,int cols) {
      zeros(rows,cols);
      ((StVector&) diag())=1.0;
    }

    // IO methods (public)

    /**
     * Loads matrix from (binary) input file stream. The current format is the
     * same as for superclass, but can also load old formats.
     *
     * @param is    Binary input file stream
     * @param noTag Do not read tag? Def.: false
     */
    ifstream& load(ifstream& is,bool noTag=false);

    // Matrix decompositions (LAPACK)

    /**
     * Cholesky decomposition (symmetric pos. def. matrix)
     * Uses LAPACK routine DPOTRF.
     * <p>
     * M must be symmetric ('upper'/'lower'). For 'upper':
     * M = U^T U, U upper triangular. For 'lower': M = L L^T, L
     * lower triangular. Upper/lower triangle of M overwritten by factor.
     * If 'clear'==true, the opposite triangle of this matrix is set to 0,
     * oth. it is not accessed.
     * NOTE: Uses structure pattern.
     *
     * @param clear S.a. Def.: false
     * @return      Return value of DPOTRF (INFO; 0 -> OK)
     */
    int cholDecomp(bool clear=false);

    /**
     * Replaces Cholesky decomp. of matrix A by inverse A^-1.
     * Uses LAPACK routine DPOTRI.
     * As opposed to 'invCholInPlace', A^-1 is written into all of this
     * matrix. This matrix must have str. patt. 'lower' / 'upper'
     * dep. on which factor is stored. The opposite triangle is overwritten
     * here. 'strpatt' is not changed here, but setting it to
     * 'normal' will result in the full inverse.
     * NOTE: Uses structure pattern (but writes into opp. triangle!)
     *
     * @return Return of DPOTRI (INFO? If so, 0 -> all OK)
     */
    int invChol() {
      int ret=invCholInPlace();
      uchar osp=strpatt;
      makeSymmAuto();
      strpatt=osp;

      return ret;
    }

    /**
     * Replaces Cholesky decomp. of matrix A by the inverse A^-1.
     * Uses LAPACK routine DPOTRI.
     * <p>
     * Replaces old version of 'invCholInPlace'. This matrix must have
     * str. patt. 'lower' / 'upper' indic. where the factor is
     * stored and whether it is upper/lower tr. As opposed to old version,
     * the diagonal is used. The opposite strict triangle of this matrix is
     * not accessed. 'strpatt' is not changed.
     * NOTE: Eventually, this should replace 'invChol'!
     * <p>
     * @return Return of DPOTRI (INFO? If so, 0 -> all OK)
     */
    int invCholInPlace();

    /**
     * Calculates log of determinant of matrix A whose Cholesky decomp.
     * is stored in this matrix (lower or upper triangle).
     *
     * @return  Log of determinant
     */
    double logDetChol() const {
      Handle<StVector> dgvec(DYNCAST(StVector,diagInt(0,RangeFull::get())));
      return 2.0*dgvec->logDet();
    }

    /** DOWNW. COMPAT.
     * Computes diagonale of the inverse of A, while the Cholesky decomp.
     * of A is stored in this matrix (in lower triangle; the diagonale is
     * given in 'p'). The matrix is not changed. This method is twice as fast
     * as using 'invChol' (if you only need the diagonal).
     * NOTE: Uses own code, not LAPACK.
     *
     * @param p     Main diagonale of Cholesky decomposition
     * @param diagV Diagonal ret. here
     */
    void diagInvChol(const StVector& p,StVector& diagV) const;

    /**
     * Computes diagonale of the inverse of A, while the Cholesky decomp.
     * of A is stored in this matrix. The matrix is not changed. This method
     * is twice as fast as using 'invChol' (if you only need the diagonal).
     * <p>
     * ATTENTION: Implemented only for the case of str. patt. 'lower'!
     * NOTE: Uses own code, not LAPACK.
     *
     * @param diagV Diagonal ret. here
     */
    void diagInvChol(StVector& diagV) const {
      if (strpatt==upper)
	throw NotImplemException(EXCEPT_MSG("StMatrix::diagInvChol implemented only for lower-triangular factors!"));
      else if (strpatt!=lower)
	throw WrongStatusException(EXCEPT_MSG("strpatt"));
      Handle<StVector> dgvec(DYNCAST(StVector,diagInt(0,RangeFull::get())));
      diagInvChol(*dgvec,diagV);
    }

    /*
     * LU decomposition done by LAPACK routine DGETRF:
     *   M = P L U
     * M overwritten by L, U (diag(L) = I, not stored).
     * P permut. matrix which is encoded in IPIV return arg. IPIV
     * describes row exchanges done on M:
     *   0 <-> IPIV[0]-1, then 1 <-> IPIV[1]-1, then ...
     * Why the -1? Because Fortran pos. indexes are 1-based, not
     * 0-based. We keep it that way, s.t. we can pass IPIV to
     * Fortran routines downstream (e.g., DGETRS in 'luSolve').
     * When solving systems, the same row exchanges have to be done
     * on the r.h.s. (left-mult. with P^T) before backsubst. This
     * means that P^T b is computed from b by:
     * - i=0,1,...,n-1:
     *   - exchange pos i <--> IPIV[i]-1
     * P b is computed from b by:
     * - i=n-1,n-2,...,0:
     *   - exchange pos i <--> IPIV[i]-1
     * The parity of the permutation is the parity of the number of
     * i s.t. IPIV[i]-1 != i.
     */

    /**
     * LU decomposition of nonsingular (square) matrix,
     * M = P L U, M==this, P permut. matrix, L lower triangular with
     * unit diag., U upper triangular.
     * <p>
     * The result overwrites M, L in lower triangle, U in upper
     * triangle and diagonal. P returned via index vector 'pind', which
     * stores the IPIV array returned by LAPACK DGETRF (see comments about
     * LU repres. above. IPIV is 1-based!).
     * <p>
     * Uses LAPACK DGETRF.
     *
     * @param pind Permut. index IPIV for P ret. here
     * @return     Return of DGETRF (INFO; 0 -> all OK)
     */
    int luDecomp(BaseVector<int>& pind);

    /**
     * Computes log|M| (log of determinant), where M = P L U (LU
     * decomp.) as a result of 'luDecomp'.
     * IPIV=='pind' encodes P, as returned by 'luDecomp'.
     *
     * @param  pind S.a.
     * @param  sign Sign of det. ret. here
     * @return log(abs(|M|))
     */
    double logDetLU(const BaseVector<int>& pind,double& sign) const;

    /**
     * Computes diagonal of the inverse of M = P^T L U (result of
     * 'luDecomp'). P given by IPIV=='pind'.
     * <p>
     * NOTE: Present implementation not efficient, M^-1 is essentially
     * computed. Replace if critical!
     *
     * @param pind  S.a.
     * @param diagV diag M^-1 ret. here
     */
    void diagInvLU(const BaseVector<int>& pind,StVector& diagV);

    /**
     * Solves M x = b, where M = P L U (LU decomp.) as a result of
     * 'luDecomp'. If 'trans'==true, M^T x = b is solved instead.
     * The permut. matrix is given by the index vector 'pind' (see
     * 'luDecomp').
     * <p>
     * Uses LAPACK DGETRS.
     * NOTE: 'b' and 'x' may be the same vector. 'b' can be matrix in
     * which case M X = B (or M^T X = B).
     *
     * @param pind  S.a.
     * @param b     Vector/matrix b
     * @param x     Solution vector/matrix
     * @param trans S.a. Def.: false
     */
    void luSolve(const BaseVector<int>& pind,const StVector& b,
		 StVector& x,bool trans=false) const {
      if (m!=n || m!=b.size() || m!=pind.size())
	throw DimMismatchException(EXCEPT_MSG(""));
      checkTS(TSARG); pind.checkTS(TSARG); b.checkTS(TSARG); x.checkTS(TSARG);
      x.ensureCapacity(n);
      StMatrix bmat,xmat;
      bmat.reassign(b); xmat.reassign(x); // sell as matrices
      luSolve(pind,bmat,xmat,trans);
    }

    void luSolve(const BaseVector<int>& pind,const StMatrix& b,
		 StMatrix& x,bool trans=false) const;

    /**
     * Compute all eigenvalues of a real symmetric matrix, but none of
     * the corr. eigenvectors.
     * Uses LAPACK driver routine DSYEVR (uses RRR driver).
     * <p>
     * A must be symmetric, pattern 'upper' / 'lower'.
     * ATTENTION: The triangle from which A is read, is overwritten by
     * undefined entries (incl. diagonal). The opposite triangle is not
     * accessed.
     * NOTE: If only a small number of eigenvalues are required, or
     * only eigenvalues in a certain interval, other LAPACK routines (or
     * DSYEVR with other arguments) should be called!
     * NOTE: The LAPACK routine needs some working space. We first perform
     * a query to find out how much, then allocate the space locally
     * (always much less than space for A).
     *
     * @param eigs Eigenvalues in ascending order
     * @return     INFO return value of DSYEVR (0: all OK)
     */
    int eigenVals(StVector& eigs);

    /**
     * Complete eigendecomposition of real symmetric matrix,
     *   A = U D U^T, U orthonormal, D diagonal.
     * Uses LAPACK driver routine DSYEVR (uses the RRR driver).
     * <p>
     * A must be symmetric, pattern 'upper' / 'lower'.
     * ATTENTION: The triangle from which A is read, is overwritten by
     * undefined entries (incl. diagonal). The opposite triangle is not
     * accessed.
     * NOTE: If only a small number of eigenvalues are required, or
     * only eigenvalues in a certain interval, other LAPACK routines (or
     * DSYEVR with other arguments) should be called!
     * NOTE: The LAPACK routine needs some working space. We first perform
     * a query to find out how much, then allocate the space locally
     * (always much less than space for A).
     * <p>
     * Guarantees:
     * - Eigenvalues in D sorted in asc. order
     * - For every col. of U, the first nonzero element is positive
     * NOTE: The cols of U are not uniquely determined if there are
     * duplicate eigenvalues.
     *
     * @param dvec Eigenvalues D in ascending order
     * @param umat U matrix
     * @return     INFO return value of DSYEVR (0: all OK)
     */
    int eigenDecomp(StVector& dvec,StMatrix& umat);

    /*
     * Update/downdate of Cholesky decomposition after rank one modification
     * or matrix extension.
     *
     * New code based on LINPACK routines dchud, dchdd, dchex. More stable
     * than old code. We do not call the LINPACK routines (except see
     * 'dchud', 'dchdd'), since they do not properly use BLAS-I.
     *
     * Dragging along matrices:
     * Dragging along means updating Z -> Z', such that
     *   Z' (L')^T = Z L^T + y v^T    [if A' = A + v v^T]
     *   Z' (L')^T = Z L^T - y v^T    [if A' = A - v v^T]
     * We generally work with columns of Z (using BLAS-I routines).
     * ATTENTION: Different from old code, which was more flexible! This is
     * in line with LINPACK dchud, dchdd, and is probably sufficient for all
     * applications.
     */

    /**
     * Update Cholesky decomp. A = L L^T after rank 1 modification:
     *   A' = A + v v^T
     * This matrix contains L in lower triangle (str. patt. 'lower') or
     * L^T in upper triangle (str. patt. 'upper'). Is overwritten by L'
     * (or L'^T), opposite triangle not accessed.
     * <p>
     * Adapted from LINPACK dchud, uses n Givens rotations with angles given
     * by c_k, s_k. These are ret. in 'cvec', 'svec' (size n each).
     * <p>
     * If Z is given in 'dragZ' (r-by-n), the r-vector y must be given as
     * well in 'dragY', and Z is overwritten by Z':
     *   Z' (L')^T = Z L^T + y v^T
     * We need a flat scratch vector of size max(n,r), where r=0 if Z is not
     * given. Passed in 'wkvec'.
     * NOTE: Can pass same ref. for 'vvec' and 'wkvec'. v is overwritten in
     * this case, but the method works.
     * <p>
     * Meaning of 'cvec', 'svec' components:
     * These are the angle values c_k, s_k as in LINPACK dchud, except that
     * if for the LIN values, we'd obtain L'_{kk} < 0, we flip signs of c_k,
     * s_k. Details are in
     *   M. Seeger
     *   Low Rank Updates for the Cholesky Decomposition
     * NOTE: The c_k, s_k def. in the TR is slightly different. See code
     * comments for this point.
     *
     * @param vvec  Vector v
     * @param cvec  Angle cosines ret. here (size n)
     * @param svec  Angle sines ret. here (size n)
     * @param wkvec Scratch vector
     * @param dragZ Matrix Z, s.a. Overwritten by Z'
     * @param dragY Vector y, s.a. Iff 'dragZ' is given
     */
    void cholUpdRk1(const StVector& vvec,StVector& cvec,StVector& svec,
		    StVector& wkvec,StMatrix* dragZ=0,const StVector* dragY=0);

    /**
     * Downdate Cholesky decomp. A = L L^T after rank 1 modification:
     *   A' = A - v v^T
     * This matrix contains L in lower triangle (str. patt. 'lower') or
     * L^T in upper triangle (str. patt. 'upper'). Is overwritten by L',
     * opposite triangle not accessed.
     * <p>
     * Adapted from LINPACK dchdd, uses n Givens rotations with angles given
     * by c_k, s_k. These are ret. in 'cvec', 'svec' (size n each).
     * We need the vector p s.t. L p = v. If 'isp'==true, p is supplied in
     * 'vvec' instead of v. Otherwise, it is computed here.
     * The routine fails if the new A' is not pos. definite, i.e. if
     *   1 - p^T p <= 0
     * In this case, we return 1, otherwise 0 (success).
     * <p>
     * If Z is given in 'dragZ' (r-by-n), the r-vector y must be given as
     * well in 'dragY', and Z is overwritten by Z':
     *   Z' (L')^T = Z L^T - y v^T
     * <p>
     * We need a flat scratch vector of size max(n,r), where r=0 if Z is not
     * given. Passed in 'wkvec'.
     * NOTE: Can pass same ref. for 'vvec' and 'wkvec'. v is overwritten in
     * this case, but the method works.
     * <p>
     * Meaning of 'cvec', 'svec' components:
     * See TR
     *   M. Seeger
     *   Low Rank Updates for the Cholesky Decomposition
     * for details. The angle values c_k, s_k parameterize the change
     * L -> L'. In LINPACK dchdd, diag(L') can have negative values, and
     * (as opposed to dchud, 'cholUpdRk1') this cannot be compensated by
     * flipping signs of angle values. In general, we have
     *  L' = LL D,
     * where LL is L' as computed by dchdd, and D is diagonal with entries
     * -1,+1, making diag(L') positive. The indices of the -1 components in
     * D are returned in 'flind' (sorted in asc. order).
     *
     * @param vvec  Vector v (or p, if 'isp'==true)
     * @param cvec  Angle cosines ret. here (size n)
     * @param svec  Angle sines ret. here (size n)
     * @param flind S.a.
     * @param wkvec Scratch vector
     * @param isp   S.a. Def.: false
     * @param dragZ Matrix Z, s.a. Overwritten by Z'
     * @param dragY Vector y, s.a. Iff 'dragZ' is given
     * @return      0: success; 1: numerical failure
     */
    int cholDndRk1(const StVector& vvec,StVector& cvec,StVector& svec,
		   BaseVector<int>& flind,StVector& wkvec,bool isp=false,
		   StMatrix* dragZ=0,const StVector* dragY=0);

    /**
     * Updates Cholesky decomp. A = L L^T (A = R^T R), if A is extended by
     * a new row/column 'avec', 'ascal'. This matrix contains L (R) which
     * must be lower (upper) triangular, str. patt. 'lower'
     * ('upper'). 'avec'==a, or 'avec'==L^-1 a (R^-T a) if 'isl'==true.
     * Otherwise, the backsubst. is done locally.
     * Returns 0 if successful, 1 if numerical error.
     * <p>
     * Dragging along: see general comments above (replace L by R^T).
     * ATTENTION: If 'isl'==true, 'avec' must NOT share memory with any of
     * the matrices managed in 'dragX<k>'!
     *
     * @param avec   New column of A (except last elem.)
     * @param ascal  Last elem. of new column of A
     * @param isl    S.a. Def.: false
     * @param dragX1 S.a. Def.: 0
     * @param dragB1 S.a. Needed iff 'dragX1' given (vectors)
     * @param dragX2 S.a. Def.: 0
     * @param dragB2 S.a. Needed iff 'dragX2' given (vectors)
     * @return       0: Success. 1: Numerical error
     */
    int cholext(const StVector& avec,double ascal,bool isl=false,
		ArrayHandle<StMatrix*>& dragX1=
		ArrayHandleZero<StMatrix*>::get(),
		ArrayHandle<StVector*>& dragB1=
		ArrayHandleZero<StVector*>::get(),
		ArrayHandle<StMatrix*>& dragX2=
		ArrayHandleZero<StMatrix*>::get(),
		ArrayHandle<StVector*>& dragB2=
		ArrayHandleZero<StVector*>::get());

    /**
     * Variant in which A grows by several rows/columns.
     *
     * @param avecs  New columns A_{.,*} of A
     * @param ascals New diag. block A_{*,*} of A. Must be symmetric
     * @param isl    S.a. Def.: false
     * @param dragX1 S.a. Def.: 0
     * @param dragB1 S.a. Needed iff 'dragX1' given (matrices)
     * @param dragX2 S.a. Def.: 0
     * @param dragB2 S.a. Needed iff 'dragX2' given (matrices)
     */
    int cholext(const StMatrix& avecs,const StMatrix& ascals,bool isl=false,
		ArrayHandle<StMatrix*>& dragX1=
		ArrayHandleZero<StMatrix*>::get(),
		ArrayHandle<StMatrix*>& dragB1=
		ArrayHandleZero<StMatrix*>::get(),
		ArrayHandle<StMatrix*>& dragX2=
		ArrayHandleZero<StMatrix*>::get(),
		ArrayHandle<StMatrix*>& dragB2=
		ArrayHandleZero<StMatrix*>::get());

    /**
     * NOTE: Direct LINPACK variant of 'cholUpdRk1'. Ours should be faster.
     * ATTENTION: Can produce negative values on diag(R)!
     * <p>
     * Direct interface to LINPACK routine DCHUD. This is basically doing
     * the job of 'cholUpdRk1' for positive 'alpha', using official code.
     * It is somewhat less general.
     * TODO: Shift own code to use this official code!
     * <p>
     * Updates Cholesky decomp. A = R^T R, R upper triangular:
     *   A' = A + x x^T
     * This matrix contains R in upper triangle, str. patt. 'upper'. Lower
     * triangle is not touched. It is replaced by R'.
     * If 'zmat' is given, it contains a n-by-r matrix Z. In this case, the
     * r-vector y is also required in 'yvec'. Z is then updated to Z', s.t.:
     *   R^T Z + x y^T = (R')^T Z'
     * This is a restricted form of "dragging along".
     * <p>
     * NOTE: DCHUD does not provide for error code return. It seems this is
     * because it is stable, given all arguments are as specified.
     * Needs two working arrays 'carr', 'sarr' of size n each. The Givens
     * rotation angles are stored and returned there. If not given, they
     * are alloc. locally.
     *
     * @param xvec Vector x
     * @param zmat Pointer to Z matrix. Def.: 0
     * @param yvec Pointer to y vector. Req. iff 'zmat' given
     * @param carr Working array. Cosines ret. there. Optional
     * @param sarr Working array. Sines ret. there. Optional
     */
    //void dchud(const StVector& xvec,StMatrix* zmat=0,const StVector* yvec=0,
    //       ArrayHandle<double>& carr=ArrayHandleZero<double>::get(),
    //       ArrayHandle<double>& sarr=ArrayHandleZero<double>::get());

    /**
     * NOTE: Direct LINPACK variant of 'cholDndRk1'. Ours should be faster.
     * ATTENTION: Can produce negative values on diag(R)!
     * <p>
     * Direct interface to LINPACK routine DCHDD. This is basically doing
     * the job of 'cholUpdRk1' for negative 'alpha', using official code.
     * It replaces 'cholUpdRk1Negative', which is own code. It is somewhat
     * less general.
     * TODO: Shift own code to use this official code!
     * <p>
     * Downdates Cholesky decomp. A = R^T R, R upper triangular:
     *   A' = A - x x^T
     * This matrix contains R in upper triangle, str. patt. 'upper'. Lower
     * triangle is not touched. It is replaced by R'.
     * If 'zmat' is given, it contains a n-by-r matrix Z. In this case, the
     * r-vector y is also required in 'yvec'. Z is then downdated to Z', s.t.:
     *   R^T Z - x y^T = (R')^T Z'
     * This is a restricted form of "dragging along".
     * <p>
     * This can go wrong, if the resulting R' is not numerically positive
     * definite. We return INFO of DCHDD: 0 for success, -1 for failure, +1
     * if the downdating of Z failed.
     *
     * @param xvec Vector x
     * @param zmat Pointer to Z matrix. Def.: 0
     * @param yvec Pointer to y vector. Req. iff 'zmat' given
     * @param carr Working array. Cosines ret. there. Optional
     * @param sarr Working array. Sines ret. there. Optional
     * @return     INFO ret. value. 0: success
     */
    //int dchdd(const StVector& xvec,StMatrix* zmat=0,const StVector* yvec=0,
    //      ArrayHandle<double>& carr=ArrayHandleZero<double>::get(),
    //      ArrayHandle<double>& sarr=ArrayHandleZero<double>::get());

    /** OLD CODE. REPLACED BY 'cholUpdRk1', 'cholDndRk1'
     * Update Cholesky decomp. A = L L^T, rank 1 modification:
     *   A' = A + 'alpha' v v^T
     * This matrix contains L in lower triangle, str. patt. 'lower'.
     * Overwritten by L', upper tri. not touched.
     * Req. p = L^-1 v. v=='vec', or p=='vec' if 'isp'==true.
     * Returns 0 if OK, 1 if numerical breakdown.
     * <p>
     * Dragging along: see general comments above.
     * ATTENTION: If 'isp'==true, 'vec' must NOT share memory with any of
     * the matrices managed in 'dragX<k>'!
     * <p>
     * Algorithms:
     * - 'alpha'>0 (positive): Direct update of L, method C1
     *   in Gill, Golub, Murray, Saunders [1]
     * - 'alpha'<0 (negative): Method C3 in [1], proposed by
     *   Stewart, part of LINPACK. This is more stable than the direct
     *   update for 'alpha'<0
     * <p>
     * Workspace:
     * Needs vector of size n (flat), two arrays of size n. Can be supplied,
     * are alloc. here otherwise.
     *
     * @param vec      See above
     * @param alpha    S.a. Def.: 1
     * @param isp      S.a. Def.: false
     * @param dragX1   S.a. Def.: 0
     * @param dragX2   S.a. Def.: 0
     * @param wkvec    S.a. Def.: 0
     * @param wkarrn1  S.a. Def.: 0
     * @param wkarrn2  S.a. Def.: 0
     * @return         0: OK, 1: Numerical error
     */
    int cholUpdRk1_old(const StVector& vec,double alpha,bool isp,
		   ArrayHandle<StMatrix*>& dragX1=
		   ArrayHandleZero<StMatrix*>::get(),
		   ArrayHandle<StMatrix*>& dragX2=
		   ArrayHandleZero<StMatrix*>::get(),StVector* wkvec=0,
		   ArrayHandle<double>& wkarrn1=ArrayHandleZero<double>::get(),
		   ArrayHandle<double>& wkarrn2=
		   ArrayHandleZero<double>::get());

  private:
    /** OLD CODE. REPLACED BY 'cholUpdRk1', 'cholDndRk1'
     * Helper for 'cholUpdRk1' for negative 'alpha', implements stable
     * downdating of CF as analyzed by Stewart (part of LINPACK).
     */
    int cholUpdRk1Negative_old(const StVector& vec,double alpha,bool isp,
			   ArrayHandle<StMatrix*>& dragX1,
			   ArrayHandle<StMatrix*>& dragX2,StVector* wkvec,
			   ArrayHandle<double>& wkarrn1,
			   ArrayHandle<double>& wkarrn2);
  public:

    // TODO: Kick this out, just call 'cholUpdRk1' multiple times!
    /** DOWNW. COMPAT.!
     * ATTENTION: Case 'alpha'<0 unstable! Use multiple calls of (negative)
     * 'cholUpdRk1'!
     * <p>
     * Variant of 'cholUpdRk1' for rank d updates. 'vmat' is n-by-d
     * matrix V, and A' = A + 'alpha' V M V^T.
     * Here M is d-by-d positive definite. If 'cholM' is given, it is the
     * Cholesky factor of M, otherwise M==I.
     * We need P = L^-1 V which is passed via 'vmat' if 'isp'==true.
     * <p>
     * If 'bMat' is given, we need a buffer matrix of size m-by-d, m the
     * largest number of rows of all 'bMat' matrices. This matrix can be
     * passed via 'buffMat' (no mask!), otherwise is alloc. locally.
     * <p>
     * ATTENTION: If 'isp'==true, 'vmat' must NOT share memory with any of
     * the matrices in 'bMat' (mutual masking, ...)!
     *
     * @param vmat    Matrix V or P
     * @param alpha   S.a.
     * @param isp     S.a. Def.: false
     * @param cholM   S.a. Def.: M==I
     * @param bMat    Optional. S.a.
     * @param buffMat S.a. Optional
     * @return        Error code. 0: OK, 1: Numerical error   
     */
    /*
    int cholUpdRkd(const StMatrix& vmat,double alpha,bool isp,
		   const StMatrix* cholM,
		   ArrayHandle<StMatrix*>& bMat,StMatrix* buffMat=0);

    int cholUpdRkd(const StMatrix& vmat,double alpha=1.0,bool isp=false,
		   const StMatrix* cholM=0,StMatrix* bMat=0,
		   StMatrix* buffMat=0) {
      ArrayHandle<StMatrix*> arrb; // empty
      if (bMat!=0) {
	arrb.changeRep(1);
	arrb[0]=bMat;
      }
      return cholUpdRkd(vmat,alpha,isp,cholM,arrb,buffMat);
    }

    int chollrup(const StMatrix& vmat,double alpha=1.0,bool isp=false,
		 ArrayHandle<StMatrix*>& mmat=
		 ArrayHandleZero<StMatrix*>::get()) {
      return cholUpdRkd(vmat,alpha,isp,0,mmat);
    }
    */

    /** DOES NOT WORK PROPERLY!
     * Update of Cholesky decomp. A = L L^T after rank 2 modification
     *   A' = (I + v u^T) A (I + u v^T)
     * This matrix must contain L in lower triangle (str. patt. 'lower').
     * It is overwritten by L' (upper triangle not accessed).
     * <p>
     * Requires w = L^T u, z = L^-1 v. 'uvec' is u or w if 'isw'==true. 'vvec'
     * is v or z if 'isz'==true.
     * Returns 0 if OK, 1 if numerical breakdown.
     * <p>
     * Dragging along: see general comments above.
     * ATTENTION: If 'isw'/'isz'==true, 'uvec'/'vvec' must NOT share memory
     * with any of the matrices managed in 'dragX<k>'!
     * <p>
     * Algorithm:
     *   Goldfarb: Factorized Variable Metric Methods for Unconstrained
     *   Optimization, Math. Comp. 30(136), 1976, 796--811
     * We implement method 1, based on Givens rotations.
     *
     * @param uvec   S.a.
     * @param vvec   S.a.
     * @param isw    S.a. Def.: false
     * @param isz    S.a. Def.: false
     * @param dragX1 S.a. Def.: 0
     * @param dragX2 S.a. Def.: 0
     * @return       0: OK, 1: Numerical error   
     */
    /*
    int cholUpdRk2Goldfarb(const StVector& uvec,const StVector& vvec,
			   bool isw=false,bool isz=false,
			   ArrayHandle<StMatrix*>& dragX1=
			   ArrayHandleZero<StMatrix*>::get(),
			   ArrayHandle<StMatrix*>& dragX2=
			   ArrayHandleZero<StMatrix*>::get());
    */

    /**
     * Update Cholesky decomp. A = R^T R, if
     *   A' = E^T A E, E special permutation matrix
     * This matrix contains R in upper triangle, str. patt. 'upper'.
     * Overwritten by R', lower tri. not touched.
     * E is defined by k=='k', l=='l', 'job'. Here, 0<=k<l<n. E reorders
     * rows/columns as follows:
     * If 'job'==1, the new ordering is ...,k-1,l,k,...,l-1,l+1,...
     * If 'job'==2, the new ordering is ...,k-1,k+1,...,l,k,l+1,...
     * The method works by computing an orthonormal U s.t. U R E = R'.
     * U is the product of l-k Givens rotations.
     * <p>
     * Dragging along: different from general comments above. If
     * X=='dragX' is given, it must have n rows. It is replaced by U X.
     * If R^T X = B, then (R')^T X' = E^T B.
     * NOTE: To move row/column i of A to the end, use k==i, l==n-1, and
     * job==2.
     * <p>
     * Uses LINPACK routine DCHEX.
     * NOTE: There is a bug in LINPACK DCHEX which is worked around here
     * (see code for details).
     *
     * @param k     S.a.
     * @param l     S.a.
     * @param job   S.a.
     * @param dragX S.a. Def.: 0
     */
    void cholUpdExch(int k,int l,int job,StMatrix* dragX=0);

    /*
     * Update Cholesky decomposition, older code
     * ATTENTION: All low-rank Cholesky methods have been modified to use
     * matrices with long columns rather than long rows.
     * ==> CHANGE OLD CODE!!
     */

    /**
     * Computes Cholesky factorization of D + A*A^T, where D is a diagonal
     * matrix with positive entries, given in 'diag', and A is this matrix.
     * Use this method only if n is subst. smaller than m. The Cholesky
     * factorization can be stored as:
     *   L^(1) ... L^(n) D^(n) L^(n)^T ... L^(1)^T.
     * Here, L^(j) is represented by two vectors p^(j) and \beta^(j), both
     * having m-1 entries, and is lower triangular with diagonal I. D^(n)
     * is diagonal with positive entries. The conventional CF is:
     *   L L^T, L = L^(1) ... L^(n) (D^(n))^{1/2}
     * Thus, the factorization can be stored in O(m n), and things like mult.
     * with L or backsubstitutions take O(m n) as well.
     * <p>
     * NOTE: In comparison with using the Sherman-Morrison formulae in this
     * context, we have a storage overhead of 2, yet using a Cholesky factor.
     * is much more stable!
     * <p>
     * The p^(j),\beta^(j) are returned in the matrix 'lFact' which is
     * 2(m-1)-by-n, namely col. j stores the concatenation \beta^(j)' p^(j)'.
     * 'diag' is overwritten by D^{(m) 1/2}.
     * NOTE: In the documentation to this routine, the p and \beta vectors
     * are n-vectors. However, the first element of each p and the last element
     * of each \beta vector are never used later and can be dropped.
     * <p>
     * NOTE: If n==1, A and the result can be repres. as vectors. To this end,
     * create a result vector of length 2(m-1), mask the vector for A and the
     * res. vector by matrices with one row. This method works with a mask
     * as long as its size is correct.
     *
     * @param diag  Diagonal matrix D. Overwritten by D^{(n) 1/2}
     * @param lFact L^(j) factors returned here (see above)
     */
    void cholDecLowRank(StVector& diag,StMatrix& lFact) const;

    /**
     * Same as 'cholDecLowRank', but here we compute the CF of D - A*A^T.
     * This works only if D - A*A^T is positive definite.
     *
     * @param diag  Diagonal matrix D. Overwritten by D^{(n) 1/2}
     * @param lFact L^(j) factors returned here
     */
    void cholDecLowRankMinus(StVector& diag,StMatrix& lFact) const;

    /**
     * Cholesky backsubstition for 'cholDecLowRank'/'cholDecLowRankMinus'
     * <p>
     * The factorization is L L', L = L^(1) ... L^(n) (D^(n))^{1/2}.
     * The L^(j) are given in this matrix (format of 'lFact', see
     * 'cholDecLowRank'), D^{(n) 1/2} in 'diag'. We solve L x = b for x, where
     * b is given in 'b', and the solution is returned via 'x'. 'b' and 'x'
     * may be the same vector variable.
     * 'b' and 'x' may also be matrices, in which case backsubstitution is
     * applied to every column of 'b', or to every row if 'trans'==true.
     *
     * @param diag  Diagonal matrix D^{(n) 1/2}
     * @param b     Vector/matrix b
     * @param x     Solution ret. here
     * @param trans S.a. Def.: false (DOWNW. COMPAT.!)
     */
    void cholLRBacksubst(const StVector& diag,const StVector& b,
			 StVector& x) const;
    void cholLRBacksubst(const StVector& diag,const StMatrix& b,
			 StMatrix& x,bool trans=false) const;

    /**
     * Cholesky backsubstition for 'cholDecLowRank'/'cholDecLowRankMinus',
     * transposed factor
     * <p>
     * The factorization is L L^T, L = L^(1) ... L^(n) (D^(n))^{1/2}.
     * The L^(j) are given in this matrix (format of 'lFact', see
     * 'cholDecLowRank'), D^{(n) 1/2} in 'diag'. We solve L' x = b for x,
     * where b is given in 'b', and the solution is returned via 'x'. 'b' and
     * 'x' may be the same vector variable.
     * 'b' and 'x' may also be matrices, in which case backsubstitution is
     * applied to every column of 'b', or to every row if 'trans'==true.
     *
     * @param diag  Diagonal matrix D^(m)
     * @param b     Vector b
     * @param x     Solution ret. here
     * @param trans S.a. Def.: false (DOWNW. COMPAT.!)
     */
    void cholLRBacksubsT(const StVector& diag,const StVector& b,
			 StVector& x) const;
    void cholLRBacksubsT(const StVector& diag,const StMatrix& b,
			 StMatrix& x,bool trans=false) const;

    /**
     * Updates conventional Chol. factor L by right-multiplication with a
     * LR Cholesky factor \tilde{L}. \tilde{L} is stored in this matrix and
     * 'diag'. L is stored in the lower triangle of 'fact' and 'dfact' which
     * are overwritten by the new factor. The upper triangle and diagonal of
     * 'fact' are not touched. This matrix remains unchanged.
     * NOTE: L and the result are stored in the conventional form as ret.
     * by 'cholDecomp', while \tilde{L} comes from 'cholDecLowRank'.
     *
     * @param diag  Diagonal matrix D for \tilde{L}
     * @param fact  Factor L (lower tri.)
     * @param dfact Factor L (diagonal)
     */
    void cholLRUpdateFact(const StVector& diag,StMatrix& fact,
			  StVector& dfact) const;


    // OLD CODE. DO NOT USE!!
    /*
    int cholDecUpdateOld(const StVector& vec,bool add,bool isp,
			 ArrayHandle<StMatrix*>& bMat);
    */

    /** OLD CODE. DO NOT USE!!!
     * Similar to 'cholDecUpdate'. Here, the decomp. is A = L D L^T,
     * D diagonal and positive, L unit lower triangular. D passed via
     * 'diag'. The update is by s*v*v^T, s=='scal', v=='vec'.
     * If 'scal'<0 and 'safe'==true, use safeguard against breakdown.
     * Namely, if 'scal' is s.t. the resulting matrix is not pos. def.,
     * it is raised by the minimum amount (plus a tolerance) to
     * render positive def.
     * The time req. is O(n^2).
     * We require the solution of L p = v. If p is already available, it
     * can be passed via 'pVec' to save time.
     *
     * @param diag    See above
     * @param vec     See above
     * @param scal    S.a.
     * @param safe    S.a. Def.: false
     * @param pVec    Optional. See above
     */
    /*
    void cholDecUpdateSpecOld(StVector& diag,const StVector& vec,
			      double scal,bool safe=false,
			      const StVector* pVec=0);
    */

    /**
     * Similar to 'cholDecUpdate'. Here, A = L D L^T, L unit lower
     * triangular (this matrix), D diagonal and positive ('diag').
     * The update is A' = A + s v v^T, s=='scal', v=='vec'.
     * We require L p = v. If 'isp'==true, p is passed in 'vec' instead
     * of v. Otherwise, p is computed here. If 'pVec' is given, p is
     * returned there. D' overwrites D in 'diag'.
     * <p>
     * Safeguarding:
     * If 'safe'==true and s is s.t. the res. matrix is not pos. def. (by
     * a pos. tolerance), the method returns with true, not changing any
     * variables except for 'pVec' (if given).
     * <p>
     * Return of L':
     * Normally, L' overwrites L in this matrix. If 'betaVec' is given,
     * this matrix is not changed, but L'' is returned s.t. L' = L L''.
     * In this case, 'pVec' must also be given or 'isp'==true. L'' is
     * represented by 'pVec' and 'betaVec'. The L -> L' update can be
     * done later by 'cholDecUpdateSpecComplete'.
     * <p>
     * Dragging along:
     * If 'mMat' is given, its rows are dragged along. This means 'mMat'
     * is M^T and L M = B for some fixed B. 'mMat' is overwritten
     * by (M')^T s.t. L' M' = B.
     *
     * @param diag    D, overwritten by D'
     * @param vec     Vector v or p (s.a.)
     * @param scal    Scalar s, non-zero
     * @param isp     S.a. Def.: false
     * @param safe    S.a. Def.: false
     * @param pVec    S.a. Def.: 0
     * @param betaVec S.a. Def.: 0
     * @param mMat    S.a. Def.: 0
     * @return        Safeguard abort?
     */
    bool cholDecUpdateSpec(StVector& diag,const StVector& vec,double scal,
			   bool isp=false,bool safe=false,StVector* pVec=0,
			   StVector* betaVec=0,StMatrix* mMat=0);

    /**
     * See 'cholDecUpdateSpec'. This method updates L to L'. L (unit lower
     * triangular) stored in this matrix, overwritten by L'.
     *
     * @param pArr    Vector p as flat array (first elem. not used)
     * @param betaArr Vector beta as flat array (last elem. not used)
     */
    void cholDecUpdateSpecComplete(const ArrayHandle<double>& pArr,
				   const ArrayHandle<double>& betaArr);

    /**
     * See 'cholDecUpdateSpec'. L is updated there to L' = L L'', where
     * L'' is defined in terms of vectors p, beta.
     * If 'mMat'==M^T, this is overwritten by M^T (L'')^-T ('trans'==false)
     * or by M^T (L'')^-1 ('trans'==true). We say that the rows of 'mMat' are
     * dragged along the update of L.
     * NOTE: This is a static method, because L is not needed.
     *
     * @param mMat    Matrix M^T
     * @param pArr    Vector p as flat array (first elem. not used)
     * @param betaArr Vector beta as flat array (last elem. not used)
     * @param trans   S.a. Def.: false
     */
    static void cholDecUpdateSpecDrag(StMatrix& mMat,
				      const ArrayHandle<double>& pArr,
				      const ArrayHandle<double>& betaArr,
				      bool trans=false);

  protected:
    // Internal methods

    /**
     * Resets 'strpatt' to 'normal' if this matrix is not square.
     * Has to be called by methods which use 'strpatt'.
     */
    void updateStrctPatt() const {
      if (m!=n) strpatt=normal;
    }

    void convertWrapped(BaseMatWrapper<double>& arg) {
      StMatrix* mat=DYNCAST(StMatrix,arg.getRep());
      uchar stct=mat->strpatt;
      BaseMatrix<double>::convertWrapped(arg);
      strpatt=stct;
    }

    BaseMatrix<double>* copy() const {
      return new StMatrix(*this);
    }

    BaseMatrix<double>* newEmpty() const {
      return new StMatrix();
    }

    BaseVector<double>* newEmptyVector() const {
      return new StVector();
    }
  };
//ENDNS

#undef TSARG

#endif
