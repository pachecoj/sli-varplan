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
 * Desc.:  Header class FastUtils
 * ------------------------------------------------------------------- */

#ifndef FASTUTILS_H
#define FASTUTILS_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/matrix/default.h"
#include "lhotse/matrix/cblas_for_cpp.h" // BLAS C interface (C linkage)
#include "lhotse/matrix/BaseLinVec.h"
#include "lhotse/matrix/BaseLinMat.h"
#include "lhotse/matrix/WriteBackMat.h"
#include <functional>
//#ifdef MATLAB_DEBUG_OLD
#include "lhotse/MatlabDebug.h"
//#endif // MATLAB_D

//BEGINNS(matrix)

#define CONVBUFFPTR(ptr,len,step) (((step)>0)?(ptr):((ptr)+((1-(len))*(step))))

#define SETTRS(flag) ((flag)?CblasTrans:CblasNoTrans)

#define SETUPLO(strct) (((strct)==WriteBackMat<double>::strctUpper || (strct)==WriteBackMat<double>::strctUppNDg)?CblasUpper:CblasLower)

#define SETDIAG(strct) (((strct)==WriteBackMat<double>::strctUpper || (strct)==WriteBackMat<double>::strctLower)?CblasNonUnit:CblasUnit)

  /**
   * Static inline methods for elementrary arithmetic matrix/vector
   * operations, using the element type T==double. BLAS functions are used
   * whenever possible. The methods are usually overloaded to also allow
   * passing position indices for the arguments (as used with mask objects).
   * <p>
   * Vector arguments:
   * In general, a vector argument is given by a buffer pointer (say: a),
   * a length (say: n), a step size (optional; say: step), a position
   * index (optional; say: aInd). If aInd is given (!=0), element i is at
   * a[aInd[i]]. If aInd is not given, element i is at *(a+i*step) for
   * step>0. For step<0, the Fortran convention is used (see below). The def.
   * value for step is 1. step is ignored if aInd is given.
   * <p>
   * BLAS support:
   * BLAS I is used for vector-vector operations whenever applicable and
   * no position indices are involved.
   * For matrix-vector and matrix-matrix (BLAS II, BLAS III), in general the
   * helper classes 'BaseLinVec', 'BaseLinMat' are used, which draw flat
   * copies for indexed mask matrices.
   * <p>
   * BLAS vector format, Fortran convention:
   * If the step size step>0, then
   *   x_i -> x[i*step]
   * If step<0, then
   *   x_i -> x[(i-n+1)*step]
   * (the step size is called incX in the CBLAS docum.). This is to avoid
   * negative indices into x[], x always points to the beginning of the
   * memory region, but if step<0, *x is the last element of the vector.
   * The methods here use this Fortran convention. We use CONVBUFFPTR to
   * convert x into our convention (where x_i -> x[i*step], even if step<0).
   * NOTE: 'BaseLinVec' objects internally use the Fortran convention as
   * well, while 'ArrayUtils' uses our convention.
   * <p>
   * BLAS matrix format:
   * We do not (in the moment) support banded matrices or triangular matrices
   * in packed format. Although CBLAS leaves the choice, we always use
   * column-major ordering (Fortran convention). A matrix is given by:
   * - buffer pointer a
   * - striding constant stride
   * - size m, n
   * We have:
   *   a_(i,j) -> a[j*stride+i], i<m,j<n
   * Here, we require stride>=m (cannot be negative as for vectors!).
   * Furthermore, there are flags which determine how matrices are to be
   * read from their rectangular array:
   * - aTrans (CBLAS_TRANSPOSE):
   *   CblasNoTrans -> use as A
   *   CblasTrans   -> use as A^T
   * - aUplo (CBLAS_UPLO): For triangular/symmetric matrices only:
   *   CblasUpper -> use upper triangle of buffer
   *   CblasLower -> use lower triangle of buffer
   * - aDiag (CBLAS_DIAG): For triangular matrices only:
   *   CblasNonUnit -> A has general diagonal
   *   CblasUnit    -> A has unit diagonal. Diagonal entries in buffer are
   *                   not used
   * <p>
   * NOTE: Man pages for most BLAS functions can be found at:
   *   http://www.mathkeisan.com/man.html  [copyright NEC]
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class FastUtils
  {
  public:
    // Nonlinear functions

    /**
     * logsumexp(a) = log sum_i exp(a_i)
     * <p>
     * This is computed in a stable way, by first determining max(a_i).
     *
     * @param a    Buffer for a
     * @param n    Length
     * @param step Step size for a
     * @param ind  Position index for a
     * @return     Function value
     */
    static double logsumexp(const double* a,int n,int step=1,const int* ind=0);

    // Vector-vector (BLAS I when applies)

    /**
     * a^T * b (inner product)
     *
     * @param a     Buffer for a
     * @param b     Buffer for b
     * @param n     Length
     * @param aStep Step size for a
     * @param bStep Step size for b
     * @param aInd  Position index for a
     * @param bInd  Position index for b
     * @return      Inner product
     */
    static double inner(const double* a,const double* b,int n,int aStep=1,
			int bStep=1,const int* aInd=0,const int* bInd=0);

    /**
     * a^T * a (inner product)
     *
     * @param a     Buffer for a
     * @param n     Length
     * @param aStep Step size for a
     * @param aInd  Position index for a. 0 -> none
     * @return      Inner product
     */
    static double inner(const double* a,int n,int aStep=1,const int* aInd=0);

    /**
     * a^T * D * b, D diagonal
     *
     * @param a     Buffer for a
     * @param b     Buffer for b
     * @param d     Buffer for d
     * @param n     Length
     * @param aStep Step size for a
     * @param bStep Step size for b
     * @param dStep Step size for d
     * @param aInd  Position index for a. 0 -> none
     * @param bInd  Position index for b. 0 -> none
     * @param dInd  Position index for b. 0 -> none
     * @return      Inner product
     */
    static double mahal(const double* a,const double* b,const double* d,int n,
			int aStep=1,int bStep=1,int dStep=1,const int* aInd=0,
			const int* bInd=0,const int* dInd=0);

    /**
     * a^T * D * a, D diagonal
     *
     * @param a     Buffer for a
     * @param d     Buffer for b
     * @param n     Length
     * @param aStep Step size for a
     * @param dStep Step size for d
     * @param aInd  Position index for a. 0 -> none
     * @param dInd  Position index for b. 0 -> none
     * @return      Inner product
     */
    static double mahal(const double* a,const double* d,int n,int aStep=1,
			int dStep=1);
    static double mahal(const double* a,const double* d,int n,int aStep,
			int dStep,const int* aInd,const int* dInd=0);

    /**
     * a += s*b, s scalar.
     * NOTE: Cases s==1,-1 treated specially.
     *
     * @param a     Vector a
     * @param b     Vector b
     * @param s     Scalar s
     * @param n     Length
     * @param aStep Step size for a
     * @param bStep Step size for b
     * @param aInd  Position index for a. 0 -> none
     * @param bInd  Position index for b. 0 -> none
     */
    static void addsmul(double* a,const double* b,double s,int n,int aStep=1,
			int bStep=1,const int* aInd=0,const int* bInd=0);

    /**
     * c = a + s*b, s a scalar.
     * NOTE: Cases s==1,-1 treated specially!
     * ATTENTION: 'addsmul' a += s*b much faster in practice!
     *
     * @param c     Target vector c
     * @param a     Vector a
     * @param b     Vector b
     * @param s     Scalar s
     * @param n     Length
     * @param cStep Step size for c
     * @param aStep Step size for a
     * @param bStep Step size for b
     * @param cInd  Position index for c. 0 -> none
     * @param aInd  Position index for a. 0 -> none
     * @param bInd  Position index for b. 0 -> none
     */
    static void addsmul(double* c,const double* a,const double* b,double s,
			int n,int cStep=1,int aStep=1,int bStep=1,
			const int* cInd=0,const int* aInd=0,const int* bInd=0);

    /**
     * a += s./b, s scalar.
     *
     * @param a     Vector a
     * @param b     Vector b
     * @param s     Scalar s
     * @param n     Length
     * @param aStep Step size for a
     * @param bStep Step size for b
     * @param aInd  Position index for a. 0 -> none
     * @param bInd  Position index for b. 0 -> none
     */
    static void addinv(double* a,const double* b,double s,int n,int aStep=1,
		       int bStep=1);
    static void addinv(double* a,const double* b,double s,int n,int aStep,
		       int bStep,const int* aInd,const int* bInd=0);

    /**
     * c = a + s./b, s a scalar.
     *
     * @param c     Target vector c
     * @param a     Vector a
     * @param b     Vector b
     * @param s     Scalar s
     * @param n     Length
     * @param cStep Step size for c
     * @param aStep Step size for a
     * @param bStep Step size for b
     * @param cInd  Position index for c. 0 -> none
     * @param aInd  Position index for a. 0 -> none
     * @param bInd  Position index for b. 0 -> none
     */
    static void addinv(double* c,const double* a,const double* b,double s,
		       int n,int cStep=1,int aStep=1,int bStep=1,
		       const int* cInd=0,const int* aInd=0,const int* bInd=0);

    /**
     * a = s*b, s scalar.
     * ATTENTION: 'smul' a *= s much faster in practice!
     *
     * @param a     Vector a
     * @param b     Vector b
     * @param s     Scalar s
     * @param n     Length
     * @param aStep Step size for a
     * @param bStep Step size for b
     * @param aInd  Position index for a. 0 -> none
     * @param bInd  Position index for b. 0 -> none
     */
    static void smul(double* a,const double* b,double s,int n,int aStep=1,
		     int bStep=1);
    static void smul(double* a,const double* b,double s,int n,int aStep,
		     int bStep,const int* aInd,const int* bInd=0);

    /**
     * a *= s, s scalar.
     *
     * @param a     Vector a
     * @param s     Scalar s
     * @param n     Length
     * @param aStep Step size for a
     * @param aInd  Position index for a. 0 -> none
     */
    static void smul(double* a,double s,int n,int aStep=1,const int* aInd=0);

    /**
     * a = b+s*1, s scalar.
     *
     * @param a     Vector a
     * @param b     Vector b
     * @param s     Scalar s
     * @param n     Length
     * @param aStep Step size for a
     * @param bStep Step size for b
     * @param aInd  Position index for a. 0 -> none
     * @param bInd  Position index for b. 0 -> none
     */
    static void addscal(double* a,const double* b,double s,int n,int aStep=1,
			int bStep=1);
    static void addscal(double* a,const double* b,double s,int n,int aStep,
			int bStep,const int* aInd,const int* bInd=0);

    /**
     * a += s*1, s scalar.
     *
     * @param a     Vector a
     * @param s     Scalar s
     * @param n     Length
     * @param aStep Step size for a
     * @param aInd  Position index for a. 0 -> none
     */
    static void addscal(double* a,double s,int n,int aStep=1,
			const int* aInd=0);

    /**
     * c = a .* b
     *
     * @param c     Target vector c
     * @param a     Vector a
     * @param b     Vector b
     * @param n     Length
     * @param cStep Step size for c
     * @param aStep Step size for a
     * @param bStep Step size for b
     * @param cInd  Position index for c. 0 -> none
     * @param aInd  Position index for a. 0 -> none
     * @param bInd  Position index for b. 0 -> none
     */
    static void prod(double* c,const double* a,const double* b,int n,
		     int cStep=1,int aStep=1,int bStep=1,const int* cInd=0,
		     const int* aInd=0,const int* bInd=0);

    /**
     * a .*= b, b vector
     *
     * @param a     Vector a
     * @param b     Vector b
     * @param n     Length
     * @param aStep Step size for a
     * @param bStep Step size for b
     * @param aInd  Position index for a. 0 -> none
     * @param bInd  Position index for b. 0 -> none
     */
    static void prod(double* a,const double* b,int n,int aStep=1,int bStep=1);
    static void prod(double* a,const double* b,int n,int aStep,int bStep,
		     const int* aInd,const int* bInd=0);

    /**
     * c = a ./ b
     *
     * @param c     Target vector c
     * @param a     Vector a
     * @param b     Vector b
     * @param n     Length
     * @param cStep Step size for c
     * @param aStep Step size for a
     * @param bStep Step size for b
     * @param cInd  Position index for c. 0 -> none
     * @param aInd  Position index for a. 0 -> none
     * @param bInd  Position index for b. 0 -> none
     */
    static void div(double* c,const double* a,const double* b,int n,
		    int cStep=1,int aStep=1,int bStep=1,const int* cInd=0,
		    const int* aInd=0,const int* bInd=0);

    /**
     * a ./= b, b vector
     *
     * @param a     Vector a
     * @param b     Vector b
     * @param n     Length
     * @param aStep Step size for a
     * @param bStep Step size for b
     * @param aInd  Position index for a. 0 -> none
     * @param bInd  Position index for b. 0 -> none
     */
    static void div(double* a,const double* b,int n,int aStep=1,int bStep=1);
    static void div(double* a,const double* b,int n,int aStep,int bStep,
		    const int* aInd,const int* bInd=0);

    /**
     * a += s*(b .* c), s scalar.
     * Cases s=+1/-1 treated sep.!
     *
     * @param a     Vector a
     * @param b     Vector b
     * @param c     Vector c
     * @param n     Length
     * @param s     Scalar s. Def.: 1
     * @param aStep Step size for a
     * @param bStep Step size for b
     * @param cStep Step size for c
     * @param aInd  Position index for a. 0 -> none
     * @param bInd  Position index for b. 0 -> none
     * @param cInd  Position index for c. 0 -> none
     */
    static void addprod(double* a,const double* b,const double* c,int n,
			double s=1.0,int aStep=1,int bStep=1,int cStep=1,
			const int* aInd=0,const int* bInd=0,const int* cInd=0);

    /**
     * a += s*(b .* b), s scalar.
     * Cases s=+1/-1 treated sep.!
     *
     * @param a     Vector a
     * @param b     Vector b
     * @param n     Length
     * @param s     Scalar s. Def.: 1
     * @param aStep Step size for a
     * @param bStep Step size for b
     * @param aInd  Position index for a. 0 -> none
     * @param bInd  Position index for b. 0 -> none
     */
    static void addprod(double* a,const double* b,int n,double s=1.0,
			int aStep=1,int bStep=1,const int* aInd=0,
			const int* bInd=0);

    /**
     * a += s*(b ./ c), s scalar.
     * Cases s=+1/-1 treated sep.!
     *
     * @param a     Vector a
     * @param b     Vector b
     * @param c     Vector c
     * @param n     Length
     * @param s     Scalar s. Def.: 1
     * @param aStep Step size for a
     * @param bStep Step size for b
     * @param cStep Step size for c
     * @param aInd  Position index for a. 0 -> none
     * @param bInd  Position index for b. 0 -> none
     * @param cInd  Position index for c. 0 -> none
     */
    static void adddiv(double* a,const double* b,const double* c,int n,
		       double s=1.0,int aStep=1,int bStep=1,int cStep=1,
		       const int* aInd=0,const int* bInd=0,const int* cInd=0);

    /**
     * a = (s*1) ./ b, b vector, s scalar
     *
     * @param a     Target vector a
     * @param b     Vector b
     * @param s     Scalar s
     * @param n     Length
     * @param aStep Step size for a
     * @param bStep Step size for b
     * @param aInd  Position index for a. 0 -> none
     * @param bInd  Position index for b. 0 -> none
     */
    static void inv(double* a,const double* b,double s,int n,int aStep=1,
		    int bStep=1);
    static void inv(double* a,const double* b,double s,int n,int aStep,
		    int bStep,const int* aInd,const int* bInd=0);

    /**
     * a = (s*1) ./ a, s scalar
     *
     * @param a     Target vector a
     * @param s     Scalar s
     * @param n     Length
     * @param aStep Step size for a
     * @param aInd  Position index for a. 0 -> none
     */
    static void inv(double* a,double s,int n,int aStep=1,const int* aInd=0);

    /**
     * Sum of elements: a^T * 1
     *
     * @param a     Vector a
     * @param n     Length
     * @param aStep Step size for a
     * @param aInd  Position index for a. 0 -> none
     */
    static double sum(const double* a,int n,int aStep=1,const int* aInd=0);

    /**
     * Stores cumulative sums of b in a, i.e.
     *   a_i = \sum_{j<=i} b_j
     *
     * @param a     Target vector a
     * @param b     Vector b
     * @param n     Length
     * @param aStep Step size for a
     * @param bStep Step size for b
     * @param aInd  Position index for a. 0 -> none
     * @param bInd  Position index for b. 0 -> none
     */
    static void cumulSum(double* a,const double* b,int n,int aStep=1,
			 int bStep=1,const int* aInd=0,const int* bInd=0);

    /*
     * Matrix-vector (BLAS II)
     * In general, we require flat matrices and vectors here ('BaseLinMat' and
     * 'BaseLinVec'). BLAS II used whenever applic.
     *  For symmetric or triangular matrices, the 'strpatt' member in
     * 'BaseLinMat' specifies where (within the buffer) the matrix is stored
     * (UPLO,DIAG in BLAS II). For the methods with triagular matrices,
     * 'strpatt' also spec. the matrix form (upper/lower triang.,
     * unit/free diagonal). Symmetric matrices can come as upper or lower
     * triangular, but not with unit diag. In spec. case, the diagonal can be
     * passed separately.
     */

    /**
     * a = alpha*op(M)*b + beta*a
     *
     * Form of M dep. on 'strpatt' in 'mmat': normal for 'strctNormal',
     * symmetric for 'strctUpper' or 'strctLower'. op(M)==M or M^T, dep. on
     * 'mtrans'. The unit diag. patterns are not allowed.
     * For triangular M, use 'trimatVec'.
     *
     * @param avec   Input/output vector
     * @param mmat   Input matrix
     * @param mtrans S.a.
     * @param bvec   Input vector
     * @param alpha  Scalar
     * @param beta   Scalar
     */
    static void matVec(BaseLinVec<double>& avec,const BaseLinMat<double>& mmat,
		       bool mtrans,const BaseLinVec<double>& bvec,double beta,
		       double alpha);

    /**
     * a = op(T)*a
     *
     * T is triangular, lower/upper, unit/free diagonal, dep. on 'strpatt' in
     * 'tmat'. op(T)==T or T^T, dep. on 'ttrans'.
     *
     * @param avec   Input/output vector
     * @param tmat   Input matrix
     * @param ttrans S.a.
     */
    static void trimatVec(BaseLinVec<double>& avec,
			  const BaseLinMat<double>& tmat,bool ttrans);

    /**
     * a = op(T^{-1}) a
     *
     * T is triangular, lower/upper, unit/free diagonal, dep. on 'strpatt' in
     * 'tmat'. op(T)==T or T^T, dep. on 'ttrans'.
     *
     * @param avec   Input/output vector
     * @param tmat   Input matrix
     * @param ttrans S.a.
     */
    static void backsubst(BaseLinVec<double>& avec,
			  const BaseLinMat<double>& tmat,bool ttrans);

    /**
     * M = M + alpha*b*c^T (rank-1 update)
     *
     * M must not have spec. structure (for symm. matrix, use 'symRankOne').
     *
     * @param mmat  Input/output matrix
     * @param bvec  Input vector
     * @param cvec  "
     * @param alpha Scalar
     */
    static void rankOne(BaseLinMat<double>& mmat,
			const BaseLinVec<double>& bvec,
			const BaseLinVec<double>& cvec,double alpha);

    /**
     * M = M + alpha*b*b^T (rank-1 update)
     *
     * M is unstructured (square) if 'strpatt'=='strctNormal', otherwise
     * symmetric (unit diag. structure not allowed).
     *
     * @param mmat  Input/output matrix
     * @param bvec  Input vector
     * @param alpha Scalar
     */
    static void symRankOne(BaseLinMat<double>& mmat,
			   const BaseLinVec<double>& bvec,double alpha);

    /**
     * M = M + alpha*(b*c^T + c*b^T) (rank-2 update)
     *
     * M must be symmetric ('strpatt' must be 'strctUpper' / 'strctLower').
     *
     * @param mmat  Input/output matrix
     * @param bvec  Input vector
     * @param cvec  "
     * @param alpha Scalar
     */
    static void symRankTwo(BaseLinMat<double>& mmat,
			   const BaseLinVec<double>& bvec,
			   const BaseLinVec<double>& cvec,double alpha);

    /**
     * M = diag(a)*B / M = B*diag(a) ('left'==true/false; 'inve'==false)
     * M = diag(a)^{-1}*B / M = B*diag(a)^{-1} ('inve'==true)
     *
     * B must not be 'strctUppNDg'/'strctLowNDg'. B,M can be the same matrix.
     * They must have same size, structure pattern. If 'incr'==true, M is
     * not overwritten by the res., but it is added to M.
     *
     * @param mmat Output matrix
     * @param avec Input vector
     * @param bmat Input matrix
     * @param inve S.a.
     * @param left S.a.
     * @param incr S.a. Def.: false
     */
    static void mulDiag(BaseLinMat<double>& mmat,
			const BaseLinVec<double>& avec,
			const BaseLinMat<double>& bmat,bool inve,bool left,
			bool incr=false);

    /* Matrix-matrix (BLAS III)
       In general, we require flat matrices (and vectors) here ('BaseLinMat'
       and 'BaseLinVec'). BLAS III used whenever applic.
       For symmetric or triangular matrices, the 'strpatt' member in
       'BaseLinMat' specifies where (within the buffer) the matrix is stored
       (UPLO,DIAG in BLAS III). For the methods with triagular matrices,
       'strpatt' also spec. the matrix form (upper/lower triang., unit/free
       diagonal). Symmetric matrices can come as upper or lower triangular,
       but not with unit diag. In spec. case, the diagonal can be passed
       separately. */

    /**
     * C = alpha*op(A)*op(B) + beta*C (matrix-matrix multiply)
     *
     * Here, op(A)==A or op(A)==A^T dep. on 'atrans'==false or true. All
     * A,B,C must have 'strctNormal' structure pattern.
     * For A or B triangular, use 'trimatMat'. For A or B symmetric, use
     * 'symmatMat'.
     * If C symmetric, use 'symRankK' / 'symRank2K'.
     *
     * @param cmat   Input/output matrix
     * @param amat   Input matrix
     * @param atrans S.a.
     * @param bmat   Input matrix
     * @param btrans S.a.
     * @param alpha  Scalar
     * @param beta   Scalar
     */
    static void matMat(BaseLinMat<double>& cmat,const BaseLinMat<double>& amat,
		       bool atrans,const BaseLinMat<double>& bmat,bool btrans,
		       double alpha,double beta);

    /**
     * C = alpha*A*B + beta*C (matrix-matrix multiply, one symmetric)
     *
     * Here, one of A,B is symmetric ('strctUpper' / 'strctLower'), but not
     * both. The other must have 'strctNormal', as must have C.
     *
     * @param cmat   Input/output matrix
     * @param amat   Input matrix
     * @param bmat   Input matrix
     * @param alpha  Scalar
     * @param beta   Scalar
     */
    static void symmatMat(BaseLinMat<double>& cmat,
			  const BaseLinMat<double>& amat,
			  const BaseLinMat<double>& bmat,double alpha,
			  double beta);

    /**
     * B = alpha*op(A)*B (matrix-matrix multiply, A triangular,
     *                    if 'trgLeft'==true)
     * B = alpha*B*op(A) (if 'trgLeft'==false)
     *
     * Here, A must be triangular, B normal ('strctNormal'). op(A)==A or A^T,
     * dep. on 'atrans'.
     *
     * @param bmat    Input/output matrix
     * @param amat    Input matrix (triagular)
     * @param atrans  S.a.
     * @param trgLeft S.a.
     * @param alpha   Scalar
     */
    static void trimatMat(BaseLinMat<double>& bmat,
			  const BaseLinMat<double>& amat,bool atrans,
			  bool trgLeft,double alpha);

    /**
     * C = alpha*A*A^T + beta*C (rank k update, 'atrans'==false)
     * C = alpha*A^T*A + beta*C (rank k update, 'atrans'==true)
     *
     * Here, C must be symmetric ('strctUpper' / 'strctLower'), A must be
     * general ('strctNormal').
     *
     * @param cmat   Input/output matrix (symmetric)
     * @param amat   Input matrix
     * @param atrans S.a.
     * @param alpha  Scalar
     * @param beta   Scalar
     */
    static void symRankK(BaseLinMat<double>& cmat,
			 const BaseLinMat<double>& amat,bool atrans,
			 double alpha,double beta);

    /**
     * C = alpha*(A*B^T + B*A^T) + beta*C (rank 2*k update,'atrans'==false)
     * C = alpha*(A^T*B + B^T*A) + beta*C (rank 2*k update,'atrans'==true)
     *
     * Here, C must be symmetric ('strctUpper' / 'strctLower'), A,B must be
     * general ('strctNormal').
     *
     * @param cmat   Input/output matrix (symmetric)
     * @param amat   Input matrix
     * @param atrans S.a.
     * @param bmat   Input matrix
     * @param alpha  Scalar
     * @param beta   Scalar
     */
    static void symRank2K(BaseLinMat<double>& cmat,
			  const BaseLinMat<double>& amat,bool atrans,
			  const BaseLinMat<double>& bmat,double alpha,
			  double beta);

    /**
     * B = alpha*op(A^{-1})*B (backsubstitution, A triangular,
     *                         if 'trgLeft'==true)
     * B = alpha*B*op(A^{-1}) (if 'trgLeft'==false)
     *
     * Here, A must be triangular, B normal ('strctNormal'). op(A)==A or A^T,
     * dep. on 'atrans'.
     *
     * @param bmat    Input/output matrix
     * @param amat    Input matrix (triagular)
     * @param atrans  S.a.
     * @param trgLeft S.a.
     * @param alpha   Scalar. Def.: 1
     */
    static void backsubstMat(BaseLinMat<double>& bmat,
			     const BaseLinMat<double>& amat,bool atrans,
			     bool trgLeft,double alpha=1.0);
  };

  // Definition of static inline methods

  inline double FastUtils::logsumexp(const double* a,int n,int step,
				     const int* ind)
  {
    int i;
    double mx,temp;

    if (n<=0) return 0.0;
    if (ind==0) {
      a=CONVBUFFPTR(a,n,step);
      const double* tempr=a;
      mx=*tempr; tempr+=step;
      for (i=1; i<n; i++,tempr+=step) {
	temp=*tempr;
	if (temp>mx) mx=temp;
      }
      for (i=0,temp=0.0; i<n; i++,a+=step)
	temp+=exp((*a)-mx);
    } else {
      const int* tempi=ind;
      mx=a[*(tempi++)];
      for (i=1; i<n; i++) {
	temp=a[*(tempi++)];
	if (temp>mx) mx=temp;
      }
      for (i=0,temp=0.0,tempi=ind; i<n; i++)
	temp+=exp((a[*(tempi++)])-mx);
    }

    return mx+log(temp);
  }

  inline double FastUtils::inner(const double* a,const double* b,int n,
				 int aStep,int bStep,const int* aInd,
				 const int* bInd)
  {
    if (n<=0) return 0.0;
    if (aInd==0) {
      if (bInd==0)
	return cblas_ddot(n,a,aStep,b,bStep);
      else {
	int i;
	double sum=0.0;
	a=CONVBUFFPTR(a,n,aStep);
	for (i=0; i<n; i++,a+=aStep)
	  sum+=((*a)*b[*(bInd++)]);
	return sum;
      }
    } else {
      int i;
      double sum=0.0;
      if (bInd==0) {
	b=CONVBUFFPTR(b,n,bStep);
	for (i=0; i<n; i++,b+=bStep)
	  sum+=((*b)*a[*(aInd++)]);
      } else {
	for (i=0; i<n; i++)
	  sum+=(a[*(aInd++)]*b[*(bInd++)]);
      }
      return sum;
    }
  }

  inline double FastUtils::inner(const double* a,int n,int aStep,
				 const int* aInd)
  {
    if (n<=0) return 0.0;
    if (aInd==0) {
      return cblas_ddot(n,a,aStep,a,aStep);
    } else {
      int i;
      double sum=0.0,temp;
      for (i=0; i<n; i++) {
	temp=a[*(aInd++)];
	sum+=(temp*temp);
      }
      return sum;
    }
  }

  inline double FastUtils::mahal(const double* a,const double* b,
				 const double* d,int n,int aStep,int bStep,
				 int dStep,const int* aInd,const int* bInd,
				 const int* dInd)
  {
    // Uahhh! Lots of cases!
    int i;
    double sum=0.0;
    if (n<=0) return 0.0;
    if (aInd==0) a=CONVBUFFPTR(a,n,aStep);
    if (bInd==0) b=CONVBUFFPTR(b,n,bStep);
    if (dInd==0) d=CONVBUFFPTR(d,n,dStep);
    if (dInd==0) {
      if (aInd==0) {
	if (bInd==0) {
	  for (i=0; i<n; i++,d+=dStep,a+=aStep,b+=bStep)
	    sum+=((*d)*(*a)*(*b));
	} else {
	  for (i=0; i<n; i++,d+=dStep,a+=aStep)
	    sum+=((*d)*(*a)*b[*(bInd++)]);
	}
      } else {
	if (bInd==0) {
	  for (i=0; i<n; i++,d+=dStep,b+=bStep)
	    sum+=((*d)*a[*(aInd++)]*(*(b++)));
	} else {
	  for (i=0; i<n; i++,d+=dStep)
	    sum+=((*d)*a[*(aInd++)]*b[*(bInd++)]);
	}
      }
    } else {
      if (aInd==0) {
	if (bInd==0) {
	  for (i=0; i<n; i++,a+=aStep,b+=bStep)
	    sum+=(d[*(dInd++)]*(*a)*(*b));
	} else {
	  for (i=0; i<n; i++,a+=aStep)
	    sum+=(d[*(dInd++)]*(*a)*b[*(bInd++)]);
	}
      } else {
	if (bInd==0) {
	  for (i=0; i<n; i++,b+=bStep)
	    sum+=(d[*(dInd++)]*a[*(aInd++)]*(*(b++)));
	} else {
	  for (i=0; i<n; i++)
	    sum+=(d[*(dInd++)]*a[*(aInd++)]*b[*(bInd++)]);
	}
      }
    }
    return sum;
  }

  inline double FastUtils::mahal(const double* a,const double* d,int n,
				 int aStep,int dStep)
  {
    if (n<=0) return 0.0;
    a=CONVBUFFPTR(a,n,aStep);
    d=CONVBUFFPTR(d,n,dStep);
    if (aStep==1 && dStep==1) {
      int i,limit=n/16;
      double sum=0.0,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,
	a14,a15;
      const double* aB=a,*dB=d;

      for (i=0; i<limit; i++) {
	a0 =*(aB++); a1 =*(aB++); a2 =*(aB++); a3 =*(aB++);
	a4 =*(aB++); a5 =*(aB++); a6 =*(aB++); a7 =*(aB++);
	a8 =*(aB++); a9 =*(aB++); a10=*(aB++); a11=*(aB++);
	a12=*(aB++); a13=*(aB++); a14=*(aB++); a15=*(aB++);
	sum+=(a0 *a0 *(*(dB++))); sum+=(a1 *a1 *(*(dB++)));
	sum+=(a2 *a2 *(*(dB++))); sum+=(a3 *a3 *(*(dB++)));
	sum+=(a4 *a4 *(*(dB++))); sum+=(a5 *a5 *(*(dB++)));
	sum+=(a6 *a6 *(*(dB++))); sum+=(a7 *a7 *(*(dB++)));
	sum+=(a8 *a8 *(*(dB++))); sum+=(a9 *a9 *(*(dB++)));
	sum+=(a10*a10*(*(dB++))); sum+=(a11*a11*(*(dB++)));
	sum+=(a12*a12*(*(dB++))); sum+=(a13*a13*(*(dB++)));
	sum+=(a14*a14*(*(dB++))); sum+=(a15*a15*(*(dB++)));
      }
      for (i=limit*16; i<n; i++) {
	a0=*(aB++);
	sum+=(a0*a0*(*(dB++)));
      }
      return sum;
    } else {
      int i=0;
      double sum=0.0,temp;
      const double* aB=a,*dB=d;

      for (; i<n; i++,aB+=aStep,dB+=dStep) {
	temp=*aB; sum+=(temp*temp*(*dB));
      }
      return sum;
    }
  }

  inline double FastUtils::mahal(const double* a,const double* d,int n,
				 int aStep,int dStep,const int* aInd,
				 const int* dInd)
  {
    if (n<=0) return 0.0;
    if (aInd==0) {
      if (dInd==0)
	return mahal(a,d,n,aStep,dStep);
      else {
	int i;
	double sum=0.0,temp;
	a=CONVBUFFPTR(a,n,aStep);
	for (i=0; i<n; i++,a+=aStep) {
	  temp=*a;
	  sum+=(temp*temp*d[*(dInd++)]);
	}
	return sum;
      }
    } else {
      int i;
      double sum=0.0,temp;
      if (dInd==0) {
	d=CONVBUFFPTR(d,n,dStep);
	for (i=0; i<n; i++,d+=dStep) {
	  temp=a[*(aInd++)];
	  sum+=(temp*temp*(*d));
	}
      } else {
	for (i=0; i<n; i++) {
	  temp=a[*(aInd++)];
	  sum+=(temp*temp*d[*(dInd++)]);
	}
      }
      return sum;
    }
  }

  inline void FastUtils::addsmul(double* a,const double* b,double s,int n,
				 int aStep,int bStep,const int* aInd,
				 const int* bInd)
  {
    if (aInd==0 && bInd==0) {
      cblas_daxpy(n,s,b,bStep,a,aStep);
    } else {
      int i;
      double scal=s;
      if (aInd==0) a=CONVBUFFPTR(a,n,aStep);
      if (bInd==0) b=CONVBUFFPTR(b,n,bStep);
      if (scal==1.0) {
	// Add b to a
	if (aInd==0) {
	  for (i=0; i<n; i++,a+=aStep)
	    (*a)+=b[*(bInd++)];
	} else if (bInd==0) {
	  for (i=0; i<n; i++,b+=bStep)
	    a[*(aInd++)]+=(*b);
	} else {
	  for (i=0; i<n; i++)
	    a[*(aInd++)]+=b[*(bInd++)];
	}
      } else if (scal==-1.0) {
	// Add b to a
	if (aInd==0) {
	  for (i=0; i<n; i++,a+=aStep)
	    (*a)-=b[*(bInd++)];
	} else if (bInd==0) {
	  for (i=0; i<n; i++,b+=bStep)
	    a[*(aInd++)]-=(*b);
	} else {
	  for (i=0; i<n; i++)
	    a[*(aInd++)]-=b[*(bInd++)];
	}
      } else {
	// General 's'
	if (aInd==0) {
	  for (i=0; i<n; i++,a+=aStep)
	    (*a)+=(scal*b[*(bInd++)]);
	} else if (bInd==0) {
	  for (i=0; i<n; i++,b+=bStep)
	    a[*(aInd++)]+=(scal*(*b));
	} else {
	  for (i=0; i<n; i++)
	    a[*(aInd++)]+=(scal*b[*(bInd++)]);
	}
      }
    }
  }

  inline void FastUtils::addsmul(double* c,const double* a,const double* b,
				 double s,int n,int cStep,int aStep,int bStep,
				 const int* cInd,const int* aInd,
				 const int* bInd)
  {
    // Uahhh! Lots of cases!
    int i;
    double scal=s;
    if (aInd==0) a=CONVBUFFPTR(a,n,aStep);
    if (bInd==0) b=CONVBUFFPTR(b,n,bStep);
    if (cInd==0) c=CONVBUFFPTR(c,n,cStep);
    if (scal==1.0) {
      if (cInd==0) {
	if (aInd==0) {
	  if (bInd==0) {
	    for (i=0; i<n; i++,c+=cStep,a+=aStep,b+=bStep)
	      (*c)=(*a)+(*b);
	  } else {
	    for (i=0; i<n; i++,c+=cStep,a+=aStep)
	      (*c)=(*a)+b[*(bInd++)];
	  }
	} else {
	  if (bInd==0) {
	    for (i=0; i<n; i++,c+=cStep,b+=bStep)
	      (*c)=a[*(aInd++)]+(*b);
	  } else {
	    for (i=0; i<n; i++,c+=cStep)
	      (*c)=a[*(aInd++)]+b[*(bInd++)];
	  }
	}
      } else {
	if (aInd==0) {
	  if (bInd==0) {
	    for (i=0; i<n; i++,a+=aStep,b+=bStep)
	      c[*(cInd++)]=(*a)+(*b);
	  } else {
	    for (i=0; i<n; i++,a+=aStep)
	      c[*(cInd++)]=(*a)+b[*(bInd++)];
	  }
	} else {
	  if (bInd==0) {
	    for (i=0; i<n; i++,b+=bStep)
	      c[*(cInd++)]=a[*(aInd++)]+(*b);
	  } else {
	    for (i=0; i<n; i++)
	      c[*(cInd++)]=a[*(aInd++)]+b[*(bInd++)];
	  }
	}
      }
    } else if (scal==-1.0) {
      if (cInd==0) {
	if (aInd==0) {
	  if (bInd==0) {
	    for (i=0; i<n; i++,c+=cStep,a+=aStep,b+=bStep)
	      (*c)=(*a)-(*b);
	  } else {
	    for (i=0; i<n; i++,c+=cStep,a+=aStep)
	      (*c)=(*a)-b[*(bInd++)];
	  }
	} else {
	  if (bInd==0) {
	    for (i=0; i<n; i++,c+=cStep,b+=bStep)
	      (*c)=a[*(aInd++)]-(*b);
	  } else {
	    for (i=0; i<n; i++,c+=cStep)
	      (*c)=a[*(aInd++)]-b[*(bInd++)];
	  }
	}
      } else {
	if (aInd==0) {
	  if (bInd==0) {
	    for (i=0; i<n; i++,a+=aStep,b+=bStep)
	      c[*(cInd++)]=(*a)-(*b);
	  } else {
	    for (i=0; i<n; i++,a+=aStep)
	      c[*(cInd++)]=(*a)-b[*(bInd++)];
	  }
	} else {
	  if (bInd==0) {
	    for (i=0; i<n; i++,b+=bStep)
	      c[*(cInd++)]=a[*(aInd++)]-(*b);
	  } else {
	    for (i=0; i<n; i++)
	      c[*(cInd++)]=a[*(aInd++)]-b[*(bInd++)];
	  }
	}
      }
    } else {
      if (cInd==0) {
	if (aInd==0) {
	  if (bInd==0) {
	    for (i=0; i<n; i++,c+=cStep,a+=aStep,b+=bStep)
	      (*c)=(*a)+scal*(*b);
	  } else {
	    for (i=0; i<n; i++,c+=cStep,a+=aStep)
	      (*c)=(*a)+scal*b[*(bInd++)];
	  }
	} else {
	  if (bInd==0) {
	    for (i=0; i<n; i++,c+=cStep,b+=bStep)
	      (*c)=a[*(aInd++)]+scal*(*b);
	  } else {
	    for (i=0; i<n; i++,c+=cStep)
	      (*c)=a[*(aInd++)]+scal*b[*(bInd++)];
	  }
	}
      } else {
	if (aInd==0) {
	  if (bInd==0) {
	    for (i=0; i<n; i++,a+=aStep,b+=bStep)
	      c[*(cInd++)]=(*a)+scal*(*b);
	  } else {
	    for (i=0; i<n; i++,a+=aStep)
	      c[*(cInd++)]=(*a)+scal*b[*(bInd++)];
	  }
	} else {
	  if (bInd==0) {
	    for (i=0; i<n; i++,b+=bStep)
	      c[*(cInd++)]=a[*(aInd++)]+scal*(*b);
	  } else {
	    for (i=0; i<n; i++)
	      c[*(cInd++)]=a[*(aInd++)]+scal*b[*(bInd++)];
	  }
	}
      }
    }
  }

  inline void FastUtils::addinv(double* a,const double* b,double s,int n,
				int aStep,int bStep)
  {
    a=CONVBUFFPTR(a,n,aStep);
    b=CONVBUFFPTR(b,n,bStep);
    if (aStep==1 && bStep==1) {
      int i,limit=n/16;
      double a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,
	scal=s;
      double* aB=a;
      const double* bB=b;

      for (i=0; i<limit; i++) {
	a0 =(*(bB++)); a1 =(*(bB++)); a2 =(*(bB++));
	a3 =(*(bB++)); a4 =(*(bB++)); a5 =(*(bB++));
	a6 =(*(bB++)); a7 =(*(bB++)); a8 =(*(bB++));
	a9 =(*(bB++)); a10=(*(bB++)); a11=(*(bB++));
	a12=(*(bB++)); a13=(*(bB++)); a14=(*(bB++));
	a15=(*(bB++));
	*(aB++)+=scal/a0;  *(aB++)+=scal/a1;  *(aB++)+=scal/a2;
	*(aB++)+=scal/a3;
	*(aB++)+=scal/a4;  *(aB++)+=scal/a5;  *(aB++)+=scal/a6;
	*(aB++)+=scal/a7;
	*(aB++)+=scal/a8;  *(aB++)+=scal/a9;  *(aB++)+=scal/a10;
	*(aB++)+=scal/a11;
	*(aB++)+=scal/a12; *(aB++)+=scal/a13; *(aB++)+=scal/a14;
	*(aB++)+=scal/a15;
      }
      for (i=limit*16; i<n; i++)
	*(aB++)+=scal/(*(bB++));
    } else {
      int i=0;
      double* aB=a;
      const double* bB=b;
      double scal=s;

      for (; i<n; i++,aB+=aStep,bB+=bStep) (*aB)+=scal/(*bB);
    }
  }

  inline void FastUtils::addinv(double* a,const double* b,double s,int n,
				int aStep,int bStep,const int* aInd,
				const int* bInd)
  {
    if (aInd==0 && bInd==0) {
      addinv(a,b,s,n,aStep,bStep);
    } else {
      int i;
      double scal=s;
      if (aInd==0) {
	a=CONVBUFFPTR(a,n,aStep);
	for (i=0; i<n; i++,a+=aStep)
	  (*a)+=(scal/b[*(bInd++)]);
      } else if (bInd==0) {
	b=CONVBUFFPTR(b,n,bStep);
	for (i=0; i<n; i++,b+=bStep)
	  a[*(aInd++)]+=(scal/(*b));
      } else {
	for (i=0; i<n; i++)
	  a[*(aInd++)]+=(scal/b[*(bInd++)]);
      }
    }
  }

  inline void FastUtils::addinv(double* c,const double* a,const double* b,
				double s,int n,int cStep,int aStep,int bStep,
				const int* cInd,const int* aInd,
				const int* bInd)
  {
    int i;
    double scal=s;
    if (cInd==0) c=CONVBUFFPTR(c,n,cStep);
    if (aInd==0) a=CONVBUFFPTR(a,n,aStep);
    if (bInd==0) b=CONVBUFFPTR(b,n,bStep);
    if (cInd==0) {
      if (aInd==0) {
	if (bInd==0) {
	  for (i=0; i<n; i++,c+=cStep,a+=aStep,b+=bStep)
	    (*c)=(*a)+scal/(*b);
	} else {
	  for (i=0; i<n; i++,c+=cStep,a+=aStep)
	    (*c)=(*a)+scal/b[*(bInd++)];
	}
      } else {
	if (bInd==0) {
	  for (i=0; i<n; i++,c+=cStep,b+=bStep)
	    (*c)=a[*(aInd++)]+scal/(*b);
	} else {
	  for (i=0; i<n; i++,c+=cStep)
	    (*c)=a[*(aInd++)]+scal/b[*(bInd++)];
	}
      }
    } else {
      if (aInd==0) {
	if (bInd==0) {
	  for (i=0; i<n; i++,a+=aStep,b+=bStep)
	    c[*(cInd++)]=(*a)+scal/(*b);
	} else {
	  for (i=0; i<n; i++,a+=aStep)
	    c[*(cInd++)]=(*a)+scal/b[*(bInd++)];
	}
      } else {
	if (bInd==0) {
	  for (i=0; i<n; i++,b+=bStep)
	    c[*(cInd++)]=a[*(aInd++)]+scal/(*b);
	} else {
	  for (i=0; i<n; i++)
	    c[*(cInd++)]=a[*(aInd++)]+scal/b[*(bInd++)];
	}
      }
    }
  }

  inline void FastUtils::smul(double* a,const double* b,double s,int n,
			      int aStep,int bStep)
  {
    a=CONVBUFFPTR(a,n,aStep);
    b=CONVBUFFPTR(b,n,bStep);
    if (aStep==1 && bStep==1) {
      int i,limit=n/16;
      double a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,
	scal=s;
      double* aB=a;
      const double* bB=b;

      for (i=0; i<limit; i++) {
	a0 =(*(bB++))*scal; a1 =(*(bB++))*scal; a2 =(*(bB++))*scal;
	a3 =(*(bB++))*scal; a4 =(*(bB++))*scal; a5 =(*(bB++))*scal;
	a6 =(*(bB++))*scal; a7 =(*(bB++))*scal; a8 =(*(bB++))*scal;
	a9 =(*(bB++))*scal; a10=(*(bB++))*scal; a11=(*(bB++))*scal;
	a12=(*(bB++))*scal; a13=(*(bB++))*scal; a14=(*(bB++))*scal;
	a15=(*(bB++))*scal;
	*(aB++)=a0;  *(aB++)=a1;  *(aB++)=a2;  *(aB++)=a3;
	*(aB++)=a4;  *(aB++)=a5;  *(aB++)=a6;  *(aB++)=a7;
	*(aB++)=a8;  *(aB++)=a9;  *(aB++)=a10; *(aB++)=a11;
	*(aB++)=a12; *(aB++)=a13; *(aB++)=a14; *(aB++)=a15;
      }
      for (i=limit*16; i<n; i++)
	*(aB++)=(scal*(*(bB++)));
    } else {
      int i=0;
      double scal=s;
      double* aB=a;
      const double* bB=b;

      for (; i<n; i++,aB+=aStep,bB+=bStep) *aB=scal*(*bB);
    }
  }

  inline void FastUtils::smul(double* a,const double* b,double s,int n,
			      int aStep,int bStep,const int* aInd,
			      const int* bInd)
  {
    if (aInd==0 && bInd==0) {
      smul(a,b,s,n,aStep,bStep);
    } else {
      int i;
      double scal=s;
      if (aInd==0) {
	a=CONVBUFFPTR(a,n,aStep);
	for (i=0; i<n; i++,a+=aStep)
	  (*a)=(scal*b[*(bInd++)]);
      } else if (bInd==0) {
	b=CONVBUFFPTR(b,n,bStep);
	for (i=0; i<n; i++,b+=bStep)
	  a[*(aInd++)]=(scal*(*b));
      } else {
	for (i=0; i<n; i++)
	  a[*(aInd++)]=(scal*b[*(bInd++)]);
      }
    }
  }

  inline void FastUtils::smul(double* a,double s,int n,int aStep,
			      const int* aInd)
  {
    int i;
    double scal=s;
    if (aInd==0) {
      cblas_dscal(n,s,a,aStep);
    } else {
      for (i=0; i<n; i++) a[*(aInd++)]*=scal;
    }
  }

  inline void FastUtils::addscal(double* a,const double* b,double s,int n,
				 int aStep,int bStep)
  {
    a=CONVBUFFPTR(a,n,aStep);
    b=CONVBUFFPTR(b,n,bStep);
    if (aStep==1 && bStep==1) {
      int i,limit=n/16;
      double a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,
	scal=s;
      double* aB=a;
      const double* bB=b;

      for (i=0; i<limit; i++) {
	a0 =(*(bB++))+scal; a1 =(*(bB++))+scal; a2 =(*(bB++))+scal;
	a3 =(*(bB++))+scal; a4 =(*(bB++))+scal; a5 =(*(bB++))+scal;
	a6 =(*(bB++))+scal; a7 =(*(bB++))+scal; a8 =(*(bB++))+scal;
	a9 =(*(bB++))+scal; a10=(*(bB++))+scal; a11=(*(bB++))+scal;
	a12=(*(bB++))+scal; a13=(*(bB++))+scal; a14=(*(bB++))+scal;
	a15=(*(bB++))+scal;
	*(aB++)=a0;  *(aB++)=a1;  *(aB++)=a2;  *(aB++)=a3;
	*(aB++)=a4;  *(aB++)=a5;  *(aB++)=a6;  *(aB++)=a7;
	*(aB++)=a8;  *(aB++)=a9;  *(aB++)=a10; *(aB++)=a11;
	*(aB++)=a12; *(aB++)=a13; *(aB++)=a14; *(aB++)=a15;
      }
      for (i=limit*16; i<n; i++)
	*(aB++)=(scal+(*(bB++)));
    } else {
      int i=0;
      double* aB=a;
      const double* bB=b;
      double scal=s;

      for (; i<n; i++,aB+=aStep,bB+=bStep)
	*aB=(scal+(*bB));
    }
  }

  inline void FastUtils::addscal(double* a,const double* b,double s,int n,
				 int aStep,int bStep,const int* aInd,
				 const int* bInd)
  {
    if (aInd==0 && bInd==0) {
      addscal(a,b,s,n,aStep,bStep);
    } else {
      int i;
      double scal=s;
      if (aInd==0) {
	a=CONVBUFFPTR(a,n,aStep);
	for (i=0; i<n; i++,a+=aStep)
	  *a=(scal+b[*(bInd++)]);
      } else if (bInd==0) {
	b=CONVBUFFPTR(b,n,bStep);
	for (i=0; i<n; i++,b+=bStep)
	  a[*(aInd++)]=(scal+(*b));
      } else {
	for (i=0; i<n; i++)
	  a[*(aInd++)]=(scal+b[*(bInd++)]);
      }
    }
  }

  inline void FastUtils::addscal(double* a,double s,int n,int aStep,
				 const int* aInd)
  {
    int i;
    double scal=s;
    if (aInd==0) {
      a=CONVBUFFPTR(a,n,aStep);
      for (i=0; i<n; i++,a+=aStep) (*a)+=scal;
    } else {
      for (i=0; i<n; i++) a[*(aInd++)]+=scal;
    }
  }

  inline void FastUtils::prod(double* c,const double* a,const double* b,int n,
			      int cStep,int aStep,int bStep,const int* cInd,
			      const int* aInd,const int* bInd)
  {
    // Uahhh! Lots of cases!
    int i;
    if (cInd==0) c=CONVBUFFPTR(c,n,cStep);
    if (aInd==0) a=CONVBUFFPTR(a,n,aStep);
    if (bInd==0) b=CONVBUFFPTR(b,n,bStep);
    if (cInd==0) {
      if (aInd==0) {
	if (bInd==0) {
	  for (i=0; i<n; i++,a+=aStep,b+=bStep,c+=cStep)
	    (*c)=(*a)*(*b);
	} else {
	  for (i=0; i<n; i++,a+=aStep,c+=cStep)
	    (*c)=(*a)*b[*(bInd++)];
	}
      } else {
	if (bInd==0) {
	  for (i=0; i<n; i++,b+=bStep,c+=cStep)
	    (*c)=a[*(aInd++)]*(*b);
	} else {
	  for (i=0; i<n; i++,c+=cStep)
	    (*c)=a[*(aInd++)]*b[*(bInd++)];
	}
      }
    } else {
      if (aInd==0) {
	if (bInd==0) {
	  for (i=0; i<n; i++,a+=aStep,b+=bStep)
	    c[*(cInd++)]=(*a)*(*b);
	} else {
	  for (i=0; i<n; i++,a+=aStep)
	    c[*(cInd++)]=(*a)*b[*(bInd++)];
	}
      } else {
	if (bInd==0) {
	  for (i=0; i<n; i++,b+=bStep)
	    c[*(cInd++)]=a[*(aInd++)]*(*b);
	} else {
	  for (i=0; i<n; i++)
	    c[*(cInd++)]=a[*(aInd++)]*b[*(bInd++)];
	}
      }
    }
  }

  inline void FastUtils::prod(double* a,const double* b,int n,int aStep,
			      int bStep)
  {
    a=CONVBUFFPTR(a,n,aStep);
    b=CONVBUFFPTR(b,n,bStep);
    if (aStep==1 && bStep==1) {
      int i,limit=n/16;
      double a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15;
      double* aB=a;
      const double* bB=b;

      for (i=0; i<limit; i++) {
	a0 =(*(bB++)); a1 =(*(bB++)); a2 =(*(bB++));
	a3 =(*(bB++)); a4 =(*(bB++)); a5 =(*(bB++));
	a6 =(*(bB++)); a7 =(*(bB++)); a8 =(*(bB++));
	a9 =(*(bB++)); a10=(*(bB++)); a11=(*(bB++));
	a12=(*(bB++)); a13=(*(bB++)); a14=(*(bB++));
	a15=(*(bB++));
	*(aB++)*=a0;  *(aB++)*=a1;  *(aB++)*=a2;  *(aB++)*=a3;
	*(aB++)*=a4;  *(aB++)*=a5;  *(aB++)*=a6;  *(aB++)*=a7;
	*(aB++)*=a8;  *(aB++)*=a9;  *(aB++)*=a10; *(aB++)*=a11;
	*(aB++)*=a12; *(aB++)*=a13; *(aB++)*=a14; *(aB++)*=a15;
      }
      for (i=limit*16; i<n; i++)
	*(aB++)*=(*(bB++));
    } else {
      int i=0;
      double* aB=a;
      const double* bB=b;

      for (; i<n; i++,aB+=aStep,bB+=bStep) (*aB)*=(*bB);
    }
  }

  inline void FastUtils::prod(double* a,const double* b,int n,int aStep,
			      int bStep,const int* aInd,const int* bInd)
  {
    if (aInd==0 && bInd==0) {
      prod(a,b,n,aStep,bStep);
    } else {
      int i;
      if (aInd==0) {
	a=CONVBUFFPTR(a,n,aStep);
	for (i=0; i<n; i++,a+=aStep)
	  (*a)*=b[*(bInd++)];
      } else if (bInd==0) {
	b=CONVBUFFPTR(b,n,bStep);
	for (i=0; i<n; i++,b+=bStep)
	  a[*(aInd++)]*=(*b);
      } else {
	for (i=0; i<n; i++)
	  a[*(aInd++)]*=b[*(bInd++)];
      }
    }
  }

  inline void FastUtils::div(double* c,const double* a,const double* b,int n,
			     int cStep,int aStep,int bStep,
			     const int* cInd,const int* aInd,const int* bInd)
  {
    // Uahhh! Lots of cases!
    int i;
    if (cInd==0) c=CONVBUFFPTR(c,n,cStep);
    if (aInd==0) a=CONVBUFFPTR(a,n,aStep);
    if (bInd==0) b=CONVBUFFPTR(b,n,bStep);
    if (cInd==0) {
      if (aInd==0) {
	if (bInd==0) {
	  for (i=0; i<n; i++,a+=aStep,b+=bStep,c+=cStep)
	    (*c)=(*a)/(*b);
	} else {
	  for (i=0; i<n; i++,a+=aStep,c+=cStep)
	    (*c)=(*a)/b[*(bInd++)];
	}
      } else {
	if (bInd==0) {
	  for (i=0; i<n; i++,b+=bStep,c+=cStep)
	    (*c)=a[*(aInd++)]/(*b);
	} else {
	  for (i=0; i<n; i++,c+=cStep)
	    (*c)=a[*(aInd++)]/b[*(bInd++)];
	}
      }
    } else {
      if (aInd==0) {
	if (bInd==0) {
	  for (i=0; i<n; i++,a+=aStep,b+=bStep)
	    c[*(cInd++)]=(*a)/(*b);
	} else {
	  for (i=0; i<n; i++,a+=aStep)
	    c[*(cInd++)]=(*a)/b[*(bInd++)];
	}
      } else {
	if (bInd==0) {
	  for (i=0; i<n; i++,b+=bStep)
	    c[*(cInd++)]=a[*(aInd++)]/(*b);
	} else {
	  for (i=0; i<n; i++)
	    c[*(cInd++)]=a[*(aInd++)]/b[*(bInd++)];
	}
      }
    }
  }

  inline void FastUtils::div(double* a,const double* b,int n,int aStep,
			     int bStep)
  {
    a=CONVBUFFPTR(a,n,aStep);
    b=CONVBUFFPTR(b,n,bStep);
    if (aStep==1 && bStep==1) {
      int i,limit=n/16;
      double a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15;
      double* aB=a;
      const double* bB=b;

      for (i=0; i<limit; i++) {
	a0 =(*(bB++)); a1 =(*(bB++)); a2 =(*(bB++));
	a3 =(*(bB++)); a4 =(*(bB++)); a5 =(*(bB++));
	a6 =(*(bB++)); a7 =(*(bB++)); a8 =(*(bB++));
	a9 =(*(bB++)); a10=(*(bB++)); a11=(*(bB++));
	a12=(*(bB++)); a13=(*(bB++)); a14=(*(bB++));
	a15=(*(bB++));
	*(aB++)/=a0;  *(aB++)/=a1;  *(aB++)/=a2;  *(aB++)/=a3;
	*(aB++)/=a4;  *(aB++)/=a5;  *(aB++)/=a6;  *(aB++)/=a7;
	*(aB++)/=a8;  *(aB++)/=a9;  *(aB++)/=a10; *(aB++)/=a11;
	*(aB++)/=a12; *(aB++)/=a13; *(aB++)/=a14; *(aB++)/=a15;
      }
      for (i=limit*16; i<n; i++)
	*(aB++)/=(*(bB++));
    } else {
      int i=0;
      double* aB=a;
      const double* bB=b;

      for (; i<n; i++,aB+=aStep,bB+=bStep) (*aB)/=(*bB);
    }
  }

  inline void FastUtils::div(double* a,const double* b,int n,int aStep,
			     int bStep,const int* aInd,const int* bInd)
  {
    if (aInd==0 && bInd==0) {
      div(a,b,n,aStep,bStep);
    } else {
      int i;
      if (aInd==0) {
	a=CONVBUFFPTR(a,n,aStep);
	for (i=0; i<n; i++,a+=aStep)
	  (*a)/=b[*(bInd++)];
      } else if (bInd==0) {
	b=CONVBUFFPTR(b,n,bStep);
	for (i=0; i<n; i++,b+=bStep)
	  a[*(aInd++)]/=(*b);
      } else {
	for (i=0; i<n; i++)
	  a[*(aInd++)]/=b[*(bInd++)];
      }
    }
  }

  inline void FastUtils::addprod(double* a,const double* b,const double* c,
				 int n,double s,int aStep,int bStep,int cStep,
				 const int* aInd,const int* bInd,
				 const int* cInd)
  {
    // Uahhh! Lots of cases!
    int i;
    double sScal=s;
    if (cInd==0) c=CONVBUFFPTR(c,n,cStep);
    if (aInd==0) a=CONVBUFFPTR(a,n,aStep);
    if (bInd==0) b=CONVBUFFPTR(b,n,bStep);
    if (s==1.0) {
      if (aInd==0) {
	if (bInd==0) {
	  if (cInd==0) {
	    for (i=0; i<n; i++,a+=aStep,b+=bStep,c+=cStep)
	      (*a)+=(*b)*(*c);
	  } else {
	    for (i=0; i<n; i++,a+=aStep,b+=bStep)
	      (*a)+=(*b)*c[*(cInd++)];
	  }
	} else {
	  if (cInd==0) {
	    for (i=0; i<n; i++,a+=aStep,c+=cStep)
	      (*a)+=b[*(bInd++)]*(*c);
	  } else {
	    for (i=0; i<n; i++,a+=aStep)
	      (*a)+=b[*(bInd++)]*c[*(cInd++)];
	  }
	}
      } else {
	if (bInd==0) {
	  if (cInd==0) {
	    for (i=0; i<n; i++,b+=bStep,c+=cStep)
	      a[*(aInd++)]+=(*b)*(*c);
	  } else {
	    for (i=0; i<n; i++,b+=bStep)
	      a[*(aInd++)]+=(*b)*c[*(cInd++)];
	  }
	} else {
	  if (cInd==0) {
	    for (i=0; i<n; i++,c+=cStep)
	      a[*(aInd++)]+=b[*(bInd++)]*(*c);
	  } else {
	    for (i=0; i<n; i++)
	      a[*(aInd++)]+=b[*(bInd++)]*c[*(cInd++)];
	  }
	}
      }
    } else if (s==-1.0) {
      if (aInd==0) {
	if (bInd==0) {
	  if (cInd==0) {
	    for (i=0; i<n; i++,a+=aStep,b+=bStep,c+=cStep)
	      (*a)-=(*b)*(*c);
	  } else {
	    for (i=0; i<n; i++,a+=aStep,b+=bStep)
	      (*a)-=(*b)*c[*(cInd++)];
	  }
	} else {
	  if (cInd==0) {
	    for (i=0; i<n; i++,a+=aStep,c+=cStep)
	      (*a)-=b[*(bInd++)]*(*c);
	  } else {
	    for (i=0; i<n; i++,a+=aStep)
	      (*a)-=b[*(bInd++)]*c[*(cInd++)];
	  }
	}
      } else {
	if (bInd==0) {
	  if (cInd==0) {
	    for (i=0; i<n; i++,b+=bStep,c+=cStep)
	      a[*(aInd++)]-=(*b)*(*c);
	  } else {
	    for (i=0; i<n; i++,b+=bStep)
	      a[*(aInd++)]-=(*b)*c[*(cInd++)];
	  }
	} else {
	  if (cInd==0) {
	    for (i=0; i<n; i++,c+=cStep)
	      a[*(aInd++)]-=b[*(bInd++)]*(*c);
	  } else {
	    for (i=0; i<n; i++)
	      a[*(aInd++)]-=b[*(bInd++)]*c[*(cInd++)];
	  }
	}
      }
    } else {
      if (aInd==0) {
	if (bInd==0) {
	  if (cInd==0) {
	    for (i=0; i<n; i++,a+=aStep,b+=bStep,c+=cStep)
	      (*a)+=sScal*(*b)*(*c);
	  } else {
	    for (i=0; i<n; i++,a+=aStep,b+=bStep)
	      (*a)+=sScal*(*b)*c[*(cInd++)];
	  }
	} else {
	  if (cInd==0) {
	    for (i=0; i<n; i++,a+=aStep,c+=cStep)
	      (*a)+=sScal*b[*(bInd++)]*(*c);
	  } else {
	    for (i=0; i<n; i++,a+=aStep)
	      (*a)+=sScal*b[*(bInd++)]*c[*(cInd++)];
	  }
	}
      } else {
	if (bInd==0) {
	  if (cInd==0) {
	    for (i=0; i<n; i++,b+=bStep,c+=cStep)
	      a[*(aInd++)]+=sScal*(*b)*(*c);
	  } else {
	    for (i=0; i<n; i++,b+=bStep)
	      a[*(aInd++)]+=sScal*(*b)*c[*(cInd++)];
	  }
	} else {
	  if (cInd==0) {
	    for (i=0; i<n; i++,c+=cStep)
	      a[*(aInd++)]+=sScal*b[*(bInd++)]*(*c);
	  } else {
	    for (i=0; i<n; i++)
	      a[*(aInd++)]+=sScal*b[*(bInd++)]*c[*(cInd++)];
	  }
	}
      }
    }
  }

  inline void FastUtils::addprod(double* a,const double* b,int n,double s,
				 int aStep,int bStep,const int* aInd,
				 const int* bInd)
  {
    int i;
    double sScal=s,temp;
    if (aInd==0) a=CONVBUFFPTR(a,n,aStep);
    if (bInd==0) b=CONVBUFFPTR(b,n,bStep);
    if (s==1.0) {
      if (aInd==0) {
	if (bInd==0) {
	  for (i=0; i<n; i++,a+=aStep,b+=bStep) {
	    temp=*b;
	    (*a)+=temp*temp;
	  }
	} else {
	  for (i=0; i<n; i++,a+=aStep) {
	    temp=b[*(bInd++)];
	    (*a)+=temp*temp;
	  }
	}
      } else {
	if (bInd==0) {
	  for (i=0; i<n; i++,b+=bStep) {
	    temp=*b;
	    a[*(aInd++)]+=temp*temp;
	  }
	} else {
	  for (i=0; i<n; i++) {
	    temp=b[*(bInd++)];
	    a[*(aInd++)]+=temp*temp;
	  }
	}	
      }
    } else if (s==-1.0) {
      if (aInd==0) {
	if (bInd==0) {
	  for (i=0; i<n; i++,a+=aStep,b+=bStep) {
	    temp=*b;
	    (*a)-=temp*temp;
	  }
	} else {
	  for (i=0; i<n; i++,a+=aStep) {
	    temp=b[*(bInd++)];
	    (*a)-=temp*temp;
	  }
	}
      } else {
	if (bInd==0) {
	  for (i=0; i<n; i++,b+=bStep) {
	    temp=*b;
	    a[*(aInd++)]-=temp*temp;
	  }
	} else {
	  for (i=0; i<n; i++) {
	    temp=b[*(bInd++)];
	    a[*(aInd++)]-=temp*temp;
	  }
	}	
      }
    } else {
      if (aInd==0) {
	if (bInd==0) {
	  for (i=0; i<n; i++,a+=aStep,b+=bStep) {
	    temp=*b;
	    (*a)+=sScal*temp*temp;
	  }
	} else {
	  for (i=0; i<n; i++,a+=aStep) {
	    temp=b[*(bInd++)];
	    (*a)+=sScal*temp*temp;
	  }
	}
      } else {
	if (bInd==0) {
	  for (i=0; i<n; i++,b+=bStep) {
	    temp=*b;
	    a[*(aInd++)]+=sScal*temp*temp;
	  }
	} else {
	  for (i=0; i<n; i++) {
	    temp=b[*(bInd++)];
	    a[*(aInd++)]+=sScal*temp*temp;
	  }
	}	
      }
    }
  }

  inline void FastUtils::adddiv(double* a,const double* b,const double* c,
				int n,double s,int aStep,int bStep,int cStep,
				const int* aInd,const int* bInd,
				const int* cInd)
  {
    // Uahhh! Lots of cases!
    int i;
    double sScal=s;
    if (cInd==0) c=CONVBUFFPTR(c,n,cStep);
    if (aInd==0) a=CONVBUFFPTR(a,n,aStep);
    if (bInd==0) b=CONVBUFFPTR(b,n,bStep);
    if (s==1.0) {
      if (aInd==0) {
	if (bInd==0) {
	  if (cInd==0) {
	    for (i=0; i<n; i++,a+=aStep,b+=bStep,c+=cStep)
	      (*a)+=(*b)/(*c);
	  } else {
	    for (i=0; i<n; i++,a+=aStep,b+=bStep)
	      (*a)+=(*b)/c[*(cInd++)];
	  }
	} else {
	  if (cInd==0) {
	    for (i=0; i<n; i++,a+=aStep,c+=cStep)
	      (*a)+=b[*(bInd++)]/(*c);
	  } else {
	    for (i=0; i<n; i++,a+=aStep)
	      (*a)+=b[*(bInd++)]/c[*(cInd++)];
	  }
	}
      } else {
	if (bInd==0) {
	  if (cInd==0) {
	    for (i=0; i<n; i++,b+=bStep,c+=cStep)
	      a[*(aInd++)]+=(*b)/(*c);
	  } else {
	    for (i=0; i<n; i++,b+=bStep)
	      a[*(aInd++)]+=(*b)/c[*(cInd++)];
	  }
	} else {
	  if (cInd==0) {
	    for (i=0; i<n; i++,c+=cStep)
	      a[*(aInd++)]+=b[*(bInd++)]/(*c);
	  } else {
	    for (i=0; i<n; i++)
	      a[*(aInd++)]+=b[*(bInd++)]/c[*(cInd++)];
	  }
	}
      }
    } else if (s==-1.0) {
      if (aInd==0) {
	if (bInd==0) {
	  if (cInd==0) {
	    for (i=0; i<n; i++,a+=aStep,b+=bStep,c+=cStep)
	      (*a)-=(*b)/(*c);
	  } else {
	    for (i=0; i<n; i++,a+=aStep,b+=bStep)
	      (*a)-=(*b)/c[*(cInd++)];
	  }
	} else {
	  if (cInd==0) {
	    for (i=0; i<n; i++,a+=aStep,c+=cStep)
	      (*a)-=b[*(bInd++)]/(*c);
	  } else {
	    for (i=0; i<n; i++,a+=aStep)
	      (*a)-=b[*(bInd++)]/c[*(cInd++)];
	  }
	}
      } else {
	if (bInd==0) {
	  if (cInd==0) {
	    for (i=0; i<n; i++,b+=bStep,c+=cStep)
	      a[*(aInd++)]-=(*b)/(*c);
	  } else {
	    for (i=0; i<n; i++,b+=bStep)
	      a[*(aInd++)]-=(*b)/c[*(cInd++)];
	  }
	} else {
	  if (cInd==0) {
	    for (i=0; i<n; i++,c+=cStep)
	      a[*(aInd++)]-=b[*(bInd++)]/(*c);
	  } else {
	    for (i=0; i<n; i++)
	      a[*(aInd++)]-=b[*(bInd++)]/c[*(cInd++)];
	  }
	}
      }
    } else {
      if (aInd==0) {
	if (bInd==0) {
	  if (cInd==0) {
	    for (i=0; i<n; i++,a+=aStep,b+=bStep,c+=cStep)
	      (*a)+=sScal*(*b)/(*c);
	  } else {
	    for (i=0; i<n; i++,a+=aStep,b+=bStep)
	      (*a)+=sScal*(*b)/c[*(cInd++)];
	  }
	} else {
	  if (cInd==0) {
	    for (i=0; i<n; i++,a+=aStep,c+=cStep)
	      (*a)+=sScal*b[*(bInd++)]/(*c);
	  } else {
	    for (i=0; i<n; i++,a+=aStep)
	      (*a)+=sScal*b[*(bInd++)]/c[*(cInd++)];
	  }
	}
      } else {
	if (bInd==0) {
	  if (cInd==0) {
	    for (i=0; i<n; i++,b+=bStep,c+=cStep)
	      a[*(aInd++)]+=sScal*(*b)/(*c);
	  } else {
	    for (i=0; i<n; i++,b+=bStep)
	      a[*(aInd++)]+=sScal*(*b)/c[*(cInd++)];
	  }
	} else {
	  if (cInd==0) {
	    for (i=0; i<n; i++,c+=cStep)
	      a[*(aInd++)]+=sScal*b[*(bInd++)]/(*c);
	  } else {
	    for (i=0; i<n; i++)
	      a[*(aInd++)]+=sScal*b[*(bInd++)]/c[*(cInd++)];
	  }
	}
      }
    }
  }

  inline void FastUtils::inv(double* a,const double* b,double s,int n,
			     int aStep,int bStep)
  {
    a=CONVBUFFPTR(a,n,aStep);
    b=CONVBUFFPTR(b,n,bStep);
    if (aStep==1 && bStep==1) {
      int i,limit=n/16;
      double a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,
	scal=s;
      double* aB=a;
      const double* bB=b;

      for (i=0; i<limit; i++) {
	a0 =(*(bB++)); a1 =(*(bB++)); a2 =(*(bB++));
	a3 =(*(bB++)); a4 =(*(bB++)); a5 =(*(bB++));
	a6 =(*(bB++)); a7 =(*(bB++)); a8 =(*(bB++));
	a9 =(*(bB++)); a10=(*(bB++)); a11=(*(bB++));
	a12=(*(bB++)); a13=(*(bB++)); a14=(*(bB++));
	a15=(*(bB++));
	*(aB++)=scal/a0;  *(aB++)=scal/a1;  *(aB++)=scal/a2;  *(aB++)=scal/a3;
	*(aB++)=scal/a4;  *(aB++)=scal/a5;  *(aB++)=scal/a6;  *(aB++)=scal/a7;
	*(aB++)=scal/a8;  *(aB++)=scal/a9;  *(aB++)=scal/a10; *(aB++)=scal/a11;
	*(aB++)=scal/a12; *(aB++)=scal/a13; *(aB++)=scal/a14; *(aB++)=scal/a15;
      }
      for (i=limit*16; i<n; i++)
	*(aB++)=scal/(*(bB++));
    } else {
      int i=0;
      double* aB=a;
      const double* bB=b;
      double scal=s;

      for (; i<n; i++,aB+=aStep,bB+=bStep) (*aB)=scal/(*bB);
    }
  }

  inline void FastUtils::inv(double* a,const double* b,double s,int n,
			     int aStep,int bStep,const int* aInd,
			     const int* bInd)
  {
    if (aInd==0 && bInd==0) {
      inv(a,b,s,n,aStep,bStep);
    } else {
      int i;
      if (aInd==0) {
	a=CONVBUFFPTR(a,n,aStep);
	for (i=0; i<n; i++,a+=aStep)
	  (*a)=s/b[*(bInd++)];
      } else if (bInd==0) {
	b=CONVBUFFPTR(b,n,bStep);
	for (i=0; i<n; i++,b+=bStep)
	  a[*(aInd++)]=s/(*b);
      } else {
	for (i=0; i<n; i++)
	  a[*(aInd++)]=s/b[*(bInd++)];
      }
    }
  }

  inline void FastUtils::inv(double* a,double s,int n,int aStep,
			     const int* aInd)
  {
    int i,pos;
    double temp;
    if (aInd==0) {
      a=CONVBUFFPTR(a,n,aStep);
      for (i=0; i<n; i++,a+=aStep) {
	temp=*a;
	(*a)=s/temp;
      }
    } else {
      for (i=0; i<n; i++) {
	pos=*(aInd++);
	a[pos]=s/a[pos];
      }
    }
  }

  inline double FastUtils::sum(const double* a,int n,int aStep,const int* aInd)
  {
    int i;
    double sum=0.0;
    if (aInd==0) {
      a=CONVBUFFPTR(a,n,aStep);
      for (i=0; i<n; i++,a+=aStep)
	sum+=(*a);
    } else {
      for (i=0; i<n; i++)
	sum+=a[*(aInd++)];
    }
    return sum;
  }

  inline void FastUtils::cumulSum(double* a,const double* b,int n,int aStep,
				  int bStep,const int* aInd,const int* bInd)
  {
    int i;
    double sum=0.0;
    if (aInd==0) {
      a=CONVBUFFPTR(a,n,aStep);
      if (bInd==0) {
	b=CONVBUFFPTR(b,n,bStep);
	for (i=0; i<n; i++,a+=aStep,b+=bStep)
	  *a=(sum+=(*b));
      } else {
	for (i=0; i<n; i++,a+=aStep)
	  *a=(sum+=b[*(bInd++)]);
      }
    } else {
      if (bInd==0) {
	b=CONVBUFFPTR(b,n,bStep);
	for (i=0; i<n; i++,b+=bStep)
	  a[*(aInd++)]=(sum+=(*b));
      } else {
	for (i=0; i<n; i++)
	  a[*(aInd++)]=(sum+=b[*(bInd++)]);
      }
    }
  }

  inline void FastUtils::matVec(BaseLinVec<double>& avec,
				const BaseLinMat<double>& mmat,bool mtrans,
				const BaseLinVec<double>& bvec,double alpha,
				double beta)
  {
    if (mmat.strpatt==WriteBackMat<double>::strctNormal) {
      // Normal matrix
      cblas_dgemv(CblasColMajor,SETTRS(mtrans),mmat.m,mmat.n,
		  alpha,mmat.buff,mmat.stride,bvec.buff,bvec.step,beta,
		  avec.buff,avec.step);
    } else {
      // Symmetric matrix ('strpatt' must be 'strctUpper' / 'strctLower')
      cblas_dsymv(CblasColMajor,SETUPLO(mmat.strpatt),mmat.n,alpha,mmat.buff,
		  mmat.stride,bvec.buff,bvec.step,beta,avec.buff,avec.step);
    }
  }

  inline void FastUtils::trimatVec(BaseLinVec<double>& avec,
				   const BaseLinMat<double>& tmat,bool ttrans)
  {
    cblas_dtrmv(CblasColMajor,SETUPLO(tmat.strpatt),SETTRS(ttrans),
		SETDIAG(tmat.strpatt),tmat.n,tmat.buff,tmat.stride,avec.buff,
		avec.step);
  }

  inline void FastUtils::backsubst(BaseLinVec<double>& avec,
				   const BaseLinMat<double>& tmat,bool ttrans)
  {
    cblas_dtrsv(CblasColMajor,SETUPLO(tmat.strpatt),SETTRS(ttrans),
		SETDIAG(tmat.strpatt),tmat.n,tmat.buff,tmat.stride,avec.buff,
		avec.step);
  }

  inline void FastUtils::rankOne(BaseLinMat<double>& mmat,
				 const BaseLinVec<double>& bvec,
				 const BaseLinVec<double>& cvec,double alpha)
  {
    cblas_dger(CblasColMajor,mmat.m,mmat.n,alpha,bvec.buff,bvec.step,
	       cvec.buff,cvec.step,mmat.buff,mmat.stride);
  }

  inline void FastUtils::symRankOne(BaseLinMat<double>& mmat,
				    const BaseLinVec<double>& bvec,
				    double alpha)
  {
    if (mmat.strpatt==WriteBackMat<double>::strctNormal) {
      // Must be square
      cblas_dger(CblasColMajor,mmat.m,mmat.n,alpha,bvec.buff,bvec.step,
		 bvec.buff,bvec.step,mmat.buff,mmat.stride);
    } else {
      // Symmetric
      cblas_dsyr(CblasColMajor,SETUPLO(mmat.strpatt),mmat.n,alpha,bvec.buff,
		 bvec.step,mmat.buff,mmat.stride);
    }
  }

  inline void FastUtils::symRankTwo(BaseLinMat<double>& mmat,
				    const BaseLinVec<double>& bvec,
				    const BaseLinVec<double>& cvec,
				    double alpha)
  {
    // M must be symmetric
    cblas_dsyr2(CblasColMajor,SETUPLO(mmat.strpatt),mmat.n,alpha,bvec.buff,
		bvec.step,cvec.buff,cvec.step,mmat.buff,mmat.stride);
  }

  inline void FastUtils::mulDiag(BaseLinMat<double>& mmat,
				 const BaseLinVec<double>& avec,
				 const BaseLinMat<double>& bmat,bool inve,
				 bool left,bool incr)
  {
    int i;
    const double* bP,*aP;
    double* mP;

    if (bmat.strpatt==WriteBackMat<double>::strctNormal) {
      // B,M normal
      if (left) {
	mP=mmat.buff; bP=bmat.buff;
	if (!inve) {
	  if (!incr) {
	    for (i=0; i<bmat.n; i++,mP+=mmat.stride,bP+=bmat.stride)
	      prod(mP,avec.buff,bP,avec.n,1,avec.step,1);
	  } else {
	    for (i=0; i<bmat.n; i++,mP+=mmat.stride,bP+=bmat.stride)
	      addprod(mP,avec.buff,bP,avec.n,1.0,1,avec.step,1);
	  }
	} else {
	  if (!incr) {
	    for (i=0; i<bmat.n; i++,mP+=mmat.stride,bP+=bmat.stride)
	      div(mP,bP,avec.buff,avec.n,1,1,avec.step);
	  } else {
	    for (i=0; i<bmat.n; i++,mP+=mmat.stride,bP+=bmat.stride)
	      adddiv(mP,bP,avec.buff,avec.n,1.0,1,1,avec.step);
	  }
	}
      } else {
	mP=mmat.buff; bP=bmat.buff;
	double temp;
	aP=CONVBUFFPTR(avec.buff,avec.n,avec.step);
	if (!incr) {
	  for (i=0; i<bmat.n; i++,mP+=mmat.stride,bP+=bmat.stride) {
	    temp=aP[i*avec.step];
	    smul(mP,bP,inve?(1/temp):temp,bmat.m);
	  }
	} else {
	  for (i=0; i<bmat.n; i++,mP+=mmat.stride,bP+=bmat.stride) {
	    temp=aP[i*avec.step];
	    addsmul(mP,bP,inve?(1/temp):temp,bmat.m);
	  }
	}
      }
    } else if (bmat.strpatt==WriteBackMat<double>::strctUpper) {
      // B,M upper triangular
      if (left) {
	mP=mmat.buff; bP=bmat.buff;
	if (!inve) {
	  if (avec.step>0) {
	    if (!incr) {
	      for (i=0; i<bmat.n; i++,mP+=mmat.stride,bP+=bmat.stride)
		prod(mP,avec.buff,bP,i+1,1,avec.step,1);
	    } else {
	      for (i=0; i<bmat.n; i++,mP+=mmat.stride,bP+=bmat.stride)
		addprod(mP,avec.buff,bP,i+1,1.0,1,avec.step,1);
	    }
	  } else {
	    aP=avec.buff+(avec.step*(1-bmat.n));
	    if (!incr) {
	      for (i=0; i<bmat.n; i++,mP+=mmat.stride,bP+=bmat.stride) {
		prod(mP,aP,bP,i+1,1,avec.step,1);
		aP+=avec.step;
	      }
	    } else {
	      for (i=0; i<bmat.n; i++,mP+=mmat.stride,bP+=bmat.stride) {
		addprod(mP,aP,bP,i+1,1.0,1,avec.step,1);
		aP+=avec.step;
	      }
	    }
	  }
	} else {
	  if (avec.step>0) {
	    if (!incr) {
	      for (i=0; i<bmat.n; i++,mP+=mmat.stride,bP+=bmat.stride)
		div(mP,bP,avec.buff,i+1,1,1,avec.step);
	    } else {
	      for (i=0; i<bmat.n; i++,mP+=mmat.stride,bP+=bmat.stride)
		adddiv(mP,bP,avec.buff,i+1,1.0,1,1,avec.step);
	    }
	  } else {
	    aP=avec.buff+(avec.step*(1-bmat.n));
	    if (!incr) {
	      for (i=0; i<bmat.n; i++,mP+=mmat.stride,bP+=bmat.stride) {
		div(mP,bP,aP,i+1,1,1,avec.step);
		aP+=avec.step;
	      }
	    } else {
	      for (i=0; i<bmat.n; i++,mP+=mmat.stride,bP+=bmat.stride) {
		adddiv(mP,bP,aP,i+1,1.0,1,1,avec.step);
		aP+=avec.step;
	      }
	    }
	  }
	}
      } else {
	mP=mmat.buff; bP=bmat.buff;
	double temp;
	aP=CONVBUFFPTR(avec.buff,avec.n,avec.step);
	if (!incr) {
	  for (i=0; i<bmat.n; i++,mP+=mmat.stride,bP+=bmat.stride) {
	    temp=aP[i*avec.step];
	    smul(mP,bP,inve?(1/temp):temp,i+1);
	  }
	} else {
	  for (i=0; i<bmat.n; i++,mP+=mmat.stride,bP+=bmat.stride) {
	    temp=aP[i*avec.step];
	    addsmul(mP,bP,inve?(1/temp):temp,i+1);
	  }
	}
      }
    } else {
      // B,M lower triangular
      if (left) {
	mP=mmat.buff; bP=bmat.buff;
	if (!inve) {
	  if (avec.step>0) {
	    if (!incr) {
	      for (i=0; i<bmat.n; i++,mP+=mmat.stride,bP+=bmat.stride)
		prod(mP+i,avec.buff+(i*avec.step),bP+i,bmat.n-i,1,
		     avec.step,1);
	    } else {
	      for (i=0; i<bmat.n; i++,mP+=mmat.stride,bP+=bmat.stride)
		addprod(mP+i,avec.buff+(i*avec.step),bP+i,bmat.n-i,1.0,1,
			avec.step,1);
	    }
	  } else {
	    if (!incr) {
	      for (i=0; i<bmat.n; i++,mP+=mmat.stride,bP+=bmat.stride)
		prod(mP+i,avec.buff,bP+i,bmat.n-i,1,avec.step,1);
	    } else {
	      for (i=0; i<bmat.n; i++,mP+=mmat.stride,bP+=bmat.stride)
		addprod(mP+i,avec.buff,bP+i,bmat.n-i,1.0,1,avec.step,1);
	    }
	  }
	} else {
	  if (avec.step>0) {
	    if (!incr) {
	      for (i=0; i<bmat.n; i++,mP+=mmat.stride,bP+=bmat.stride)
		div(mP+i,bP+i,avec.buff+(i*avec.step),bmat.n-i,1,1,avec.step);
	    } else {
	      for (i=0; i<bmat.n; i++,mP+=mmat.stride,bP+=bmat.stride)
		addprod(mP+i,bP+i,avec.buff+(i*avec.step),bmat.n-i,1.0,1,1,
			avec.step);
	    }
	  } else {
	    if (!incr) {
	      for (i=0; i<bmat.n; i++,mP+=mmat.stride,bP+=bmat.stride)
		div(mP+i,bP+i,avec.buff,bmat.m-i,1,1,avec.step);
	    } else {
	      for (i=0; i<bmat.n; i++,mP+=mmat.stride,bP+=bmat.stride)
		addprod(mP+i,bP+i,avec.buff,bmat.m-i,1.0,1,1,avec.step);
	    }
	  }
	}
      } else {
	mP=mmat.buff; bP=bmat.buff;
	double temp;
	aP=CONVBUFFPTR(avec.buff,avec.n,avec.step);
	if (!incr) {
	  for (i=0; i<bmat.n; i++,mP+=mmat.stride,bP+=bmat.stride) {
	    temp=aP[i*avec.step];
	    smul(mP+i,bP+i,inve?(1/temp):temp,bmat.m-i);
	  }
	} else {
	  for (i=0; i<bmat.n; i++,mP+=mmat.stride,bP+=bmat.stride) {
	    temp=aP[i*avec.step];
	    addsmul(mP+i,bP+i,inve?(1/temp):temp,bmat.m-i);
	  }
	}
      }
    }
  }

  inline void FastUtils::matMat(BaseLinMat<double>& cmat,
				const BaseLinMat<double>& amat,bool atrans,
				const BaseLinMat<double>& bmat,bool btrans,
				double alpha,double beta)
  {
    cblas_dgemm(CblasColMajor,SETTRS(atrans),SETTRS(btrans),cmat.m,cmat.n,
		atrans?amat.m:amat.n,alpha,amat.buff,amat.stride,bmat.buff,
		bmat.stride,beta,cmat.buff,cmat.stride);
  }

  inline void FastUtils::symmatMat(BaseLinMat<double>& cmat,
				   const BaseLinMat<double>& amat,
				   const BaseLinMat<double>& bmat,
				   double alpha,double beta)
  {
    if (amat.strpatt==WriteBackMat<double>::strctNormal) {
      // B symmetric
      cblas_dsymm(CblasColMajor,CblasRight,SETUPLO(bmat.strpatt),cmat.m,
		  cmat.n,alpha,amat.buff,amat.stride,bmat.buff,bmat.stride,
		  beta,cmat.buff,cmat.stride);
    } else {
      // A symmetric
      cblas_dsymm(CblasColMajor,CblasLeft,SETUPLO(amat.strpatt),cmat.m,cmat.n,
		  alpha,amat.buff,amat.stride,bmat.buff,bmat.stride,beta,
		  cmat.buff,cmat.stride);
    }
  }

  inline void FastUtils::trimatMat(BaseLinMat<double>& bmat,
				   const BaseLinMat<double>& amat,bool atrans,
				   bool trgLeft,double alpha)
  {
    cblas_dtrmm(CblasColMajor,trgLeft?CblasLeft:CblasRight,
		SETUPLO(amat.strpatt),SETTRS(atrans),SETDIAG(amat.strpatt),
		bmat.m,bmat.n,alpha,amat.buff,amat.stride,bmat.buff,
		bmat.stride);
  }

  inline void FastUtils::symRankK(BaseLinMat<double>& cmat,
				  const BaseLinMat<double>& amat,bool atrans,
				  double alpha,double beta)
  {
    cblas_dsyrk(CblasColMajor,SETUPLO(cmat.strpatt),SETTRS(atrans),cmat.n,
		atrans?amat.m:amat.n,alpha,amat.buff,amat.stride,beta,
		cmat.buff,cmat.stride);
  }

  inline void FastUtils::symRank2K(BaseLinMat<double>& cmat,
				  const BaseLinMat<double>& amat,bool atrans,
				  const BaseLinMat<double>& bmat,double alpha,
				  double beta)
  {
    cblas_dsyr2k(CblasColMajor,SETUPLO(cmat.strpatt),SETTRS(atrans),cmat.n,
		 atrans?amat.m:amat.n,alpha,amat.buff,amat.stride,bmat.buff,
		 bmat.stride,beta,cmat.buff,cmat.stride);
  }

  inline void FastUtils::backsubstMat(BaseLinMat<double>& bmat,
				      const BaseLinMat<double>& amat,
				      bool atrans,bool trgLeft,double alpha)
  {
    cblas_dtrsm(CblasColMajor,trgLeft?CblasLeft:CblasRight,
		SETUPLO(amat.strpatt),SETTRS(atrans),SETDIAG(amat.strpatt),
		bmat.m,bmat.n,alpha,amat.buff,amat.stride,bmat.buff,
		bmat.stride);
  }
#undef CONVBUFFPTR
#undef SETDIAG
#undef SETUPLO
#undef SETTRS

//ENDNS

#endif
