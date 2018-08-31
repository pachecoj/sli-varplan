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
 * Module: quad
 * Desc.:  Header class GaussProdQuadrature
 * ------------------------------------------------------------------- */

#ifndef GAUSSPRODQUADRATURE_H
#define GAUSSPRODQUADRATURE_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/quad/default.h"
#include "lhotse/quad/MultiQuadrature.h"
#include "lhotse/matrix/StVector.h"
#include "lhotse/matrix/StMatrix.h"
#include "lhotse/matrix/ArrayUtils.h"

//BEGINNS(quad)
  /**
   * Uses the product rule based on 3 point Gauss-Hermite quadrature
   * (in 1-d). All polynomials whose 1-d projections are all of degree
   * <= 5 are integrated exactly.
   * The 3-point Gauss-Hermite rule is:
   *   E[f(x)] \approx (2/3) f(0) + (1/6) (f(sqrt(3)) + f(-sqrt(3)))
   * <p>
   * NOTE: Requires 3^d function evaluations. Suitable only for rather
   * small d.
   * <p>
   * NOTE: Complexity of 'compMoments' / 'compLogMoments' is
   * 3^(d+1) - 4 d - 3 additions without the computation of the f vector.
   * If an evaluation of f(.) is O(d), then the dominating cost is in
   * computing the f vector which is O(3^d d).
   * <p>
   * NOTE: If this is time-critical, some further tricks should be done
   * to allow easy inlining and loop unrolling by the compiler:
   * - fix dimensionality 'd' to a constant (separate class for each req.
   *   'd'). Allows the compiler to "unroll" the recursive method calls
   * - if this does not lead to unrolling, recode recursive bit as unrolled
   *   loops
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class GaussProdQuadrature : public virtual MultiQuadrature
  {
  protected:
    // Constants

    static const double off,logw0,logw1; // parameters for Gauss-Hermite

    // Members

    int d;                        // dimension of u
    StVector fvec;                // transfer information
    int fvecFlag;                 /* 0: invalid; 1: 'compMoments';
				     2: compLogMoments */
    ArrayHandle<StVector> uvec;   // Needed by 'compLogPartHelper'
    StVector tempVec;             // "
    ArrayHandle<StVector> colMsk; // "
    ArrayHandle<double> logwval;  // "
    int depth;                    // "
    StMatrix cfactCopy;           // "
    int posFvec;                  // "
    StVector expFvec;             // Needed by 'compMoments'
    StVector expTilf;             // Needed by 'compMoments'
    StMatrix logvecs;             // Needed by 'compLogPartMultiHelper'
    StVector valvec;              // "
    StVector rowMsk;              // "

  public:
    // Public methods

    /**
     * Constructor.
     *
     * @param dim Dimension of u
     */
    GaussProdQuadrature(int dim) : d(dim),fvecFlag(0) {
      int i,j;

      if (dim<1) throw InvalidParameterException("'dim' must be positive");
      uvec.changeRep(dim); logwval.changeRep(dim);
      colMsk.changeRep(dim);
      for (i=0,j=1; i<dim; i++,j*=3) uvec[i].zeros(dim);
      tempVec.zeros(dim);
      fvec.zeros(j);
      expFvec.zeros(j); expTilf.zeros(j/3);
      cfactCopy.zeros(dim);
      char msg[80];
      sprintf(msg,"GaussProdQuadrature: Rule needs %d evals.",j);
      printMsgStdout(msg);
    }

    int getDim() const {
      return d;
    }

    bool hasNonnegWeights() const {
      return true;
    }

    /**
     * We store the vector (f(u) + log(w_i) - log Z) in 'fvec'. See code for
     * the order of u. w_i the prod. of corr. weights.
     */
    double compLogPart(const StVector& mu,const StMatrix& cfact);

    /**
     * ATTENTION: Replace by more efficient implementation for special
     * log likelihood!
     */
    double compLogPartDiag(const StVector& mu,const StVector& cfact) {
      StMatrix cfmat; cfmat.zeros(getDim());
      ((StVector&) cfmat.diag())=cfact;
      return compLogPart(mu,cfmat);
    }

    /**
     * We need the vector 'fvec' computed by 'compLogPart' (for the same
     * parameters).
     */
    void compMoments(const StVector& mu,const StMatrix& cfact,
		     StVector& tilmu,StMatrix& tilcov,bool raw=false);

    void compMomentsDiag(const StVector& mu,const StVector& cfact,
			 StVector& tilmu,StVector& tildiag,bool raw=false);

    void compLogPartMulti(const StVector& mu,const StMatrix& cfact,
			  StVector& logz);

    /**
     * ATTENTION: Replace by more efficient implementation for special
     * log likelihood!
     */
    void compLogPartMultiDiag(const StVector& mu,const StVector& cfact,
			      StVector& logz) {
      StMatrix cfmat; cfmat.zeros(getDim());
      ((StVector&) cfmat.diag())=cfact;
      compLogPartMulti(mu,cfmat,logz);
    }

    /**
     * ATTENTION: Not efficient, just for prototyping! Replace by more
     * efficient implementation for special log likelihood!
     */
    void compMultiSpecDiag(const StVector& mu,const StVector& cfact,
			   StVector& logz,double logpart);

    /**
     * The vector (w_i f(u)) is stored in 'fvec', required by
     * 'compLogMoments'. Same order of u as for 'compLogPart', same def.
     * of w_i.
     */
    double compExpectLog(const StVector& mu,const StMatrix& cfact);

    /**
     * ATTENTION: Replace by more efficient implementation for special
     * log likelihood!
     */
    double compExpectLogDiag(const StVector& mu,const StVector& cfact) {
      StMatrix cfmat; cfmat.zeros(getDim());
      ((StVector&) cfmat.diag())=cfact;
      return compExpectLog(mu,cfmat);
    }

    /**
     * Requires 'fvec' computed by 'compExpectLog' for the same pars.
     * NOTE: Given 'fvec', we do not access 'mu', 'cfact' here.
     */
    void compLogMoments(const StVector& mu,const StMatrix& cfact,
			StVector& cvec,StMatrix& gmat);

    /**
     * ATTENTION: Only log moments procedure 0 is implemented!
     */
    void compLogMomentsDiag(const StVector& mu,const StVector& cfact,
			    StVector& cvec,StVector& gvec,double* elog=0);

  protected:
    // Internal methods

    void resetTransfer() {
      fvecFlag=0;
    }

    /**
     * Helper for 'compLogPart'. Called in a recursive way. See code.
     * Uses members 'uvec','logwval','depth'. uvec[0] has to be init.
     * to the mean, logwval[0] to 0 (depth -1 variables). 'colMsk[i]'
     * has to mask the i-th CF column.
     */
    void compLogPartHelper();

    /**
     * Much the same as 'compLogPartHelper', but a matrix is accum.
     * in 'logvecs'.
     */
    void compLogPartMultiHelper();

    void compMultiSpecDiagHelper();

    /**
     * Much the same as 'compLogPartHelper', but (w_i f(u)) is built.
     * Uses same members as 'compLogPartHelper'.
     */
    void compExpectLogHelper();

    /**
     * General helper. The size of 'vec' is 3^k. We sum positions
     * j, 3^(k-1)+j, 2*3^(k-1)+j and write this into pos. j,
     * j=0,...,3^(k-1)-1.
     *
     * @param vec Input/output vector
     * @param sz  Size of 'vec' / 3, must be 3^(k-1), k>=1. Returned
     *            'vec' has size 'sz'
     */
    static void marginalize(double* vec,int sz);
  };

  // Inline methods

  /*
   * The 1-d Gauss-Hermite rule is
   *   w_0 g(0) + w_1 (g(off) + g(-off)),
   * where g is the integrand relative to N(0,I). The product rule is just
   * the product of that and requires 3^d evaluations. We traverse in the
   * order
   *   (0,0,0,0,0), (0,0,0,0,1), ..., (0,0,0,1,0), ...
   */

  /*
   * For the log partition function,
   *   g(x) = exp( f(L*x + mu) ).
   * We traverse in the natural order and compute the components of 'fvec'
   * as f(L*x + mu) + log(w_i), w_i the appropr. product of weights.
   * We then use stable summation ('stableLogSum') to compute log Z,
   * subtract that off 'fvec'.
   * All the work is done by the recursive method 'compLogPartHelper'
   */
  inline double GaussProdQuadrature::compLogPart(const StVector& mu,
						 const StMatrix& cfact)
  {
    if (mu.size()!=d || cfact.rows()!=d || cfact.cols()!=d)
      throw WrongDimensionException("mu, cfact have wrong size");

    posFvec=0; // pos. into 'fvec'
    uvec[0]=mu; // depth -1
    logwval[0]=0.0; // depth -1
    cfactCopy=cfact; // 'cfactCopy' is not a mask, 'cfact' may be one!
    for (int i=0; i<d; i++)
      colMsk[i].reassign((StVector&) cfactCopy(RangeFull::get(),i));
    depth=0;
    compLogPartHelper(); // Does all the work (compute 'fvec')
    if (posFvec!=fvec.size())
      throw InternalException(EXCEPT_MSG("")); // Sanity check
    double logz=fvec.stableLogSum();
    fvec.addscal(-logz); // complete 'fvec'
    fvecFlag=1;

    return logz;
  }

  /*
   * Recursive helper method for 'compLogPart'. Information is transferred
   * via 'uvec','logwval'. 'depth' is the current depth of computation
   * (0,...,d-1).
   * For 'depth'<d-1, we call ourselves 3 times, once for the different
   * signs 0,+1,-1. Before each call, we compute the current 'uvec[depth+1]'
   * based on 'uvec[depth]', the sign and col. 'depth' of 'cfactCopy'.
   * ATTENTION: Both 'uvec' and 'logwval' are indexed differently,
   * 'uvec[depth+1]' for depth 'depth' (like base 1 arrays).
   *
   * The real work is done for 'depth'==d-1 (no more rec. calls), by doing
   * function eval for the 'uvec[depth+1]'.
   * NOTE: No parameters to speed up recursive calls and allow compiler
   * inlining.
   */
  inline void GaussProdQuadrature::compLogPartHelper()
  {
    double temp;
    if (depth<d-1) {
      // Recursive calls
      depth++; // next level
      // Sign 0
      uvec[depth]=uvec[depth-1];
      temp=logwval[depth-1];
      logwval[depth]=temp+logw0;
      compLogPartHelper();
      // Sign +1
      uvec[depth].addprod(uvec[depth-1],off,colMsk[depth-1]);
      logwval[depth]=temp+logw1;
      compLogPartHelper();
      // Sign -1
      uvec[depth].addprod(uvec[depth-1],-off,colMsk[depth-1]);
      compLogPartHelper();
      depth--; // back to this level
    } else {
      // The real work: evaluations
      // Sign 0
      temp=logwval[depth];
      fvec[posFvec++]=compF(uvec[depth])+logw0+temp;
      // Sign +1
      temp+=logw1;
      tempVec.addprod(uvec[depth],off,colMsk[depth]);
      fvec[posFvec++]=compF(tempVec)+temp;
      // Sign -1
      tempVec.addprod(uvec[depth],-off,colMsk[depth]);
      fvec[posFvec++]=compF(tempVec)+temp;
    }
  }

  /*
   * For the mean,
   *   g(x) = exp( f(u) - log Z) * u
   * for the covariance matrix
   *   g(x) = exp( f(u) - log Z) * u*u^T,
   * where u = L*x + mu.
   * If
   *   alpha(x) = exp( f(u) - log Z),
   * we accumulate a matrix alpha(x)*x*x^T and a vector alpha(x)*x (the
   * "raw" statistics). The x have simple structure which is used by our
   * method. 'fvec' contains values of f(u) - log Z + log(w_i) computed by
   * 'compLogPart'.
   *
   * If t_i = 0,1,2, s(t_i) = 0,+1,-1, t = (t_i), then:
   * \sigma_i     = \sum_t g(t) s(t_i)
   * \sigma_{i,i} = \sum_t g(t) s(t_i)^2
   * \sigma_{j,i} = \sum_t g(t) s(t_j) s(t_i)
   * Here, \sigma_i comp. of 'tilmu', \sigma_{j,i} comp. of 'tilcov'.
   * To compute them we need "marginals" of the g = exp(f) vector,
   * summed out over all but one and all but two comp. and these comp.
   * taking signs +1 / -1.
   */
  inline void GaussProdQuadrature::compMoments(const StVector& mu,
					       const StMatrix& cfact,
					       StVector& tilmu,
					       StMatrix& tilcov,bool raw)
  {
    register int l;
    register double sum;
    register const double* fP;
    int i,j,k1,k2;
    double temp;

    if (mu.size()!=d || cfact.rows()!=d || cfact.cols()!=d)
      throw WrongDimensionException(EXCEPT_MSG("mu,cfact"));
    if (fvecFlag!=1)
      throw WrongStatusException(EXCEPT_MSG("'compLogPart' has to be called"));
    StVector colMsk;
    tilmu.zeros(d); tilcov.zeros(d); // accus
    tilcov.setStrctPatt(WriteBackMat<double>::strctNormal);
    ArrayHandle<int> pow3(d);
    for (i=d-1,j=1; i>=0; i--,j*=3) pow3[i]=j; // 3^(d-i-1)
    expFvec.apply1(fvec,ptr_fun(exp)); // exp(f)
    ArrayHandle<double> expFA=expFvec.getFlatBuff();
    expTilf.zeros(fvec.size()/3);
    ArrayHandle<double> expTilfA=expTilf.getFlatBuff();

    for (i=0; i<d; i++) {
      /* Here, 'expFA' has size 3^(d-i), obtained from orig. one by
	 summing over dim. < i */
      colMsk.reassign((StVector&) tilcov(RangeFull::get(),i));
      // Sign i: +1
      k1=pow3[i]; // 3^(d-i-1)
      // Copy +1 part to 'expTilfA'
      ArrayUtils<double>::copy(expTilfA,expFA+k1,k1);
      for (j=i+1; j<d; j++) {
	/* Here, 'expTilfA' has size 3^(d-j), obtained from orig. one by
	   summing over dim. < j, except for i: +1 */
	// Sign j: +1
	k2=pow3[j]; // 3^(d-j-1)
	for (l=0,fP=expTilfA.p()+k2,sum=0.0; l<k2; l++) sum+=*(fP++);
	colMsk[j]+=sum; // +,+
	// Sign j: -1
	for (l=0,sum=0.0; l<k2; l++) sum+=*(fP++);
	colMsk[j]-=sum; // +,-
	// Marginalize 'expTilfA'
	marginalize(expTilfA,k2);
      }
      temp=*(expTilfA.p());
      tilmu[i]+=temp; colMsk[i]+=temp; // +
      // Sign i: -1
      // Copy -1 part to 'expTilfA'
      ArrayUtils<double>::copy(expTilfA,expFA+2*k1,k1);
      for (j=i+1; j<d; j++) {
	/* Here, 'expTilfA' has size 3^(d-j), obtained from orig. one by
	   summing over dim. < j, except for i: -1 */
	// Sign j: +1
	k2=pow3[j]; // 3^(d-j-1)
	for (l=0,fP=expTilfA.p()+k2,sum=0.0; l<k2; l++) sum+=*(fP++);
	colMsk[j]-=sum; // -,+
	// Sign j: -1
	for (l=0,sum=0.0; l<k2; l++) sum+=*(fP++);
	colMsk[j]+=sum; // -,-
	// Marginalize 'expTilfA'
	marginalize(expTilfA,k2);
      }
      temp=*(expTilfA.p());
      tilmu[i]-=temp; colMsk[i]+=temp; // -
      // Marginalize 'expFA'
      marginalize(expFA,k1);
    }
    // As a byproduct, 'expFA[0]' contains sum over orig. 'expFvec'
    double scalacc=*(expFA.p());

    tilmu.prod(off);
    tilcov.smul(off*off);
    if (!raw) {
      cfact.setStrctPatt(WriteBackMat<double>::strctLower);
      cfact.triMulVec(tilmu);
      tilcov.makeSymm(true);
      cfact.triMul(tilcov,false,true);
      cfact.triMul(tilcov,true,false);
      tilcov.setStrctPatt(WriteBackMat<double>::strctLower); // is symm.!
      tilcov.symRankTwo(tilmu,mu);
      tilcov.rankOne(mu,scalacc);
      tilmu.addprod(scalacc,mu); // mean complete
      tilcov.rankOne(tilmu,-1.0);
    }
    tilcov.makeSymm(true); // cov. complete
    tilcov.setStrctPatt(WriteBackMat<double>::strctNormal);
  }

  /*
   * As comp. to 'compMoments', do not need inner loop over j or 'expTilf'
   * copy, instead 'expF' is marginalized directly to obtain "one-site"
   * marginals.
   */
  inline void GaussProdQuadrature::compMomentsDiag(const StVector& mu,
						   const StVector& cfact,
						   StVector& tilmu,
						   StVector& tildiag,
						   bool raw)
  {
    register int l;
    register double sum;
    register const double* fP;
    int i,j,k1;

    if (mu.size()!=d || cfact.size()!=d)
      throw WrongDimensionException(EXCEPT_MSG("mu,cfact"));
    if (fvecFlag!=1)
      throw WrongStatusException(EXCEPT_MSG("'compLogPartDiag' has to be called"));
    tilmu.zeros(d); tildiag.zeros(d); // accus
    ArrayHandle<int> pow3(d);
    for (i=d-1,j=1; i>=0; i--,j*=3) pow3[i]=j; // 3^(d-i-1)
    expFvec.apply1(fvec,ptr_fun(exp)); // exp(f)
    ArrayHandle<double> expFA=expFvec.getFlatBuff();

    for (i=0; i<d; i++) {
      /* Here, 'expFA' has size 3^(d-i), obtained from orig. one by
	 summing over dim. < i */
      // Sign i: +1
      k1=pow3[i]; // 3^(d-i-1)
      for (l=0,fP=expFA.p()+k1,sum=0.0; l<k1; l++) sum+=*(fP++);
      tilmu[i]+=sum; tildiag[i]+=sum; // +
      // Sign i: -1
      for (l=0,sum=0.0; l<k1; l++) sum+=*(fP++);
      tilmu[i]-=sum; tildiag[i]+=sum; // -
      // Marginalize 'expFA'
      marginalize(expFA,k1);
    }
    // As a byproduct, 'expFA[0]' contains sum over orig. 'expFvec'
    double scalacc=*(expFA.p());

    tilmu.prod(off);
    tildiag.prod(off*off);
    if (!raw) {
      tilmu.prod(cfact);
      tildiag.prod(cfact); tildiag.prod(cfact);
      tildiag.addsprod(tildiag,tilmu,mu,2.0);
      tildiag.addsprod(tildiag,mu,mu,scalacc);
      tilmu.addprod(scalacc,mu); // mean complete
      tildiag.addsprod(tildiag,tilmu,tilmu,-1.0); // cov. diag. complete
    }
  }

  inline void GaussProdQuadrature::marginalize(double *vec,int sz)
  {
    register double* trgP=vec;
    register const double* src1P=vec+sz,*src2P=vec+2*sz;
    register int i;

    for (i=0; i<sz; i++)
      *(trgP++)+=(*(src1P++)+*(src2P++));
  }

  /*
   * Vectorized version of 'compLogPart' above.
   */
  inline void GaussProdQuadrature::compLogPartMulti(const StVector& mu,
						    const StMatrix& cfact,
						    StVector& logz)
  {
    int i,numcomp=getFNum();

    if (mu.size()!=d || cfact.rows()!=d || cfact.cols()!=d)
      throw WrongDimensionException(EXCEPT_MSG(""));
    logvecs.zeros(fvec.size(),numcomp);
    posFvec=0; // row pos. into 'logvecs'
    uvec[0]=mu; // depth -1
    logwval[0]=0.0; // depth -1
    cfactCopy=cfact; // 'cfactCopy' is not a mask, 'cfact' may be one!
    for (i=0; i<d; i++)
      colMsk[i].reassign((StVector&) cfactCopy(RangeFull::get(),i));
    depth=0;
    compLogPartMultiHelper(); // Does all the work (compute 'logvecs')
    if (posFvec!=logvecs.rows())
      throw InternalException(EXCEPT_MSG("")); // Sanity check

    logz.zeros(numcomp);
    for (i=0; i<numcomp; i++)
      logz[i]=logvecs(RangeFull::get(),i)->stableLogSum();
  }

  /*
   * Same as 'compLogPartHelper', but a matrix is computed here in
   * 'logvecs'.
   */
  inline void GaussProdQuadrature::compLogPartMultiHelper()
  {
    double temp;
    if (depth<d-1) {
      // Recursive calls
      depth++;
      // Sign 0
      uvec[depth]=uvec[depth-1];
      temp=logwval[depth-1];
      logwval[depth]=temp+logw0;
      compLogPartMultiHelper();
      // Sign +1
      uvec[depth].addprod(uvec[depth-1],off,colMsk[depth-1]);
      logwval[depth]=temp+logw1;
      compLogPartMultiHelper();
      // Sign -1
      uvec[depth].addprod(uvec[depth-1],-off,colMsk[depth-1]);
      compLogPartMultiHelper();
      depth--;
    } else {
      // The real work: evaluations
      // Sign 0
      temp=logwval[depth];
      rowMsk.reassign((StVector&) logvecs(posFvec++,RangeFull::get()));
      compFMulti(uvec[depth],valvec);
      valvec.addscal(logw0+temp); rowMsk.addprod(1.0,valvec);
      // Sign +1
      temp+=logw1;
      tempVec.addprod(uvec[depth],off,colMsk[depth]);
      rowMsk.reassign((StVector&) logvecs(posFvec++,RangeFull::get()));
      compFMulti(tempVec,valvec);
      valvec.addscal(temp); rowMsk.addprod(1.0,valvec);
      // Sign -1
      tempVec.addprod(uvec[depth],-off,colMsk[depth]);
      rowMsk.reassign((StVector&) logvecs(posFvec++,RangeFull::get()));
      compFMulti(tempVec,valvec);
      valvec.addscal(temp); rowMsk.addprod(1.0,valvec);
    }
  }

  /*
   * Not efficient, just for prototypes!
   */
  inline void GaussProdQuadrature::compMultiSpecDiag(const StVector& mu,
						     const StVector& cfact,
						     StVector& logz,
						     double logpart)
  {
    int i,numcomp=getFNum();

    if (mu.size()!=d || cfact.size()!=d)
      throw WrongDimensionException(EXCEPT_MSG(""));
    logvecs.zeros(fvec.size(),numcomp);
    posFvec=0; // row pos. into 'logvecs'
    uvec[0]=mu; // depth -1
    logwval[0]=0.0; // depth -1
    cfactCopy.zeros(d);
    ((StVector&) cfactCopy.diag())=cfact; // diagonal
    for (i=0; i<d; i++)
      colMsk[i].reassign((StVector&) cfactCopy(RangeFull::get(),i));
    depth=0;
    compMultiSpecDiagHelper(); // Does all the work (compute 'logvecs')
    if (posFvec!=logvecs.rows())
      throw InternalException(EXCEPT_MSG("")); // Sanity check

    logz.zeros(numcomp);
    StVector lvMsk;
    for (i=0; i<numcomp; i++) {
      lvMsk.reassign((StVector&) logvecs(RangeFull::get(),i));
      lvMsk.addscal(-logpart);
      logz[i]=lvMsk.stableLogSum();
    }
  }

  /*
   * Not efficient, just for prototypes!
   */
  inline void GaussProdQuadrature::compMultiSpecDiagHelper()
  {
    double temp;
    if (depth<d-1) {
      // Recursive calls
      depth++;
      // Sign 0
      uvec[depth]=uvec[depth-1];
      temp=logwval[depth-1];
      logwval[depth]=temp+logw0;
      compMultiSpecDiagHelper();
      // Sign +1
      uvec[depth].addprod(uvec[depth-1],off,colMsk[depth-1]);
      logwval[depth]=temp+logw1;
      compMultiSpecDiagHelper();
      // Sign -1
      uvec[depth].addprod(uvec[depth-1],-off,colMsk[depth-1]);
      compMultiSpecDiagHelper();
      depth--;
    } else {
      // The real work: evaluations
      // Sign 0
      temp=logwval[depth];
      rowMsk.reassign((StVector&) logvecs(posFvec++,RangeFull::get()));
      compFMulti(uvec[depth],valvec);
      valvec.addscal(logw0+temp+compF(uvec[depth]));
      rowMsk.addprod(1.0,valvec);
      // Sign +1
      temp+=logw1;
      tempVec.addprod(uvec[depth],off,colMsk[depth]);
      rowMsk.reassign((StVector&) logvecs(posFvec++,RangeFull::get()));
      compFMulti(tempVec,valvec);
      valvec.addscal(temp+compF(tempVec));
      rowMsk.addprod(1.0,valvec);
      // Sign -1
      tempVec.addprod(uvec[depth],-off,colMsk[depth]);
      rowMsk.reassign((StVector&) logvecs(posFvec++,RangeFull::get()));
      compFMulti(tempVec,valvec);
      valvec.addscal(temp+compF(tempVec));
      rowMsk.addprod(1.0,valvec);
    }
  }

  /*
   * Similar to 'compLogPart', but integrand is f(u) here. 'fvec' returns
   * values w_i f(u).
   */
  inline double GaussProdQuadrature::compExpectLog(const StVector& mu,
						   const StMatrix& cfact)
  {
    if (mu.size()!=d || cfact.rows()!=d || cfact.cols()!=d)
      throw WrongDimensionException("mu, cfact have wrong size");

    posFvec=0; // pos. into 'fvec'
    uvec[0]=mu; // depth -1
    logwval[0]=0.0; // depth -1
    cfactCopy=cfact; // 'cfactCopy' is not a mask, 'cfact' may be one!
    for (int i=0; i<d; i++)
      colMsk[i].reassign((StVector&) cfactCopy(RangeFull::get(),i));
    depth=0;
    compExpectLogHelper(); // Does all the work (compute 'fvec')
    if (posFvec!=fvec.size())
      throw InternalException(EXCEPT_MSG("")); // Sanity check

    fvecFlag=2;
    return fvec.sum();
  }

  /*
   * Same as 'compLogPartHelper', but (w_i f(u)) is built in 'fvec'.
   */
  inline void GaussProdQuadrature::compExpectLogHelper()
  {
    double temp;
    if (depth<d-1) {
      // Recursive calls
      depth++;
      // Sign 0
      uvec[depth]=uvec[depth-1];
      temp=logwval[depth-1];
      logwval[depth]=temp+logw0;
      compExpectLogHelper();
      // Sign +1
      uvec[depth].addprod(uvec[depth-1],off,colMsk[depth-1]);
      logwval[depth]=temp+logw1;
      compExpectLogHelper();
      // Sign -1
      uvec[depth].addprod(uvec[depth-1],-off,colMsk[depth-1]);
      compExpectLogHelper();
      depth--;
    } else {
      // The real work: evaluations
      // Sign 0
      temp=logwval[depth];
      fvec[posFvec++]=compF(uvec[depth])*exp(logw0+temp);
      // Sign +1
      temp=exp(temp+logw1);
      tempVec.addprod(uvec[depth],off,colMsk[depth]);
      fvec[posFvec++]=compF(tempVec)*temp;
      // Sign -1
      tempVec.addprod(uvec[depth],-off,colMsk[depth]);
      fvec[posFvec++]=compF(tempVec)*temp;
    }
  }

  /*
   * Same as first part of 'compMoments', except that vectors are not
   * exponentiated here. 'fvec' contains values w_i f(u)
   */
  inline void GaussProdQuadrature::compLogMoments(const StVector& mu,
						  const StMatrix& cfact,
						  StVector& cvec,
						  StMatrix& gmat)
  {
    register int l;
    register double sum;
    register const double* fP;
    int i,j,k1,k2;
    double temp;

    if (mu.size()!=d || cfact.rows()!=d || cfact.cols()!=d)
      throw WrongDimensionException(EXCEPT_MSG("mu,cfact"));
    if (fvecFlag!=2)
      throw WrongStatusException(EXCEPT_MSG("'compExpectLog' has to be called"));
    StVector colMsk;
    cvec.zeros(d); gmat.zeros(d); // accus
    gmat.setStrctPatt(WriteBackMat<double>::strctNormal);
    ArrayHandle<int> pow3(d);
    for (i=d-1,j=1; i>=0; i--,j*=3) pow3[i]=j; // 3^(d-i-1)
    expFvec=fvec; // copy of 'fvec' (which is not changed here)
    ArrayHandle<double> expFA=expFvec.getFlatBuff();
    expTilf.zeros(fvec.size()/3);
    ArrayHandle<double> expTilfA=expTilf.getFlatBuff();

    for (i=0; i<d; i++) {
      /* Here, 'expFA' has size 3^(d-i), obtained from orig. one by
	 summing over dim. < i */
      colMsk.reassign((StVector&) gmat(RangeFull::get(),i));
      // Sign i: +1
      k1=pow3[i];
      // Copy +1 part to 'expTilfA'
      ArrayUtils<double>::copy(expTilfA,expFA+k1,k1);
      for (j=i+1; j<d; j++) {
	/* Here, 'expTilfA' has size 3^(d-j), obtained from orig. one by
	   summing over dim. < j, except for i: +1 */
	// Sign j: +1
	k2=pow3[j];
	for (l=0,fP=expTilfA.p()+k2,sum=0.0; l<k2; l++) sum+=*(fP++);
	colMsk[j]+=sum; // +,+
	// Sign j: -1
	for (l=0,sum=0.0; l<k2; l++) sum+=*(fP++);
	colMsk[j]-=sum; // +,-
	// Marginalize 'expTilfA'
	marginalize(expTilfA,k2);
      }
      temp=*(expTilfA.p());
      cvec[i]+=temp; colMsk[i]+=temp; // +
      // Sign i: -1
      // Copy -1 part to 'expTilfA'
      ArrayUtils<double>::copy(expTilfA,expFA+2*k1,k1);
      for (j=i+1; j<d; j++) {
	/* Here, 'expTilfA' has size 3^(d-j), obtained from orig. one by
	   summing over dim. < j, except for i: -1 */
	// Sign j: +1
	k2=pow3[j];
	for (l=0,fP=expTilfA.p()+k2,sum=0.0; l<k2; l++) sum+=*(fP++);
	colMsk[j]-=sum; // -,+
	// Sign j: -1
	for (l=0,sum=0.0; l<k2; l++) sum+=*(fP++);
	colMsk[j]+=sum; // -,-
	// Marginalize 'expTilfA'
	marginalize(expTilfA,k2);
      }
      temp=*(expTilfA.p());
      cvec[i]-=temp; colMsk[i]+=temp; // -
      // Marginalize 'expFA'
      marginalize(expFA,k1);
    }

    cvec.prod(off);
    gmat.smul(off*off);
    gmat.makeSymm(true);
  }

  /*
   * See 'compMomentsDiag'.
   */
  inline void GaussProdQuadrature::compLogMomentsDiag(const StVector& mu,
						      const StVector& cfact,
						      StVector& cvec,
						      StVector& gvec,
						      double* elog)
  {
    register int l;
    register double sum;
    register const double* fP;
    int i,j,k1;

    if (elog!=0 && logMomentsProc()!=0)
      throw NotImplemException(EXCEPT_MSG(""));
    if (mu.size()!=d || cfact.size()!=d)
      throw WrongDimensionException(EXCEPT_MSG("mu,cfact"));
    if (fvecFlag!=2)
      throw WrongStatusException(EXCEPT_MSG("'compExpectLog' has to be called"));
    cvec.zeros(d); gvec.zeros(d); // accus
    ArrayHandle<int> pow3(d);
    for (i=d-1,j=1; i>=0; i--,j*=3) pow3[i]=j; // 3^(d-i-1)
    expFvec=fvec; // copy of 'fvec' (which is not changed here)
    ArrayHandle<double> expFA=expFvec.getFlatBuff();

    for (i=0; i<d; i++) {
      /* Here, 'expFA' has size 3^(d-i), obtained from orig. one by
	 summing over dim. < i */
      // Sign i: +1
      k1=pow3[i];
      for (l=0,fP=expFA.p()+k1,sum=0.0; l<k1; l++) sum+=*(fP++);
      cvec[i]+=sum; gvec[i]+=sum; // +
      // Sign i: -1
      for (l=0,sum=0.0; l<k1; l++) sum+=*(fP++);
      cvec[i]-=sum; gvec[i]+=sum; // -
      // Marginalize 'expFA'
      marginalize(expFA,k1);
    }

    cvec.prod(off);
    gvec.prod(off*off);
  }
//ENDNS

#endif
