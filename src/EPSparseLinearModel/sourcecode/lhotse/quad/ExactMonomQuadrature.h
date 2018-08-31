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
 * Desc.:  Header class ExactMonomQuadrature
 * ------------------------------------------------------------------- */

#ifndef EXACTMONOMQUADRATURE_H
#define EXACTMONOMQUADRATURE_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/quad/default.h"
#include "lhotse/quad/MultiQuadrature.h"
#include "lhotse/matrix/StVector.h"
#include "lhotse/matrix/StMatrix.h"

//BEGINNS(quad)
  /**
   * Uses method of exact monomials (see Davis/Rabinowitz: Methods of
   * Numerical Integration, p. 369) to approximate mean and covariance
   * matrix of a tilted Gaussian prop. to
   *   (exp f(u)) N(u | mu,sigma),
   * where u is fairly low-dimensional (for 1-dim. u, Gaussian quadrature
   * should be used!) and f, mu, sigma can be specified. We use a rule of
   * precision 5, requiring 2*d^2 + 1 function evaluations, d the dimension
   * of u. f has to be smooth. It has to be provided in 'compF'.
   * <p>
   * NOTE: This is the precision 5 rule of McNamee/Stenger for N(0,I)
   * as fully symmetric weight function. It is exact for monomials up to
   * total degree 5. For d>=5, weight w_1 is negative which is considered
   * undesirable (see Davis/Rabinowitz). The accuracy for this rule and
   * esp. rules of even higher precision is unclear.
   * NOTE: A more accurate alternative are product rules which however need
   * a number of evaluations growing exponential in d.
   * <p>
   * Possible alternatives:
   * Genz (SIAM J. Numer. Anal. 23, p. 1273ff) describes a rule for which
   * the weights can be computed analytically given the generators, the
   * latter are chosen as "Patterson" sequence. This rule is reported to
   * behave better ==> TODO: UNDERSTAND AND IMPLEMENT!
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class ExactMonomQuadrature : public virtual MultiQuadrature
  {
  protected:
    // Members

    int d;                        // dimension of u
    double off,logw0,logw1,logw2; // parameters for exact monomials rule
    int w1Stat;                   // w1: 0: <0; 1: >0; 2: ==0
    StVector fvec;                // transfer information
    int fvecFlag;                 /* 0: invalid; 1: 'compMoments';
				     2: 'compLogMoments' */

  public:
    // Public methods

    /**
     * Constructor. Initializes the parameters for the quadrature routine.
     * How to obtain the weights? The rule is automatically exact for all
     * monomials which have at least one odd degree. There are essentially
     * 4 different monomials with all degrees even: 1, x_1^2, x_1^4,
     * x_1^2 x_2^2 (with total order <= 5). Gives 4 (nonlinear!) equations
     * for u, w_0, w_1, w_2, the system has a unique solution for N(0,I)
     * as weight function.
     *
     * @param dim Dimension of u
     */
    ExactMonomQuadrature(int dim) : d(dim),fvecFlag(0) {
      //mexPrintf("ExactM., this=%d\n",this);
      if (dim<1) throw InvalidParameterException("'dim' must be positive");
      off=sqrt(3.0);
      double dd=(double) dim;
      logw0=log(1.0+dd*(dd-7.0)/18.0);
      logw2=-log(36.0);
      if (dim>4) {
	logw1=log((dd-4.0)/18.0);
	w1Stat=0; // <0
      } else if (dim<4) {
	logw1=log((4.0-dd)/18.0);
	w1Stat=1; // >0
      } else
	w1Stat=2; // ==0
    }

    int getDim() const {
      return d;
    }

    /**
     * This rule has negative weights for dimension >= 5.
     */
    bool hasNonnegWeights() const {
      return (d<=4);
    }

    /**
     * We store the vector
     *   (f(u) + log(|w_i|) - log Z)
     * in 'fvec' (req. by 'compMoments').
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
     * The vector (w_i f(u)) is stored in 'fvec', required by
     * 'compLogMoments'.
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
     * ATTENTION: Log moments procedure 1 not implemented!
     */
    void compLogMomentsDiag(const StVector& mu,const StVector& cfact,
			    StVector& cvec,StVector& gvec,double* elog=0);

    /**
     * ATTENTION: Not implemented!
     */
    double compExpectLogEVecsDiag(const StVector& mu,const StVector& cfact,
				  StVector& e1vec,StVector& e2vec) {
      throw NotImplemException(EXCEPT_MSG(""));
    }

  protected:
    // Internal methods

    void resetTransfer() {
      fvecFlag=0;
    }
  };

  // Inline methods

  /*
   * The precision 5 rule has the form
   *   w_0 g(0) + w_1 g[off] + w_2 g[off,off],
   * where g is the integrand relative to N(0,I). Here,
   *   g[...] = \sum_{x\in [...]} g(x)
   * and [a1,...,ar] means the set of all points obtained from
   *   (a1,...,ar,0,..,0) in \R^d
   * by arbitrary permutations and sign flips. g[off] requires 2*d, g[off,off]
   * 2*d*(d-1) evalutions of g.
   */

  /*
   * For the log partition function,
   *   g(x) = exp( f(L*x + mu) ).
   * We use stable summation for the different terms. The eval. points x
   * are sparse, i.e. 'off' times 0/1 vectors with <=2 1s. Thus just have
   * to add <=2 columns of L to mu.
   * In 'fvec', we return (f(u) + log(|w_i|) - log Z), u the eval. points used
   * by the rule (the one for w0: mu, the ones for w2 (2*d*(d-1) points), the
   * ones for w1 (2*(d-1) points, only if w1!=0).
   */
  inline double ExactMonomQuadrature::compLogPart(const StVector& mu,
						  const StMatrix& cfact)
  {
    int i,j,pos1,pos2;

    if (mu.size()!=d || cfact.rows()!=d || cfact.cols()!=d)
      throw WrongDimensionException("mu, cfact have wrong size");
    StVector u1vec(d),u2vec(d),mask1,mask2;
    fvec.zeros((w1Stat!=2)?(2*d*d+1):(2*d*(d-1)+1));

    pos1=0; pos2=2*d*(d-1)+1; // terms corr. to w1 start at 'pos2' (if w1!=0)
    fvec[pos1++]=compF(mu)+logw0;
    for (i=0; i<d; i++) {
      mask1.reassign(cfact(RangeFull::get(),i)); // column of L
      u1vec.addprod(mu,off,mask1);
      if (w1Stat!=2)
	fvec[pos2++]=compF(u1vec)+logw1;
      for (j=i+1; j<d; j++) {
	mask2.reassign(cfact(RangeFull::get(),j)); // column of L
	u2vec.addprod(u1vec,off,mask2);
	fvec[pos1++]=compF(u2vec)+logw2;
	u2vec.addprod(u1vec,-off,mask2);
        fvec[pos1++]=compF(u2vec)+logw2;
      }
      u1vec.addprod(mu,-off,mask1);
      if (w1Stat!=2)
	fvec[pos2++]=compF(u1vec)+logw1;
      for (j=i+1; j<d; j++) {
	mask2.reassign(cfact(RangeFull::get(),j)); // column of L
	u2vec.addprod(u1vec,off,mask2);
	fvec[pos1++]=compF(u2vec)+logw2;
	u2vec.addprod(u1vec,-off,mask2);
	fvec[pos1++]=compF(u2vec)+logw2;
      }
    }

    double logz;
    if (w1Stat==0) {
      // w1 negative
      double posPart=fvec(Range(0,2*d*(d-1)))->stableLogSum();
      double negPart=fvec(Range(2*d*(d-1)+1))->stableLogSum(); // w1 part
      if (posPart<=negPart)
	throw NumericalException("ExactMonomQuadrature::compLogPart: Integral not positive!");
      logz=posPart+log(1.0-exp(negPart-posPart)); // stable
    } else {
      // w1 positive or w1==0 (no w1 part in latter case)
      logz=fvec.stableLogSum();
    }
    fvec.addscal(-logz);
    fvecFlag=1;

    return logz;
  }

  /*
   * For the mean,
   *   g(x) = exp( f(u) - log Z) * u
   * for the covariance matrix
   *   g(x) = exp( f(u) - log Z) * u*u^T,
   * where u = L*x + mu.
   * If
   *   alpha(x) = exp( f(u) - log Z),
   * we accumulate a matrix alpha(x)*x*x^T, a vector alpha(x)*x and
   * a scalar alpha(x).
   * x has simple structure, =0 for w0, off*(+- e_i) for w1 (!=0),
   * off*(+- e_i +- e_j), i<j, for w2. We accumulate the e_i terms.
   *
   * 'fvec' contains values of f(u) - log Z + log(|w_i|).
   */
  inline void ExactMonomQuadrature::compMoments(const StVector& mu,
						const StMatrix& cfact,
						StVector& tilmu,
						StMatrix& tilcov,bool raw)
  {
    int i,j,pos1,pos2;
    double alphap,alpham,scalacc,muacc,temp;

    if (mu.size()!=d || cfact.rows()!=d || cfact.cols()!=d)
      throw WrongDimensionException(EXCEPT_MSG("mu,cfact"));
    if (fvecFlag!=1)
      throw WrongStatusException(EXCEPT_MSG("'compLogPart' has to be called"));
    StVector colMsk;
    tilmu.zeros(d); tilcov.zeros(d); // accus
    tilcov.setStrctPatt(WriteBackMat<double>::strctNormal);

    pos1=0; pos2=2*d*(d-1)+1; // pos. into 'fvec'
    scalacc=exp(fvec[pos1++]);
    for (i=0; i<d; i++) {
      colMsk.reassign((StVector&) tilcov(RangeFull::get(),i));
      // First comp. +
      if (w1Stat!=2) {
	// x = +e_i
	alphap=-exp(fvec[pos2++]);
	if (w1Stat==1) alphap=-alphap;
	scalacc+=alphap;
	tilmu[i]+=alphap;
	tilcov.set(i,i,tilcov.get(i,i)+alphap); // (i,i)
      }
      muacc=0.0; // accu for comp. i of 'tilmu'
      for (j=i+1; j<d; j++) {
	// Second comp. +: x = +e_i+e_j
	alphap=exp(fvec[pos1++]);
	// Second comp. -: x = +e_i-e_j
	alpham=exp(fvec[pos1++]);
	temp=alphap+alpham;
	muacc+=temp;
	tilcov.set(j,j,tilcov.get(j,j)+temp); // (j,j)
	temp=alphap-alpham;
	tilmu[j]+=temp;
	colMsk[j]+=temp; // (j,i)
      }
      tilmu[i]+=muacc;
      colMsk[i]+=muacc; // (i,i)
      scalacc+=muacc;
      // First comp. -
      if (w1Stat!=2) {
	// x = -e_i
	alphap=-exp(fvec[pos2++]);
	if (w1Stat==1) alphap=-alphap;
	scalacc+=alphap;
	tilmu[i]-=alphap;
	tilcov.set(i,i,tilcov.get(i,i)+alphap); // (i,i)
      }
      muacc=0.0; // accu for comp. i of 'tilmu'
      for (j=i+1; j<d; j++) {
	// Second comp. +: x = -e_i+e_j
	alphap=exp(fvec[pos1++]);
	// Second comp. -: x = -e_i-e_j
	alpham=exp(fvec[pos1++]);
	temp=alphap+alpham;
	muacc+=temp;
	tilcov.set(j,j,tilcov.get(j,j)+temp); // (j,j)
	temp=alphap-alpham;
	tilmu[j]+=temp;
	colMsk[j]-=temp; // (j,i)
      }
      tilmu[i]-=muacc;
      colMsk[i]+=muacc; // (i,i)
      scalacc+=muacc;
    }
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

  inline void ExactMonomQuadrature::compMomentsDiag(const StVector& mu,
						    const StVector& cfact,
						    StVector& tilmu,
						    StVector& tildiag,
						    bool raw)
  {
    int i,j,pos1,pos2;
    double alphap,alpham,scalacc,muacc,temp;

    if (mu.size()!=d || cfact.size()!=d)
      throw WrongDimensionException(EXCEPT_MSG("mu,cfact"));
    if (fvecFlag!=1)
      throw WrongStatusException(EXCEPT_MSG("'compLogPart' has to be called"));
    tilmu.zeros(d); tildiag.zeros(d); // accus

    pos1=0; pos2=2*d*(d-1)+1; // pos. into 'fvec'
    scalacc=exp(fvec[pos1++]);
    for (i=0; i<d; i++) {
      // First comp. +
      if (w1Stat!=2) {
	// x = +e_i
	alphap=-exp(fvec[pos2++]);
	if (w1Stat==1) alphap=-alphap;
	scalacc+=alphap;
	tilmu[i]+=alphap;
	tildiag[i]+=alphap;
      }
      muacc=0.0; // accu for comp. i of 'tilmu'
      for (j=i+1; j<d; j++) {
	// Second comp. +: x = +e_i+e_j
	alphap=exp(fvec[pos1++]);
	// Second comp. -: x = +e_i-e_j
	alpham=exp(fvec[pos1++]);
	temp=alphap+alpham;
	muacc+=temp;
	tildiag[j]+=temp; // (j,j)
	temp=alphap-alpham;
	tilmu[j]+=temp; // (j,i)
      }
      tilmu[i]+=muacc;
      tildiag[i]+=muacc; // (i,i)
      scalacc+=muacc;
      // First comp. -
      if (w1Stat!=2) {
	// x = -e_i
	alphap=-exp(fvec[pos2++]);
	if (w1Stat==1) alphap=-alphap;
	scalacc+=alphap;
	tilmu[i]-=alphap;
	tildiag[i]+=alphap; // (i,i)
      }
      muacc=0.0; // accu for comp. i of 'tilmu'
      for (j=i+1; j<d; j++) {
	// Second comp. +: x = -e_i+e_j
	alphap=exp(fvec[pos1++]);
	// Second comp. -: x = -e_i-e_j
	alpham=exp(fvec[pos1++]);
	temp=alphap+alpham;
	muacc+=temp;
	tildiag[j]+=temp; // (j,j)
	temp=alphap-alpham;
	tilmu[j]+=temp; // (j,i)
      }
      tilmu[i]-=muacc;
      tildiag[i]+=muacc; // (i,i)
      scalacc+=muacc;
    }
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

  /*
   * Vectorized version of 'compLogPart' above.
   */
  inline void ExactMonomQuadrature::compLogPartMulti(const StVector& mu,
						     const StMatrix& cfact,
						     StVector& logz)
  {
    int i,j,pos1,pos2,numcomp=getFNum();

    if (mu.size()!=d || cfact.rows()!=d || cfact.cols()!=d)
      throw WrongDimensionException("mu, cfact have wrong size");
    StVector u1vec(d),u2vec(d),valvec(numcomp),mask1,mask2,rowMsk;
    i=(w1Stat!=2)?(2*d*d+1):(2*d*(d-1)+1);
    StMatrix logvecs(i,numcomp); logvecs.zeros();
    logz.zeros(numcomp);

    pos1=0; pos2=2*d*(d-1)+1;
    rowMsk.reassign((StVector&) logvecs(pos1++,RangeFull::get()));
    compFMulti(mu,valvec);
    valvec.addscal(logw0); rowMsk.addprod(1.0,valvec);
    for (i=0; i<d; i++) {
      mask1.reassign(cfact(RangeFull::get(),i)); // column of L
      u1vec.addprod(mu,off,mask1);
      if (w1Stat!=2) {
	rowMsk.reassign((StVector&) logvecs(pos2++,RangeFull::get()));
	compFMulti(u1vec,valvec);
	valvec.addscal(logw1); rowMsk.addprod(1.0,valvec);
      }
      for (j=i+1; j<d; j++) {
	mask2.reassign(cfact(RangeFull::get(),j)); // column of L
	u2vec.addprod(u1vec,off,mask2);
	rowMsk.reassign((StVector&) logvecs(pos1++,RangeFull::get()));
	compFMulti(u2vec,valvec);
	valvec.addscal(logw2); rowMsk.addprod(1.0,valvec);
	u2vec.addprod(u1vec,-off,mask2);
	rowMsk.reassign((StVector&) logvecs(pos1++,RangeFull::get()));
	compFMulti(u2vec,valvec);
	valvec.addscal(logw2); rowMsk.addprod(1.0,valvec);
      }
      u1vec.addprod(mu,-off,mask1);
      if (w1Stat!=2) {
	rowMsk.reassign((StVector&) logvecs(pos2++,RangeFull::get()));
	compFMulti(u1vec,valvec);
	valvec.addscal(logw1); rowMsk.addprod(1.0,valvec);
      }
      for (j=i+1; j<d; j++) {
	mask2.reassign(cfact(RangeFull::get(),j)); // column of L
	u2vec.addprod(u1vec,off,mask2);
	rowMsk.reassign((StVector&) logvecs(pos1++,RangeFull::get()));
	compFMulti(u2vec,valvec);
	valvec.addscal(logw2); rowMsk.addprod(1.0,valvec);
	u2vec.addprod(u1vec,-off,mask2);
	rowMsk.reassign((StVector&) logvecs(pos1++,RangeFull::get()));
	compFMulti(u2vec,valvec);
	valvec.addscal(logw2); rowMsk.addprod(1.0,valvec);
      }
    }

    if (w1Stat==0) {
      // w1 negative
      double posPart,negPart;
      for (i=0; i<numcomp; i++) {
	posPart=logvecs(Range(0,2*d*(d-1)),i)->stableLogSum();
	negPart=logvecs(Range(2*d*(d-1)+1),i)->stableLogSum();
	if (posPart<=negPart)
	  throw NumericalException("ExactMonomQuadrature::compLogPart: Integral not positive!");
	logz[i]=posPart+log(1.0-exp(negPart-posPart));
      }
    } else {
      // w1 positive or w1==0
      for (i=0; i<numcomp; i++)
	logz[i]=logvecs(RangeFull::get(),i)->stableLogSum();
    }
  }

  /*
   * Similar to 'compLogPart', but integrand is f(u) here. 'fvec' returns
   * values w_i f(u).
   */
  inline double ExactMonomQuadrature::compExpectLog(const StVector& mu,
						    const StMatrix& cfact)
  {
    int i,j,pos1,pos2;

    if (mu.size()!=d || cfact.rows()!=d || cfact.cols()!=d)
      throw WrongDimensionException("mu, cfact have wrong size");
    StVector u1vec(d),u2vec(d),mask1,mask2;
    fvec.zeros((w1Stat!=2)?(2*d*d+1):(2*d*(d-1)+1));

    pos1=0; pos2=2*d*(d-1)+1; // terms corr. to w1 start at 'pos2' (if w1!=0)
    fvec[pos1++]=compF(mu)*exp(logw0); // w_0 f(u)
    double w2=exp(logw2),w1=0.0;
    if (w1Stat==0) w1=-exp(logw1);
    else if (w1Stat==1) w1=exp(logw1);
    for (i=0; i<d; i++) {
      mask1.reassign(cfact(RangeFull::get(),i)); // column of L
      u1vec.addprod(mu,off,mask1);
      if (w1Stat!=2)
	fvec[pos2++]=compF(u1vec)*w1; // w_1 f(u)
      for (j=i+1; j<d; j++) {
	mask2.reassign(cfact(RangeFull::get(),j)); // column of L
	u2vec.addprod(u1vec,off,mask2);
	fvec[pos1++]=compF(u2vec)*w2; // w_2 f(u)
	u2vec.addprod(u1vec,-off,mask2);
	fvec[pos1++]=compF(u2vec)*w2; // w_2 f(u)
      }
      u1vec.addprod(mu,-off,mask1);
      if (w1Stat!=2)
	fvec[pos2++]=compF(u1vec)*w1; // w_1 f(u)
      for (j=i+1; j<d; j++) {
	mask2.reassign(cfact(RangeFull::get(),j)); // column of L
	u2vec.addprod(u1vec,off,mask2);
	fvec[pos1++]=compF(u2vec)*w2; // w_2 f(u)
	u2vec.addprod(u1vec,-off,mask2);
	fvec[pos1++]=compF(u2vec)*w2; // w_2 f(u)
      }
    }
    if (pos1!=2*d*(d-1)+1 || (w1Stat!=2 && pos2!=2*d*d+1))
      throw InternalException(EXCEPT_MSG("")); // Sanity check
    if (fvec.size()!=pos2)
      throw InternalException(EXCEPT_MSG("")); // Sanity check

    fvecFlag=2;
    return fvec.sum();
  }

  /*
   * Like 'compMoments', but simpler. For c, we accumulate w_i*f(u)*x.
   * For G, we accum. w_i*f(u)*x*x^T. f(u) values given by 'fvec', simple
   * structure of x eval. points is used as in 'compMoments'.
   *
   * 'fvec' contains values w_i f(u).
   */
  inline void ExactMonomQuadrature::compLogMoments(const StVector& mu,
						   const StMatrix& cfact,
						   StVector& cvec,
						   StMatrix& gmat)
  {
    int i,j,pos1,pos2;
    double alphap,alpham,muacc,temp;

    if (fvecFlag!=2)
      throw WrongStatusException(EXCEPT_MSG("'compExpectLog' has to be called"));
    StVector colMsk;
    cvec.zeros(d); gmat.zeros(d); // reset accus
    gmat.setStrctPatt(WriteBackMat<double>::strctNormal);

    /* ATTENTION: Skip first elem. in 'fvec', not needed here (as opposed
       to 'compMoments')! */
    pos1=1; pos2=2*d*(d-1)+1; // pos. into 'fvec'
    for (i=0; i<d; i++) {
      colMsk.reassign((StVector&) gmat(RangeFull::get(),i));
      // First comp. +
      if (w1Stat!=2) {
	// x = +e_i
	alphap=fvec[pos2++]; // w_1 f(u)
	cvec[i]+=alphap;
	colMsk[i]+=alphap; // (i,i)
      }
      muacc=0.0; // accu for comp. i of 'cvec'
      for (j=i+1; j<d; j++) {
	// Second comp. +: x = +e_i+e_j
	alphap=fvec[pos1++]; // w_2 f(u)
	// Second comp. -: x = +e_i-e_j
	alpham=fvec[pos1++]; // w_2 f(u)
	temp=alphap+alpham;
	muacc+=temp;
	gmat.set(j,j,gmat.get(j,j)+temp); // (j,j)
	temp=alphap-alpham;
	cvec[j]+=temp;
	colMsk[j]+=temp; // (j,i)
      }
      cvec[i]+=muacc;
      colMsk[i]+=muacc; // (i,i)
      // First comp. -
      if (w1Stat!=2) {
	// x = -e_i
	alphap=fvec[pos2++]; // w_1 f(u)
	cvec[i]-=alphap;
	colMsk[i]+=alphap; // (i,i)
      }
      muacc=0.0; // accu for comp. i of 'cvec'
      for (j=i+1; j<d; j++) {
	// Second comp. +: x = -e_i+e_j
	alphap=fvec[pos1++]; // w_2 f(u)
	// Second comp. -: x = -e_i-e_j
	alpham=fvec[pos1++]; // w_2 f(u)
	temp=alphap+alpham;
	muacc+=temp;
	gmat.set(j,j,gmat.get(j,j)+temp); // (j,j)
	temp=alphap-alpham;
	cvec[j]+=temp;
	colMsk[j]-=temp; // (j,i)
      }
      cvec[i]-=muacc;
      colMsk[i]+=muacc; // (i,i)
    }
    if (pos1!=2*d*(d-1)+1 || (w1Stat!=2 && pos2!=2*d*d+1))
      throw InternalException(EXCEPT_MSG("")); // Sanity check
    cvec.prod(off);
    gmat.smul(off*off);
    gmat.makeSymm(true);
  }

  inline void ExactMonomQuadrature::compLogMomentsDiag(const StVector& mu,
						       const StVector& cfact,
						       StVector& cvec,
						       StVector& gvec,
						       double* elog)
  {
    int i,j,pos1,pos2;
    double alphap,alpham,muacc,temp;

    if (elog!=0 && logMomentsProc()!=0)
      throw NotImplemException(EXCEPT_MSG(""));
    if (fvecFlag!=2)
      throw WrongStatusException(EXCEPT_MSG("'compExpectLog' has to be called"));
    cvec.zeros(d); gvec.zeros(d); // reset accus

    /* ATTENTION: Skip first elem. in 'fvec', not needed here (as opposed
       to 'compMoments')! */
    pos1=1; pos2=2*d*(d-1)+1; // pos. into 'fvec'
    for (i=0; i<d; i++) {
      // First comp. +
      if (w1Stat!=2) {
	// x = +e_i
	alphap=fvec[pos2++]; // w_1 f(u)
	cvec[i]+=alphap;
	gvec[i]+=alphap; // (i,i)
      }
      muacc=0.0; // accu for comp. i of 'cvec'
      for (j=i+1; j<d; j++) {
	// Second comp. +: x = +e_i+e_j
	alphap=fvec[pos1++]; // w_2 f(u)
	// Second comp. -: x = +e_i-e_j
	alpham=fvec[pos1++]; // w_2 f(u)
	temp=alphap+alpham;
	muacc+=temp;
	gvec[j]+=temp; // (j,j)
	temp=alphap-alpham;
	cvec[j]+=temp; // (j,i)
      }
      cvec[i]+=muacc;
      gvec[i]+=muacc; // (i,i)
      // First comp. -
      if (w1Stat!=2) {
	// x = -e_i
	alphap=fvec[pos2++]; // w_1 f(u)
	cvec[i]-=alphap;
	gvec[i]+=alphap; // (i,i)
      }
      muacc=0.0; // accu for comp. i of 'cvec'
      for (j=i+1; j<d; j++) {
	// Second comp. +: x = -e_i+e_j
	alphap=fvec[pos1++]; // w_2 f(u)
	// Second comp. -: x = -e_i-e_j
	alpham=fvec[pos1++]; // w_2 f(u)
	temp=alphap+alpham;
	muacc+=temp;
	gvec[j]+=temp; // (j,j)
	temp=alphap-alpham;
	cvec[j]+=temp; // (j,i)
      }
      cvec[i]-=muacc;
      gvec[i]+=muacc; // (i,i)
    }
    if (pos1!=2*d*(d-1)+1 || (w1Stat!=2 && pos2!=2*d*d+1))
      throw InternalException(EXCEPT_MSG("")); // Sanity check
    cvec.prod(off);
    gvec.prod(off*off);
  }
//ENDNS

#endif
