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
 * Desc.:  Header class SoftmaxGaussQuadrature
 * ------------------------------------------------------------------- */

#ifndef SOFTMAXGAUSSQUADRATURE_H
#define SOFTMAXGAUSSQUADRATURE_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/quad/default.h"
#include "lhotse/quad/GaussProdQuadrature.h"
#include "lhotse/quad/SoftmaxADFMoments.h"
#include "lhotse/matrix/FastUtils.h"
#include "lhotse/matrix/ArrayUtils.h"

//BEGINNS(quad)
  /**
   * Combines 'GaussProdQuadrature' and 'SoftmaxADFMoments'.
   * The methods 'compLogPart', 'compLogPartDiag', 'compLogPartMulti',
   * 'compLogPartMultiDiag', 'compExpectLog', 'compExpectLogDiag' have
   * special implementations here which are much faster than the generic
   * code (up to a factor of C in complexity, and they can be inlined more
   * efficiently).
   * ==> Optimizing compilation very important here!
   * <p>
   * NOTE: The C function 'log1p' (math.h) computes log(1+x) in a way
   * which is stable even if x close to 0. Apparently, 'log' may not be
   * stable in this case.
   * ==> ATTENTION: 'log1p' is not in the ANSI standard!
   * <p>
   * 'fvecFlag' has the following stati:
   * - 0: nothing
   * - 1: Precomp. from 'compLogPart' in 'fvec', used by 'compMoments'
   * - 2: Precomp. from 'compExpectLog' in 'fvec', used by 'compLogMoments'
   * - 3: Precomp. from 'compExpectLogMultiLogPartMultiDiag' in 'logvecs2'
   *      (ExpectLog), used by 'compLogMomentsMultiDiag'
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class SoftmaxGaussQuadrature : public GaussProdQuadrature,
				 public SoftmaxADFMoments
  {
  protected:
    // Additional members

    ArrayHandle<ArrayHandle<double> > vvecP; // For 'compLogPart'
    ArrayHandle<ArrayHandle<double> > lvecP; // "
    ArrayHandle<double> vzero,vplus,vminus;  // For 'compLogPartDiag'
    ArrayHandle<double> logsumexp;           // "
    double v_cls;                            // "
    StVector tempVec2;                       // "
    ArrayHandle<double> vcurr;               // For 'compLogPartMultiDiag'
    ArrayHandle<double> logvecs2;            // "
    double* posLV2;                          // "
    ArrayHandle<double> logvecs3;            // For 'compExpectLogMultiLogPartMultiDiag'
    double* posLV3;                          // "

  public:
    // Public methods

    /**
     * Constructor.
     *
     * @param dim Dimension of u
     */
    SoftmaxGaussQuadrature(int dim) : GaussProdQuadrature(dim),
				      SoftmaxADFMoments(dim) {
      int i;

      vzero.changeRep(dim);
      vplus.changeRep(dim);
      vminus.changeRep(dim);
      vcurr.changeRep(dim);
      logsumexp.changeRep(dim);
      tempVec2.zeros(dim);
      vvecP.changeRep(dim);
      lvecP.changeRep(dim);
      for (i=0; i<dim; i++) {
	vvecP[i]=uvec[i].getFlatBuff();
	lvecP[i]=cfactCopy(Range(i,dim-1),i)->getFlatBuff();
      }
    }

    double compLogPart(const StVector& mu,const StMatrix& cfact);

    double compLogPartDiag(const StVector& mu,const StVector& cfact);

    void compLogPartMulti(const StVector& mu,const StMatrix& cfact,
			  StVector& logz);

    void compLogPartMultiDiag(const StVector& mu,const StVector& cfact,
			      StVector& logz);

    /**
     * ATTENTION: As opposed to 'compLogPartMultiDiag', the current class
     * has to be set here ('setCls').
     */
    void compMultiSpecDiag(const StVector& mu,const StVector& cfact,
			   StVector& logz,double logpart);

    double compExpectLog(const StVector& mu,const StMatrix& cfact);

    double compExpectLogDiag(const StVector& mu,const StVector& cfact);

    void compExpectLogMultiDiag(const StVector& mu,const StVector& cfact,
				StVector& rvec);

    double compExpectLogLogPartMultiDiag(const StVector& mu,
					 const StVector& cfact,StVector& logz);

    void compExpectLogMultiLogPartMultiDiag(const StVector& mu,
					    const StVector& cfact,
					    StVector& rvec,StVector& logz);

    /**
     * See 'GaussProdQuadrature' for comments about computing log moments.
     * This is 'compLogMomentsDiag' in parallel for the different components
     * of the softmax log likelihood.
     */
    void compLogMomentsMultiDiag(const StVector& mu,const StVector& cfact,
				 StMatrix& cmat,StMatrix& gmat);

    double compExpectLogEVecsDiag(const StVector& mu,const StVector& cfact,
				  StVector& e1vec,StVector& e2vec);

    // Internal methods

    // Helper methods: 'xxxHelper' is helper method for 'xxx'

    void compLogPartHelper(double lse,double lwv);

    void compLogPartDiagHelper(double lse,double lwv);

    void compLogPartMultiHelper(double lse,double lwv);

    void compLogPartMultiDiagHelper(double lse,double lwv);

    void compMultiSpecDiagHelper();

    void compExpectLogHelper(double lse,double lwv);

    void compExpectLogDiagHelper(double lse,double lwv);

    void compExpectLogMultiDiagHelper(double lse,double lwv);

    void compExpectLogLogPartMultiDiagHelper(double lse,double lwv);

    void compExpectLogMultiLogPartMultiDiagHelper(double lse,double lwv);

    /**
     * Helper for 'compExpectLogEVecsDiag'. Let j=='level'. 'logvecs2' is
     * a matrix of size d-by-3^(d-j). The first j+1 rows are left
     * constant, for the rem. d-j-1 rows we do 3^(d-j-1) small logsumexp
     * inv. elements i, i+3^(d-j-1), i+2*3^(d-j-1) each, writing the
     * result into element i.
     *
     * @param level S.a.
     */
    void accumLogSumExp(int level);
  };

  // Inline methods

  /*
   * Technique similar to 'compLogPartDiagHelper', see comments there.
   * Here, v has the form
   *   v = v_0 + \sum_i s_i l_i,
   * l_i columns of a lower triangular matrix. That means that comp.
   * 0,...,i of v depend on the signs s_0,...,s_i only. We can build
   * 'logsumexp' sequentially as in 'compLogPartDiagHelper'. At level i,
   * we know that in deeper levels only comp. >i of v are required, so
   * these have to be propagated as well. This is done in 'uvec' (we use
   * 'vvecP' here pointing into 'uvec'). The l_i are kept in 'lvecP' and
   * are actually cols of 'cfactCopy'.
   */
  inline void SoftmaxGaussQuadrature::compLogPartHelper(double lse,double lwv)
  {
    const double* vprevP=vvecP[depth-1].p()+1,*lP=lvecP[depth].p();
    if (depth<d-1) {
      // Recursive calls (note that 'depth'>0)
      double* vactP=vvecP[depth].p();
      int len=d-depth;
      depth++; // Next level
      if (depth-1!=cls) {
	// Sign 0
	ArrayUtils<double>::copy(vactP,vprevP,len);
	compLogPartHelper(lse-log1p(exp((*vactP)+lse)),lwv+logw0);
	// Sign +1
	FastUtils::addsmul(vactP,vprevP,lP,1.0,len);
	lwv+=logw1;
	compLogPartHelper(lse-log1p(exp((*vactP)+lse)),lwv);
	// Sign -1
	FastUtils::addsmul(vactP,vprevP,lP,-1.0,len);
	compLogPartHelper(lse-log1p(exp((*vactP)+lse)),lwv);
      } else {
	// Need to set 'v_cls' in addition
	// Sign 0
	ArrayUtils<double>::copy(vactP,vprevP,len); v_cls=*vactP;
	compLogPartHelper(lse-log1p(exp(v_cls+lse)),lwv+logw0);
	// Sign +1
	FastUtils::addsmul(vactP,vprevP,lP,1.0,len); v_cls=*vactP;
	lwv+=logw1;
	compLogPartHelper(lse-log1p(exp(v_cls+lse)),lwv);
	// Sign -1
	FastUtils::addsmul(vactP,vprevP,lP,-1.0,len); v_cls=*vactP;
	compLogPartHelper(lse-log1p(exp(v_cls+lse)),lwv);
      }
      depth--; // Back to this level
    } else {
      // Final depth d-1
      if (cls!=depth) {
	// Sign 0
	double temp=*vprevP;
	fvec[posFvec++]=lse-log1p(exp(temp+lse))+v_cls+logw0+lwv;
	// Sign +1
	temp+=(*lP);
	lwv+=logw1;
	fvec[posFvec++]=lse-log1p(exp(temp+lse))+v_cls+lwv;
	// Sign -1
	temp=(*vprevP)-(*lP);
	fvec[posFvec++]=lse-log1p(exp(temp+lse))+v_cls+lwv;
      } else {
	// Sign 0
	v_cls=*vprevP;
	fvec[posFvec++]=lse-log1p(exp(v_cls+lse))+v_cls+logw0+lwv;
	// Sign +1
	v_cls+=(*lP);
	lwv+=logw1;
	fvec[posFvec++]=lse-log1p(exp(v_cls+lse))+v_cls+lwv;
	// Sign -1
	v_cls=(*vprevP)-(*lP);
	fvec[posFvec++]=lse-log1p(exp(v_cls+lse))+v_cls+lwv;
      }
    }
  }

  /*
   * See 'compLogPartHelper'. We have
   *   v = (gamma h + beta) + off gamma L s.
   * We store off gamma L in 'cfactCopy', whose lower tr. columns are masked
   * by 'lvecP'. We store gamma h + beta in 'tempVec'.
   * NOTE: 'uvec[i]' is masked by 'vvecP[i]', but we use only the first
   * d-i entries.
   * The calls here to 'compLogPartHelper' are at level 1, not 0.
   */
  inline double SoftmaxGaussQuadrature::compLogPart(const StVector& mu,
						    const StMatrix& cfact)
  {
    if (mu.size()!=d || cfact.rows()!=d || cfact.cols()!=d)
      throw WrongDimensionException("mu, cfact have wrong size");

    // Precomputations
    tempVec.addprod(bias,gamma,mu);
    uchar spatt=cfact.getStrctPatt();
    cfact.setStrctPatt(WriteBackMat<double>::strctNormal);
    cfactCopy.smul(cfact,off*gamma);
    cfact.setStrctPatt(spatt);
    posFvec=0; // pos. into 'fvec'
    double* vactP=vvecP[0].p();
    ArrayHandle<double> tP=tempVec.getFlatBuff();
    if (cls!=0) {
      // Sign 0
      uvec[0]=tempVec;
      depth=1;
      compLogPartHelper(-(*vactP),logw0);
      // Sign +1
      FastUtils::addsmul(vactP,tP,lvecP[0].p(),1.0,d);
      depth=1;
      compLogPartHelper(-(*vactP),logw1);
      // Sign -1
      FastUtils::addsmul(vactP,tP,lvecP[0].p(),-1.0,d);
      depth=1;
      compLogPartHelper(-(*vactP),logw1);
    } else {
      // Also have to set 'v_cls'
      // Sign 0
      uvec[0]=tempVec; v_cls=*vactP;
      depth=1;
      compLogPartHelper(-v_cls,logw0);
      // Sign +1
      FastUtils::addsmul(vactP,tP,lvecP[0].p(),1.0,d); v_cls=*vactP;
      depth=1;
      compLogPartHelper(-v_cls,logw1);
      // Sign -1
      FastUtils::addsmul(vactP,tP,lvecP[0].p(),-1.0,d); v_cls=*vactP;
      depth=1;
      compLogPartHelper(-v_cls,logw1);
    }
    if (posFvec!=fvec.size())
      throw InternalException(EXCEPT_MSG("")); // Sanity check
    double logz=fvec.stableLogSum();
    fvec.addscal(-logz); // complete 'fvec'
    fvecFlag=1;

    return logz;
  }

  /*
   * Same principle as recursive method 'compLogPartHelper', see comments
   * there. We use the members 'logwval', 'depth', 'posFvec' with the
   * same semantics.
   * Here, f(u) is evaluated in pieces, allowing to construct the value
   * at depth t from the value at depth t-1. For this, regard f as function
   * of v = gamma*u+beta (see 'SoftmaxADFMoments'):
   *   f(v) = v_y - log 1^T exp(v)
   * The second part lse_C(v) = -log 1^T exp(v) is built in 'logsumexp'. Its
   * "eval. at depth t" is defined as
   *   lse_t(v) = -log 1^T exp((v_0,...,v_{t-1}))
   * and constant 0 at depth 0. Now for t>0:
   *   lse_t(v) = lse_{t-1}(v) - log(1 + exp(v_t + lse_{t-1}(v))),
   * thus can be computed in O(1) from lse_{t-1}.
   * The component v_t depends on the current sign there (0,+1,-1), the
   * component values are precomputed in 'vzero', 'vplus', 'vminus'.
   * We also make sure that 'v_cls' contains v_y once depth 'd'-1 is
   * attained. Then, f(v) is assembled from 'logsumexp[depth-1]' and
   * 'v_cls'.
   * NOTE: A difference to 'compLogPartHelper': 'logwval', 'logsumexp' are
   * indexed normally here by depth, because the first call is at level 1
   * (level -1 is not used).
   * <p>
   * If the depth is <d-1: For any sign we compute the 'logsumexp' value at
   * depth t from its value at depth t-1, then do a recursive call for depth
   * t+1. At the final depth, we compute the f values (for each sign) from
   * the 'logsumexp' value at level t-1 and 'v_cls' and store them in 'fvec'.
   */
  inline void SoftmaxGaussQuadrature::compLogPartDiagHelper(double lse,
							    double lwv)
  {
    if (depth<d-1) {
      // Recursive calls (note that 'depth'>0)
      depth++; // Next level
      if (depth-1!=cls) {
	// Sign 0
	compLogPartDiagHelper(lse-log1p(exp(vzero[depth-1]+lse)),lwv+logw0);
	// Sign +1
	lwv+=logw1;
	compLogPartDiagHelper(lse-log1p(exp(vplus[depth-1]+lse)),lwv);
	// Sign -1
	compLogPartDiagHelper(lse-log1p(exp(vminus[depth-1]+lse)),lwv);
      } else {
	// Need to set 'v_cls' in addition
	// Sign 0
	v_cls=vzero[cls];
	compLogPartDiagHelper(lse-log1p(exp(v_cls+lse)),lwv+logw0);
	// Sign +1
	v_cls=vplus[cls];
	lwv+=logw1;
	compLogPartDiagHelper(lse-log1p(exp(v_cls+lse)),lwv);
	// Sign -1
	v_cls=vminus[cls];
	compLogPartDiagHelper(lse-log1p(exp(v_cls+lse)),lwv);
      }
      depth--; // Back to this level
    } else {
      // Final depth d-1
      if (cls!=depth) {
	// Sign 0
	fvec[posFvec++]=lse-log1p(exp(vzero[depth]+lse))+v_cls+logw0+lwv;
	// Sign +1
	lwv+=logw1;
	fvec[posFvec++]=lse-log1p(exp(vplus[depth]+lse))+v_cls+lwv;
	// Sign -1
	fvec[posFvec++]=lse-log1p(exp(vminus[depth]+lse))+v_cls+lwv;
      } else {
	// Sign 0
	v_cls=vzero[depth];
	fvec[posFvec++]=lse-log1p(exp(v_cls+lse))+v_cls+logw0+lwv;
	// Sign +1
	v_cls=vplus[depth];
	lwv+=logw1;
	fvec[posFvec++]=lse-log1p(exp(v_cls+lse))+v_cls+lwv;
	// Sign -1
	v_cls=vminus[depth];
	fvec[posFvec++]=lse-log1p(exp(v_cls+lse))+v_cls+lwv;
      }
    }
  }

  /*
   * Similar to 'compLogPart', but calling rec. method
   * 'compLogPartDiagHelper'. Also, we call the helper method at level 1,
   * not level 0.
   *
   * If s is the vector of current signs, the current eval. point is
   *   v = gamma (h + off L s) + beta = (gamma h + beta) + off gamma L s
   * Because L is diagonal, v_t depends on s_t only. The values for the
   * different signs are precomp. in 'vzero', 'vplus', 'vminus'.
   */
  inline double SoftmaxGaussQuadrature::compLogPartDiag(const StVector& mu,
							const StVector& cfact)
  {
    if (mu.size()!=d || cfact.size()!=d)
      throw WrongDimensionException("mu, cfact have wrong size");

    // Pre-compute v_t values (using members of 'SoftmaxADFMoments')
    tempVec.addprod(bias,gamma,mu);
    vzero.copy(tempVec.getFlatBuff());
    tempVec2.addprod(tempVec,off*gamma,cfact);
    vplus.copy(tempVec2.getFlatBuff());
    tempVec2.addprod(tempVec,-off*gamma,cfact);
    vminus.copy(tempVec2.getFlatBuff());
    posFvec=0; // pos. into 'fvec'
    if (cls!=0) {
      // Call 'compLogPartDiagHelper', sign 0
      depth=1;
      compLogPartDiagHelper(-vzero[0],logw0);
      // Call 'compLogPartDiagHelper', sign +1
      depth=1;
      compLogPartDiagHelper(-vplus[0],logw1);
      // Call 'compLogPartDiagHelper', sign -1
      depth=1;
      compLogPartDiagHelper(-vminus[0],logw1);
    } else {
      // Also have to set 'v_cls'
      // Call 'compLogPartDiagHelper', sign 0
      v_cls=vzero[0];
      depth=1;
      compLogPartDiagHelper(-v_cls,logw0);
      // Call 'compLogPartDiagHelper', sign +1
      v_cls=vplus[0];
      depth=1;
      compLogPartDiagHelper(-v_cls,logw1);
      // Call 'compLogPartDiagHelper', sign -1
      v_cls=vminus[0];
      depth=1;
      compLogPartDiagHelper(-v_cls,logw1);
    }
    if (posFvec!=fvec.size())
      throw InternalException(EXCEPT_MSG("")); // Sanity check
    double logz=fvec.stableLogSum();
    fvec.addscal(-logz); // complete 'fvec'
    fvecFlag=1;

    return logz;
  }

  /*
   * Similar to 'compLogPartMultiDiagHelper', but using the buildup of
   * v from 'compLogPartHelper'. The current v is built in 'vcurr'.
   */
  inline void SoftmaxGaussQuadrature::compLogPartMultiHelper(double lse,
							     double lwv)
  {
    double* vcu=vcurr.p()+depth;
    if (depth<d-1) {
      // Recursive calls (note that 'depth'>0)
      const double* vprevP=vvecP[depth-1].p()+1,*lP=lvecP[depth].p();
      double* vactP=vvecP[depth].p();
      int len=d-depth;
      depth++; // Next level
      // Sign 0
      ArrayUtils<double>::copy(vactP,vprevP,len);
      compLogPartMultiHelper(lse-log1p(exp(((*vcu)=(*vactP))+lse)),lwv+logw0);
      // Sign +1
      FastUtils::addsmul(vactP,vprevP,lP,1.0,len);
      lwv+=logw1;
      compLogPartMultiHelper(lse-log1p(exp(((*vcu)=(*vactP))+lse)),lwv);
      // Sign -1
      FastUtils::addsmul(vactP,vprevP,lP,-1.0,len);
      compLogPartMultiHelper(lse-log1p(exp(((*vcu)=(*vactP))+lse)),lwv);
      depth--; // Back to this level
    } else {
      // Final depth d-1
      const double* vprevP=vvecP[depth-1].p()+1,*lP=lvecP[depth].p();
      // Sign 0
      v_cls=(*vcu)=*vprevP;
      FastUtils::addscal(posLV2,vcurr.p(),
			 lse-log1p(exp(v_cls+lse))+logw0+lwv,d);
      posLV2+=d;
      // Sign +1
      (*vcu)=(v_cls+=(*lP));
      lwv+=logw1;
      FastUtils::addscal(posLV2,vcurr.p(),
			 lse-log1p(exp(v_cls+lse))+lwv,d);
      posLV2+=d;
      // Sign -1
      (*vcu)=v_cls=(*vprevP)-(*lP);
      FastUtils::addscal(posLV2,vcurr.p(),
			 lse-log1p(exp(v_cls+lse))+lwv,d);
      posLV2+=d;
    }
  }

  /*
   * Similar to 'compLogPart' and 'compLogPartMultiDiag'.
   */
  inline void SoftmaxGaussQuadrature::compLogPartMulti(const StVector& mu,
						       const StMatrix& cfact,
						       StVector& logz)
  {
    double temp;

    if (mu.size()!=d || cfact.rows()!=d || cfact.cols()!=d)
      throw WrongDimensionException("mu, cfact have wrong size");
    // Precomputations
    tempVec.addprod(bias,gamma,mu);
    uchar spatt=cfact.getStrctPatt();
    cfact.setStrctPatt(WriteBackMat<double>::strctNormal);
    cfactCopy.smul(cfact,off*gamma);
    cfact.setStrctPatt(spatt);
    if (logvecs2.isZero())
      logvecs2.changeRep(fvec.size()*d);
    posLV2=logvecs2.p(); // pos. into 'logvecs2'
    double* vactP=vvecP[0].p();
    ArrayHandle<double> tP=tempVec.getFlatBuff();
    // Sign 0
    uvec[0]=tempVec; temp=vcurr[0]=*vactP;
    depth=1;
    compLogPartMultiHelper(-temp,logw0);
    // Sign +1
    FastUtils::addsmul(vactP,tP,lvecP[0].p(),1.0,d);
    temp=vcurr[0]=*vactP;
    depth=1;
    compLogPartMultiHelper(-temp,logw1);
    // Sign -1
    FastUtils::addsmul(vactP,tP,lvecP[0].p(),-1.0,d);
    temp=vcurr[0]=*vactP;
    depth=1;
    compLogPartMultiHelper(-temp,logw1);
    if ((posLV2-logvecs2.p())!=logvecs2.size())
      throw InternalException(EXCEPT_MSG("")); // Sanity check
    /* Stable summations. The "rows" of the "matrix" 'logvecs2' are picked
       out using mask vector 'msk' */
    logz.zeros(d);
    StVector msk;
    int sz=(fvec.size()-1)*d;
    for (int i=0; i<d; i++) {
      msk.reassign(logvecs2.p(),logvecs2.size(),Range(i,sz+i,d));
      logz[i]=msk.stableLogSum();
    }
  }

  /*
   * Almost the same as 'compLogPartDiagHelper', except that the current
   * v vector is built in 'vcurr' (instead only comp. 'cls' in 'v_cls')
   * and the comp. in the final depth is a vector-scalar addition.
   * NOTE: We use 'logvecs2' instead of 'logvecs' and avoid using
   * 'StVector', calling 'FastUtils' methods directly. 'logvecs' is never
   * allocated.
   *
   * 'logvecs2' is like a matrix containing the 'd'-vector values as cols
   * (stored in column order). 'posLV2' points to the beginning of the next
   * free column.
   */

  inline void SoftmaxGaussQuadrature::compLogPartMultiDiagHelper(double lse,
								 double lwv)
  {
    double* vcu=vcurr.p()+depth;
    if (depth<d-1) {
      // Recursive calls (note that 'depth'>0)
      depth++; // Next level
      // Sign 0
      compLogPartMultiDiagHelper(lse-log1p(exp(((*vcu)=vzero[depth-1])+lse)),
				 lwv+logw0);
      // Sign +1
      lwv+=logw1;
      compLogPartMultiDiagHelper(lse-log1p(exp(((*vcu)=vplus[depth-1])+lse)),
				 lwv);
      // Sign -1
      compLogPartMultiDiagHelper(lse-log1p(exp(((*vcu)=vminus[depth-1])+lse)),
				 lwv);
      depth--; // Back to this level
    } else {
      // Final depth d-1
      // Sign 0
      v_cls=(*vcu)=vzero[depth];
      FastUtils::addscal(posLV2,vcurr.p(),
			 lse-log1p(exp(v_cls+lse))+logw0+lwv,d);
      posLV2+=d;
      // Sign +1
      lwv+=logw1;
      v_cls=(*vcu)=vplus[depth];
      FastUtils::addscal(posLV2,vcurr.p(),
			 lse-log1p(exp(v_cls+lse))+lwv,d);
      posLV2+=d;
      // Sign -1
      v_cls=(*vcu)=vminus[depth];
      FastUtils::addscal(posLV2,vcurr.p(),
			 lse-log1p(exp(v_cls+lse))+lwv,d);
      posLV2+=d;
    }
  }

  /*
   * Union of 'compLogPartDiag' and 'compLogPartMulti'.
   */
  inline void
  SoftmaxGaussQuadrature::compLogPartMultiDiag(const StVector& mu,
					       const StVector& cfact,
					       StVector& logz)
  {
    double temp;

    if (mu.size()!=d || cfact.size()!=d)
      throw WrongDimensionException("mu, cfact have wrong size");
    // Pre-compute v_t values (using members of 'SoftmaxADFMoments')
    tempVec.addprod(bias,gamma,mu);
    vzero.copy(tempVec.getFlatBuff());
    tempVec2.addprod(tempVec,off*gamma,cfact);
    vplus.copy(tempVec2.getFlatBuff());
    tempVec2.addprod(tempVec,-off*gamma,cfact);
    vminus.copy(tempVec2.getFlatBuff());
    if (logvecs2.isZero())
      logvecs2.changeRep(fvec.size()*d);
    posLV2=logvecs2.p(); // pos. into 'logvecs2'
    // Sign 0
    temp=vcurr[0]=vzero[0];
    depth=1;
    compLogPartMultiDiagHelper(-temp,logw0);
    // Sign +1
    temp=vcurr[0]=vplus[0];
    depth=1;
    compLogPartMultiDiagHelper(-temp,logw1);
    // Sign -1
    temp=vcurr[0]=vminus[0];
    depth=1;
    compLogPartMultiDiagHelper(-temp,logw1);
    if ((posLV2-logvecs2.p())!=logvecs2.size())
      throw InternalException(EXCEPT_MSG("")); // Sanity check
    /* Stable summations. The "rows" of the "matrix" 'logvecs2' are picked
       out using mask vector 'msk' */
    logz.zeros(d);
    StVector msk;
    int sz=(fvec.size()-1)*d;
    for (int i=0; i<d; i++) {
      msk.reassign(logvecs2.p(),logvecs2.size(),Range(i,sz+i,d));
      logz[i]=msk.stableLogSum();
    }
  }

  /*
   * Union of 'compLogPartMultiDiagHelper' and 'compLogPartDiagHelper'.
   * The function to evaluate is now
   *   v_{y_i} + v_y - 2 log 1^T exp(v) for all y,
   * where y_i is 'cls'. We build -log 1^T exp(v) in 'logsumexp', v in
   * 'vcurr' and v_{y_i} in 'v_cls'.
   */
  inline void SoftmaxGaussQuadrature::compMultiSpecDiagHelper()
  {
    double temp,temp2;
    double* vcu=vcurr.p()+depth;
    if (depth<d-1) {
      // Recursive calls (note that 'depth'>0)
      double* lse=logsumexp.p()+(depth-1),*lwv=logwval.p()+(depth-1);
      temp=*(lse++);
      temp2=*(lwv++);
      depth++; // Next level
      if (depth-1!=cls) {
	// Sign 0
	*lse=temp-log1p(exp(((*vcu)=vzero[depth-1])+temp));
	*lwv=temp2+logw0;
	compMultiSpecDiagHelper();
	// Sign +1
	*lse=temp-log1p(exp(((*vcu)=vplus[depth-1])+temp));
	*lwv=temp2+logw1;
	compMultiSpecDiagHelper();
	// Sign -1
	*lse=temp-log1p(exp(((*vcu)=vminus[depth-1])+temp));
	compMultiSpecDiagHelper();
      } else {
	// Need to set 'v_cls' in addition
	// Sign 0
	*lse=temp-log1p(exp((v_cls=(*vcu)=vzero[depth-1])+temp));
	*lwv=temp2+logw0;
	compMultiSpecDiagHelper();
	// Sign +1
	*lse=temp-log1p(exp((v_cls=(*vcu)=vplus[depth-1])+temp));
	*lwv=temp2+logw1;
	compMultiSpecDiagHelper();
	// Sign -1
	*lse=temp-log1p(exp((v_cls=(*vcu)=vminus[depth-1])+temp));
	compMultiSpecDiagHelper();
      }
      depth--; // Back to this level
    } else {
      // Final depth d-1
      temp=logsumexp[depth-1];
      temp2=logwval[depth-1];
      if (cls!=depth) {
	// Sign 0
	(*vcu)=vzero[depth];
	FastUtils::addscal(posLV2,vcurr.p(),v_cls+
			   2.0*(temp-log1p(exp((*vcu)+temp)))+logw0+temp2,d);
	posLV2+=d;
	// Sign +1
	temp2+=logw1;
	(*vcu)=vplus[depth];
	FastUtils::addscal(posLV2,vcurr.p(),v_cls+
			   2.0*(temp-log1p(exp((*vcu)+temp)))+temp2,d);
	posLV2+=d;
	// Sign -1
	(*vcu)=vminus[depth];
	FastUtils::addscal(posLV2,vcurr.p(),v_cls+
			   2.0*(temp-log1p(exp((*vcu)+temp)))+temp2,d);
	posLV2+=d;
      } else {
	// Sign 0
	v_cls=(*vcu)=vzero[depth];
	FastUtils::addscal(posLV2,vcurr.p(),v_cls+
			   2.0*(temp-log1p(exp(v_cls+temp)))+logw0+temp2,d);
	posLV2+=d;
	// Sign +1
	temp2+=logw1;
	v_cls=(*vcu)=vplus[depth];
	FastUtils::addscal(posLV2,vcurr.p(),v_cls+
			   2.0*(temp-log1p(exp(v_cls+temp)))+temp2,d);
	posLV2+=d;
	// Sign -1
	v_cls=(*vcu)=vminus[depth];
	FastUtils::addscal(posLV2,vcurr.p(),v_cls+
			   2.0*(temp-log1p(exp(v_cls+temp)))+temp2,d);
	posLV2+=d;
      }
    }
  }

  /*
   * Union of 'compLogPartMultiDiag' and 'compLogPartDiag'. In addition
   * to 'compLogPartMultiDiag', we have to be careful about 'v_cls' and
   * subtract 'logpart' before the final log-sum-exp.
   */
  inline void
  SoftmaxGaussQuadrature::compMultiSpecDiag(const StVector& mu,
					    const StVector& cfact,
					    StVector& logz,double logpart)
  {
    double temp;

    if (mu.size()!=d || cfact.size()!=d)
      throw WrongDimensionException("mu, cfact have wrong size");
    // Pre-compute v_t values (using members of 'SoftmaxADFMoments')
    tempVec.addprod(bias,gamma,mu);
    vzero.copy(tempVec.getFlatBuff());
    tempVec2.addprod(tempVec,off*gamma,cfact);
    vplus.copy(tempVec2.getFlatBuff());
    tempVec2.addprod(tempVec,-off*gamma,cfact);
    vminus.copy(tempVec2.getFlatBuff());
    if (logvecs2.isZero())
      logvecs2.changeRep(fvec.size()*d);
    posLV2=logvecs2.p(); // pos. into 'logvecs2'
    if (cls!=0) {
      // Sign 0
      temp=vcurr[0]=vzero[0];
      logsumexp[0]=-temp;
      logwval[0]=logw0;
      depth=1;
      compMultiSpecDiagHelper();
      // Sign +1
      temp=vcurr[0]=vplus[0];
      logsumexp[0]=-temp;
      logwval[0]=logw1;
      depth=1;
      compMultiSpecDiagHelper();
      // Sign -1
      temp=vcurr[0]=vminus[0];
      logsumexp[0]=-temp;
      depth=1;
      compMultiSpecDiagHelper();
    } else {
      // Also have to set 'v_cls'
      // Sign 0
      v_cls=vcurr[0]=vzero[0];
      logsumexp[0]=-v_cls;
      logwval[0]=logw0;
      depth=1;
      compMultiSpecDiagHelper();
      // Sign +1
      v_cls=vcurr[0]=vplus[0];
      logsumexp[0]=-v_cls;
      logwval[0]=logw1;
      depth=1;
      compMultiSpecDiagHelper();
      // Sign -1
      v_cls=vcurr[0]=vminus[0];
      logsumexp[0]=-v_cls;
      depth=1;
      compMultiSpecDiagHelper();
    }
    if ((posLV2-logvecs2.p())!=logvecs2.size())
      throw InternalException(EXCEPT_MSG("")); // Sanity check
    /* Stable summations. The "rows" of the "matrix" 'logvecs2' are picked
       out using mask vector 'msk'. Also subtract 'logpart' from vectors
       before log-sum-exp. */
    logz.zeros(d);
    StVector msk;
    int sz=(fvec.size()-1)*d;
    FastUtils::addscal(logvecs2.p(),-logpart,logvecs2.size());
    for (int i=0; i<d; i++) {
      msk.reassign(logvecs2.p(),logvecs2.size(),Range(i,sz+i,d));
      logz[i]=msk.stableLogSum();
    }
  }

  /*
   * Similar to 'compExpectLogHelper'.
   */
  inline void SoftmaxGaussQuadrature::compExpectLogHelper(double lse,
							  double lwv)
  {
    const double* vprevP=vvecP[depth-1].p()+1,*lP=lvecP[depth].p();
    if (depth<d-1) {
      // Recursive calls (note that 'depth'>0)
      double* vactP=vvecP[depth].p();
      int len=d-depth;
      depth++; // Next level
      if (depth-1!=cls) {
	// Sign 0
	ArrayUtils<double>::copy(vactP,vprevP,len);
	compExpectLogHelper(lse-log1p(exp((*vactP)+lse)),lwv+logw0);
	// Sign +1
	FastUtils::addsmul(vactP,vprevP,lP,1.0,len);
	lwv+=logw1;
	compExpectLogHelper(lse-log1p(exp((*vactP)+lse)),lwv);
	// Sign -1
	FastUtils::addsmul(vactP,vprevP,lP,-1.0,len);
	compExpectLogHelper(lse-log1p(exp((*vactP)+lse)),lwv);
      } else {
	// Need to set 'v_cls' in addition
	// Sign 0
	ArrayUtils<double>::copy(vactP,vprevP,len); v_cls=*vactP;
	compExpectLogHelper(lse-log1p(exp(v_cls+lse)),lwv+logw0);
	// Sign +1
	FastUtils::addsmul(vactP,vprevP,lP,1.0,len); v_cls=*vactP;
	lwv+=logw1;
	compExpectLogHelper(lse-log1p(exp(v_cls+lse)),lwv);
	// Sign -1
	FastUtils::addsmul(vactP,vprevP,lP,-1.0,len); v_cls=*vactP;
	compExpectLogHelper(lse-log1p(exp(v_cls+lse)),lwv);
      }
      depth--; // Back to this level
    } else {
      // Final depth d-1
      if (cls!=depth) {
	// Sign 0
	double temp=*vprevP;
	fvec[posFvec++]=(lse-log1p(exp(temp+lse))+v_cls)*exp(logw0+lwv);
	// Sign +1
	lwv=exp(lwv+logw1);
	temp+=(*lP);
	fvec[posFvec++]=(lse-log1p(exp(temp+lse))+v_cls)*lwv;
	// Sign -1
	temp=(*vprevP)-(*lP);
	fvec[posFvec++]=(lse-log1p(exp(temp+lse))+v_cls)*lwv;
      } else {
	// Sign 0
	v_cls=*vprevP;
	fvec[posFvec++]=(lse-log1p(exp(v_cls+lse))+v_cls)*exp(logw0+lwv);
	// Sign +1
	v_cls+=(*lP);
	lwv=exp(lwv+logw1);
	fvec[posFvec++]=(lse-log1p(exp(v_cls+lse))+v_cls)*lwv;
	// Sign -1
	v_cls=(*vprevP)-(*lP);
	fvec[posFvec++]=(lse-log1p(exp(v_cls+lse))+v_cls)*lwv;
      }
    }
  }

  /*
   * Similar to 'compLogPart'.
   */
  inline double SoftmaxGaussQuadrature::compExpectLog(const StVector& mu,
						      const StMatrix& cfact)
  {
    if (mu.size()!=d || cfact.rows()!=d || cfact.cols()!=d)
      throw WrongDimensionException("mu, cfact have wrong size");

    // Precomputations
    tempVec.addprod(bias,gamma,mu);
    uchar spatt=cfact.getStrctPatt();
    cfact.setStrctPatt(WriteBackMat<double>::strctNormal);
    cfactCopy.smul(cfact,off*gamma);
    cfact.setStrctPatt(spatt);
    posFvec=0; // pos. into 'fvec'
    double* vactP=vvecP[0].p();
    ArrayHandle<double> tP=tempVec.getFlatBuff();
    if (cls!=0) {
      // Sign 0
      uvec[0]=tempVec;
      depth=1;
      compExpectLogHelper(-(*vactP),logw0);
      // Sign +1
      FastUtils::addsmul(vactP,tP,lvecP[0].p(),1.0,d);
      depth=1;
      compExpectLogHelper(-(*vactP),logw1);
      // Sign -1
      FastUtils::addsmul(vactP,tP,lvecP[0].p(),-1.0,d);
      depth=1;
      compExpectLogHelper(-(*vactP),logw1);
    } else {
      // Also have to set 'v_cls'
      // Sign 0
      uvec[0]=tempVec; v_cls=*vactP;
      depth=1;
      compExpectLogHelper(-v_cls,logw0);
      // Sign +1
      FastUtils::addsmul(vactP,tP,lvecP[0].p(),1.0,d); v_cls=*vactP;
      depth=1;
      compExpectLogHelper(-v_cls,logw1);
      // Sign -1
      FastUtils::addsmul(vactP,tP,lvecP[0].p(),-1.0,d); v_cls=*vactP;
      depth=1;
      compExpectLogHelper(-v_cls,logw1);
    }
    if (posFvec!=fvec.size())
      throw InternalException(EXCEPT_MSG("")); // Sanity check

    fvecFlag=2;
    return fvec.sum();
  }

  /*
   * Similar to 'compLogPartDiagHelper'.
   */
  inline void SoftmaxGaussQuadrature::compExpectLogDiagHelper(double lse,
							      double lwv)
  {
    if (depth<d-1) {
      // Recursive calls (note that 'depth'>0)
      depth++; // Next level
      if (depth-1!=cls) {
	// Sign 0
	compExpectLogDiagHelper(lse-log1p(exp(vzero[depth-1]+lse)),logw0+lwv);
	// Sign +1
	lwv+=logw1;
	compExpectLogDiagHelper(lse-log1p(exp(vplus[depth-1]+lse)),lwv);
	// Sign -1
	compExpectLogDiagHelper(lse-log1p(exp(vminus[depth-1]+lse)),lwv);
      } else {
	// Need to set 'v_cls' in addition
	// Sign 0
	v_cls=vzero[cls];
	compExpectLogDiagHelper(lse-log1p(exp(v_cls+lse)),lwv+logw0);
	// Sign +1
	v_cls=vplus[cls];
	lwv+=logw1;
	compExpectLogDiagHelper(lse-log1p(exp(v_cls+lse)),lwv);
	// Sign -1
	v_cls=vminus[cls];
	compExpectLogDiagHelper(lse-log1p(exp(v_cls+lse)),lwv);
      }
      depth--; // Back to this level
    } else {
      // Final depth d-1
      if (cls!=depth) {
	// Sign 0
	fvec[posFvec++]=(lse-log1p(exp(vzero[depth]+lse))+v_cls)*
	  exp(logw0+lwv);
	// Sign +1
	lwv=exp(lwv+logw1);
	fvec[posFvec++]=(lse-log1p(exp(vplus[depth]+lse))+v_cls)*lwv;
	// Sign -1
	fvec[posFvec++]=(lse-log1p(exp(vminus[depth]+lse))+v_cls)*lwv;
      } else {
	// Sign 0
	v_cls=vzero[depth];
	fvec[posFvec++]=(lse-log1p(exp(v_cls+lse))+v_cls)*exp(logw0+lwv);
	// Sign +1
	v_cls=vplus[depth];
	lwv=exp(lwv+logw1);
	fvec[posFvec++]=(lse-log1p(exp(v_cls+lse))+v_cls)*lwv;
	// Sign -1
	v_cls=vminus[depth];
	fvec[posFvec++]=(lse-log1p(exp(v_cls+lse))+v_cls)*lwv;
      }
    }
  }

  /*
   * Similar to 'compLogPartDiag'.
   */
  inline double
  SoftmaxGaussQuadrature::compExpectLogDiag(const StVector& mu,
					    const StVector& cfact)
  {
    if (mu.size()!=d || cfact.size()!=d)
      throw WrongDimensionException("mu, cfact have wrong size");

    // Pre-compute v_t values (using members of 'SoftmaxADFMoments')
    tempVec.addprod(bias,gamma,mu);
    vzero.copy(tempVec.getFlatBuff());
    tempVec2.addprod(tempVec,off*gamma,cfact);
    vplus.copy(tempVec2.getFlatBuff());
    tempVec2.addprod(tempVec,-off*gamma,cfact);
    vminus.copy(tempVec2.getFlatBuff());
    posFvec=0; // pos. into 'fvec'
    if (cls!=0) {
      // Sign 0
      depth=1;
      compExpectLogDiagHelper(-vzero[0],logw0);
      // Sign +1
      depth=1;
      compExpectLogDiagHelper(-vplus[0],logw1);
      // Sign -1
      depth=1;
      compExpectLogDiagHelper(-vminus[0],logw1);
    } else {
      // Also have to set 'v_cls'
      // Sign 0
      v_cls=vzero[0];
      depth=1;
      compExpectLogDiagHelper(-v_cls,logw0);
      // Call 'compExpectLogDiagHelper', sign +1
      v_cls=vplus[0];
      depth=1;
      compExpectLogDiagHelper(-v_cls,logw1);
      // Call 'compExpectLogDiagHelper', sign -1
      v_cls=vminus[0];
      depth=1;
      compExpectLogDiagHelper(-v_cls,logw1);
    }
    if (posFvec!=fvec.size())
      throw InternalException(EXCEPT_MSG("")); // Sanity check

    fvecFlag=2;
    return fvec.sum();
  }

  /*
   * Similar to 'compLogPartMultiDiagHelper', only that 'ExpectLog' is
   * computed rather than 'LogPart'.
   */
  inline void SoftmaxGaussQuadrature::compExpectLogMultiDiagHelper(double lse,
								   double lwv)
  {
    double* vcu=vcurr.p()+depth;
    if (depth<d-1) {
      // Recursive calls (note that 'depth'>0)
      depth++; // Next level
      // Sign 0
      compExpectLogMultiDiagHelper(lse-log1p(exp(((*vcu)=vzero[depth-1])+lse)),
				   lwv+logw0);
      // Sign +1
      lwv+=logw1;
      compExpectLogMultiDiagHelper(lse-log1p(exp(((*vcu)=vplus[depth-1])+lse)),
				   lwv);
      // Sign -1
      compExpectLogMultiDiagHelper(lse-log1p(exp(((*vcu)=vminus[depth-1])+
						 lse)),lwv);
      depth--; // Back to this level
    } else {
      // Final depth d-1
      // Sign 0
      v_cls=(*vcu)=vzero[depth];
      FastUtils::addscal(posLV2,vcurr.p(),lse-log1p(exp(v_cls+lse)),d);
      FastUtils::smul(posLV2,exp(logw0+lwv),d);
      posLV2+=d;
      // Sign +1
      lwv=exp(lwv+logw1);
      v_cls=(*vcu)=vplus[depth];
      FastUtils::addscal(posLV2,vcurr.p(),lse-log1p(exp(v_cls+lse)),d);
      FastUtils::smul(posLV2,lwv,d);
      posLV2+=d;
      // Sign -1
      v_cls=(*vcu)=vminus[depth];
      FastUtils::addscal(posLV2,vcurr.p(),lse-log1p(exp(v_cls+lse)),d);
      FastUtils::smul(posLV2,lwv,d);
      posLV2+=d;
    }
  }

  /*
   * Similar to 'compLogPartMultiDiag'.
   */
  inline void
  SoftmaxGaussQuadrature::compExpectLogMultiDiag(const StVector& mu,
						 const StVector& cfact,
						 StVector& rvec)
  {
    double temp;

    if (mu.size()!=d || cfact.size()!=d)
      throw WrongDimensionException("mu, cfact have wrong size");
    // Pre-compute v_t values (using members of 'SoftmaxADFMoments')
    tempVec.addprod(bias,gamma,mu);
    vzero.copy(tempVec.getFlatBuff());
    tempVec2.addprod(tempVec,off*gamma,cfact);
    vplus.copy(tempVec2.getFlatBuff());
    tempVec2.addprod(tempVec,-off*gamma,cfact);
    vminus.copy(tempVec2.getFlatBuff());
    if (logvecs2.isZero())
      logvecs2.changeRep(fvec.size()*d);
    posLV2=logvecs2.p(); // pos. into 'logvecs2'
    // Sign 0
    temp=vcurr[0]=vzero[0];
    depth=1;
    compExpectLogMultiDiagHelper(-temp,logw0);
    // Sign +1
    temp=vcurr[0]=vplus[0];
    depth=1;
    compExpectLogMultiDiagHelper(-temp,logw1);
    // Sign -1
    temp=vcurr[0]=vminus[0];
    depth=1;
    compExpectLogMultiDiagHelper(-temp,logw1);
    if ((posLV2-logvecs2.p())!=logvecs2.size())
      throw InternalException(EXCEPT_MSG("")); // Sanity check
    // Sum up "rows" of the "matrix" 'logvecs2'
    rvec.zeros(d);
    StVector msk;
    int sz=(fvec.size()-1)*d;
    for (int i=0; i<d; i++) {
      msk.reassign(logvecs2.p(),logvecs2.size(),Range(i,sz+i,d));
      rvec[i]=msk.sum();
    }
  }

  /*
   * Union of 'compLogPartMultiDiagHelper' and 'compExpectLogDiagHelper'.
   * We build -log 1^T exp(v) in 'logsumexp', v in 'vcurr' and v_{y_i} in
   * 'v_cls'.
   */
  inline void
  SoftmaxGaussQuadrature::compExpectLogLogPartMultiDiagHelper(double lse,
							      double lwv)
  {
    double* vcu=vcurr.p()+depth;
    if (depth<d-1) {
      // Recursive calls (note that 'depth'>0)
      depth++; // Next level
      if (depth-1!=cls) {
	// Sign 0
	compExpectLogLogPartMultiDiagHelper(lse-log1p(exp(((*vcu)=vzero[depth-1])+lse)),logw0+lwv);
	// Sign +1
	lwv+=logw1;
	compExpectLogLogPartMultiDiagHelper(lse-log1p(exp(((*vcu)=vplus[depth-1])+lse)),lwv);
	// Sign -1
	compExpectLogLogPartMultiDiagHelper(lse-log1p(exp(((*vcu)=vminus[depth-1])+lse)),lwv);
      } else {
	// Need to set 'v_cls' in addition
	// Sign 0
	v_cls=(*vcu)=vzero[cls];
	compExpectLogLogPartMultiDiagHelper(lse-log1p(exp(v_cls+lse)),
					    lwv+logw0);
	// Sign +1
	v_cls=(*vcu)=vplus[cls];
	lwv+=logw1;
	compExpectLogLogPartMultiDiagHelper(lse-log1p(exp(v_cls+lse)),lwv);
	// Sign -1
	v_cls=(*vcu)=vminus[cls];
	compExpectLogLogPartMultiDiagHelper(lse-log1p(exp(v_cls+lse)),lwv);
      }
      depth--; // Back to this level
    } else {
      // Final depth d-1
      double temp;
      if (cls!=depth) {
	// Sign 0
	temp=lse-log1p(exp(((*vcu)=vzero[depth])+lse));
	fvec[posFvec++]=(temp+v_cls)*exp(logw0+lwv);
	FastUtils::addscal(posLV2,vcurr.p(),temp+logw0+lwv,d);
	posLV2+=d;
	// Sign +1
	lwv+=logw1;
	temp=lse-log1p(exp(((*vcu)=vplus[depth])+lse));
	fvec[posFvec++]=(temp+v_cls)*exp(lwv);
	FastUtils::addscal(posLV2,vcurr.p(),temp+lwv,d);
	posLV2+=d;
	// Sign -1
	temp=lse-log1p(exp(((*vcu)=vminus[depth])+lse));
	fvec[posFvec++]=(temp+v_cls)*exp(lwv);
	FastUtils::addscal(posLV2,vcurr.p(),temp+lwv,d);
	posLV2+=d;
      } else {
	// Sign 0
	v_cls=(*vcu)=vzero[depth];
	temp=lse-log1p(exp(v_cls+lse));
	fvec[posFvec++]=(temp+v_cls)*exp(logw0+lwv);
	FastUtils::addscal(posLV2,vcurr.p(),temp+logw0+lwv,d);
	posLV2+=d;
	// Sign +1
	lwv+=logw1;
	v_cls=(*vcu)=vplus[depth];
	temp=lse-log1p(exp(v_cls+lse));
	fvec[posFvec++]=(temp+v_cls)*exp(lwv);
	FastUtils::addscal(posLV2,vcurr.p(),temp+lwv,d);
	posLV2+=d;
	// Sign -1
	v_cls=(*vcu)=vminus[depth];
	temp=lse-log1p(exp(v_cls+lse));
	fvec[posFvec++]=(temp+v_cls)*exp(lwv);
	FastUtils::addscal(posLV2,vcurr.p(),temp+lwv,d);
	posLV2+=d;
      }
    }
  }

  /*
   * Union of 'compExpectLogDiag' and 'compLogPartMultiDiag'.
   */
  inline double
  SoftmaxGaussQuadrature::compExpectLogLogPartMultiDiag(const StVector& mu,
							const StVector& cfact,
							StVector& logz)
  {
    double temp;

    if (mu.size()!=d || cfact.size()!=d)
      throw WrongDimensionException("mu, cfact have wrong size");
    // Pre-compute v_t values (using members of 'SoftmaxADFMoments')
    tempVec.addprod(bias,gamma,mu);
    vzero.copy(tempVec.getFlatBuff());
    tempVec2.addprod(tempVec,off*gamma,cfact);
    vplus.copy(tempVec2.getFlatBuff());
    tempVec2.addprod(tempVec,-off*gamma,cfact);
    vminus.copy(tempVec2.getFlatBuff());
    posFvec=0; // pos. into 'fvec'
    if (logvecs2.isZero())
      logvecs2.changeRep(fvec.size()*d);
    posLV2=logvecs2.p(); // pos. into 'logvecs2'
    if (cls!=0) {
      // Sign 0
      depth=1;
      temp=vcurr[0]=vzero[0];
      compExpectLogLogPartMultiDiagHelper(-temp,logw0);
      // Sign +1
      depth=1;
      temp=vcurr[0]=vplus[0];
      compExpectLogLogPartMultiDiagHelper(-temp,logw1);
      // Sign -1
      depth=1;
      temp=vcurr[0]=vminus[0];
      compExpectLogLogPartMultiDiagHelper(-temp,logw1);
    } else {
      // Also have to set 'v_cls'
      // Sign 0
      v_cls=vcurr[0]=vzero[0];
      depth=1;
      compExpectLogLogPartMultiDiagHelper(-v_cls,logw0);
      // Call 'compExpectLogDiagHelper', sign +1
      v_cls=vcurr[0]=vplus[0];
      depth=1;
      compExpectLogLogPartMultiDiagHelper(-v_cls,logw1);
      // Call 'compExpectLogDiagHelper', sign -1
      v_cls=vcurr[0]=vminus[0];
      depth=1;
      compExpectLogLogPartMultiDiagHelper(-v_cls,logw1);
    }
    if (posFvec!=fvec.size())
      throw InternalException(EXCEPT_MSG("")); // Sanity check
    if ((posLV2-logvecs2.p())!=logvecs2.size())
      throw InternalException(EXCEPT_MSG("")); // Sanity check
    /* Stable summations. The "rows" of the "matrix" 'logvecs2' are picked
       out using mask vector 'msk' */
    logz.zeros(d);
    StVector msk;
    int sz=(fvec.size()-1)*d;
    for (int i=0; i<d; i++) {
      msk.reassign(logvecs2.p(),logvecs2.size(),Range(i,sz+i,d));
      logz[i]=msk.stableLogSum();
    }

    fvecFlag=2;
    return fvec.sum();
  }

  /*
   * Union of 'compLogPartMultiDiagHelper' and 'compExpectLogMultiDiagHelper'.
   * We build -log 1^T exp(v) in 'logsumexp', v in 'vcurr'.
   */
  inline void
  SoftmaxGaussQuadrature::compExpectLogMultiLogPartMultiDiagHelper(double lse,
								   double lwv)
  {
    double* vcu=vcurr.p()+depth;
    if (depth<d-1) {
      // Recursive calls (note that 'depth'>0)
      depth++; // Next level
      // Sign 0
      compExpectLogMultiLogPartMultiDiagHelper(lse-log1p(exp(((*vcu)=vzero[depth-1])+lse)),logw0+lwv);
      // Sign +1
      lwv+=logw1;
      compExpectLogMultiLogPartMultiDiagHelper(lse-log1p(exp(((*vcu)=vplus[depth-1])+lse)),lwv);
      // Sign -1
      compExpectLogMultiLogPartMultiDiagHelper(lse-log1p(exp(((*vcu)=vminus[depth-1])+lse)),lwv);
      depth--; // Back to this level
    } else {
      // Final depth d-1
      double temp;
      // Sign 0
      temp=lse-log1p(exp(((*vcu)=vzero[depth])+lse));
      FastUtils::addscal(posLV2,vcurr.p(),temp,d);
      FastUtils::smul(posLV2,exp(logw0+lwv),d);
      posLV2+=d;
      FastUtils::addscal(posLV3,vcurr.p(),temp+logw0+lwv,d);
      posLV3+=d;
      // Sign +1
      lwv+=logw1;
      temp=lse-log1p(exp(((*vcu)=vplus[depth])+lse));
      FastUtils::addscal(posLV2,vcurr.p(),temp,d);
      FastUtils::smul(posLV2,exp(lwv),d);
      posLV2+=d;
      FastUtils::addscal(posLV2,vcurr.p(),temp+lwv,d);
      posLV2+=d;
      // Sign -1
      temp=lse-log1p(exp(((*vcu)=vminus[depth])+lse));
      FastUtils::addscal(posLV2,vcurr.p(),temp,d);
      FastUtils::smul(posLV2,exp(lwv),d);
      posLV2+=d;
      FastUtils::addscal(posLV2,vcurr.p(),temp+lwv,d);
      posLV2+=d;
    }
  }

  /*
   * Union of 'compExpectLogMultiDiag' and 'compLogPartMultiDiag'. The
   * 'ExpectLogMulti' stuff is kept in 'logvecs2' and 'fvecFlag' is set to 3.
   */
  inline void
  SoftmaxGaussQuadrature::compExpectLogMultiLogPartMultiDiag(const StVector& mu,
							     const StVector& cfact,
							     StVector& rvec,
							     StVector& logz)
  {
    double temp;

    if (mu.size()!=d || cfact.size()!=d)
      throw WrongDimensionException("mu, cfact have wrong size");
    // Pre-compute v_t values (using members of 'SoftmaxADFMoments')
    tempVec.addprod(bias,gamma,mu);
    vzero.copy(tempVec.getFlatBuff());
    tempVec2.addprod(tempVec,off*gamma,cfact);
    vplus.copy(tempVec2.getFlatBuff());
    tempVec2.addprod(tempVec,-off*gamma,cfact);
    vminus.copy(tempVec2.getFlatBuff());
    if (logvecs2.isZero())
      logvecs2.changeRep(fvec.size()*d);
    posLV2=logvecs2.p(); // pos. into 'logvecs2'
    if (logvecs3.isZero())
      logvecs3.changeRep(fvec.size()*d);
    posLV3=logvecs3.p(); // pos. into 'logvecs3'
    // Sign 0
    depth=1;
    temp=vcurr[0]=vzero[0];
    compExpectLogMultiLogPartMultiDiagHelper(-temp,logw0);
    // Sign +1
    depth=1;
    temp=vcurr[0]=vplus[0];
    compExpectLogMultiLogPartMultiDiagHelper(-temp,logw1);
    // Sign -1
    depth=1;
    temp=vcurr[0]=vminus[0];
    compExpectLogMultiLogPartMultiDiagHelper(-temp,logw1);
    if ((posLV2-logvecs2.p())!=logvecs2.size())
      throw InternalException(EXCEPT_MSG("")); // Sanity check
    if ((posLV3-logvecs3.p())!=logvecs3.size())
      throw InternalException(EXCEPT_MSG("")); // Sanity check
    /* Stable summations. The "rows" of the "matrix" 'logvecs3' are picked
       out using mask vector 'msk' */
    logz.zeros(d);
    StVector msk;
    int sz=(fvec.size()-1)*d,i;
    for (i=0; i<d; i++) {
      msk.reassign(logvecs3.p(),logvecs3.size(),Range(i,sz+i,d));
      logz[i]=msk.stableLogSum();
    }
    // Sum up "rows" of the "matrix" 'logvecs2'
    rvec.zeros(d);
    for (i=0; i<d; i++) {
      msk.reassign(logvecs2.p(),logvecs2.size(),Range(i,sz+i,d));
      rvec[i]=msk.sum();
    }

    fvecFlag=3;
  }

  /*
   * Requires 'fvecFlag'==3. We copy and transpose 'logvecs2' into 'logvecs3',
   * then do every component like 'GaussProdQuadrature::compLogMomentsDiag'.
   * The "column" of 'logvecs3' has the role of 'expFVec' there.
   */
  inline void
  SoftmaxGaussQuadrature::compLogMomentsMultiDiag(const StVector& mu,
						  const StVector& cfact,
						  StMatrix& cmat,
						  StMatrix& gmat)
  {
    register int l;
    register double sum;
    register const double* fP;
    int i,j,k1;

    if (mu.size()!=d || cfact.size()!=d)
      throw WrongDimensionException(EXCEPT_MSG("mu,cfact"));
    if (fvecFlag!=3)
      throw WrongStatusException(EXCEPT_MSG("'compExpectLogMultiLogPartMultiDiag' has to be called"));
    cmat.zeros(d); gmat.zeros(d); // accus
    ArrayHandle<int> pow3(d);
    for (i=d-1,j=1; i>=0; i--,j*=3) pow3[i]=j; // 3^(d-i-1)
    int fsz=j; // 3^d
    if (logvecs2.size()!=d*fsz)
      throw InvalidParameterException(EXCEPT_MSG("")); // sanity check
    if (logvecs3.isZero() || logvecs3.size()!=d*fsz)
      logvecs3.changeRep(d*fsz);
    posLV3=logvecs3.p(); // pos. into 'logvecs3'
    // Copy transpose of 'logvecs2' into 'logvecs3'
    {
      StMatrix msk2; msk2.reassign(logvecs2,d,fsz,d);
      StMatrix msk3; msk3.reassign(logvecs3,fsz,d,fsz);
      msk3.trans(msk2);
    }
    // Loop over components
    StVector cvec,gvec;
    for (j=0; j<d; j++,posLV3+=fsz) {
      cvec.reassign((StVector&) cmat(RangeFull::get(),j));
      gvec.reassign((StVector&) gmat(RangeFull::get(),j));
      for (i=0; i<d; i++) {
	/* Here, 'posLV3' has size 3^(d-i), obtained from orig. one by
	   summing over dim. < i */
	// Sign i: +1
	k1=pow3[i];
	for (l=0,fP=posLV3+k1,sum=0.0; l<k1; l++) sum+=*(fP++);
	cvec[i]+=sum; gvec[i]+=sum; // +
	// Sign i: -1
	for (l=0,sum=0.0; l<k1; l++) sum+=*(fP++);
	cvec[i]-=sum; gvec[i]+=sum; // -
	// Marginalize 'expFA'
	marginalize(posLV3,k1);
      }
    }
    cmat.smul(off);
    gmat.smul(off*off);
  }

  inline void SoftmaxGaussQuadrature::accumLogSumExp(int level)
  {
    int i,j,fact,mxi;
    double e1,e2,e3,mx;
    double* p1;
    const double* p2,*p3;

    for (i=d-1,fact=1; i>level; i--,fact*=3); // 3^(d-1-level)
    p1=logvecs2.p();
    for (i=0; i<fact; i++) {
      p1+=(level+1);
      p2=p1+(fact*d); p3=p2+(fact*d);
      for (j=level+1; j<d; j++) {
	mx=e1=*p1; mxi=1;
	if ((e2=*(p2++))>mx) { mx=e2; mxi=2; }
	if ((e3=*(p3++))>mx) mxi=3;
	switch (mxi) {
	case 1:
	  *(p1++)=e1+log1p(exp(e2-e1)+exp(e3-e1));
	  break;
	case 2:
	  *(p1++)=e2+log1p(exp(e1-e2)+exp(e3-e2));
	  break;
	default:
	  *(p1++)=e3+log1p(exp(e1-e3)+exp(e2-e3));
	}
      }
    }
  }

  /*
   * Apart from post-processing, this is 'compExpectLogLogPartMultiDiag'.
   * The accus for 'LogPartMulti' are stored in 'logvecs2', and they are
   * used to compute (E[s_y P(y | u_i)])_y as well (for 'e2vec'). For the
   * latter, we loop over y. For each y, we use 2 logsumexp to accum.
   * elements for s_y=+1, s_y=-1. We then use 'accumLogSumExp' to sum out
   * level y in 'logvecs2'.
   */
  inline double
  SoftmaxGaussQuadrature::compExpectLogEVecsDiag(const StVector& mu,
						 const StVector& cfact,
						 StVector& e1vec,
						 StVector& e2vec)
  {
    double temp;

    if (mu.size()!=d || cfact.size()!=d)
      throw WrongDimensionException("mu, cfact have wrong size");
    // Pre-compute v_t values (using members of 'SoftmaxADFMoments')
    tempVec.addprod(bias,gamma,mu);
    vzero.copy(tempVec.getFlatBuff());
    tempVec2.addprod(tempVec,off*gamma,cfact);
    vplus.copy(tempVec2.getFlatBuff());
    tempVec2.addprod(tempVec,-off*gamma,cfact);
    vminus.copy(tempVec2.getFlatBuff());
    posFvec=0; // pos. into 'fvec'
    if (logvecs2.isZero())
      logvecs2.changeRep(fvec.size()*d);
    posLV2=logvecs2.p(); // pos. into 'logvecs2'
    if (cls!=0) {
      // Sign 0
      depth=1;
      temp=vcurr[0]=vzero[0];
      compExpectLogLogPartMultiDiagHelper(-temp,logw0);
      // Sign +1
      depth=1;
      temp=vcurr[0]=vplus[0];
      compExpectLogLogPartMultiDiagHelper(-temp,logw1);
      // Sign -1
      depth=1;
      temp=vcurr[0]=vminus[0];
      compExpectLogLogPartMultiDiagHelper(-temp,logw1);
    } else {
      // Also have to set 'v_cls'
      // Sign 0
      v_cls=vcurr[0]=vzero[0];
      depth=1;
      compExpectLogLogPartMultiDiagHelper(-v_cls,logw0);
      // Call 'compExpectLogDiagHelper', sign +1
      v_cls=vcurr[0]=vplus[0];
      depth=1;
      compExpectLogLogPartMultiDiagHelper(-v_cls,logw1);
      // Call 'compExpectLogDiagHelper', sign -1
      v_cls=vcurr[0]=vminus[0];
      depth=1;
      compExpectLogLogPartMultiDiagHelper(-v_cls,logw1);
    }
    if (posFvec!=fvec.size())
      throw InternalException(EXCEPT_MSG("")); // Sanity check
    if ((posLV2-logvecs2.p())!=logvecs2.size())
      throw InternalException(EXCEPT_MSG("")); // Sanity check
    // Vector e1
    // Stable summations. The "rows" of the "matrix" 'logvecs2' are picked
    // out using mask vector 'msk'
    e1vec.zeros(d);
    StVector msk;
    int i,sz=(fvec.size()-1)*d;
    for (i=0; i<d; i++) {
      msk.reassign(logvecs2,Range(i,sz+i,d));
      e1vec[i]=exp(msk.stableLogSum());
    }
    e1vec[cls]-=1.0;
    // Vector e2
    // See comment above
    e2vec.zeros(d);
    int fact;
    for (i=0,fact=1; i<d-1; i++,fact*=3); // 3^(d-1)
    for (i=0; i<d; i++) {
      // Elements 0,...,3*fact-1 of rows >= i in 'logvecs2' are valid
      // 'e2vec[i]' is obtained from row i. 'fact' is 3^(d-i-1)
      // Sign +
      msk.reassign(logvecs2,Range(fact*d+i,(2*fact-1)*d+i,d));
      e2vec[i]=exp(msk.stableLogSum());
      // Sign -
      msk.reassign(logvecs2,Range(2*fact*d+i,(3*fact-1)*d+i,d));
      e2vec[i]=(e2vec[i]-exp(msk.stableLogSum()))*0.5*off/cfact[i];
      // Sum out level i
      accumLogSumExp(i);
      fact/=3;
    }

    fvecFlag=2;
    return fvec.sum();
  }

//ENDNS

#endif
