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
 * Desc.:  Header class MultiProbitQuadrature
 * ------------------------------------------------------------------- */

#ifndef MULTIPROBITQUADRATURE_H
#define MULTIPROBITQUADRATURE_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/quad/default.h"
#include "lhotse/quad/MultiClassQuadrature.h"
#include "lhotse/matrix/StVector.h"
#include "lhotse/matrix/StMatrix.h"
#include "lhotse/matrix/ArrayUtils.h"
#include "lhotse/matrix/FastUtils.h"
#include "lhotse/quad/GaussQuad.h"
#include "lhotse/rando/TestUtils.h"

//BEGINNS(quad)

  /**
   * Implementation of 'MultiClassQuadrature' for the multivariate probit
   * likelihood defined as
   *   P(y | u) = Pr\{ v_y >= v_c for all c \},
   * where v = u+b+n and n is N(0,I). The Pr is over n.
   * <p>
   * With this likelihood, the methods here require most one-dimensional
   * and some two-dim. Gaussian quadratures (we use a product rule). The
   * one-dim. rule is provided by 'GaussQuad'.
   * The rule is setup upon construction.
   * <p>
   * NOTE: Only the 'Diag' methods are implemented!
   * NOTE: Expect. over log likelihood deactivated! See code in
   * 'MultiProbitQuadrature_withlog'.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class MultiProbitQuadrature : public MultiClassQuadrature
  {
  protected:
    // Members

    int d;                        // Dimension of u
    int cls;                      // Current target (0,...,'getDim()'-1)
    StVector bias;                // Intercept parameters
    int smooMode;                 // 0: (mu,sigma^2); 1: (b1,b2,sigma^2)
    double smooMu,smooSigsq;      // Parameters smoothing distrib. R
    double smooB1,smooB2;         // Mode 1
    double smooKL;                // D[R || N(0,1)] (mode 0 only)
    int expLogMode;               // See constructor, 'emode'

    int quadN;                    // 1-d quad.: number of points
    ArrayHandle<double> quadW;    // 1-d quad.: log weights
    ArrayHandle<double> equadW;   // 1-d quad.: weights
    ArrayHandle<double> quadS;    // 1-d quad.: abscissas
    ArrayHandle<double> fvec;     // Info transfer phase 1 -> phase 2
    double fLogZ;                 // "
    int fvecFlag;                 // ". 0: none; 1: 'compMomentsDiag';
                                  // 2: 'compLogMomentsDiag'
    ArrayHandle<double> muDiff,sigma;

    ArrayHandle<double> a0mat,mu0vec;
    ArrayHandle<double> aypr,muypr;
    ArrayHandle<double> i1part,i2part,zpart;
    ArrayHandle<double> accu3,accu6;

  public:
    // Public methods

    /**
     * Constructor. Number of points for 1-d quadrature (Gauss-Hermite)
     * must be passed in 'n'. The 2-d quadratures are done using the
     * product rule with n^2 evaluations.
     *
     * @param dim   Dimension of u
     * @param n     Number of q points
     */
    MultiProbitQuadrature(int dim,int n);

    int getDim() const {
      return d;
    }

    bool hasParams() const {
      return false;
    }

    void setCls(int cl) {
      if (cl<0 || cl>=getDim())
	throw InvalidParameterException("Invalid class target 'cl'");
      if (cl!=cls) resetTransfer(); // reset
      cls=cl;
    }

    int getCls() const {
      return cls;
    }

    void setBias(const StVector& bvec) {
      if (bvec.size()!=getDim())
	throw WrongDimensionException("'bvec' has wrong size");
      bias=bvec;
      resetTransfer(); // reset
    }

    const StVector& getBias() const {
      return bias;
    }

    void resetBias() {
      bias.zeros(getDim());
      resetTransfer(); // reset
    }

    bool hasNonnegWeights() const {
      return true;
    }

    double compLogPart(const StVector& mu,const StMatrix& cfact) {
      throw NotImplemException(EXCEPT_MSG(""));
    }

    double compLogPartDiag(const StVector& mu,const StVector& cfact);

    void compMoments(const StVector& mu,const StMatrix& cfact,StVector& tilmu,
		     StMatrix& tilcov,bool raw) {
      throw NotImplemException(EXCEPT_MSG(""));
    }

    /**
     * ATTENTION: 'raw'==true not supported! Could be done if required.
     */
    void compMomentsDiag(const StVector& mu,const StVector& cfact,
			 StVector& tilmu,StVector& tildiag,bool raw=false);

    void compLogPartMulti(const StVector& mu,const StMatrix& cfact,
			  StVector& logz) {
      throw NotImplemException(EXCEPT_MSG(""));
    }

    void compLogPartMultiDiag(const StVector& mu,const StVector& cfact,
			      StVector& logz);

    double compExpectLog(const StVector& mu,const StMatrix& cfact) {
      throw NotImplemException(EXCEPT_MSG(""));
    }

    /**
     * Depends on smoothing distribution R ('smooXXX').
     */
    double compExpectLogDiag(const StVector& mu,const StVector& cfact) {
      switch (expLogMode) {
      case 0:
	return compExpectLogDiag0(mu,cfact);
      case 1:
      case 2:
	return compExpectLogDiag1(mu,cfact);
      }
    }

    double compExpectLogDiag0(const StVector& mu,const StVector& cfact);
    double compExpectLogDiag1(const StVector& mu,const StVector& cfact);

    /**
     * The log moments procedure is 0 for 'expLogMode'==0, 1 for
     * 'expLogMode'==1,2.
     */
    int logMomentsProc() const {
      return (expLogMode==0)?0:1;
    }

    void compLogMoments(const StVector& mu,const StMatrix& cfact,
			StVector& cvec,StMatrix& gmat) {
      throw NotImplemException(EXCEPT_MSG(""));
    }

    /**
     * Depends on smoothing distribution R ('smooXXX'). Requires
     * 'compExpectLogDiag' to be called before ('fvecFlag'==2) for log
     * moments proc. 0 ('expLogMode').
     */
    void compLogMomentsDiag(const StVector& mu,const StVector& cfact,
			    StVector& cvec,StVector& gvec,double* elog) {
      switch (expLogMode) {
      case 0:
	compLogMomentsDiag0(mu,cfact,cvec,gvec);
	break;
      case 1:
	compLogMomentsDiag1(mu,cfact,cvec,gvec,elog);
	break;
      case 2:
	compLogMomentsDiag2(mu,cfact,cvec,gvec,elog);
	break;
      }
    }

    void compLogMomentsDiag0(const StVector& mu,const StVector& cfact,
			     StVector& cvec,StVector& gvec);
    void compLogMomentsDiag1(const StVector& mu,const StVector& cfact,
			     StVector& cvec,StVector& gvec,double* elog=0);
    void compLogMomentsDiag2(const StVector& mu,const StVector& cfact,
			     StVector& cvec,StVector& gvec,double* elog=0);

    /**
     * Depends on smoothing distribution R ('smooXXX').
     */
    double compExpectLogEVecsDiag(const StVector& mu,const StVector& cfact,
				  StVector& e1vec,StVector& e2vec) {
      switch (expLogMode) {
      case 0:
	return compExpectLogEVecsDiag0(mu,cfact,e1vec,e2vec);
      case 1:
	throw NotImplemException(EXCEPT_MSG(""));
      case 2:
	return compExpectLogEVecsDiag2(mu,cfact,e1vec,e2vec);
      }
    }

    double compExpectLogEVecsDiag0(const StVector& mu,const StVector& cfact,
				   StVector& e1vec,StVector& e2vec);

    double compExpectLogEVecsDiag2(const StVector& mu,const StVector& cfact,
				   StVector& e1vec,StVector& e2vec);

    // Internal methods

    /**
     * Not a general implementation, so 'compF', 'compFMulti' not provided.
     */
    double compF(const StVector& u) const {
      throw NotImplemException(EXCEPT_MSG(""));
    }

    void compFMulti(const StVector& u,StVector& val) const {
      throw NotImplemException(EXCEPT_MSG(""));
    }

    int getFNum() const {
      return d;
    }

    void resetTransfer() {
      fvecFlag=0;
    }

    /**
     * Helper for 'compLogMomentsDiag1'. See documentation.
     * Matrix 'a0mat' is lower triangular, passed as column-scanned vector:
     * (0,0), (1,0), ..., (0,1), (1,1), ...
     * The number of y' looped over here dep. on length of 'aypr', 'muypr'.
     * The res. 'i1part', 'i2part' have the same size. If 'zbnd' is given
     * it has the same size.
     *
     * @param a0mat  4-by-4 cov. matrix A_0
     * @param mu0vec 4 mean vector mu_0
     * @param aypr   Variances a_{y'}, y'\ne y,j
     * @param muypr  Means h_{y'} + \beta_{y'}, y'\ne y,j
     * @param i1part Parts for I_1 (mean) ret. here. Sum them up
     * @param i2part Parts for I_2 (var.) ret. here. Sum them up
     * @param zbnd   Parts for z(y) lower bound ret. here. Sum up. Optional
     */
    void helperLogMomentsExt(const ArrayHandle<double>& a0mat,
			     const ArrayHandle<double>& mu0vec,
			     const ArrayHandle<double>& aypr,
			     const ArrayHandle<double>& muypr,
			     ArrayHandle<double>& i1part,
			     ArrayHandle<double>& i2part,
			     ArrayHandle<double>* zbnd=0);

    /**
     * The same as 'helperLogMomentsExt', but only for the upper left
     * 3-by-3 part of A. Method returns 'zbnd' value (summed up).
     *
     * @param a0mat  3-by-3 cov. matrix A_0
     * @param mu0vec 3 mean vector mu_0
     * @param aypr   Variances a_{y'}, y'\ne y,j
     * @param muypr  Means h_{y'} + \beta_{y'}, y'\ne y,j
     */
    double helperLogMomentsExtSmall(const ArrayHandle<double>& a0mat,
				    const ArrayHandle<double>& mu0vec,
				    const ArrayHandle<double>& aypr,
				    const ArrayHandle<double>& muypr);

    /**
     * In log moments procedures 1,2, we need to determine which Q marg.
     * is "closest" to Q(u_y), y=='cls'. This is done here using a
     * symmetrized rel. entropy.
     *
     * @param mu    Mean Q
     * @param cfact Variances Q
     * @return      Index for marg. closest to Q(u_y)
     */
    int findClosest(const StVector& mu,const StVector& cfact) const;
  };

  /*
   * The 1-d quadrature for log Z can be written as logsumexp where the
   * args are sums over y'\ne y, y the current target 'cls'. 'fvec' is a
   * 'quadN'-by-'d' matrix. In the cols of 'fvec' we store these sums over
   * y'\ne y,j, j the col. number (including the quad. weights). This is for
   * j\ne y, the column y remains undefined (it is used as summation accu
   * here). 'fvecFlag' is set to 1. We also add log N(z_j) - log Z to the
   * elements of col. j.
   */
  inline double MultiProbitQuadrature::compLogPartDiag(const StVector& mu,
						       const StVector& cfact)
  {
    int i,j;
    double temp,mudf,sig,sigy,fact,off,z;
    double* fvecP,*fsumP;

    // Initialization
    if (mu.size()!=d || cfact.size()!=d)
      throw WrongDimensionException(EXCEPT_MSG(""));
    ArrayUtils<double>::fill(fvec,0.0,quadN*d);

    // Main loop
    fvecP=fvec.p();
    temp=cfact[cls]; sigy=sqrt(1.0+temp*temp);
    for (j=0; j<d; j++) {
      if (j==cls) {
	// Skip y
	fvecP+=quadN;
	continue;
      }
      fsumP=fvec.p()+(cls*quadN); // use column y for summation accu
      mudf=mu[cls]+bias[cls]-mu[j]-bias[j]; // diff of mu's
      temp=cfact[j]; sig=sqrt(1.0+temp*temp);
      fact=sigy/sig; off=mudf/sig;
      for (i=0; i<quadN; i++) {
	z=quadS[i]*fact+off;
	*(fsumP++)+=(temp=TestUtils::logCdfNormal(z));
	*(fvecP++)=-temp+TestUtils::logPdfNormal(z);
      }
    }
    fsumP=fvec.p()+(cls*quadN);
    FastUtils::addsmul(fsumP,quadW.p(),1.0,quadN);
    for (j=0,fvecP=fvec.p(); j<d; j++,fvecP+=quadN) {
      if (j==cls) continue;
      FastUtils::addsmul(fvecP,fsumP,1.0,quadN);
    }

    // Final logsumexp
    StVector msk; msk.reassign(fvec,Range(quadN*cls,quadN*(cls+1)-1));
    fLogZ=msk.stableLogSum(); // log Z
    FastUtils::addscal(fvec.p(),-fLogZ,quadN*d); // subtract log Z
    fvecFlag=1; // info in 'fvec', 'fLogZ'

    return fLogZ;
  }

  inline void MultiProbitQuadrature::compMomentsDiag(const StVector& mu,
						     const StVector& cfact,
						     StVector& tilmu,
						     StVector& tildiag,
						     bool raw)
  {
    int i,i2,j;
    double temp,temp2,fact,off,sum,z,muacc,varacc,sig,a,n0;
    double* fvecP;

    if (raw)
      throw NotImplemException(EXCEPT_MSG("raw moments not supported"));
    if (mu.size()!=d || cfact.size()!=d)
      throw WrongDimensionException(EXCEPT_MSG("mu,cfact"));
    if (fvecFlag!=1)
      throw WrongStatusException(EXCEPT_MSG("'compLogPartDiag' has to be called"));

    // Initialization
    temp=mu[cls]+bias[cls];
    for (j=0; j<d; j++) {
      muDiff[j]=temp-mu[j]-bias[j];
      temp2=cfact[j];
      sigma[j]=sqrt(1.0+temp2*temp2);
    }
    tilmu.zeros(d); tildiag.zeros(d);

    // Components other than 'cls'
    // Means
    StVector msk;
    for (j=0; j<d; j++) {
      if (j==cls) continue;
      msk.reassign(fvec,Range(j*quadN,(j+1)*quadN-1));
      temp=cfact[j];
      tilmu[j]=mu[j]-temp*temp*exp(msk.stableLogSum())/sigma[j];
    }
    // Variances
    fvecP=fvec.p();
    for (j=0; j<d; j++) {
      if (j==cls) {
	fvecP+=quadN; // skip this column
	continue;
      }
      sig=sigma[j];
      temp=cfact[j]; a=temp*temp;
      fact=sigma[cls]/sig;
      off=muDiff[j]/sig;
      for (i=0,varacc=0.0; i<quadN; i++) {
	z=quadS[i]*fact+off;
	varacc+=z*exp(*(fvecP++));
      }
      temp=a/sig;
      tildiag[j]=-varacc*temp*temp;
      temp=tilmu[j]-mu[j];
      tildiag[j]+=(a-temp*temp);
    }

    // Component 'cls'. This needs a 2-d quadrature (nested)
    // Outer loop (n1)
    muacc=varacc=0.0;
    msk.reassign(fvec,Range(cls*quadN,(cls+1)*quadN-1));
    for (i=0; i<quadN; i++) {
      fvecP=fvec.p()+(cls*quadN); // use as temp. storage (for logsumexp)
      fact=cfact[cls]*quadS[i]; // u_y - h_y
      // Inner loop (n0)
      for (i2=0; i2<quadN; i2++) {
	n0=fact+quadS[i2];
	for (j=0,sum=0.0; j<cls; j++)
	  sum+=TestUtils::logCdfNormal((n0+muDiff[j])/sigma[j]);
	for (j=cls+1; j<d; j++)
	  sum+=TestUtils::logCdfNormal((n0+muDiff[j])/sigma[j]);
	*(fvecP++)=sum-fLogZ+quadW[i2];
      }
      temp=exp(msk.stableLogSum()+quadW[i]); // inner quad.
      fact+=mu[cls]; // u_y
      muacc+=(temp*=fact); // E[u_y ...]
      varacc+=temp*fact; // E[u_y^2 ...]
    }
    tilmu[cls]=muacc;
    tildiag[cls]=varacc-muacc*muacc;
  }

  /*
   * Similar to 'compLogPartDiag', with outer loop over y. We use col.
   * 'cls' of 'fvec'. If 'fvecFlag'==1, this does not invalidate 'fvec',
   * but if 'fvecFlag'==2, it is reset to 0.
   */
  inline void
  MultiProbitQuadrature::compLogPartMultiDiag(const StVector& mu,
					      const StVector& cfact,
					      StVector& logz)
  {
    int i,j,k;
    double temp,sum,n0;
    double* fvecP;

    // Initialization
    if (mu.size()!=d || cfact.size()!=d)
      throw WrongDimensionException(EXCEPT_MSG(""));
    logz.zeros(d);
    for (i=0; i<d; i++) {
      muDiff[i]=mu[i]+bias[i];
      temp=cfact[i];
      sigma[i]=sqrt(1.0+temp*temp);
    }

    // Main loop
    // We use the 'cls' column of 'fvec' as temp. storage. Even if 'fvecFlag'
    // ==1, this column is not used in 'compMomentsDiag'
    if (fvecFlag==2) fvecFlag=0; // reset
    StVector msk; msk.reassign(fvec,Range(cls*quadN,(cls+1)*quadN-1));
    for (j=0; j<d; j++) {
      fvecP=fvec.p()+(cls*quadN);
      for (i=0; i<quadN; i++) {
	n0=quadS[i]*sigma[j]+muDiff[j];
	for (k=0,sum=0.0; k<j; k++)
	  sum+=TestUtils::logCdfNormal((n0-muDiff[k])/sigma[k]);
	for (k=j+1; k<d; k++)
	  sum+=TestUtils::logCdfNormal((n0-muDiff[k])/sigma[k]);
	*(fvecP++)=sum+quadW[i];
      }
      logz[j]=msk.stableLogSum();
    }
  }

  /*
   * We use first d-1 elem. of 'fvec' to transfer info to 'compLogMomentsDiag'
   * ('fvecFlag' is set to 2). These comp. have the same form as the z(y)
   * lower bound (the value returned), but for the sum running over y'\ne y,j,
   * where j\ne y.
   */
  inline double
  MultiProbitQuadrature::compExpectLogDiag0(const StVector& mu,
					    const StVector& cfact)
  {
    int i,j;
    double temp,mudf,fact,sum,crit;
    double* fvecP;

    if (mu.size()!=d || cfact.size()!=d)
      throw WrongDimensionException(EXCEPT_MSG(""));

    fvecP=fvec.p();
    crit=0.0;
    for (j=0; j<d; j++) {
      if (j==cls) continue; // skip
      if (smooMode==0) {
	mudf=smooMu+mu[cls]+bias[cls]-mu[j]-bias[j];
	temp=cfact[cls]; fact=smooSigsq+temp*temp;
      } else {
	temp=smooB1+1.0;
	mudf=temp*mu[cls]+bias[cls]+smooB2-mu[j]-bias[j];
	temp*=cfact[cls];
	fact=smooSigsq+temp*temp;
      }
      temp=cfact[j];
      fact=sqrt(fact+temp*temp);
      for (i=0,sum=0.0; i<quadN; i++)
	sum+=equadW[i]*TestUtils::logCdfNormal(fact*quadS[i]+mudf);
      crit+=sum;
      *(fvecP++)=-sum;
    }
    crit-=smooKL;
    if (smooMode==1) {
      temp=smooB1*mu[cls]+smooB2;
      fact=smooB1*cfact[cls];
      crit-=0.5*(fact*fact+temp*temp);
    }
    FastUtils::addscal(fvec.p(),crit,d-1);
    fvecFlag=2;

    return crit;
  }

  inline void
  MultiProbitQuadrature::compLogMomentsDiag0(const StVector& mu,
					     const StVector& cfact,
					     StVector& cvec,StVector& gvec)
  {
    int i,i2,j;
    double temp,temp2,sum1,sum2,off1,off2,mudf,sig,sigy,cacc1,gacc1,cacc2,
      gacc2,fact,add;

    if (mu.size()!=d || cfact.size()!=d)
      throw WrongDimensionException(EXCEPT_MSG("mu,cfact"));
    if (fvecFlag!=2)
      throw WrongStatusException(EXCEPT_MSG("'compExpectLogDiag' has to be called"));
    cvec.zeros(d); gvec.zeros(d);

    // Main loop (over comp. != 'cls')
    for (j=0; j<cls; j++) gvec[j]=fvec[j];
    for (j=cls+1; j<d; j++) gvec[j]=fvec[j-1];
    if (smooMode==0) {
      fact=1.0; add=smooMu;
    } else {
      fact=1.0+smooB1; add=smooB2;
    }
    temp=cfact[cls]*fact;
    sigy=sqrt(smooSigsq+temp*temp);
    cacc2=gacc2=0.0; // accus for 'cls' comp.
    for (j=0; j<d; j++) {
      if (j==cls) continue; // skip
      mudf=fact*mu[cls]+add+bias[cls]-mu[j]-bias[j];
      temp=cfact[j]; sig=sqrt(smooSigsq+temp*temp);
      // 'cacc1','gacc1' accus for 'j' comp.
      for (i=0,cacc1=gacc1=0.0; i<quadN; i++) {
	off1=mudf-cfact[j]*quadS[i];
	off2=mudf+cfact[cls]*fact*quadS[i];
	for (i2=0,sum1=sum2=0.0; i2<quadN; i2++) {
	  temp=quadS[i2]; temp2=equadW[i2];
	  sum1+=temp2*TestUtils::logCdfNormal(sigy*temp+off1);
 	  sum2+=temp2*TestUtils::logCdfNormal(sig*temp+off2);
	}
	temp=quadS[i]*equadW[i];
	cacc1+=(sum1*temp);
	cacc2+=(sum2*temp);
	temp*=quadS[i];
	gacc1+=(sum1*temp);
	gacc2+=(sum2*temp);
      }
      cvec[j]+=cacc1; gvec[j]+=gacc1;
    }
    cvec[cls]+=cacc2; gvec[cls]+=gacc2;
  }

  inline double
  MultiProbitQuadrature::compExpectLogEVecsDiag0(const StVector& mu,
						 const StVector& cfact,
						 StVector& e1vec,
						 StVector& e2vec)
  {
    int i,j;
    double temp,temp2,mudf,fact,crit,eps1,eps2,eps1acc,eps2acc;

    // ATTENTION: Mode 1 not implemented!
    if (smooMode==1)
      throw NotImplemException(EXCEPT_MSG(""));
    if (mu.size()!=d || cfact.size()!=d)
      throw WrongDimensionException(EXCEPT_MSG(""));
    e1vec.zeros(d); e2vec.zeros(d);

    crit=eps1acc=eps2acc=0.0;
    for (j=0; j<d; j++) {
      if (j==cls) continue; // skip
      mudf=smooMu+mu[cls]+bias[cls]-mu[j]-bias[j];
      temp=cfact[cls]; fact=smooSigsq+temp*temp;
      temp=cfact[j];
      fact=sqrt(fact+temp*temp);
      for (i=0,eps1=eps2=0.0; i<quadN; i++) {
	temp=fact*quadS[i]+mudf;
	temp2=equadW[i];
	crit+=temp2*TestUtils::logCdfNormal(temp);
	temp2*=TestUtils::derivLogCdfNormal(temp);
	eps1+=temp2;
	eps2+=temp2*quadS[i];
      }
      eps1acc+=eps1;
      e1vec[j]=eps1;
      eps2*=(0.5/fact);
      eps2acc+=eps2;
      e2vec[j]=-eps2;
    }
    e1vec[cls]=-eps1acc;
    e2vec[cls]=-eps2acc;
    // DEBUG:
    //cout << "EVecsDiag: z=" << crit-smooKL << ",mu=" << smooMu << ",sigsq=" << smooSigsq << endl;

    return crit-smooKL;
  }

  /*
   * See documentation. We do not do precomp. here, 'fvecFlag' is not changed
   * and 'fvec' is not used.
   */
  inline double
  MultiProbitQuadrature::compExpectLogDiag1(const StVector& mu,
					    const StVector& cfact)
  {
    int i1,cind,k,l;
    double temp,ay,ac,muc,b1,b2,sigc,arg,del,s;
    double zbnd=0.0;

    if (smooMode!=0)
      throw WrongStatusException(EXCEPT_MSG("Need smoothing mode 0"));
    if (mu.size()!=d || cfact.size()!=d)
      throw WrongDimensionException(EXCEPT_MSG(""));
    // Select marginal "closest" to Q(u_y)
    cind=findClosest(mu,cfact);
    s=sqrt(1.0+smooSigsq);
    ay=cfact[cls]; ay*=ay;
    ac=cfact[cind]; ac*=ac;
    del=mu[cls]+bias[cls]-mu[cind]-bias[cind];
    // Tough part first
    mu0vec[0]=(smooMu+del)/s;
    mu0vec[1]=smooMu+mu[cls]+bias[cls];
    mu0vec[2]=smooMu+del;
    a0mat[0]=(ay+ac)/(1.0+smooSigsq);
    a0mat[1]=ay/s;
    a0mat[2]=(ay+ac)/s;
    a0mat[3]=a0mat[4]=smooSigsq+ay;
    a0mat[5]=smooSigsq+ay+ac;
    for (k=l=0; l<d; l++) {
      if (l==cls || l==cind) continue;
      temp=cfact[l];
      aypr[k]=temp*temp;
      muypr[k++]=mu[l]+bias[l];
    }
    zbnd+=helperLogMomentsExtSmall(a0mat,mu0vec,aypr,muypr);
    //cout << "zbnd=" << zbnd << endl; // DEBUG!

    // Expect. over log Z_c part and rel. ent. part (1-d quad.)
    muc=(smooMu+del)/s; // mu_1
    sigc=sqrt(ay+ac)/s; // l_{1,1}
    if (smooMu==0.0 && smooSigsq==1.0) {
      for (i1=0; i1<quadN; i1++)
	zbnd+=equadW[i1]*TestUtils::logCdfNormal(sigc*quadS[i1]+muc);
    } else {
      b1=smooSigsq*smooMu/s;
      b2=-0.5*smooSigsq*(smooSigsq-1.0)/(smooSigsq+1.0);
      zbnd-=smooKL;
      for (i1=0; i1<quadN; i1++) {
	arg=sigc*quadS[i1]+muc;
	zbnd+=equadW[i1]*(TestUtils::logCdfNormal(arg)-
			  (b1+b2*arg)*TestUtils::derivLogCdfNormal(arg));
      }
    }

    return zbnd;
  }

  inline void
  MultiProbitQuadrature::compLogMomentsDiag1(const StVector& mu,
					     const StVector& cfact,
					     StVector& cvec,StVector& gvec,
					     double* elog)
  {
    int i1,i2,cind,j,k,l;
    double temp,temp2,ay,ac,s,del,muacc,varacc,muc,off,fact1,fact2,
      n1,sum,b1,b2,arg,zbnd,sigc;

    if (smooMode!=0)
      throw WrongStatusException(EXCEPT_MSG("Need smoothing mode 0"));
    if (mu.size()!=d || cfact.size()!=d)
      throw WrongDimensionException(EXCEPT_MSG(""));
    cvec.zeros(d); gvec.zeros(d);
    ay=cfact[cls]; ay*=ay;
    s=sqrt(1.0+smooSigsq);

    // 1) Components other than 'cls'
    // Main loop over components other than 'cls'
    for (j=0; j<d; j++) {
      if (j==cls) continue; // skip
      del=mu[cls]+bias[cls]-mu[j]-bias[j];
      // Tough part first
      mu0vec[0]=(smooMu+del)/s;
      mu0vec[1]=smooMu+mu[cls]+bias[cls];
      mu0vec[2]=smooMu+del;
      mu0vec[3]=0.0;
      ac=cfact[j]; ac*=ac;
      a0mat[0]=(ay+ac)/(1.0+smooSigsq);
      a0mat[1]=ay/s;
      a0mat[2]=(ay+ac)/s;
      a0mat[3]=-cfact[j]/s;
      a0mat[4]=a0mat[5]=smooSigsq+ay;
      a0mat[6]=0.0;
      a0mat[7]=smooSigsq+ay+ac;
      a0mat[8]=-cfact[j];
      a0mat[9]=1.0;
      for (k=l=0; l<d; l++) {
	if (l==cls || l==j) continue;
	temp=cfact[l];
	aypr[k]=temp*temp;
	muypr[k++]=mu[l]+bias[l];
      }
      helperLogMomentsExt(a0mat,mu0vec,aypr,muypr,i1part,i2part);
      cvec[j]+=FastUtils::sum(i1part,d-2);
      gvec[j]+=FastUtils::sum(i2part,d-2);
      // Expect. over log Z_c part (2-d quad.)
      muc=(smooMu+del)/s;
      fact1=cfact[cls]/s; fact2=cfact[j]/s;
      muacc=varacc=0.0;
      for (i1=0; i1<quadN; i1++) {
	n1=quadS[i1];
	off=muc-fact2*n1;
	for (i2=0,sum=0.0; i2<quadN; i2++)
	  sum+=equadW[i2]*TestUtils::logCdfNormal(fact1*quadS[i2]+off);
	muacc+=(temp=equadW[i1]*n1*sum);
	varacc+=(n1*temp);
      }
      cvec[j]+=muacc; gvec[j]+=varacc;
      // Rel. ent. part (2-d quad.)
      if (smooMu!=0.0 || smooSigsq!=1.0) {
	b1=smooSigsq*smooMu/s;
	b2=-0.5*smooSigsq*(smooSigsq-1.0)/(smooSigsq+1.0);
	muacc=0.0;
	varacc=-smooKL;
	for (i1=0; i1<quadN; i1++) {
	  n1=quadS[i1];
	  off=muc-fact2*n1;
	  for (i2=0,sum=0.0; i2<quadN; i2++) {
	    arg=fact1*quadS[i2]+off;
	    sum+=equadW[i2]*(b1+b2*arg)*TestUtils::derivLogCdfNormal(arg);
	  }
	  muacc-=(temp=equadW[i1]*n1*sum);
	  varacc-=(n1*temp);
	}
	cvec[j]+=muacc; gvec[j]+=varacc;
      }
    }

    // 2) Component 'cls'
    // Select marginal "closest" to Q(u_y)
    cind=findClosest(mu,cfact);
    del=mu[cls]+bias[cls]-mu[cind]-bias[cind];
    zbnd=0.0; // for 'elog'
    // Tough part first
    mu0vec[0]=(smooMu+del)/s;
    mu0vec[1]=smooMu+mu[cls]+bias[cls];
    mu0vec[2]=smooMu+del;
    mu0vec[3]=0.0;
    ac=cfact[cind]; ac*=ac;
    a0mat[0]=(ay+ac)/(1.0+smooSigsq);
    a0mat[1]=ay/s;
    a0mat[2]=(ay+ac)/s;
    a0mat[3]=cfact[cls]/s;
    a0mat[4]=a0mat[5]=smooSigsq+ay;
    a0mat[6]=cfact[cls];
    a0mat[7]=smooSigsq+ay+ac;
    a0mat[8]=cfact[cls];
    a0mat[9]=1.0;
    for (k=l=0; l<d; l++) {
      if (l==cls || l==cind) continue;
      temp=cfact[l];
      aypr[k]=temp*temp;
      muypr[k++]=mu[l]+bias[l];
    }
    helperLogMomentsExt(a0mat,mu0vec,aypr,muypr,i1part,i2part,&zpart);
    cvec[cls]+=FastUtils::sum(i1part,d-2);
    gvec[cls]+=FastUtils::sum(i2part,d-2);
    zbnd+=FastUtils::sum(zpart,d-2);
    // Expect. over log Z_c part (2-d quad.). 1-d quad. for 'zbnd'
    muc=(smooMu+del)/s;
    fact1=cfact[cls]/s; fact2=cfact[cind]/s;
    sigc=sqrt(ay+ac)/s;
    muacc=varacc=0.0;
    for (i1=0; i1<quadN; i1++) {
      n1=quadS[i1];
      off=muc+fact1*n1;
      for (i2=0,sum=0.0; i2<quadN; i2++)
	sum+=equadW[i2]*TestUtils::logCdfNormal(off-fact2*quadS[i2]);
      muacc+=(temp=equadW[i1]*n1*sum);
      varacc+=(n1*temp);
      zbnd+=equadW[i1]*TestUtils::logCdfNormal(sigc*quadS[i1]+muc);
    }
    cvec[cls]+=muacc; gvec[cls]+=varacc;
    // Rel. ent. part (2-d quad.). 1-d quad. for 'zbnd'
    if (smooMu!=0.0 || smooSigsq!=1.0) {
      b1=smooSigsq*smooMu/s;
      b2=-0.5*smooSigsq*(smooSigsq-1.0)/(smooSigsq+1.0);
      muacc=0.0;
      varacc=-smooKL;
      for (i1=0; i1<quadN; i1++) {
	n1=quadS[i1];
	off=muc+fact1*n1;
	for (i2=0,sum=0.0; i2<quadN; i2++) {
	  arg=off-fact2*quadS[i2];
	  sum+=equadW[i2]*(b1+b2*arg)*TestUtils::derivLogCdfNormal(arg);
	}
	muacc-=(temp=equadW[i1]*n1*sum);
	varacc-=(n1*temp);
	arg=sigc*quadS[i1]+muc;
	zbnd-=equadW[i1]*(b1+b2*arg)*TestUtils::derivLogCdfNormal(arg);
      }
      cvec[cls]+=muacc; gvec[cls]+=varacc;
    }

    if (elog!=0) *elog=zbnd;
  }

  inline void
  MultiProbitQuadrature::compLogMomentsDiag2(const StVector& mu,
					     const StVector& cfact,
					     StVector& cvec,StVector& gvec,
					     double* elog)
  {
    int i1,i2,cind,j,k,l;
    double temp,temp2,ay,ac,s,del,muacc,varacc,muc,off,fact1,fact2,
      n1,sum,b1,b2,arg,zbnd,sigc;

    if (smooMode!=0)
      throw WrongStatusException(EXCEPT_MSG("Need smoothing mode 0"));
    if (mu.size()!=d || cfact.size()!=d)
      throw WrongDimensionException(EXCEPT_MSG(""));
    cvec.zeros(d); gvec.zeros(d);
    ay=cfact[cls]; ay*=ay;
    s=sqrt(1.0+smooSigsq);

    // 1) Component 'cls'
    // Select marginal "closest" to Q(u_y) -> c
    cind=findClosest(mu,cfact);
    del=mu[cls]+bias[cls]-mu[cind]-bias[cind];
    ac=cfact[cind]; ac*=ac;
    zbnd=0.0; // for 'elog'
    // Tough part first
    mu0vec[0]=(smooMu+del)/s;
    mu0vec[1]=smooMu+mu[cls]+bias[cls];
    mu0vec[2]=smooMu+del;
    mu0vec[3]=0.0;
    a0mat[0]=(ay+ac)/(1.0+smooSigsq);
    a0mat[1]=ay/s;
    a0mat[2]=(ay+ac)/s;
    a0mat[3]=cfact[cls]/s;
    a0mat[4]=a0mat[5]=smooSigsq+ay;
    a0mat[6]=cfact[cls];
    a0mat[7]=smooSigsq+ay+ac;
    a0mat[8]=cfact[cls];
    a0mat[9]=1.0;
    for (k=l=0; l<d; l++) {
      if (l==cls || l==cind) continue;
      temp=cfact[l];
      aypr[k]=temp*temp;
      muypr[k++]=mu[l]+bias[l];
    }
    // We require 'zpart' again below
    helperLogMomentsExt(a0mat,mu0vec,aypr,muypr,i1part,i2part,&zpart);
    cvec[cls]+=FastUtils::sum(i1part,d-2);
    gvec[cls]+=FastUtils::sum(i2part,d-2);
    zbnd+=FastUtils::sum(zpart,d-2);
    FastUtils::smul(zpart,-1.0,d-2); // completed below
    // Expect. over log Z_c part (2-d quad.). 1-d quad. for 'zbnd'
    muc=(smooMu+del)/s; // mu_1
    fact1=cfact[cls]/s; fact2=cfact[cind]/s;
    sigc=sqrt(ay+ac)/s; // l_{1,1}
    muacc=varacc=0.0;
    for (i1=0; i1<quadN; i1++) {
      n1=quadS[i1];
      off=muc+fact1*n1;
      for (i2=0,sum=0.0; i2<quadN; i2++)
	sum+=equadW[i2]*TestUtils::logCdfNormal(off-fact2*quadS[i2]);
      muacc+=(temp=equadW[i1]*n1*sum);
      varacc+=(n1*temp);
      zbnd+=equadW[i1]*TestUtils::logCdfNormal(sigc*quadS[i1]+muc);
    }
    cvec[cls]+=muacc; gvec[cls]+=varacc;
    // Rel. ent. part (2-d quad.). 1-d quad. for 'zbnd'
    if (smooMu!=0.0 || smooSigsq!=1.0) {
      b1=smooSigsq*smooMu/s;
      b2=-0.5*smooSigsq*(smooSigsq-1.0)/(smooSigsq+1.0);
      muacc=0.0;
      varacc=-smooKL;
      zbnd-=smooKL;
      for (i1=0; i1<quadN; i1++) {
	n1=quadS[i1];
	off=muc+fact1*n1;
	for (i2=0,sum=0.0; i2<quadN; i2++) {
	  arg=off-fact2*quadS[i2];
	  sum+=equadW[i2]*(b1+b2*arg)*TestUtils::derivLogCdfNormal(arg);
	}
	muacc-=(temp=equadW[i1]*n1*sum);
	varacc-=(n1*temp);
	arg=sigc*quadS[i1]+muc;
	zbnd-=equadW[i1]*(b1+b2*arg)*TestUtils::derivLogCdfNormal(arg);
      }
      cvec[cls]+=muacc; gvec[cls]+=varacc;
    }
    if (elog!=0) *elog=zbnd; // return z(y) bound
    FastUtils::addscal(zpart,zbnd,d-2); // 'zpart' complete (req. below)

    // 2) Component 'cind'
    // Tough part first
    mu0vec[0]=(smooMu+del)/s;
    mu0vec[1]=smooMu+mu[cls]+bias[cls];
    mu0vec[2]=smooMu+del;
    mu0vec[3]=0.0;
    a0mat[0]=(ay+ac)/(1.0+smooSigsq);
    a0mat[1]=ay/s;
    a0mat[2]=(ay+ac)/s;
    a0mat[3]=-cfact[cind]/s;
    a0mat[4]=a0mat[5]=smooSigsq+ay;
    a0mat[6]=0.0;
    a0mat[7]=smooSigsq+ay+ac;
    a0mat[8]=-cfact[cind];
    a0mat[9]=1.0;
    for (k=l=0; l<d; l++) {
      if (l==cls || l==cind) continue;
      temp=cfact[l];
      aypr[k]=temp*temp;
      muypr[k++]=mu[l]+bias[l];
    }
    helperLogMomentsExt(a0mat,mu0vec,aypr,muypr,i1part,i2part);
    cvec[cind]+=FastUtils::sum(i1part,d-2);
    gvec[cind]+=FastUtils::sum(i2part,d-2);
    // Expect. over log Z_c part (2-d quad.)
    muc=(smooMu+del)/s;
    fact1=cfact[cls]/s; fact2=cfact[cind]/s;
    muacc=varacc=0.0;
    for (i1=0; i1<quadN; i1++) {
      n1=quadS[i1];
      off=muc-fact2*n1;
      for (i2=0,sum=0.0; i2<quadN; i2++)
	sum+=equadW[i2]*TestUtils::logCdfNormal(fact1*quadS[i2]+off);
      muacc+=(temp=equadW[i1]*n1*sum);
      varacc+=(n1*temp);
    }
    cvec[cind]+=muacc; gvec[cind]+=varacc;
    // Rel. ent. part (2-d quad.)
    if (smooMu!=0.0 || smooSigsq!=1.0) {
      b1=smooSigsq*smooMu/s;
      b2=-0.5*smooSigsq*(smooSigsq-1.0)/(smooSigsq+1.0);
      muacc=0.0;
      varacc=-smooKL;
      for (i1=0; i1<quadN; i1++) {
	n1=quadS[i1];
	off=muc-fact2*n1;
	for (i2=0,sum=0.0; i2<quadN; i2++) {
	  arg=fact1*quadS[i2]+off;
	  sum+=equadW[i2]*(b1+b2*arg)*TestUtils::derivLogCdfNormal(arg);
	}
	muacc-=(temp=equadW[i1]*n1*sum);
	varacc-=(n1*temp);
      }
      cvec[cind]+=muacc; gvec[cind]+=varacc;
    }

    // 3) Components other than 'cls', 'cind'
    int truej=0;
    ArrayHandle<double> aj(1),muj(1);
    for (j=0; j<d; j++) {
      if (j==cls || j==cind) continue; // skip
      cvec[j]=0.0; gvec[j]=zpart[truej++]; // Init. with 'zpart' comp. above
      // Tough part first
      mu0vec[0]=(smooMu+del)/s;
      mu0vec[1]=smooMu+mu[cls]+bias[cls];
      mu0vec[2]=smooMu+del;
      mu0vec[3]=0.0;
      a0mat[0]=(ay+ac)/(1.0+smooSigsq);
      a0mat[1]=ay/s;
      a0mat[2]=(ay+ac)/s;
      a0mat[3]=0.0;
      a0mat[4]=a0mat[5]=smooSigsq+ay;
      a0mat[6]=-cfact[j];
      a0mat[7]=smooSigsq+ay+ac;
      a0mat[8]=0.0;
      a0mat[9]=1.0;
      temp=cfact[j];
      aj[0]=temp*temp;
      muj[0]=mu[j]+bias[j];
      // Only for single y' = j:
      helperLogMomentsExt(a0mat,mu0vec,aj,muj,i1part,i2part);
      cvec[j]+=i1part[0];
      gvec[j]+=i2part[0];
    }
  }

  /*
   * Complicated! See external notes.
   */
  inline double
  MultiProbitQuadrature::compExpectLogEVecsDiag2(const StVector& mu,
						 const StVector& cfact,
						 StVector& e1vec,
						 StVector& e2vec)
  {
    int i1,i2,ypr,cind,j;
    double mlphir11,n1,n2,mu1,mu2,mu3,lphir22,r32,lf0,arg,zbnd,temp,temp2,
      r11,r22,lhr11,lphir11;
    double s,del,ay,ac;
    double l11_1,l11_2;
    double l21_1,l21_2;
    double l31_1,l31_2;
    double l22_1,l22_2,l22_3;
    double l32_1,l32_2,l32_3;
    double l33_1,l33_2,l33_3;
    double r11_1,r11_2;
    double dr22,dr32;
    double accu1,accu2,accu3,accu4,accu5,accu6;
    double fact1,fact2,fact3,fact4,isqrt;
    double b1,b2;
    bool isKL=(smooMu!=0.0 || smooSigsq!=1.0);

    // ATTENTION: Mode 1 not implemented!
    if (smooMode==1)
      throw NotImplemException(EXCEPT_MSG(""));
    if (mu.size()!=d || cfact.size()!=d)
      throw WrongDimensionException(EXCEPT_MSG(""));
    e1vec.zeros(d); e2vec.zeros(d);
    ay=cfact[cls]; ay*=ay;
    s=sqrt(1.0+smooSigsq);

    // Select marginal "closest" to Q(u_y) -> c
    cind=findClosest(mu,cfact);
    del=mu[cls]+bias[cls]-mu[cind]-bias[cind];
    ac=cfact[cind]; ac*=ac;
    zbnd=-smooKL;
    // Precomp.
    for (j=ypr=0; j<d; j++) {
      if (j==cls || j==cind) continue;
      temp=cfact[j]; aypr[ypr]=temp*temp;
      muypr[ypr++]=mu[j]+bias[j];
    }
    a0mat[0]=(ay+ac)/(1.0+smooSigsq); a0mat[1]=ay/s; a0mat[2]=(ay+ac)/s;
    a0mat[3]=a0mat[4]=smooSigsq+ay; a0mat[5]=smooSigsq+ay+ac;
    b1=smooSigsq*smooMu/s; b2=0.5*smooSigsq*(smooSigsq-1.0)/(smooSigsq+1.0);

    // n1 loop
    register double l11=sqrt(a0mat[0]);
    mu3=smooMu+del; mu1=mu3/s;
    l11_1=l11_2=(0.5/(1.0+smooSigsq))/l11;
    accu1=accu2=accu4=accu5=0.0;
    for (i1=0; i1<quadN; i1++) {
      n1=quadS[i1];
      // Quad. weights for n1 merged into 'mlphir11'
      r11=l11*n1+mu1;
      mlphir11=-(lphir11=TestUtils::logCdfNormal(r11))+quadW[i1];
      lhr11=TestUtils::logPdfNormal(r11)-lphir11;
      r11_1=l11_1*n1;
      r11_2=l11_2*n1;
      // y' loop
      for (j=0,ypr=-1; j<d; j++) {
	if (j==cls || j==cind) continue; // skip
	ypr++;
	accu3=accu6=0.0;
	// Compute Cholesky factor
	register double l21,l22,l31,l32,l33;
	l21=a0mat[1]/l11; l31=a0mat[2]/l11;
	l22=(a0mat[3]+aypr[ypr])-l21*l21;
	if (l22<=0.0) throw NumericalException(EXCEPT_MSG(""));
	l22=sqrt(l22);
	l32=(a0mat[4]-l21*l31)/l22;
	l33=a0mat[5]-l31*l31-l32*l32;
	if (l33<=0.0) throw NumericalException(EXCEPT_MSG(""));
	l33=sqrt(l33);
	mu2=smooMu+mu[cls]+bias[cls]-muypr[ypr];
	// Compute derivative factors
	l21_1=(1.0/s-l21*l11_1)/l11;
	l31_1=(1.0/s-l31*l11_1)/l11;
	l22_1=(0.5-l21*l21_1)/l22;
	l32_1=(1.0-l32*l22_1-l21*l31_1-l31*l21_1)/l22;
	l33_1=(0.5-l31*l31_1-l32*l32_1)/l33;
	l21_2=(0.0-l21*l11_2)/l11;
	l31_2=(1.0/s-l31*l11_2)/l11;
	l22_2=(0.0-l21*l21_2)/l22;
	l32_2=(0.0-l32*l22_2-l21*l31_2-l31*l21_2)/l22;
	l33_2=(0.5-l31*l31_2-l32*l32_2)/l33;
	//l21_3=(0.0-l21*l11_3)/l11; // --> 0
	//l31_3=(0.0-l31*l11_3)/l11; // --> 0
	//l22_3=(0.5-l21*l21_3)/l22;
	l22_3=0.5/l22;
	//l32_3=(0.0-l32*l22_3-l21*l31_3-l31*l21_3)/l22;
	l32_3=(-l32*l22_3)/l22;
	//l33_3=(0.0-l31*l31_3-l32*l32_3)/l33;
	l33_3=(-l32*l32_3)/l33;
	// n_2 loop
	for (i2=0; i2<quadN; i2++) {
	  n2=quadS[i2];
	  // Quad. weights for n2 merged into 'lphir22'
	  r22=l22*n2+l21*n1+mu2;
	  lphir22=equadW[i2]*(temp=TestUtils::logCdfNormal(r22));
	  r32=l32*n2+l31*n1+mu3;
	  isqrt=1.0/sqrt(1.0+l33*l33);
	  arg=r32*isqrt;
	  lf0=TestUtils::logCdfNormal(arg);
	  zbnd+=lphir22*exp(lf0+mlphir11);
	  fact1=isqrt*lphir22*exp(mlphir11+TestUtils::logPdfNormal(arg));
	  fact2=exp(quadW[i2]+lf0+mlphir11+TestUtils::logPdfNormal(r22)-temp);
	  fact3=-lphir22*exp(lf0+mlphir11+lhr11);
	  fact4=-arg*l33*isqrt;
	  // xx_1 factors (mu1_1 = 1/s, mu2_1 = 1, mu3_1 = 1)
	  dr22=l22_1*n2+l21_1*n1;
	  dr32=l32_1*n2+l31_1*n1;
	  accu1+=(fact1*(dr32+fact4*l33_1)+fact2*dr22+fact3*r11_1);
	  accu4+=(fact1+fact2+fact3/s);
	  // xx_2 factors (mu1_2 = -1/s, mu2_2 = 0, mu3_2 = -1)
	  dr22=l22_2*n2+l21_2*n1;
	  dr32=l32_2*n2+l31_2*n1;
	  accu2+=(fact1*(dr32+fact4*l33_2)+fact2*dr22+fact3*r11_2);
	  accu5-=(fact1+fact3/s);
	  // xx_3 factors (mu1_3 = 0, mu2_3 = -1, mu3_3 = 0)
	  //dr22=l22_3*n2+l21_3*n1;
	  dr22=l22_3*n2;
	  //dr32=l32_3*n2+l31_3*n1;
	  dr32=l32_3*n2;
	  //accu3+=(fact1*(dr32+fact4*l33_3)+fact2*dr22+fact3*r11_3);
	  accu3+=(fact1*(dr32+fact4*l33_3)+fact2*dr22);
	  accu6-=fact2;
	}
	e1vec[j]-=accu6; e2vec[j]-=accu3;
      }
      // log Z_c, rel. ent. parts (variable w = r11)
      if (isKL) {
	temp=b1-r11*b2;
	temp2=exp(lhr11); // H(r11)
	zbnd+=equadW[i1]*(lphir11-temp*temp2);
	temp2=exp(quadW[i1]+lhr11)*(1.0+b2+(r11+temp2)*temp);
      } else {
	zbnd+=equadW[i1]*lphir11;
	temp2=exp(quadW[i1]+lhr11);
      }
      // mu1_1 = 1/s, mu1_2 = -1/s, l11_1 = l11_2
      temp=temp2*n1*l11_1;
      temp2/=s;
      accu1+=temp; accu2+=temp;
      accu4+=temp2; accu5-=temp2;
    }
    e1vec[cls]=-accu4; e1vec[cind]=-accu5;
    e2vec[cls]=-accu1; e2vec[cind]=-accu2;

    return zbnd;
  }

  // Internal methods

  inline int MultiProbitQuadrature::findClosest(const StVector& mu,
						const StVector& cfact) const
  {
    int j,cind;
    double temp,var,muc,temp2,crit;

    // Select marginal "closest" to Q(u_y) -> c
    temp=cfact[cls]; var=temp*temp+1.0;
    muc=mu[cls]+bias[cls];
    j=(cls>0)?0:1;
    temp2=cfact[j]; temp=var+temp2*temp2;
    temp2=mu[j]+bias[j]-muc;
    crit=log(temp)+temp2*temp2/temp;
    cind=j++;
    for (; j<d; j++) {
      if (j==cls) continue;
      temp2=cfact[j]; temp=var+temp2*temp2;
      temp2=mu[j]+bias[j]-muc;
      temp=log(temp)+temp2*temp2/temp;
      if (temp<crit) {
	cind=j; crit=temp;
      }
    }

    return cind;
  }

  /*
   * See documentation for details. A is A_0 + a_{y'} \delta_2 \delta_2^T.
   * \mu is \mu_0 - (h_{y'}+\beta_{y'}) \delta_2. We compute the Cholesky
   * factor. for A here and do 2-d quadrature based on the coefficients.
   */
  inline void
  MultiProbitQuadrature::helperLogMomentsExt(const ArrayHandle<double>& a0mat,
					     const ArrayHandle<double>& mu0vec,
					     const ArrayHandle<double>& aypr,
					     const ArrayHandle<double>& muypr,
					     ArrayHandle<double>& i1part,
					     ArrayHandle<double>& i2part,
					     ArrayHandle<double>* zbnd)
  {
    int i1,i2,ypr,nypr;
    double mlphir11,n1,n2,mu1,mu2,mu3,mu4,lphir22,r32,lf0,lf1,arg,r42,
      i1acc,i2acc,zacc;

    if (a0mat.size()<10 || mu0vec.size()<4)
      throw WrongDimensionException(EXCEPT_MSG(""));
    nypr=aypr.size();
    if (muypr.size()!=nypr || i1part.size()<nypr ||
	i2part.size()<nypr || (zbnd!=0 && zbnd->size()<nypr))
      throw WrongDimensionException(EXCEPT_MSG(""));
    // n_1 loop
    ArrayUtils<double>::fill(i1part,0.0,nypr);
    ArrayUtils<double>::fill(i2part,0.0,nypr);
    if (zbnd!=0)
      ArrayUtils<double>::fill(zbnd->p(),0.0,nypr);
    if (a0mat[0]<=0.0)
      throw NumericalException(EXCEPT_MSG(""));
    register double l11=sqrt(a0mat[0]);
    mu1=mu0vec[0]; mu3=mu0vec[2]; mu4=mu0vec[3];
    for (i1=0; i1<quadN; i1++) {
      n1=quadS[i1];
      // Quad. weights for n1 merged into 'mlphir11'
      mlphir11=-TestUtils::logCdfNormal(l11*n1+mu1)+quadW[i1];
      // y' loop
      for (ypr=0; ypr<nypr; ypr++) {
	// Compute Cholesky factor
	register double l21,l22,l31,l32,l33,l41,l42,l43,l44,temp;
	l21=a0mat[1]/l11; l31=a0mat[2]/l11; l41=a0mat[3]/l11;
	l22=(a0mat[4]+aypr[ypr])-l21*l21;
	if (l22<=0.0) throw NumericalException(EXCEPT_MSG(""));
	l22=sqrt(l22);
	l32=(a0mat[5]-l21*l31)/l22;
	l42=(a0mat[6]-l21*l41)/l22;
	l33=a0mat[7]-l31*l31-l32*l32;
	if (l33<=0.0) throw NumericalException(EXCEPT_MSG(""));
	l33=sqrt(l33);
	l43=(a0mat[8]-l31*l41-l32*l42)/l33;
	l44=a0mat[9]-l41*l41-l42*l42-l43*l43;
	if (l44<=0.0) throw NumericalException(EXCEPT_MSG(""));
	l44=sqrt(l44);
	mu2=mu0vec[1]-muypr[ypr];
	i1acc=i2acc=zacc=0.0;
	// n_2 loop
	for (i2=0; i2<quadN; i2++) {
	  n2=quadS[i2];
	  // Quad. weights for n2 merged into 'lphir22'
	  lphir22=equadW[i2]*TestUtils::logCdfNormal(l22*n2+l21*n1+mu2);
	  r32=l32*n2+l31*n1+mu3;
	  r42=l42*n2+l41*n1+mu4;
	  arg=r32/sqrt(1.0+l33*l33);
	  lf0=TestUtils::logCdfNormal(arg);
	  zacc+=lphir22*exp(lf0+mlphir11);
	  lf1=TestUtils::logPdfNormal(arg)+log(l33)-0.5*log1p(l33*l33);
	  i1acc+=lphir22*(r42*exp(lf0+mlphir11)+l43*exp(lf1+mlphir11));
	  temp=l33*l33;
	  i2acc+=lphir22*((l44*l44+r42*r42+
			   l43*l43*(1.0-(temp/(1.0+temp))*arg*
				    TestUtils::derivLogCdfNormal(arg)))
			  *exp(lf0+mlphir11)+2.0*l43*r42*exp(lf1+mlphir11));
	}
	i1part[ypr]+=i1acc;
	i2part[ypr]+=i2acc;
	if (zbnd!=0)
	  (*zbnd)[ypr]+=zacc;
      }
    }
  }

  inline double
  MultiProbitQuadrature::helperLogMomentsExtSmall(const ArrayHandle<double>& a0mat,const ArrayHandle<double>& mu0vec,const ArrayHandle<double>& aypr,const ArrayHandle<double>& muypr)
  {
    int i1,i2,ypr;
    double mlphir11,n1,n2,mu1,mu2,mu3,lphir22,r32,lf0,arg;

    if (a0mat.size()<6 || mu0vec.size()<3)
      throw WrongDimensionException(EXCEPT_MSG(""));
    if (aypr.size()!=d-2 || muypr.size()!=d-2)
      throw WrongDimensionException(EXCEPT_MSG(""));
    // n_1 loop
    double zbnd=0.0;
    if (a0mat[0]<=0.0)
      throw NumericalException(EXCEPT_MSG(""));
    register double l11=sqrt(a0mat[0]);
    mu1=mu0vec[0]; mu3=mu0vec[2];
    for (i1=0; i1<quadN; i1++) {
      n1=quadS[i1];
      // Quad. weights for n1 merged into 'mlphir11'
      mlphir11=-TestUtils::logCdfNormal(l11*n1+mu1)+quadW[i1];
      // y' loop
      for (ypr=0; ypr<d-2; ypr++) {
	// Compute Cholesky factor
	register double l21,l22,l31,l32,l33,temp;
	l21=a0mat[1]/l11; l31=a0mat[2]/l11;
	l22=(a0mat[3]+aypr[ypr])-l21*l21;
	if (l22<=0.0) throw NumericalException(EXCEPT_MSG(""));
	l22=sqrt(l22);
	l32=(a0mat[4]-l21*l31)/l22;
	l33=a0mat[5]-l31*l31-l32*l32;
	if (l33<=0.0) throw NumericalException(EXCEPT_MSG(""));
	l33=sqrt(l33);
	mu2=mu0vec[1]-muypr[ypr];
	// n_2 loop
	for (i2=0; i2<quadN; i2++) {
	  n2=quadS[i2];
	  // Quad. weights for n2 merged into 'lphir22'
	  lphir22=equadW[i2]*TestUtils::logCdfNormal(l22*n2+l21*n1+mu2);
	  r32=l32*n2+l31*n1+mu3;
	  arg=r32/sqrt(1.0+l33*l33);
	  lf0=TestUtils::logCdfNormal(arg);
	  zbnd+=lphir22*exp(lf0+mlphir11);
	}
      }
    }

    return zbnd;
  }
//ENDNS

#endif
