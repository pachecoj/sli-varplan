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
 * Project source file
 * Module: eplin
 * Desc.:  Header class EPGaussLaguerreSite
 * ------------------------------------------------------------------- */

#ifndef EPLIN_EPGAUSSLAGUERRESITE_H
#define EPLIN_EPGAUSSLAGUERRESITE_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eplin/default.h"
#include "src/eplin/EPSingleSite.h"
#include "lhotse/quad/GaussQuad.h"
#include "lhotse/specfun/Specfun.h"
#include "lhotse/rando/Generator.h"
#include "lhotse/rando/Random.h"
#include "lhotse/matrix/StVector.h"

//#define DO_DEBCOMP 1
#ifdef DO_DEBCOMP
#include "lhotse/DebugMacros.h"
#endif

//BEGINNS(eplin)
  /**
   * Abstract class. Implements 'compJointMoments' using Gauss-Laguerre
   * quadrature, calling 'compLogI0' as subroutine (not implemented
   * here).
   * <p>
   * The present implementation works for sites of the form
   *   C tt exp[ -tt f(a) ],
   * where tt = tau/sigma, C is a constant, and C, f(.) do not depend on
   * tau, sigma. See 'EPLaplaceSite', 'EPExponentialSite' for examples.
   * <p>
   * Both fractional and joint moments are supported. For the joint moments,
   * the outer expectations over sigma^2 are done using Gauss-Laguerre
   * quadrature (see 'GaussQuad').
   * <p>
   * Gauss-Laguerre needs abscissas and weights, which depend on 'alpha' passed
   * to 'compJointMoments'. We maintain a set of these for different 'alpha'
   * values. For a new 'alpha', we find the entry with the closest alpha value,
   * and use it if the abs. diff. is <= 'alphaDiff'. Otherwise, a new entry
   * is computed for 'alpha', it replaces a randomly drawn one.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class EPGaussLaguerreSite : public EPSingleSite
  {
  protected:
    // Members

    double tau;
    int numev;
    mutable MAP_TYPE(double,int) alphaVals;
    mutable ArrayHandle<double> storage;
    int maxEntries;
    double alphaDiff;
    mutable Handle<Generator> prnGen;
    mutable ArrayHandle<double> wkbuff;
    int verbose;
#ifdef DO_DEBCOMP
    DCOMP_DEFSTOREVARS; // DEBUG
#endif

  public:
    mutable bool otherRule;

    // Public methods

    /**
     * Constructor. If this one is used, joint moments are not supported.
     *
     * @param itau Init. value for tau
     */
    EPGaussLaguerreSite(double itau) : tau(itau),numev(0) {
      if (itau<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
    }

    /**
     * Constructor. Joint moments are supported. 'nev' is the number of
     * function eval. in the Gauss-Laguerre quadrature.
     * <p>
     * A set of at most 'lmax' entries for abscissas/weights, one for each
     * value of 'alpha', is maintained. Requires a buffer of size
     * 2*'nev'*'lmax'. Can be passed in 'buff', or is alloc. here.
     * 'delAlph' is the value for 'alphaDiff', see header comment.
     *
     * @param itau    Init. value for tau
     * @param nev     S.a.
     * @param lmax    S.a.
     * @param delAlph S.a.
     * @param gen     PRN generator
     * @param buff    S.a. Def.: 0
     */
    EPGaussLaguerreSite(double itau,int nev,int lmax,double delAlph,
		      const Handle<Generator>& gen,
		      const ArrayHandle<double>& buff=
		      ArrayHandleZero<double>::get()) :
      tau(itau),numev(nev),maxEntries(lmax),storage(buff),prnGen(gen),
      alphaDiff(delAlph),verbose(0),otherRule(false) {
      if (itau<=0.0 || nev<2) throw InvalidParameterException(EXCEPT_MSG(""));
      if (lmax<1 || delAlph<=0.0 || (!(buff==0) && buff.size()<2*nev*lmax))
	throw InvalidParameterException(EXCEPT_MSG(""));
      if (buff==0) storage.changeRep(2*nev*lmax);
      wkbuff.changeRep(5*nev);
    }

    bool suppJointMoments() const {
      return (numev>0);
    }

    /**
     * @return Current tau value
     */
    double getTau() const {
      return tau;
    }

    /**
     * @param tv New tau value
     */
    void setTau(double tv) {
      if (tv<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
      tau=tv;
    }

    bool suppFractional() const {
      return true;
    }

    /**
     * @return Verbosity level
     */
    int getVerbose() const {
      return verbose;
    }

    /**
     * Set verbosity level:
     * 0: No output
     * 1: Some output
     * 2: All output
     *
     * @param verb New verb. level
     */
    void setVerbose(int verb) {
      if (verb<0) throw InvalidParameterException(EXCEPT_MSG(""));
      verbose=verb;
    }

    /**
     * The computation of 'addVec' (if 'addVec'!=0) is done in 'compAddVec'.
     */
    double compJointMoments(double h,double rho,double alpha,double zeta,
			    double& hp,double& rhop,double& alphap,
			    double& zetap,double frac=1.0,StVector* addVec=0)
      const;

    double compMeanSigsq(double h,double rho,double alpha,double zeta,
			 double frac) const;

  protected:
    // Internal methods

    /**
     * Computes log I_0, essentially as in 'compCondMoments', but tau =1,
     * sigma^2 =1, frac = 1.
     * Also, if t_i(.) = C tt exp( -tt f(a) ), then I_0 is based on
     * tt exp( -tt f(a) ) only, the leading constant C is not used!
     *
     * @param h    Mean
     * @param rho  Variance
     * @param beta Beta ret. here
     * @param nu   Nu ret. here
     * @return    log I_0
     */
    virtual double compLogI0(double h,double rho,double* beta,double* nu)
      const = 0;

    /**
     * If t_i(.) = C tt exp( -tt f(a) ) for a leading constant C, this
     * method returns log(C).
     *
     * @return log(C)
     */
    virtual double logLeadConst() const = 0;

    /**
     * Does 'addVec' computation, based on all (intermediate) results comp.
     * in 'compJointMoments'.
     * Has to be implemented by subclasses.
     */
    virtual void compAddVec(StVector& addVec,double h,double rho,double alpha,
			    double zeta,double hp,double rhop,double alphap,
			    double zetap,double frac,double logz,double cval,
			    double dela,const double* logi0,const double* logx,
			    const double* lwgts,const double* absc,
			    const double* beta,const double* nu) const {
      throw NotImplemException(EXCEPT_MSG(""));
    }

    /**
     * Helper for 'compJointMoments', 'compMeanSigsq'.
     */
    void updateWgtsAbsc(double aval,double& al,int& lpos) const;
  };

  // Inline methods

  inline double
  EPGaussLaguerreSite::compJointMoments(double h,double rho,double alpha,
					double zeta,double& hp,double& rhop,
					double& alphap,double& zetap,
					double frac,StVector* addVec)
    const
  {
    int i,lpos;
    double m1,m2,m3,m4,m6,temp,al,temp2,temp3,temp4,temp5;
    double* logi0,*logx,*sbuff,*beta,*nu;
    const double* lwgts,*absc;
    double ftau,aval,cval,lzeta,dela,xval,logz;
    StVector vecmsk;
#ifdef DO_DEBCOMP
    StVector debMsk; // DEBUG
#endif

    if (rho<=0.0 || alpha<= 0.0 || zeta<=0.0 || frac<=0.0 || frac>1.0)
      throw InvalidParameterException(EXCEPT_MSG(""));
    if (numev==0)
      throw WrongStatusException(EXCEPT_MSG("Joint moments not supported"));
    ftau=frac*tau; aval=alpha+0.5*frac-1.0;
    if (otherRule) aval-=2.0; // DEBUG!
#ifdef DO_DEBCOMP
    DCOMP_OPENSTOREFILE("debug_file");
    DCOMP_STORESCAL(h); DCOMP_STORESCAL(rho); DCOMP_STORESCAL(alpha);
    DCOMP_STORESCAL(zeta); DCOMP_STORESCAL(frac); DCOMP_STORESCAL(tau);
#endif

    // Make sure there are "close-by" weights/abscissas
    updateWgtsAbsc(aval,al,lpos);
    // First loop, to compute normalization constant log Z_i
    i=2*lpos;
    lwgts=storage.p()+(numev*i); absc=lwgts+numev;
    //debLwgts=lwgts; debAbsc=absc; // DEBUG!!
    logi0=wkbuff.p(); logx=logi0+numev; sbuff=logx+numev;
    beta=sbuff+numev; nu=beta+numev;
    dela=aval-al; lzeta=log(zeta);
    if (otherRule) dela+=2.0; // DEBUG
    cval=frac*(log(tau)+logLeadConst())-
      0.5*frac*lzeta-Specfun::logGamma(alpha);
    temp2=ftau*ftau*rho; temp3=ftau*h/sqrt(zeta);
    for (i=0; i<numev; i++) {
      xval=absc[i];
      // sigma^2 = tau = 1
      logi0[i]=temp=compLogI0(temp3*sqrt(xval),temp2,&beta[i],&nu[i])+lwgts[i];
      sbuff[i]=temp+cval+dela*(logx[i]=log(xval));
    }
    vecmsk.reassign(wkbuff,Range(2*numev,3*numev-1));
    logz=vecmsk.stableLogSum();
#ifdef DO_DEBCOMP
    debMsk.reassign(lwgts,numev); DCOMP_STOREVEC(debMsk);
    debMsk.reassign(absc,numev); DCOMP_STOREVEC(debMsk);
    debMsk.reassign(logi0,numev); DCOMP_STOREVEC(debMsk);
    debMsk.reassign(logx,numev); DCOMP_STOREVEC(debMsk);
    debMsk.reassign(beta,numev); DCOMP_STOREVEC(debMsk);
    debMsk.reassign(nu,numev); DCOMP_STOREVEC(debMsk);
    debMsk.reassign(sbuff,numev); DCOMP_STOREVEC(debMsk);
    DCOMP_STORESCAL(aval); DCOMP_STORESCAL(al);
    DCOMP_STORESCAL(lzeta); DCOMP_STORESCAL(cval);
    DCOMP_STORESCAL(logz);
#endif
    if (verbose>1)
      cout << "compJointMoments: log M_0=" << logz << endl;
    cval-=logz;
    // M_1, M_2, also using stable log summation
    temp=lzeta-logz;
    for (i=0; i<numev; i++)
      sbuff[i]+=(temp-logx[i]);
    m1=exp(temp2=vecmsk.stableLogSum());
#ifdef DO_DEBCOMP
    debMsk.reassign(sbuff,numev); DCOMP_STOREVEC(debMsk);
    DCOMP_STORESCAL(temp2);
#endif
    if (verbose>1)
      cout << "compJointMoments: M_1=" << m1 << endl;
    for (i=0; i<numev; i++)
      sbuff[i]+=(lzeta-logx[i]);
    m2=exp(temp3=vecmsk.stableLogSum());
    #ifdef DO_DEBCOMP
    debMsk.reassign(sbuff,numev); DCOMP_STOREVEC(debMsk);
    DCOMP_STORESCAL(temp3);
    #endif
    if (verbose>1)
      cout << "compJointMoments: M_2=" << m2 << endl;
    temp=2.0*temp2-temp3;
    if (temp>=0.0)
      throw NumericalException(EXCEPT_MSG("Variance must be positive"));
    alphap=1.0+1.0/(1.0-(temp3=exp(temp)));
    zetap=exp(temp2-log1p(-temp3));
    #ifdef DO_DEBCOMP
    DCOMP_STORESCAL(alphap); DCOMP_STORESCAL(zetap);
    #endif
    if (verbose>1)
      cout << "compJointMoments: alphap=" << alphap << ",zetap=" << zetap
	   << endl;
    // Second loop for remaining moments
    m3=m4=m6=0.0;
    temp2=ftau*rho*sqrt(zeta); temp5=rho*zeta;
    for (i=0; i<numev; i++) {
      xval=absc[i];
      temp=exp(logi0[i]+cval+(dela+1.0)*logx[i]-lzeta);
      m6+=temp;
      temp3=temp2/sqrt(xval); // kappa = eta tau sigma rho
      m3+=(temp*(temp4=h+temp3*beta[i]));
      m4+=(temp*(temp4*temp4+temp5/xval-temp3*temp3*nu[i]));
      #ifdef DO_DEBCOMP
      DCOMP_STORESCAL(temp); DCOMP_STORESCAL(temp3);
      DCOMP_STORESCAL(temp4);
      #endif
    }
    hp=m3*zetap/alphap;
    if (verbose>1)
      cout << "compJointMoments: M_3=" << m3 << ",hp=" << hp << endl;
    rhop=m4+hp*(hp*m6-2.0*m3);
    #ifdef DO_DEBCOMP
    DCOMP_STORESCAL(hp); DCOMP_STORESCAL(rhop);
    #endif
    if (verbose>1)
      cout << "compJointMoments: M_4=" << m4 << ",rhop=" << rhop << endl;
    // Optional: 'addVec'
    if (addVec!=0) compAddVec(*addVec,h,rho,alpha,zeta,hp,rhop,alphap,zetap,
			      frac,logz,cval,dela,logi0,logx,lwgts,absc,
			      beta,nu);
    #ifdef DO_DEBCOMP
    if (addVec!=0)
      DCOMP_STORESCAL((*addVec)[0]);
    DCOMP_CLOSESTOREFILE;
    #endif

    return logz;
  }

  inline double
  EPGaussLaguerreSite::compMeanSigsq(double h,double rho,double alpha,
				     double zeta,double frac) const
  {
    int i,lpos;
    double m1,temp,al,temp2,temp3;
    double* logx,*sbuff;
    const double* lwgts,*absc;
    double ftau,aval,cval,lzeta,dela,xval,logz;
    StVector vecmsk;

    if (rho<=0.0 || alpha<= 0.0 || zeta<=0.0 || frac<=0.0 || frac>1.0)
      throw InvalidParameterException(EXCEPT_MSG(""));
    if (numev==0)
      throw WrongStatusException(EXCEPT_MSG("Joint moments not supported"));
    ftau=frac*tau; aval=alpha+0.5*frac-1.0;
    if (otherRule) aval-=2.0; // DEBUG!

    // Make sure there are "close-by" weights/abscissas
    updateWgtsAbsc(aval,al,lpos);
    // First loop, to compute normalization constant log Z_i
    i=2*lpos;
    lwgts=storage.p()+(numev*i); absc=lwgts+numev;
    logx=wkbuff.p(); sbuff=logx+numev;
    dela=aval-al; lzeta=log(zeta);
    if (otherRule) dela+=2.0; // DEBUG
    cval=frac*(log(tau)+logLeadConst())-
      0.5*frac*lzeta-Specfun::logGamma(alpha);
    temp2=ftau*ftau*rho; temp3=ftau*h/sqrt(zeta);
    for (i=0; i<numev; i++) {
      xval=absc[i];
      // sigma^2 = tau = 1
      temp=compLogI0(temp3*sqrt(xval),temp2,0,0)+lwgts[i];
      sbuff[i]=temp+cval+dela*(logx[i]=log(xval));
    }
    vecmsk.reassign(wkbuff,Range(numev,2*numev-1));
    logz=vecmsk.stableLogSum();
    if (verbose>1)
      cout << "compMeanSigsq: log M_0=" << logz << endl;
    cval-=logz;
    // M_1, also using stable log summation
    temp=lzeta-logz;
    for (i=0; i<numev; i++)
      sbuff[i]+=(temp-logx[i]);
    m1=exp(vecmsk.stableLogSum());
    if (verbose>1)
      cout << "compMeanSigsq: M_1=" << m1 << endl;

    return m1;
  }

  inline void
  EPGaussLaguerreSite::updateWgtsAbsc(double aval,double& al,int& lpos) const
  {
    int i,sz;
    double temp;
    bool compThem=false;

    // Make sure "close-by" weights and abscissas are given at pos. 'lpos' in
    // 'storage'. Namely, the corr. alpha value is 'al', which is at most
    // 'alphaDiff' from 'aval'
    //cout << "a=" << aval << endl; // DEBUG!
    if ((sz=alphaVals.size())<maxEntries) {
      // There is still space, but must not insert the same twice
      MAP_CONSTITER(double,int) it=alphaVals.find(aval);
      if (it==alphaVals.end()) {
	lpos=sz; compThem=true;
      } else {
	lpos=it->second; al=aval;
      }
    } else {
      // Find closest entry
      // First element with key >= 'aval'
      MAP_CONSTITER(double,int) it=alphaVals.lower_bound(aval);
      if (it==alphaVals.end()) {
	// All have key < 'aval': pick last one
	--it;
	al=it->first; lpos=it->second;
      } else {
	al=it->first; lpos=it->second;
	if (it!=alphaVals.begin() && al!=aval) {
	  // Look at prev. one as well, whose key < 'aval'
	  --it; temp=it->first;
	  if (fabs(temp-aval)<fabs(al-aval)) {
	    // Prev. one is closer
	    al=temp; lpos=it->second;
	  }
	}
      }
      //cout << "  Closest: " << al << endl; // DEBUG
      if (fabs(al-aval)>alphaDiff) {
	// Not close enough -> recompute
	// Replace a random entry
	compThem=true;
	i=Random::devUniInt(0,maxEntries-1,prnGen);
	MAP_ITER(double,int) it2=alphaVals.begin();
	for (; i>0; ++it2,i--); // Go to pos. i (this is not random access!)
	lpos=it2->second;
	alphaVals.erase(it2);
      }
    }
    // DEBUG!!
    //if (!compThem && fabs(al-aval)>1e-6) compThem=true;
    // END DEBUG
    if (compThem) {
      //cout << "  Compute new" << endl; // DEBUG
      // Compute and store log weights/abscissas for 'alpha'
      i=2*lpos;
      GaussQuad::gaussLaguerre(numev,aval,storage.p()+(numev*i),
			       storage.p()+(numev*(i+1)),wkbuff.p());
      alphaVals[aval]=lpos;
      al=aval;
    } //else
    //cout << "  Use al=" << al << endl;
  }
//ENDNS

#ifdef DO_DEBCOMP
#include "lhotse/DebugMacros_undef.h"
#undef DO_DEBCOMP
#endif

#endif
