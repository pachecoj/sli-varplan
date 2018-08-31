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
 * Desc.:  Header class EPLaplaceSite
 * ------------------------------------------------------------------- */

#ifndef EPLIN_EPLAPLACESITE_H
#define EPLIN_EPLAPLACESITE_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eplin/default.h"
#include "src/eplin/EPGaussLaguerreSite.h"
#include "lhotse/quad/LaplaceSingleLikehood.h"

//BEGINNS(eplin)
  /**
   * Represents Laplace distribution site, in context of EP for linear or
   * generalized linear model:
   *   t(a | sigma^2) = tt/2 exp[ -tt |a| ], tt = tau/sigma > 0.
   * Typically, another class manages this object and maintains these
   * parameters.
   * <p>
   * Both fractional and joint moments are supported. The cond. moments
   * are computed by calling 'LaplaceSingleLikehood::compLogPartStatic' in
   * module 'quad'. For the joint moments, the outer expectations over
   * sigma^2 are done using Gauss-Laguerre quadrature (see
   * 'EPGaussLaguerreSite').
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class EPLaplaceSite : public EPGaussLaguerreSite
  {
  public:
    // Public methods

    /**
     * Constructor. If this one is used, joint moments are not supported.
     *
     * @param itau Init. value for tau
     */
    EPLaplaceSite(double itau) : EPGaussLaguerreSite(itau) {}

    /**
     * Constructor. Joint moments are supported. See 'EPGaussLaguerreSite'.
     *
     * @param itau    Init. value for tau
     * @param nev     S.a.
     * @param lmax    S.a.
     * @param delAlph S.a.
     * @param gen     PRN generator
     * @param buff    S.a. Def.: 0
     */
    EPLaplaceSite(double itau,int nev,int lmax,double delAlph,
		  const Handle<Generator>& gen,
		  const ArrayHandle<double>& buff=
		  ArrayHandleZero<double>::get()) :
      EPGaussLaguerreSite(itau,nev,lmax,delAlph,gen,buff) {}

    /**
     * If 'addVec' is given (!=0), E[-|a|] is ret. in the first component,
     * where E[.] over the tilted distribution.
     * NOTE: This is correct for 'frac'==1 only!
     */
    double compCondMoments(double sigsq,double h,double rho,double* beta=0,
			   double* nu=0,double frac=1.0,StVector* addVec=0)
      const {
      bool add=(addVec!=0);
      double tiltau=tau/sqrt(sigsq),logz,temp;

      if (add) {
	if (beta==0) throw InvalidParameterException(EXCEPT_MSG(""));
	if (frac!=1.0)
	  throw NotImplemException(EXCEPT_MSG("Not implemented for fractional"));
	addVec->zeros(1);
      }
      temp=frac*tau;
      logz=LaplaceSingleLikehood::compLogPartStatic(h*tiltau*frac,
						    rho*temp*temp,1.0,beta,nu,
						    1.0,add?&(*addVec)[0]:0);
      if (beta!=0) {
	temp=frac*tiltau;
	(*beta)*=temp; (*nu)*=temp; (*nu)*=temp;
	if (add) (*addVec)[0]/=tiltau;
      }

      return logz+frac*log(0.5*tiltau);
    }

  protected:

    double compLogI0(double h,double rho,double* beta,double* nu) const {
      return LaplaceSingleLikehood::compLogPartStatic(h,rho,1.0,beta,nu);
    }

    double logLeadConst() const {
      return log(0.5);
    }

    /**
     * The 'addVec' computation in 'compCondMoments':
     * E[|a|/sigma] is ret. in the first component, where E[.] over the
     * tilted distribution.
     * NOTE: This is correct for 'frac'==1 only!
     */
    void compAddVec(StVector& addVec,double h,double rho,double alpha,
		    double zeta,double hp,double rhop,double alphap,
		    double zetap,double frac,double logz,double cval,
		    double dela,const double* logi0,const double* logx,
		    const double* lwgts,const double* absc,const double* beta,
		    const double* nu) const;
  };

  // Inline methods

  inline void
  EPLaplaceSite::compAddVec(StVector& addVec,double h,double rho,double alpha,
			    double zeta,double hp,double rhop,double alphap,
			    double zetap,double frac,double logz,double cval,
			    double dela,const double* logi0,const double* logx,
			    const double* lwgts,const double* absc,
			    const double* beta,const double* nu) const
  {
    int i;
    double m5,temp,temp6,temp7,temp8,temp9,xval,ftau;

    m5=0.0;
    ftau=frac*tau;
    temp6=h*h/rho/zeta; temp7=sqrt(rho*M_2_PI);
    temp8=h/sqrt(zeta); temp9=ftau*rho;
    for (i=0; i<numev; i++) {
      xval=absc[i];
      temp=dela*logx[i]+cval;
      m5+=(exp(lwgts[i]+temp-0.5*temp6*xval)*temp7-
	   exp(logi0[i]+temp)*(temp8*sqrt(xval)*beta[i]+temp9));
    }
    addVec.zeros(1); addVec[0]=m5;
    if (verbose>1)
      cout << "compJointMoments: M_5=" << m5 << endl;
  }
//ENDNS

#endif
