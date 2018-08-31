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
 * Desc.:  Header class EPExponentialSite
 * ------------------------------------------------------------------- */

#ifndef EPLIN_EPEXPONENTIALSITE_H
#define EPLIN_EPEXPONENTIALSITE_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eplin/default.h"
#include "src/eplin/EPGaussLaguerreSite.h"

//BEGINNS(eplin)

#define TAIL_APPROX_THRES 5.0 // See 'compLogPart'

  /**
   * Represents exponential distribution site, in context of EP for linear or
   * generalized linear model:
   *   t(a | sigma^2) = tt exp[ -tt a ] I{a>=0}, tt=tau/sigma.
   * This has the form supported by 'EPGaussLaguerreSite', where the outer
   * integration over sigma^2 (for joint moments) is done using Gauss-
   * Laguerre quadrature.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class EPExponentialSite : public EPGaussLaguerreSite
  {
  public:
    // Public methods

    /**
     * Constructor. If this one is used, joint moments are not supported.
     *
     * @param itau  Init. value for tau
     */
    EPExponentialSite(double itau) : EPGaussLaguerreSite(itau) {}

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
    EPExponentialSite(double itau,int nev,int lmax,double delAlph,
		      const Handle<Generator>& gen,
		      const ArrayHandle<double>& buff=
		      ArrayHandleZero<double>::get()) :
      EPGaussLaguerreSite(itau,nev,lmax,delAlph,gen,buff) {}

    /**
     * If 'addVec'!=0, E[-a] is ret. in '*addVec[0]', where E[.] is over the
     * tilted distribution. In this case, 'beta' and 'nu' must be given.
     * NOTE: This is correct for 'frac'==1 only!
     */
    double compCondMoments(double sigsq,double h,double rho,double* beta=0,
			   double* nu=0,double frac=1.0,StVector* addVec=0)
      const {
      double tiltau=tau/sqrt(sigsq),logz;
      bool add=(addVec!=0);

      if (add && beta==0) throw InvalidParameterException(EXCEPT_MSG(""));
      logz=compLogPartStatic(h,sigsq*rho,tiltau,beta,nu,frac)+frac*log(tiltau);
      if (add)
	addVec->fill(1,-h-sigsq*rho*(*beta));

      return logz;
    }

  protected:

    double compLogI0(double h,double rho,double* beta,double* nu) const {
      return compLogPartStatic(h,rho,1.0,beta,nu);
    }

    double logLeadConst() const {
      return 0.0;
    }

    /**
     * If 'addVec'!=0, E[a/sigma] is returned in '*addVec[0]', where E[.]
     * is over the tilted distribution.
     */
    void compAddVec(StVector& addVec,double h,double rho,double alpha,
		    double zeta,double hp,double rhop,double alphap,
		    double zetap,double frac,double logz,double cval,
		    double dela,const double* logi0,const double* logx,
		    const double* lwgts,const double* absc,const double* beta,
		    const double* nu) const;

    static double compLogPartStatic(double mean,double var,double tau,
				    double* beta,double* nu,double frac=1.0);
  };

  // Inline methods

  /*
   * See 'LaplaceSingleLikehood::compLogPartStatic' for how this works
   */
  inline double
  EPExponentialSite::compLogPartStatic(double mean,double var,double tau,
				       double* beta,double* nu,double frac)
  {
    double temp,a,sqa,h,argm,logi0,tailm;

    if (frac<=0.0) throw InvalidParameterException(EXCEPT_MSG("frac"));
    if (tau<=0.0) throw InvalidParameterException(EXCEPT_MSG("tau"));
    tau*=frac;
    h=mean*tau; a=var*tau*tau; sqa=tau*sqrt(var); argm=sqa-h/sqa;
    // log(1-Phi(x_-)) = F(x_-)
    if (argm>=TAIL_APPROX_THRES) {
      // Use tail approximation
      temp=1.0/argm/argm;
      // g(x_-) = 1 + 'tailm'
      tailm=-temp*(1.0-3.0*temp*(1.0-5.0*temp*(1.0-7.0*temp)));
      logi0=-0.5*argm*argm-log(argm)+log1p(tailm)-0.5*(M_LN2+M_LNPI);
    } else
      logi0=log1p(-Specfun::cdfNormal(argm));
    logi0+=(0.5*a-h); // log I_0 complete
    if (beta!=0) {
      if (nu==0) throw InvalidParameterException(EXCEPT_MSG(""));
      *beta=temp=tau*(exp(-0.5*mean*mean/var-logi0)/M_SQRT2/M_SQRTPI/sqa-1.0);
      *nu=(temp+tau)*(temp+mean/var);
    }

    return logi0;
  }

  inline void
  EPExponentialSite::compAddVec(StVector& addVec,double h,double rho,
				double alpha,double zeta,double hp,
				double rhop,double alphap,double zetap,
				double frac,double logz,double cval,
				double dela,const double* logi0,
				const double* logx,const double* lwgts,
				const double* absc,const double* beta,
				const double* nu) const
  {
    int i;
    double m5,temp2,temp3,xval;

    // Same as M_3 in 'EPGaussLaguerreSite', but a/sigma here instead of
    // a/sigma^2 there
    m5=0.0;
    temp2=frac*tau*rho*sqrt(zeta);
    temp3=cval-0.5*log(zeta); dela+=0.5;
    for (i=0; i<numev; i++) {
      xval=absc[i];
      m5+=exp(logi0[i]+temp3+dela*logx[i])*(h+temp2*beta[i]/sqrt(xval));
    }
    addVec.zeros(1); addVec[0]=m5;
    if (verbose>1)
      cout << "compJointMoments: M_5=" << m5 << endl;
  }
//ENDNS

#endif
