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
 * Desc.:  Header class EPSingleLaplaceManager
 * ------------------------------------------------------------------- */

/*
 * Workaround for dynamic_cast bug in MATLAB_MEX mode:
 * If MATLAB_MEX is defined, the inheritance from 'EPSingleSiteManager'
 * is not virtual. This avoids bugs due to the fact that we cannot use
 * dynamic_cast in MEX files.
 * NOTE: The diamond-shaped hierarchy does not work in MEX mode then!
 */

/*
 * TODO:
 * If needed: allow for diff. tau per site.
 */

#ifndef EPLIN_EPSINGLELAPLACEMANAGER
#define EPLIN_EPSINGLELAPLACEMANAGER

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eplin/default.h"
#include "src/eplin/EPSingleSiteManager.h"
#include "src/eplin/EPLaplaceSite.h"
#include "src/eplin/EPExponentialSite.h"
//#include "src/eplin/EPDebugSite.h" // DEBUG!

//BEGINNS(eplin)
  /**
   * Manages n identical Laplace single sites, for fixed noise variance
   * sigma^2. The sites share a single tau parameter (see 'EPLaplaceSite').
   * <p>
   * Can also manage n identical exponential joint sites (see
   * 'EPExponentialSite').
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class EPSingleLaplaceManager :
#ifndef MATLAB_MEX
  public virtual EPSingleSiteManager
#else
  public EPSingleSiteManager
#endif
  {
  protected:
    // Members

    int n;                           // Number of sites
    StVector siteB,sitePi;           // Site parameters
    Handle<EPLaplaceSite> siteL;     // Laplace site
    Handle<EPExponentialSite> siteE; // Exponential site
    //Handle<EPDebugSite> siteE;       // DEBUG!
    double tau;                      // Site hyperparameter (positive)

  public:
    // Public methods

    /**
     * Constructor.
     * Site pars. are init. to all 0, have to be set to appropriate values!
     *
     * @param np     Number of sites
     * @param taup   Init. value for tau
     * @param expon  Use exponential sites instead of Laplace? Def.: false
     */
    EPSingleLaplaceManager(int np,double taup,bool expon=false) :
      tau(taup),n(np) {
      if (np<1 || taup<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
      if (!expon)
	siteL.changeRep(new EPLaplaceSite(taup));
      else
	siteE.changeRep(new EPExponentialSite(taup));
      //siteE.changeRep(new EPDebugSite(taup)); // DEBUG!
      siteB.zeros(np); sitePi.zeros(np);
    }

    /**
     * Special constructor.
     * Site pars passed in 'msiteb', 'msitepi'. These vectors are masked
     * here, no copies are drawn!
     *
     * @param msiteb  b_i vector. Masked here
     * @param msitepi pi_i vector. Masked here
     * @param taup    Init. value for tau
     * @param expon   Use exponential sites instead of Laplace? Def.: false
     */
    EPSingleLaplaceManager(const StVector& msiteb,const StVector& msitepi,
			   double taup,bool expon=false) :
      tau(taup),n(msiteb.size()) {
      if (n<1 || msitepi.size()!=n || taup<=0.0)
	throw InvalidParameterException(EXCEPT_MSG(""));
      if (!expon)
	siteL.changeRep(new EPLaplaceSite(taup));
      else
	siteE.changeRep(new EPExponentialSite(taup));
      //siteE.changeRep(new EPDebugSite(taup)); // DEBUG!
      siteB.reassign(msiteb); sitePi.reassign(msitepi);
    }

    int size() const {
      return n;
    }

    const StVector& getSiteB() const {
      return siteB;
    }

    const StVector& getSitePi() const {
      return sitePi;
    }

    void setSiteB(const StVector& newb) {
      if (newb.size()!=n) throw WrongDimensionException(EXCEPT_MSG(""));
      siteB=newb;
    }

    void setSiteB(double newb) {
      siteB.fill(n,newb);
    }

    void setSitePi(const StVector& newpi) {
      if (newpi.size()!=n) throw WrongDimensionException(EXCEPT_MSG(""));
      sitePi=newpi;
    }

    void setSitePi(double newpi) {
      sitePi.fill(n,newpi);
    }

    void setSiteB(int i,double newb) {
      if (i<0 || i>=n) throw InvalidParameterException(EXCEPT_MSG(""));
      siteB[i]=newb;
    }

    void setSitePi(int i,double newpi) {
      if (i<0 || i>=n) throw InvalidParameterException(EXCEPT_MSG(""));
      sitePi[i]=newpi;
    }

    /**
     * @return Site hyperparameter tau
     */
    double getTau() const {
      return tau;
    }

    /**
     * @param taun New value for tau
     */
    void setTau(double taun) {
      if (taun<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
      tau=taun;
      if (siteE==0) siteL->setTau(taun);
      else siteE->setTau(taun);
    }

    bool suppFractional(int i) const {
      return true;
    }

    double compMoments(int i,double sigsq,double h,double rho,double* beta=0,
		       double* nu=0,double frac=1.0,StVector* addVec=0) const {
      if (i<0 || i>=n) throw OutOfRangeException(EXCEPT_MSG(""));

      return (siteE==0)?
	siteL->compCondMoments(sigsq,h,rho,beta,nu,frac,addVec):
	siteE->compCondMoments(sigsq,h,rho,beta,nu,frac,addVec);
    }
  };
//ENDNS

#endif
