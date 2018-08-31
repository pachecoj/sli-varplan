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
 * Desc.:  Header abstract class EPSingleSiteManager
 * ------------------------------------------------------------------- */

/*
 * Workaround for dynamic_cast bug in MATLAB_MEX mode:
 * If MATLAB_MEX is defined, the inheritance from 'EPSiteManager' is
 * not virtual. This avoids bugs due to the fact that we cannot use
 * dynamic_cast in MEX files.
 * NOTE: The diamond-shaped hierarchy does not work in MEX mode then!
 */

#ifndef EPLIN_EPSINGLESITEMANAGER_H
#define EPLIN_EPSINGLESITEMANAGER_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eplin/default.h"
#include "src/eplin/EPSiteManager.h"
#include "lhotse/matrix/predecl.h"

//BEGINNS(eplin)
  /**
   * Manages complete set of sites to be used for EP on the linear or
   * generalized linear model. Relays moment services.
   * <p>
   * For example, in the linear model there is one (prior) site for each
   * parameter. In GLM/GP applications, there is one (likelihood) site for
   * each datapoint. Sites are numbered 0,...,n-1, n ret. by 'size'.
   * <p>
   * This class is for conditional moments (over single variable a, for
   * fixed sigma^2). There may be site hyperpars., but this is specific to
   * subclasses.
   * <p>
   * Site parameters:
   * Site parameters are b, pi for each site. The site approximation is
   *   N^U(a | sigma^-2 b, sigma^-2 pi).
   * If there are constraints on these parameters, they are not checked
   * here.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class EPSingleSiteManager :
#ifndef MATLAB_MEX
  public virtual EPSiteManager
#else
  public EPSiteManager
#endif
  {
  public:
    // Public methods

    virtual ~EPSingleSiteManager() {}

    /**
     * @return Vector of b site pars.
     */
    virtual const StVector& getSiteB() const = 0;

    /**
     * @return Vector of pi site pars.
     */
    virtual const StVector& getSitePi() const = 0;

    /**
     * @param newb New vector of b site pars.
     */
    virtual void setSiteB(const StVector& newb) = 0;
    virtual void setSiteB(double newb) = 0;

    /**
     * @param newpi New vector of pi site pars.
     */
    virtual void setSitePi(const StVector& newpi) = 0;
    virtual void setSitePi(double newpi) = 0;

    /**
     * @param i    Site index
     * @param newb New value of b site par.
     */
    virtual void setSiteB(int i,double newb) = 0;

    /**
     * @param i     Site index
     * @param newpi New vector of pi site pars.
     */
    virtual void setSitePi(int i,double newpi) = 0;

    /**
     * @param i Site index
     * @return  Supports frac. moments?
     */
    virtual bool suppFractional(int i) const = 0;

    /**
     * Relays 'EPSingleSite::compCondMoments' to site with index 'i'.
     *
     * @param i      Site index
     * @param sigsq  Value sigma^2
     * @param h      Mean
     * @param rho    Variance / sigma^2
     * @param beta   S.a. Def.: 0
     * @param nu     S.a. Def.: 0
     * @param frac   S.a. Def.: 1
     * @param addVec S.a. Def.: 0
     * @return       log Z
     */
    virtual double compMoments(int i,double sigsq,double h,double rho,
			       double* beta=0,double* nu=0,double frac=1.0,
			       StVector* addVec=0) const = 0;
  };
//ENDNS

#endif
