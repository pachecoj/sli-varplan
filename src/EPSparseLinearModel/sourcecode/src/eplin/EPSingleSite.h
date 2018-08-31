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
 * Desc.:  Header abstract class EPSingleSite
 * ------------------------------------------------------------------- */

#ifndef EPLIN_EPSINGLESITE_H
#define EPLIN_EPSINGLESITE_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eplin/default.h"

//BEGINNS(eplin)
  /**
   * Represents site t(a | sigma^2) to be used in the context of EP on the
   * linear or a generalized linear model.
   * <p>
   * Here, a is an unconstrained scalar, and sigma^2 is a positive noise
   * variance parameter. For example, this may be a prior site in a linear
   * model with weights a_i and noise variance sigma^2.
   * In most applications, t(a | sigma^2) is a log-concave function of a,
   * for all fixed sigma^2.
   * <p>
   * Conditional moments:
   * The main service method is 'compCondMoments', which returns moments
   * of the tilted distribution (over a) prop. to
   *   t(a | sigma^2) N(a | h,sigma^2 rho)
   * This computation may require numerical quadrature.
   * <p>
   * Fractional moments:
   * Class supports this iff 'suppFractional' returns true. In this case,
   * a fraction f in (0,1] can be passed to the moment comp. methods, and
   * they then use t(a | sigma^2)^f.
   * <p>
   * Joint moments:
   * Class supports this iff 'suppJointMoments' returns true. In this case,
   * moments of a bivariate tilted distribution are computed upon call of
   * 'compJointMoments'. The distrib. is over (a,sigma^2), prop. to
   *   t(a | sigma^2) N(a | h,sigma^2 rho) IG(sigma^2 | alpha,zeta),
   * where IG denotes inverse-Gamma. This computation may require numerical
   * quadrature.
   * <p>
   * Dependence on additional parameters:
   * This is specific to subclass implementations. Typically, there is a
   * class managing all sites together, and parameters are obtained from
   * that one.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class EPSingleSite
  {
  public:
    // Public methods

    virtual ~EPSingleSite() {}

    /**
     * See header comment.
     * NOTE: If fractional and joint moments are supported, then fractional
     * joint moments have to be supp. as well.
     *
     * @return Do we support fractional moments?
     */
    virtual bool suppFractional() const {
      return false; // not supp. by def.
    }

    /**
     * See header comment.
     *
     * @return Do we support joint moments?
     */
    virtual bool suppJointMoments() const {
      return false; // not supp. by def.
    }

    /**
     * The tilted distribution is prop. to
     *   t(a | sigma^2) N(a | h,sigma^2 rho),
     * where sigma^2=='sigsq', and has normalization constant
     *   Z = E[t(a | sigma^2)], E[.] over N(a | h,sigma^2 rho)
     * The method returns log Z.
     * If 'beta', 'nu' are used (both or none), mean hhat and variance
     * sigma^2 rhohat of the tilted distribution are also computed.
     * Namely,
     *   hhat   = h + sigma^2 rho beta
     *   rhohat = rho (1 - sigma^2 rho nu)
     * Note that
     *   beta   = (d/dh) log Z
     *   nu     = -(d^2/dh^2) log Z
     * For a log-concave site, nu must be nonnegative.
     * If fractional moments are supported (see header), f=='frac' can
     * be used, in that f in (0,1). t(a | sigma^2)^f is used then.
     * <p>
     * Further results can be returned in comp. of the vector 'addRes',
     * if given. Format depends on subclass.
     *
     * @param sigsq  Value sigma^2
     * @param h      Mean
     * @param rho    Variance / sigma^2
     * @param beta   S.a. Def.: 0
     * @param nu     S.a. Def.: 0
     * @param frac   S.a. Def.: 1
     * @param addVec S.a. Def.: 0
     * @return       log Z
     */
    virtual double compCondMoments(double sigsq,double h,double rho,
				   double* beta=0,double* nu=0,
				   double frac=1.0,StVector* addVec=0)
      const = 0;

    /**
     * For joint moments over (a,sigma^2), the tilted distribution is prop.
     * to
     *   t(a | sigma^2) N(a | h,sigma^2 rho) IG(sigma^2 | alpha,zeta),
     * where sigma^2=='sigsq', and the normal. constant is Z (so that
     * Z^-1 ... is normalized). The inverse-Gamma distribution is
     *   IG(x | alpha,zeta) propto x^(-alpha-1) exp(-zeta/x), x>0.
     * The moments we compute here are the parameters of
     *   N(a | hp,sigma^2 rhop) IG(sigma^2 | alphap,zetap),
     * which minimizes D[ Tilted | . ], and log Z. For example, alphap,
     * zetap determine mean, variance of the tilted distribution over
     * sigma^2.
     * <p>
     * If fractional moments are supported (see header), f=='frac' can
     * be used, in that f in (0,1). t(a | sigma^2)^f is used then.
     * <p>
     * Further results can be returned in comp. of the vector 'addRes',
     * if given. Format depends on subclass.
     *
     * @param h      S.a.
     * @param rho    S.a.
     * @param alpha  S.a.
     * @param zeta   S.a.
     * @param hp     S.a. Res. ret. here
     * @param rhop   S.a. Res. ret. here
     * @param alphap S.a. Res. ret. here
     * @param zetap  S.a. Res. ret. here
     * @param frac   S.a. Def.: 1
     * @param addRes S.a. Def.: 0
     * @return       log Z (s.a.)
     */
    virtual double compJointMoments(double h,double rho,double alpha,
				    double zeta,double& hp,double& rhop,
				    double& alphap,double& zetap,
				    double frac=1.0,StVector* addRes=0)
      const {
      throw NotImplemException(EXCEPT_MSG(""));
    }

    /**
     * Does part of 'compJointMoments', in that only the mean of sigma^2 is
     * computed under the tilted distribution.
     *
     * @param h      See 'compCondMoments'
     * @param rho    "
     * @param alpha  "
     * @param zeta   "
     * @param frac   "
     * @return       E[sigma^2]
     */
    virtual double compMeanSigsq(double h,double rho,double alpha,
				 double zeta,double frac) const {
      throw NotImplemException(EXCEPT_MSG(""));
    }
  };
//ENDNS

#endif
