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
 * Desc.:  Header abstract class SingleLikehoodFactor
 * ------------------------------------------------------------------- */

#ifndef SINGLELIKEHOODFACTOR_H
#define SINGLELIKEHOODFACTOR_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/quad/default.h"
#include "lhotse/matrix/StVector.h"

class DataIndex;

//BEGINNS(quad)
  /**
   * Abstract base class for single (natural) parameter likelihood
   * factors in the context of EP approximate inference (ADF projections).
   * The factors are t_i(u), i the data point index. An implementation has
   * to maintain the responses y_i. The factors may depend on hyperparameters
   * (jointly for all factors), in which case 'compExpDeriv' has to be
   * implemented (and 'getNumHyperpars' has to return the number of pars.).
   * The hyperpars. are uncontrained double values, with 'getHyperpars',
   * 'setHyperpars' access methods.
   * NOTE: The factors may depend on additional parameters, but these
   * cannot be optimized over generically.
   * <p>
   * The 'SingleQuadrature' hierarchy provides implementations for the
   * quadrature required here.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class SingleLikehoodFactor
  {
  public:
    // Public methods

    virtual ~SingleLikehoodFactor() {}

    /**
     * @return Number of sites
     */
    virtual int size() const = 0;

    /**
     * Define the partition function
     *   Z_i = E[t_i(u)], E[.] over N(h,a), h=='mean',a=='var', i=='ind'.
     * The method returns log Z_i.
     * <p>
     * ADF terms (tilted moments):
     * If 'alpha', 'nu' are given (either both or none), the method also
     * returns
     *   alpha = (d/dh) log Z_i
     *   nu = -(d^2/dh^2) log Z_i
     * The moments hhat, ahat of the tilted distribution prop. to
     * t_i(u) N(u | h,a) are then computed as
     *   hhat = h + a alpha
     *   ahat = a (1 - nu a)
     * <p>
     * Fractional sites:
     * If fractional sites are supported ('suppFracSites' returns true),
     * a posit. value != 1 can be passed in 'frac'. Within this method
     * (only!), Z_i is replaced by Z_i^f then, where f=='frac'.
     *
     * @param ind   Site index i
     * @param mean  S.a., h
     * @param var   S.a., a
     * @param alpha S.a. Optional
     * @param nu    S.a. Optional
     * @param frac  S.a. Def.: 1
     * @return      log Z_i
     */
    virtual double compLogPart(int ind,double mean,double var,
			       double* alpha=0,double* nu=0,
			       double frac=1.0) const = 0;

    /**
     * See 'compLogPart'.
     *
     * @return Are fractional sites supported (only in 'compLogPart')?
     */
    virtual bool suppFracSites() const {
      return false; // not supp. by def.
    }

    /**
     * Computes
     *   E[log t_i(u)], E[.] over N(h,a), h=='mean', a=='var', i=='ind'.
     * If 'lmean', 'lvar' are given (either both or none), we also compute
     * the standardized moments
     *   lmean = E[log t_i(u) v], lvar = E[log t_i(u) v^2],
     *   v = (u-h)/sqrt(a)
     * NOTE: An implementation is not required by some users of this class.
     * This is why the def. implementation is given here, and the method is
     * not abstract. The def. implementation throws an exception.
     *
     * @param ind   Site index i
     * @param mean  S.a.
     * @param var   S.a.
     * @param lmean S.a. Optional
     * @param lvar  "
     * @return      S.a.
     */
    virtual double compExpectLog(int ind,double mean,double var,
				 double* lmean=0,double* lvar=0) const {
      throw NotImplemException(EXCEPT_MSG("Not implemented"));
    }

    /**
     * @param data Dataset
     * @return     Is the target type of 'data' compatible with this
     *             likelihood factor?
     */
    virtual bool checkDataset(const DataIndex& data) const = 0;

    /**
     * If the likelihood factor has no hyperpars., this method returns 0
     * (def. implement.). Used to test whether there are hyperpars. to be
     * optimized.
     *
     * @return Number of hyperpars.
     */
    virtual int getNumHyperpars() const {
      return 0;
    }

    /**
     * @return Hyperparam. vector
     */
    virtual const StVector& getHyperpars() const {
      throw NotImplemException("Likelihood has no hyperparameters");
    }

    /**
     * @param vec New hyperparam. vector
     */
    virtual void setHyperpars(const StVector& vec) {
      throw NotImplemException("Likelihood has no hyperparameters");
    }

    /**
     * Computes
     *   Z_i^{-1} E[(d t_i(u))/(d p)], E[.] over N(h,a), h=='mean',
     * a=='var', i=='ind' for hyperparameter p, which has index 'comp'.
     * NOTE: An implementation is not required by some users of this class.
     * This is why the def. implementation is given here, and the method is
     * not abstract. The def. implementation throws an exception.
     * Method has to be implemented iff 'getNumHyperpars' returns >0.
     *
     * @param ind   Site index i
     * @param comp  Index of hyperpar.
     * @param mean  S.a.
     * @param var   S.a.
     * @return      S.a.
     */
    virtual double compExpDeriv(int ind,int comp,double mean,double var)
      const {
      throw NotImplemException(EXCEPT_MSG("Likelihood has no hyperparameters"));
    }

    /**
     * Only if there are hyperpars.
     *
     * @return Can hyperpars. be printed?
     */
    virtual bool canPrint() const {
      return false;
    }

    virtual void print(ostream& os) const {
      throw NotImplemException(EXCEPT_MSG("Not implemented for this likelihood function"));
    }

  protected:
    // Internal methods

    /**
     * Helper for 'compExpDeriv'. Computes Z_i^{-1} E[(d t_i(u))/(d p)],
     * E[.] over N(h,a), h=='mean', a=='var', i=='ind', for the case that
     * p is a location parameter: t(u,p) = t(u+p). In this case, integration
     * by parts gives
     *   Z^{-1} E[(d t(u))/(d p)] = a^-1 Z^{-1} E[t(u) (u-h)] = alpha
     *
     * @param ind   Site index i
     * @param mean  S.a.
     * @param var   S.a.
     * @return      S.a.
     */
    double compExpDerivLocation(int ind,double mean,double var) const {
      double alpha,nu;
      compLogPart(ind,mean,var,&alpha,&nu);
      return alpha;
    }

    /**
     * Helper for 'compExpDeriv'. Computes Z_i^{-1} E[(d t_i(u))/(d p)],
     * E[.] over N(h,a), h=='mean', a=='var', i=='ind', for the case that
     * exp(p) is a scale parameter: t(u,p) = t(exp(p) u). In this case,
     * integration by parts gives
     *   Z^{-1} E[(d t(u))/(d p)] = -Z^{-1} E[ t(u) (1 - u a^-1 (u-h)) ]
     *   = h*alpha + a (alpha^2 - nu)
     *
     * @param ind   Site index i
     * @param mean  S.a.
     * @param var   S.a.
     * @return      S.a.
     */
    double compExpDerivScale(int ind,double mean,double var) const {
      double alpha,nu;
      compLogPart(ind,mean,var,&alpha,&nu);
      return mean*alpha+var*(alpha*alpha-nu);
    }
  };
//ENDNS

#endif
