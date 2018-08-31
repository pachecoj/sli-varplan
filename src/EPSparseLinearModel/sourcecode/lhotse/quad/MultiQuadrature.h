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
 * Desc.:  Header abstract class MultiQuadrature
 * ------------------------------------------------------------------- */

#ifndef MULTIQUADRATURE_H
#define MULTIQUADRATURE_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/quad/default.h"
#include "lhotse/matrix/predecl.h"

//BEGINNS(quad)
  /**
   * Some rules are configured by additional parameters which can be
   * collected as members in subclasses of this one.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class MultiQuadParams
  {
  public:
    virtual ~MultiQuadParams() {}

    /**
     * Resets members to def. values (same as for default constructor).
     */
    virtual void reset() = 0;
  };

  /**
   * Interface for multidimensional quadrature rules to approximate
   * statistics for GP Bayesian inference which can be expressed as
   * multivariate Gaussian expectations of a smooth function (typically
   * the likelihood, log likelihood).
   * <p>
   * Concrete subclasses inherit in diamond shape from 2 subclasses of
   * this one, one implementing the concrete rule, the other implementing
   * the integrand function f.
   * <p>
   * Vector-valued f:
   * 'compLogPartMulti' computes log partition functions for the tilted
   * Gaussians prop. to
   *   (exp f_i(u)) N(u | mu,sigma),
   * where now f=(f_i) is vector-valued (needs to be implemented in
   * 'compFMulti'). This makes sense if the comp. of f require shared
   * computations, so 'compFMulti' is sign. faster than calling 'compF'
   * separ. for each component. 'getFNum' has to return the number of comp.
   * of f.
   * NOTE: If this is not the case, implement 'compFMulti' by calling
   * 'compF' for each dimension.
   * <p>
   * Naming conventions:
   * Roughly this is
   *   'comp'<func><multi><diag>
   * where:
   * - <multi>=='Multi': Vector-valued f, rather than real-val. f
   * - <diag>=='Diag': Cov. matrix sigma is diagonal
   * If two kinds of results are computed together, the <func><multi><diag>
   * part is repeated.
   * <p>
   * Partial implementations:
   * In some cases, methods are not implemented in subclasses, for example
   * the <diag>=='' methods may not be provided if the subclass is specific
   * to diagonal cov. matrices. The subclass must provide an implementation
   * which throws a 'NotImplemException'.
   * NOTE: If methods here are only implem. in special subclasses, the def.
   * implementation here is such a dummy.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class MultiQuadrature
  {
  public:
    // Public methods

    MultiQuadrature() {}

    virtual ~MultiQuadrature() {}

    /**
     * @return Dimensionality
     */
    virtual int getDim() const = 0;

    /**
     * Returns true if the quadrature rule has nonnegative weights only,
     * which implies that it realizes an expectation under a discrete
     * probability measure, so moment characteristics (constraints) are
     * preserved.
     *
     * @return S.a.
     */
    virtual bool hasNonnegWeights() const = 0;

    /**
     * Some rules come with additional parameters which influence some of
     * the 'compXXX' methods. The parameters are members of a corr. subclass
     * of 'MultiQuadParams', this method is used to set them.
     * NOTE: The def. implementation does nothing.
     *
     * @param pars New parameter values
     */
    virtual void setParams(const MultiQuadParams* pars) {}

    /**
     * @return Does the rule come with add. params. which can be set via
     *         'setParams'?
     */
    virtual bool hasParams() const {
      return false;
    }

    /**
     * Creates and returns object of appropriate 'MultiQuadParams' subclass
     * (the type accepted by 'setParams').
     * NOTE: Wrap into a handle.
     *
     * @return New 'MultiQuadParams' object
     */
    virtual MultiQuadParams* getParamObj() const {
      throw WrongStatusException(EXCEPT_MSG("No additional parameters"));
    }

    /**
     * Computes log partition function
     *   log Z = log E[exp(f(u))],
     * E[] over N(mu,sigma), sigma = L*L^T (Cholesky).
     * <p>
     * NOTE: For many rules, it is useful to store some information here
     * (e.g. a vector) which is used by 'compMoments'. If so, store it in
     * a member.
     *
     * @param mu    Mean mu
     * @param cfact Cholesky factor L (lower tr.)
     * @return      log Z
     */
    virtual double compLogPart(const StVector& mu,const StMatrix& cfact) = 0;

    /**
     * Same as 'compLogPart', but Cholesky factor L is diagonal. For
     * certain functions f, this can save up to a factor of C (dimension).
     * Otherwise, just call 'compLogPart'.
     * <p>
     * NOTE: For many rules, it is useful to store some information here
     * (e.g. a vector) which is used by 'compMoments'. If so, store it in
     * a member.
     *
     * @param mu    Mean mu
     * @param cfact Cholesky factor L (diagonal)
     * @return      log Z
     */
    virtual double compLogPartDiag(const StVector& mu,
				   const StVector& cfact) = 0;

    /**
     * Computes mean 'tilmu' and covariance matrix 'tilcov' for the tilted
     * distribution prop. to
     *   exp(f(u)) N(u | mu,sigma), sigma = L*L^T.
     * <p>
     * NOTE: For many rules, information is passed from 'compLogPart' to
     * 'compMoments'. Do this via a member.
     * If g(u) = exp(f(u) - log Z), the moments are computed using
     *   E[g(u) x], E[g(u) x x^T],
     * where u = L x + mu so that x is N(0,I). If 'raw'==true, these raw
     * statistics are returned in 'tilmu', 'tilcov' rather than mean and cov.
     * matrix.
     *
     * @param mu     Mean mu
     * @param cfact  Cholesky factor L (lower tr.)
     * @param tilmu  Mean of tilted ret. here
     * @param tilcov Cov. matrix of tilted ret. here
     * @param raw    S.a. Def.: false
     */
    virtual void compMoments(const StVector& mu,const StMatrix& cfact,
			     StVector& tilmu,StMatrix& tilcov,
			     bool raw=false) = 0;

    /**
     * Same as 'compMoments', but we only compute the diagonal of the
     * tilted cov. matrix. Also, the Cholesky factor L is diagonal (since
     * the cov. matrix of Q is diagonal).
     *
     * @param mu      Mean mu
     * @param cfact   Cholesky factor L, diagonal
     * @param tilmu   Mean of tilted ret. here
     * @param tildiag Diag. of tilted cov. matrix ret. here
     * @param raw     See 'compMoments'. Def.: false
     */
    virtual void compMomentsDiag(const StVector& mu,const StVector& cfact,
				 StVector& tilmu,StVector& tildiag,
				 bool raw=false) = 0;

    /**
     * Computes log partition functions
     *   log Z_i = log E[exp(f_i(u))],
     * E[] over N(mu,sigma), sigma = L*L^T (Cholesky). Here, f=(f_i)
     * is multi-valued (number of comp. ret. by 'getFNum'), implem. by
     * 'compFMulti'.
     *
     * @param mu    Mean mu
     * @param cfact Cholesky factor L (lower tr.)
     * @param logz  Vector (log Z_i) ret. here
     */
    virtual void compLogPartMulti(const StVector& mu,const StMatrix& cfact,
				  StVector& logz) = 0;

    /**
     * Same as 'compLogPartMulti', but for diagonal L.
     */
    virtual void compLogPartMultiDiag(const StVector& mu,const StVector& cfact,
				      StVector& logz) = 0;

    /**
     * Computes expressions
     *   log E[exp(f(u) + f_i(u) - log Z)],
     * E[] over N(mu,sigma), sigma = L*L^T (L diag.). Here, f given by
     * 'compF', (f_i) by 'compFMulti'. log Z (constant) passed in 'logpart'.
     * <p>
     * ATTENTION: No general implementation in the moment!
     */
    virtual void compMultiSpecDiag(const StVector& mu,const StVector& cfact,
				   StVector& logz,double logpart) {
      throw NotImplemException(EXCEPT_MSG(""));
    }

    /**
     * Computes expected value
     *   E[f(u)], E[] over N(mu,sigma), sigma = L*L^T (Cholesky).
     * <p>
     * NOTE: For many rules, information is passed from 'compExpectLog'
     * to 'compLogMoments'. Use a member for that.
     *
     * @param mu    Mean mu
     * @param cfact Cholesky factor L (lower tr.)
     * @return      S.a.
     */
    virtual double compExpectLog(const StVector& mu,const StMatrix& cfact) = 0;

    /**
     * Same as 'compExpectLog', but for diagonal L.
     */
    virtual double compExpectLogDiag(const StVector& mu,
				     const StVector& cfact) = 0;

    /**
     * Computes expected values
     *   E[f_i(u)], E[] over N(mu,sigma), sigma = L*L^T (Cholesky).
     * Here, (f_i) given by 'compFMulti'.
     * <p>
     * ATTENTION: No general implementation in the moment!
     *
     * @param mu    Mean mu
     * @param cfact Cholesky factor L (lower tr.)
     * @param rvec  Vector with expect. values ret. here
     */
    virtual void compExpectLogMulti(const StVector& mu,const StMatrix& cfact,
				    StVector& rvec) {
      throw NotImplemException(EXCEPT_MSG(""));
    }

    /**
     * Same as 'compExpectLogMulti', but for diagonal L.
     * <p>
     * ATTENTION: No general implementation in the moment!
     */
    virtual void compExpectLogMultiDiag(const StVector& mu,
					const StVector& cfact,StVector& rvec) {
      throw NotImplemException(EXCEPT_MSG(""));
    }

    /**
     * The same as 'compExpectLogDiag', but in parallel 'compLogPartMultiDiag'
     * is computed (return via 'logz').
     * NOTE: 'compLogMomentsDiag' can be run afterwards just as well.
     * NOTE: The def. implementation here simply calls 'compLogPartMultiDiag',
     * followed by 'compExpectLogDiag'. Overwrite in subclasses to save
     * double computations!
     */
    virtual double compExpectLogLogPartMultiDiag(const StVector& mu,
						 const StVector& cfact,
						 StVector& logz) {
      // Default implementation: Overwrite in subclasses!
      compLogPartMultiDiag(mu,cfact,logz);
      return compExpectLogDiag(mu,cfact);
    }

    /**
     * Same as 'compExpectLogMultiDiag', but also 'compLogPartMultiDiag'
     * is computed (return via 'logz').
     * The 'ExpectLogMulti' precomp. is kept in 'logvecs2' and 'fvecFlag'
     * is set to 3, s.t. 'compLogMomentsMultiDiag' can be called.
     * The 'LogPartMulti' computation uses 'logvecs3' of the same size as
     * 'logvecs2'.
     * <p>
     * ATTENTION: No general implementation in the moment!
     *
     * @param mu      Mean mu
     * @param cfact   Cholesky factor L, diagonal
     * @param rvec    ExpectLogMulti vector ret. here (see
     *                'compExpectLogMultiDiag')
     * @param logz    LogPartMulti vector ret. here
     */
    virtual void compExpectLogMultiLogPartMultiDiag(const StVector& mu,
						    const StVector& cfact,
						    StVector& rvec,
						    StVector& logz) {
      throw NotImplemException(EXCEPT_MSG(""));
    }

    /**
     * If both the results of 'compExpectLogXXX' and 'compLogMomentsXXX'
     * are required, there are the following procedures:
     * - 0: 'compExpectLogXXX' is called first. Some precomp. is done there
     *      which allows 'compLogMomentsXXX' to be called
     * - 1: 'compLogMomentsXXX' returns both results together
     * The def. procedure is 0. Other procedures can be defined in subclasses
     * and assigned new numbers.
     *
     * @return Procedure code
     */
    virtual int logMomentsProc() const {
      return 0;
    }

    /**
     * Computes moment expressions
     *   c = L^-1 E[f(u) (u-mu)],
     *   G = L^-1 E[f(u) (u-mu) (u-mu)^T] L^-T,
     * E[] over N(u | mu,sigma), sigma = L*L^T.
     * <p>
     * NOTE: As opposed to 'compMoments', the linear and quadratic
     * parts are already centered here, and the covariance matrix is
     * divided out.
     * <p>
     * ATTENTION: Supports log moments procedure 0 only
     * (see 'logMomentsProc')!
     *
     * @param mu    Mean mu
     * @param cfact Cholesky factor L (lower tr.)
     * @param cvec  Vector c ret. here
     * @param gmat  Matrix G ret. here
     */
    virtual void compLogMoments(const StVector& mu,const StMatrix& cfact,
				StVector& cvec,StMatrix& gmat) = 0;

    /**
     * Same as 'compLogMoments', but here L, sigma are diagonal and we
     * only compute the diagonal of G.
     * <p>
     * Log moments procedure:
     * For 0 (def.), 'compExpectLogDiag' has to be called before. For 1,
     * the result of 'compExpectLogDiag' is ret. here via 'elog'.
     *
     * @param mu    Mean mu
     * @param cfact Cholesky factor L (diagonal)
     * @param cvec  Vector c ret. here
     * @param gvec  Vector diag(G) ret. here
     * @param elog  S.a. For procedure 1 only
     */
    virtual void compLogMomentsDiag(const StVector& mu,const StVector& cfact,
				    StVector& cvec,StVector& gvec,
				    double* elog=0) = 0;

    /**
     * Does same as 'compLogMomentsDiag', but for several f_i in parallel
     * (ret. by 'compFMulti'). The vectors for the indiv. components are ret.
     * in the cols of 'cmat' and 'gmat'.
     * <p>
     * ATTENTION: No general implementation in the moment!
     *
     * @param mu    Mean mu
     * @param cfact Cholesky factor L (diagonal)
     * @param cmat  c vectors ret. in cols
     * @param gmat  diag(G) vectors ret. in cols
     */
    virtual void compLogMomentsMultiDiag(const StVector& mu,
					 const StVector& cfact,
					 StMatrix& cmat,StMatrix& gmat) {
      throw NotImplemException(EXCEPT_MSG(""));
    }

    /**
     * Method required to compute gradient of learning criterion in the
     * 'Diag' case. The first part of the criterion is G_1 = -\sum_i z_i. The
     * vectors e1, e2 are
     * such that
     *   d G_1 = \sum_{c,i} e1_{i,c} d h_{i,c} + e2_{i,c} d a_{i,c}.
     * Here, i is fixed, the Q moments h_i, a_i^{1/2} are given in 'mu',
     * 'cfact', the i-parts of e1, e2 are ret. via 'e1vec', 'e2vec', the
     * method returns z_i.
     *
     * @param mu    S.a.
     * @param fact  S.a.
     * @param e1vec S.a.
     * @param e2vec S.a.
     * @return      S.a.
     */
    virtual double compExpectLogEVecsDiag(const StVector& mu,
					  const StVector& cfact,
					  StVector& e1vec,StVector& e2vec) = 0;

  protected:
    // Internal methods

    /**
     * Has to be provided by subclasses!
     *
     * @param u Eval. point
     * @return  Function value f(u)
     */
    virtual double compF(const StVector& u) const = 0;

    /**
     * Has to be provided by subclasses!
     *
     * @param u   Eval. point
     * @param val Function value f(u) ret. here
     */
    virtual void compFMulti(const StVector& u,StVector& val) const = 0;

    /**
     * Has to be provided by subclasses!
     *
     * @return Number of comp. of f (ret. by 'compFMulti')
     */
    virtual int getFNum() const = 0;

    /**
     * Some methods may precompute information required by others which
     * is stored in members. If this method is called, this "transfer"
     * information is marked as invalid.
     */
    virtual void resetTransfer() = 0;
  };
//ENDNS

#endif
