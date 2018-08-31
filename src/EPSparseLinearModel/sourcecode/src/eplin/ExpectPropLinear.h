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
 * Desc.:  Header class ExpectPropLinear
 * ------------------------------------------------------------------- */

/*
 * TODO:
 * This is quite a mess! Simplify the following:
 * -
 */

#ifndef EPLIN_EXPECTPROPLINEAR_H
#define EPLIN_EXPECTPROPLINEAR_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eplin/default.h"
#include "src/eplin/EPPrior.h"
#include "src/eplin/EPSingleSiteManager.h"
#include "lhotse/matrix/StVector.h"
#include "lhotse/matrix/StMatrix.h"
#include "lhotse/matrix/IndexVector.h"
#include "lhotse/rando/Generator.h"
//#include "lhotse/DebugMacros.h" // DEBUG

//BEGINNS(eplin)
#define CHECKCONSISTENT do { if (prior->cols()!=sites->size()) \
  throw InvalidParameterException(EXCEPT_MSG("")); } while (0)

#define UPDATE_REPRES do { if (!up2date) refresh(); } while (0)

  /**
   * Maintains and updates expectation propagation (EP) approximate inference
   * posterior representation for model with degenerate Gaussian prior
   * (flat along subspace, not normalizable) of particular form and product
   * of univariate sites (type 'EPSingleSite', all sites maintained by
   * 'EPSingleSiteManager'). The degen. Gaussian part is proportional to
   *   N^U(u | sigma^-2 b^(0), sigma^-2 X^T), b^(0) = X^T v,
   * and is represented by a 'EPPrior' object (which represents X, v).
   * <p>
   * Example is linear model on n variables, m datapoints, where likelihood
   * (Gaussian) gives Gaussian factor, sigma^2 is noise variance. Sites
   * come from factorizing prior, f.ex. Laplace for sparse linear model.
   * <p>
   * The posterior approximation is Gaussian, obtained by replacing the
   * likelihood by the site approximation
   *   N^U(u | sigma^-2 b, sigma^-2 Pi),
   * Pi diagonal with positive entries. Because of the degeneracy of the
   * prior, all entries pi_i must be positive, we cannot start EP as usual
   * from Pi = 0.
   * EP is in fact started from Pi = eps I, eps>0 passed at construction,
   * and b = 0, which corresponds to the posterior for a Gaussian likelihood
   * N(u | 0,sigma^2 eps^-1 I), instead of the one given by 'sites'. Updates
   * which would result in very small positive or even non-pos. pi_i values,
   * are rejected (see PI_LOWER_BOUND constant).
   * <p>
   * X is m-by-n. We distinguish two cases: degenerate (n > m), and
   * non-degenerate. The repres. in the degen. case scales as O(m^2), the
   * non-degen. as O(n^2). The non-degen. repres. is used if m >= 'nonDegenM',
   * the degen. repres. is used otherwise. We may switch in 'priorChange'.
   * 'nonDegenM' must be <= n. If 'nonDegenM'==0, the non-degen. repres. is
   * always used.
   * NOTE: The non-degen. repres. leads to more stable computations. The
   * degen. is more efficient if m << n, but the advantage vanishes at about
   * m = n/2. 'nonDegenM' > n/2 may lead to numerical instabilities.
   *
   * ATTENTION: Distinguish between degen./non-degen. case and
   * degen./non-degen. repres.! Degen. case iff n > m. In the non-degen.
   * case, have to use non-degen. repres., but can use non-degen. repres.
   * also in the degen. case (for stability). The constraints on Pi depend
   * on the case, not the repres. We need pi_i >= PI_LOWER_BOUND if n > m,
   * even if the non-degen. repres. is used.
   * <p>
   * Degenerate EP representation (m < 'nonDegenM'):
   * The posterior cov. matrix is
   *   A = sigma^2 [ X^T X + Pi ]^-1
   *     = sigma^2 [ Pi^-1 - Pi^-1 X^T (I + X Pi^-1 X^T)^-1 X Pi^-1 ]
   * It is represented by L L^T = I + X Pi^-1 X^T (Cholesky), L in
   * 'factL' (lower tri.). Furthermore,
   *   gamma = L^-1 X Pi^-1 (b^(0) + b) in 'gamma',
   * so that the posterior mean is
   *   h = Pi^-1 (b^(0) + b - X^T L^-T gamma).
   * <p>
   * Non-degenerate EP representation (m >= 'nonDegenM'):
   * The posterior cov. matrix is
   *   A = sigma^2 [ X^T X + Pi ]^-1
   * It is represented by L L^T = sigma^2 A^-1 (Cholesky), L in 'factL'
   * (lower tri.). Furthermore,
   *   gamma = L^-1 (b^(0) + b),
   * so that the posterior mean is
   *   h = L^-T gamma.
   * <p>
   * Change in prior:
   * The EP repres. (L, gamma) must be consistent with b, Pi, and the prior
   * at all times. In general, 'refresh' is recomputing the repres. from
   * scratch, which must be done if the prior changes in general.
   * In many applications, the prior changes as follows: X' = [X; x^T] (new
   * row x), and v' = [v; v_*]^T (new element v_*). In this case,
   * 'priorChange' can be used instead of 'refresh'. The EP repres. is
   * updated using low rank formulae in this case, which is much faster
   * than 'refresh'. 'priorChange' must be called AFTER the 'prior' object
   * has been modified in the said way.
   * If 'priorChange' is called and m=='nonDegenM', this means we switch
   * from the degenerate repres. to the non-degenerate one. In this case,
   * 'refresh' is called to compute the repres. from scratch: it cannot be
   * derived from the degen. repres.
   * <p>
   * EP updates:
   * This is done locally for a single site using 'updateSite'. From the
   * current site marginal, the cavity marginal is computed, then the ADF
   * projection is done through the site manager 'sites'. Finally, the EP
   * representation is updated by low rank formulae. Note that the site
   * marginals are computed on demand, they are not maintained explicitly
   * between updates.
   * <p>
   * Refreshing:
   * The EP repres. L, gamma should be recomp. from scratch now and then,
   * maybe after each complete sweep over all sites. Done by 'refresh'.
   * Also has to be used if the prior (see 'priorChange') or the sites change.
   * <p>
   * Temp. buffers:
   * Can be passed (handles) upon construction, oth. they are alloc. here.
   * They are used within methods here, but they can be changed between calls
   * methods in any way. If several 'ExpectPropLinear' objects are used in a
   * serial way, they can share their buffers that way.
   * <p>
   * Keep marginals up-to-date at all times:
   * 
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class ExpectPropLinear
  {
  protected:
    // Members

    double eps;
    Handle<EPSingleSiteManager> sites; // Sites and their pars.
    StMatrix factL;                    // Posterior repres.
    StVector gamma;
    bool up2date;                      // Is repres. up-2-date?
    Handle<EPPrior> prior;             // Maintains X, v
    int nonDegenM;
    int verbose;                       // Verbosity level
    mutable Handle<StVector> tempvec;
    mutable Handle<StMatrix> tempmat,tempmat2,tempmat3;
    mutable StVector tempvec2,tempvec3;
    mutable IndexVector tempind;
    //public:
    //DCOMP_DEFVARS; // DEBUG
    //protected:
    // DEBUG:
    //public:
    //mutable  double debScal;
    //mutable StVector debVec,debVec2;
    //mutable StMatrix debMat,debMat2;
    //mutable ifstream* deb_is;
    // DEBUG:
    //mutable StVector deb_vvec,deb_gamma;
    //mutable StMatrix deb_factL,deb_mat;
    //mutable int deb_cnt;

  public:
    // Public methods

    /**
     * Constructor. Starts from b = 0, Pi = 'peps' I. In the degenerate case
     * (n>m), 'peps' must be larger than PI_LOWER_BOUND.
     * 'nonDegM' is the value for 'nonDegenM' (see header comment), must be
     * <= n. If 'nonDegM'==0, the non-degen. repres. is always used.
     * The default is n/2.
     * We need some temp. matrices, which can be passed via 'buffmatXXX',
     * otherwise they are alloc. here.
     * <p>
     * Initial representation and site pars. are set acc. to 'istat':
     * - 0: Reset site pars. to 0, 'peps' (def.)
     * - 1: Do not change site pars. in 'psites'
     * The representation is not computed here, we start with
     * 'up2date'==false. It is computed (calling 'refresh') once first
     * needed.
     *
     * @param peps     Value for 'eps'. S.a.
     * @param pprior   Prior object
     * @param psites   Site manager object
     * @param nonDegM  S.a. Def.: n/2
     * @param istat    S.a. Def.: 0
     * @param buffmat1 S.a. Optional
     * @param buffmat2 S.a. Optional
     * @param buffmat3 S.a. Optional
     * @param buffvec  S.a. Optional
     */
    ExpectPropLinear(double peps,const Handle<EPPrior>& pprior,
		     const Handle<EPSingleSiteManager>& psites,int nonDegM=-1,
		     int istat=0,const Handle<StMatrix>& buffmat1=
		     HandleZero<StMatrix>::get(),
		     const Handle<StMatrix>& buffmat2=
		     HandleZero<StMatrix>::get(),
		     const Handle<StMatrix>& buffmat3=
		     HandleZero<StMatrix>::get(),
		     const Handle<StVector>& buffvec=
		     HandleZero<StVector>::get());

    virtual ~ExpectPropLinear() {}

#ifdef MATLAB_MEX
    /**
     * Specialized method. Do not use!
     *
     * Used to run EP within a Matlab MEX file. In this case, all buffers
     * should come from Matlab, there should not be allocations here.
     * The buffers required here are 'factL', 'gamma', 'tempvec2',
     * 'tempvec3' ('*tempvec1' is passed at construction already).
     * NOTE: This method has to be called directly after construction,
     * when these variables are still empty and can become masks.
     *
     * ATTENTION: In MEX mode, 'tempvecXXX' and 'tempmat' are masks, the
     * vectors of size n. Their size cannot be changed, but some code
     * requires shorter vectors. Has to be dealt with in code which is
     * required by the MEX file!
     */
    virtual void provideBuffers(const StMatrix& mfactL,const StVector& mgamma,
				const StVector& mtempvec2,
				const StVector& mtempvec3) {
      int n=prior->cols();
      int sz=isDegenerate()?(prior->rows()):n;
      if (!factL.isMaskOrEmpty() || !gamma.isMaskOrEmpty() ||
	  !tempvec2.isMaskOrEmpty() || !tempvec3.isMaskOrEmpty())
	throw WrongStatusException(EXCEPT_MSG(""));
      if (mfactL.rows()!=sz || mfactL.cols()!=sz || mgamma.size()!=sz ||
	  mtempvec2.size()!=n || mtempvec3.size()!=n)
	throw WrongDimensionException(EXCEPT_MSG(""));
      factL.reassign(mfactL); gamma.reassign(mgamma);
      tempvec2.reassign(mtempvec2); tempvec3.reassign(mtempvec3);
    }
#endif

    /**
     * The non-degen. representation is used iff m >= 'nonDegenM', where
     * m number of datapoints ('prior->rows()').
     *
     * @return Do we use the degenerate representation?
     */
    virtual bool isDegenerate() const {
      CHECKCONSISTENT;
      return (prior->rows()<nonDegenM);
    }

    /**
     * @return Prior object
     */
    virtual Handle<EPPrior>& getPrior() {
      return prior;
    }

    /**
     * Change prior object. 'up2date' is set to false, so the repres. will be
     * recomp. when needed next.
     * NOTE: This is different from 'priorChange', which is about extending
     * the prior, and the repres. need not be recomp. from scratch there.
     *
     * @param pr New prior object
     */
    virtual void setPrior(const Handle<EPPrior>& pr) {
      if (pr->cols()!=sites->size())
	throw InvalidParameterException(EXCEPT_MSG(""));
      prior=pr;
      up2date=false;
    }

    /**
     * @return Sites manager object
     */
    virtual Handle<EPSingleSiteManager>& getSites() {
      return sites;
    }

    /**
     * Change sites object. 'up2date' is set to false, so repres. will be
     * recomp. from scratch when next needed.
     *
     * @param si New sites manager object
     */
    virtual void setSites(const Handle<EPSingleSiteManager>& si) {
      if (prior->cols()!=si->size())
	throw InvalidParameterException(EXCEPT_MSG(""));
      sites=si;
      up2date=false;
    }

    /**
     * @return Value of 'eps' constant
     */
    double getEps() const {
      return eps;
    }

    /**
     * Verbosity level:
     * - 0: No messages
     * - 1: Some messages
     * - 2: Many messages
     *
     * @param verb New level
     */
    void setVerbose(int verb) {
      if (verb<0) throw InvalidParameterException(EXCEPT_MSG(""));
      verbose=verb;
    }

    /**
     * Recomputes EP representation from scratch. Prior and/or sites may
     * have been changed.
     */
    virtual void refresh();

    /**
     * ADF update local at site 'i', which changes the EP repres. globally.
     * For the degenerate repres., the update is not done if the change in
     * 1/pi_i would be smaller than 'thres'. For the non-degen. repres., it
     * is not done if the change in pi_i would be smaller than 'thres'.
     * <p>
     * NOTE: In the degenerate case, if the update of pi_i would lead to a
     * value < PI_LOWER_BOUND, it is fixed to that value. In the
     * non-degen. case, we do not allow pi_i to become negative (fix it to
     * 0 if that is proposed).
     * If 'delta' is given, then |h-h'| and |v-v'| are computed, where
     * the old marginal is N(h,v), and the max. of these is ret. there.
     * <p>
     * Fractional update:
     * If 'frac'!=1, it must be in (0,1). A fractional ADF update is done in
     * this case, with fraction 'frac'. Requires 'sites' to support fractional
     * sites.
     * <p>
     * Update type:
     * Depending on 'utype', other types of update can be done (they differ
     * only in the comput. of the new pi_i, b_i):
     * - 0: Normal ADF update (def.)
     * - 1: ADF update with b_i constrained to 0
     * - 2: Tipping SBL update (b_i=0)
     * - 3: Girolami variational update (b_i=0)
     *
     * @param i     Site index
     * @param thres S.a. Positive
     * @param frac  S.a. Def.: 1
     * @param delta S.a. Def.: 0
     * @param utype S.a. Def.: 0
     * @return      Has update been done?
     */
    virtual bool updateSite(int i,double thres,double frac=1.0,
			    double* delta=0,int utype=0);

    /**
     * Computes log Z values for all sites, given marg. means, variances
     * in 'means', 'vars'.
     * If 'exabs' is given, we also return the E[-|a|] values there.
     * This requires Laplace sites, 'sites' must be of type
     * 'EPSingleLaplaceManager'.
     *
     * @param means Marg. means
     * @param vars  Marg. vars.
     * @param logz  Log Z values ret. here
     * @param exabs S.a. Def.: 0
     */
    virtual void compLogZ(const StVector& means,const StVector& vars,
			  StVector& logz,StVector* exabs=0);

    /**
     * Computes approximation to -log P(D) (neg. log marginal likelihood),
     * based on the current representation. This makes sense only (at
     * least for the gradient) if the site parameters have converged.
     * 'meth' selects the approx. inf. method, see 'updateSite'.
     *
     * @param gradX    Gradient w.r.t. code matrix X ret. here
     * @param gradSig  Deriv. w.r.t. sigma^-2 ret. here
     * @param gradTau  Deriv. w.r.t. tau ret. here
     * @param meth     S.a. Def.: 0 (EP)
     * @param wkvecm   Scratch vector (size m). Def.: 0 (alloc. here)
     * @param wkvecn   Scratch vector (size n). Req. only for grads.
     *                 Def.: 0 (alloc. here)
     */
    double compCritGrad(StMatrix& gradX,double& gradSig,double& gradTau,
			int meth=0,StVector* wkvecm=0,StVector* wkvecn=0,
			int gradStat=2);

    double compCritGrad(double& gradSig,double& gradTau,int meth=0,
			StVector* wkvecm=0,StVector* wkvecn=0) {
      StMatrix dummy;

      return compCritGrad(dummy,gradSig,gradTau,meth,wkvecm,wkvecn,1);
    }

    double compCritGrad(int meth=0,StVector* wkvecm=0) {
      StMatrix dummy;
      double temp,temp2;

      return compCritGrad(dummy,temp,temp2,meth,wkvecm,0,0);
    }

    /**
     * NOTE: Means diff. things in degenerate/non-degen. repres.!
     *
     * @return Matrix L (EP repres.)
     */
    virtual const StMatrix& getFactL() {
      CHECKCONSISTENT;
      UPDATE_REPRES;

      return factL;
    }

    /**
     * NOTE: Means diff. things in degenerate/non-degen. repres.!
     *
     * @return Vector gamma (EP repres.)
     */
    virtual const StVector& getGamma() {
      CHECKCONSISTENT;
      UPDATE_REPRES;

      return gamma;
    }

    /**
     * This method must be called if the prior in 'prior' just has been
     * changed, in that X' = [X; x^T] (new row x), and a new entry v_* has
     * been appended to the v vector.
     * NOTE: 'prior' must have been modified in this way BEFORE this method
     * is called, so m (# rows) is already at its new value.
     * If the change results in m=='nonDegenM', we switch to the non-degen.
     * repres., by simply calling 'refresh'.
     */
    virtual void priorChange();

    /**
     * Resets site parameters to b = 0, Pi = eps I. The representation is
     * recomputed.
     */
    virtual void reset();

    /**
     * Initializes the EP repres. from the arguments. L, gamma must be passed,
     * they are not checked. If they are not passed, 'refresh' is called.
     * ATTENTION: L, gamma depend on repres. type (degen./non-degen.).
     *
     * @param newb   New value for b vector
     * @param newpi  New value for pi vector
     * @param newl   New value for L factor. Optional
     * @param newgam New value for gamma vector. Optional
     */
    virtual void initRepres(const StVector& newb,const StVector& newpi,
			    const StMatrix* newl=0,const StVector* newgam=0);

    /**
     * @param means Marginal means h ret. here
     * @param vars  Marg. variances ret. here. Optional
     */
    virtual void getMargMoments(StVector& means);

    virtual void getMargMoments(StVector& means,StVector& vars);

    /**
     * Generates sample of 'num' i.i.d. vectors from the currently repres.
     * posterior.
     *
     * @param num  Number of vectors to be sampled
     * @param smat Samples ret. in cols. of this matrix
     * @param gen  PRN generator
     */
    virtual void sample(int num,StMatrix& smat,Handle<Generator>& gen);

    /**
     * Used to run MCMC method of Park, Casella within our framework.
     * See 'ExpectPropJointLinear'. Init. the repres. and site pars., where
     * (pi_i) are passed, b_i are set to 0.
     *
     * @param newpi New values for (pi_i)
     */
    virtual void mcmcInitRepres(const StVector& newpi) {
      int n=prior->cols();
      StVector d1;

      if (newpi.size()!=n || !newpi.checkBounds(DefIVal<double>::posit())) {
	//newpi.print(cout); cout << endl; // DEBUG
	throw InvalidParameterException(EXCEPT_MSG("newpi"));
      }
      d1.zeros(n);
      initRepres(d1,newpi);
    }

    /**
     * Used to run MCMC method of Park, Casella within our framework.
     * See 'ExpectPropJointLinear'. Samples a, given the current pi.
     *
     * @param newa     New a vector ret. here
     * @param gen      PRN generator
     */
    virtual void mcmcSampleA(StVector& newa,Handle<Generator>& gen)
    {
      newa.zeros(prior->cols());
      sample(1,StMatrix::mask(newa),gen);
    }

    /**
     * Used to run MCMC method of Park, Casella within our framework.
     * See 'ExpectPropJointLinear'. Samples pi, given a. These are written
     * into 'sites', and the repres. is updated.
     *
     * @param curra     Current a vector
     * @param gen       PRN generator
     */
    virtual void mcmcSamplePi(const StVector& curra,Handle<Generator>& gen);

    /**
     * Used to estimate log marginal likelihood by Chib's method alongside
     * an MCMC run. Evaluates the log of the full conditional density
     * of a, given the current pi. The normalization constant has
     * be correct here. This is the density we sample from in
     * 'mcmcSampleA'.
     *
     * @param a     Vector a to eval. at
     * @return      log full cond. density
     */
    virtual double mcmcEvalFullCond(const StVector& a);

    /**
     * Computes a number of scores, dep. on the current posterior. 'exl1'
     * is the expected L_1 norm |u-v|, where the expect. is over u from the
     * posterior.
     *
     * @param v       Input vector
     * @param ent     Posterior entropy ret. here
     * @param logprob Log-posterior density value of v ret. here
     * @param exl1    S.a.
     */
    virtual void compScores(const StVector& v,double& ent,double& logprob,
			    double& exl1);

    /**
     * DEBUG method
     *
     * @param amat Poster. cov. (or its inv.) ret. here
     * @param inv  Get inverse? Def.: false
     */
    virtual void debugGetPostCov(StMatrix& amat,bool inv=false);

    // FLORIAN
    virtual void addMargVarRank(BaseVector<int>& ranks);
    // END FLORIAN
  };
#undef CHECKCONSISTENT
#undef UPDATE_REPRES
//ENDNS
//#include "lhotse/DebugMacros_undef.h" // DEBUG

#endif
