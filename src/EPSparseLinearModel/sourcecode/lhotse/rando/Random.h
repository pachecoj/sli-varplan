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
 * Module: rando
 * Desc.:  Header class Random
 * ------------------------------------------------------------------- */

/*
 * TODO:
 * - In general: wrap one of the free PRN packages!
 * - Implement 'devPoisson', 'devBinomial' !!
 * - The "quick" accept./reject. tests are possibly outdated on modern
 *   processors which can computed exp,log and trigonometric functions
 *   fast. Test if there's a real need for them, or if they even slow
 *   things down!
 */

#ifndef RANDOM_H
#define RANDOM_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/rando/default.h"
#include "lhotse/rando/Generator.h"
#include "lhotse/matrix/predecl.h"
#include "lhotse/matrix/BaseVector.h"

USING(matrix);

//BEGINNS(rando)
  /**
   * Provides static methods to sample from common parametric distributions
   * and other random-related utility methods
   * <p>
   * We use the following source for algorithms generating non-uniform
   * random variates:
   *   Devroye, L.
   *   Non-Uniform Random Variate Generation
   *   Online at cgm.cs.mcgill.ca/~luc/rnbookindex.html
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class Random
  {
  protected:
    // Internal static members

    static Generator* defGen; // default generator

  public:
    // Public static methods

    /**
     * Sets default generator. The old generator is destroyed
     */
    static void setDefGen(Generator* gen) {
      if (defGen!=0) delete defGen;
      defGen=gen;
    }

    /**
     * Returns default generator
     */
    static Generator* getDefGen() {
      return defGen;
    }

    // Deviates with univariate continuous densities

    /**
     * Generates uniform [a,b] deviate.
     */
    static double devUniform(double a,double b,Generator* gen=defGen) {
      if (a>b) throw InvalidParameterException(EXCEPT_MSG(""));

      return (b-a)*gen->get()+a;
    }

    /**
     * Generates normal deviate using Box-Muller method (polar method). The
     * default is a zero-mean, unit-variance variate
     *   P(x|mu,rho) = (2\pi \rho^2)^{-1} exp(-0.5*((x-mu)/rho)^2)
     * <p>
     * NOTE: In each active turn, we generate TWO normal variates. Every
     * second turn is non-active and simply returns the second of the two
     * variates generated in the last turn. Therefore, 'gen' is only really
     * used every second turn!
     */
    static double devNormal(Generator* gen=defGen);
    static double devNormal(double mean,double stddev,Generator* gen=defGen);

    /**
     * Generates normal deviate using Ratio-of-Uniforms method. See
     * Devroye, p.194ff. This method is usually faster than Box-Muller.
     * The default is a zero-mean, unit-variance variate
     */
    static double devNormalRatUnif(Generator* gen=defGen);
    static double devNormalRatUnif(double mean,double stddev,
				   Generator* gen=defGen);

    /**
     * Generates exponential deviate using simple inversion method. When no
     * mean is given, unit mean is assumed.
     * <p>
     * Note: Faster methods exist, but only if machine-code implementations
     * are considered
     *   P(x|mu) = mu^{-1} exp(-x/mu)
     */
    static double devExp(Generator* gen=defGen);
    static double devExp(double mean,Generator* gen=defGen);

    /**
     * Generates gamma deviate using Best's rejection algorithm based
     * upon rejection from the t distribution with 2 degrees of freedom
     * (for a>=1) or a transformation of an EPD variate (obt. by rejection
     * method) (for a<1). The default scale parameter is b=1.
     * Mean and variance are a*b and a*b^2.
     *   P(x|a,b) = (b^a \Gamma(a))^{-1} x^{a-1} exp(-x/b)
     * <p>
     * Another standard parameterization of the Gamma distribution is in
     * terms of mean mu and shape parameter alpha (see Neal). The variance
     * in this case is 2*mu^2/alpha. Conversion in our format is done by:
     *   alpha=2*a, mu=a*b or a=alpha/2,b=2*mu/alpha.
     * <p>
     * Another standard parameterization uses a shape parameter alpha more
     * compatible with the Chi-Square distribution, and an inverse scale
     * parameter beta. Conversion is done by:
     *   alpha=2*a, beta=2/b or a=alpha/2, b=2/beta
     * <p>
     * NOTE: Best's algorithm is chosen for simplicity. If speed is more of
     * an issue, the recommendations in Devroye, p.401ff, should be
     * considered. Maybe an almost-exact-inversion method should be
     * implemented.
     */
    static double devGamma(double a,Generator* gen=defGen);
    static double devGamma(double a,double b,Generator* gen=defGen);

    /**
     * Generates beta deviate using Cheng's rejection algorithm (if
     * min(a,b) > 1 or a+b > 1.5) or Johnk's method (otherwise). There
     * are special algorithms for a<1, b>1 that outperform Cheng's method,
     * but they are curr. not implemented.
     *   P(x|a,b) = (Gamma(a+b)/(Gamma(a)*Gamma(b)))*x^{a-1}*(1-x)^{b-1}
     */
    static double devBeta(double a,double b,Generator* gen=defGen);

    /**
     * Generates t deviate (a shape parameter, "degrees of freedom"). We
     * use a ratio-of-uniforms method for a>=1. For the rare case a<1,
     * we use a transformation method which might be slow.
     * The mean exists for a>1 (default is 0). The variance exists for a>2
     * (default is 1).
     * <p>
     * Note: If special cases like Cauchy are frequently required, it
     * would be better to implement special generators!
     *   P(x|mu,sigma,a) = (Gamma((a+1)/2)/(Gamma(a/2)*sqrt(\pi*a)*rho))*
     *                     [1 + rho^{-1}*((x-mu)/rho)^2]^{-(a+1)/2}
     */
    static double devT(double a,Generator* gen=defGen);
    static double devT(double a,double mean,double stddev,
		       Generator* gen=defGen);

    // Deviates with univariate discrete densities

    /**
     * Generates Bernoulli variate with mean p. The default mean is 1/2
     */
    static int devBernoulli(double p=0.5,Generator* gen=defGen);

    /**
     * Generates geometric variate with success prob. p. The default p is 1/2
     */
    static int devGeometric(double p=0.5,Generator* gen=defGen);

    /**
     * Generates Poisson variate with rate gamma>0.
     * NOT IMPLEMENTED IN THE MOMENT!
     */
    static int devPoisson(double gamma,Generator* gen=defGen) {
      throw NotImplemException(EXCEPT_MSG("NOT IMPLEMENTED"));
    }

    /**
     * Generates binomial variate with parameters n>=1 and pp.
     * NOT IMPLEMENTED IN THE MOMENT!
     */
    static int devBinomial(int n,double pp=0.5,Generator* gen=defGen) {
      throw NotImplemException(EXCEPT_MSG("NOT IMPLEMENTED"));
    }

    /**
     * Returns variate drawn uniformly from a,...,b, where a<b, a,b are
     * integers.
     *
     * @param a   Limit a
     * @param b   Limit b. Must be >a
     * @param gen Optional. PRN gen.
     * @return    Variate
     */
    static int devUniInt(int a,int b,Generator* gen=defGen);

    /**
     * Generates variate from discrete distribution p_0,...,p_{n-1}, given in
     * 'distr' (contains values log p_i). The variate lies in {0,...,n-1}.
     * <p>
     * NOTE: Needs O(n), inefficient if many variates are drawn from the
     * same distribution, in which case after some O(n) precomputation a
     * variate can be drawn in O(log n) (binary search on the cumul.
     * distrib. function), or even in O(1) (see Knuth: Art of Computer
     * Programming).
     *
     * @param distr Log probabilities (must be normalized)
     * @param gen   PRN generator to be used
     * @return      Variate
     */
    static unsigned devDiscrete(const BaseVector<double>& distr,
				Generator* gen=defGen);

    /**
     * Generates variate from inverse Gaussian distribution, using a
     * Chi^2 (1 dof) and a uniform variate. This is a simple method
     * described in Chhikara and Folks (1989), Sect. 4.5, and implemented
     * in the STATMOD R package as INVGAUSS.
     * <p>
     * The inverse Gaussian density is
     *   P(x | mu,lam) = sqrt(lam/(2*pi)) x^{-3/2}
     *     exp[ -lam (x-mu)^2 / (2 mu^2 x) ], x>0
     * It has mean mu>0 and scale parameter lam>0. Note that y = a x, a>0,
     * has inverse Gaussian distribution with mean mu a, scale lam a.
     *
     * @param mu  Mean (positive)
     * @param lam Scale (positive). Def.: 1
     * @param gen PRNG. Def.: defGen
     * @return    Variate
     */
    static double devInvGauss(double mu,Generator* gen=defGen) {
      return devInvGauss(mu,1.0,defGen);
    }

    static double devInvGauss(double mu,double lam,Generator* gen=defGen);

    /**
     * Generates variate from distribution with log-concave density P(x)
     * and mode at m. Rejection method, where each trial requires two uniform
     * variates (three for types 1,2) and a single evaluation of log P.
     * There are several types, selected by 'dtype':
     * - 0: P supported on [m,infty)
     * - 1: P supported on both sides of m, symmetric around m
     * - 2: P supported on both sides of m
     * The accept probability in each trial is >= 1/2 for types 0,1, >= 1/4
     * for type 2. Taken from Devroye, p. 291.
     * <p>
     * Let P(x) = c f(x), c>0. 'lfunc' is a pointer to log f(x). 'lfm' is
     * log f(m), 'logc' is log(c).
     *
     * @param lfunc Pointer to log f(x) function
     * @param m     Value of m
     * @param lfm   Value of log f(m)
     * @param logc  Value of log(c)
     * @param dtype Type of f. Def.: 0
     * @param gen   PRNG. Def.: defGen
     * @return      Variate
     */
    static double devLogConcave(double (*lfunc) (double),double m,double lfm,
				double logc,int dtype=0,Generator* gen=defGen);

    // Deviates with multivariate continuous densities

    /**
     * Generates multivariate normal vector, given mean and Cholesky decomp.
     * of covariance matrix. Here, 'chCov' is the lower-tr. Chol. factor.
     * <p>
     * NOTE: In any case, 'mean' and 'targ' may be the same vector. In this
     * case, 'mean' is overwritten by the sample.
     *
     * @param targ   Random vector returned here. Must have correct size
     * @param mean   Optional. Default is 0
     * @param chCov  Optional. Cholesky decomp. of cov. matrix. Default cov.
     *               is I.
     * @param gen    Optional. PRN generator
     */
    static void devMultiNormal(BaseVector<double>& targ,Generator* gen=defGen);
    static void devMultiNormal(BaseVector<double>& targ,
			       const BaseVector<double>& mean,
			       Generator* gen=defGen);
    static void devMultiNormal(StVector& targ,
			       const StVector& mean,
			       const StMatrix& chCov,
			       Generator* gen=defGen);

    /**
     * Frequently useful special case of 'devMultiNormal'. Here, we produce
     * a sample from N(\mu,\Sigma), where \mu = \Sigma r. The vector r and
     * the Cholesky decomp. of the inverse \Sigma^{-1} are given.
     * <p>
     * NOTE: 'r' and 'targ' may be the same vector. In this case, 'r' is
     * overwritten by the sample.
     *
     * @param targ    Random vector returned here. Must have correct size
     * @param r       r vector. See above
     * @param chICov  Cholesky factor for \Sigma^{-1}
     * @param gen     Optional. PRN generator
     */
    static void devMultiNormalInv(StVector& targ,
				  const StVector& r,
				  const StMatrix& chICov,
				  Generator* gen=defGen);

    /**
     * Generates random vector from the Dirichlet(\alpha_1,\dots,\alpha_k)
     * distribution. This involves sampling of k Gamma variates.
     * <p>
     * NOTE: For k=2, 'devBeta' is more efficient!
     *
     * @param alpha Vector with Dirichlet parameters
     * @param rv    Random vector returned here
     * @param gen   PRN generator to be used
     */
    static void devDirichlet(const BaseVector<double>& alpha,
			     BaseVector<double>& rv,
			     Generator* gen=defGen);

    // Other random objects

    /**
     * Generates random permutation of 0,...,n-1 and stores it into given
     * integer vector. The swapping method is used.
     * NOTE: If 'jumble' is true (def.: false), we assume that 'iv' already
     * contains 'n' entries and simply apply a random permutation on this
     * vector.
     *
     * @param n      Size for random perm.
     * @param iv     Permutation vector ret. here. If 'jumble'==true: See above
     * @param gen    Optional. PRN generator
     * @param jumble Optional. Def.: false. See above
     */
    static void randomPerm(int n,BaseVector<int>& iv,Generator* gen=defGen,
			   bool jumble=false);

    /**
     * Picks d elements from 0,...,n-1 at random, without repetition.
     * <p>
     * NOTE: If 'shift' is true (def.: false), we do the following: we assume
     * that 'iv' already contains exactly n elements. We pick d of these at
     * random and move them to the first d positions within 'iv' by swapping
     * them with other further down. More precisely, we do d random
     * transpositions on the vector 'iv', where the j-th transposition exchanges
     * elements j-1 and k, where k is chosen at random from j-1,...,n-1. Note
     * that for 'shift'==false, we do the same on a temp. vector init. with
     * 0,...,n-1, then copy its first d elements to 'iv'.
     * Note also that 'randomPerm' uses the same technique with d==n.
     *
     * @param n     Range
     * @param d     Number of entries to select
     * @param iv    Result returned here. If 'shift'==true: See above!
     * @param gen   Optional. PRN gen.
     * @param shift Optional. Def.: false. See above
     */
    static void randomSet(int n,int d,BaseVector<int>& iv,
			  Generator* gen=defGen,bool shift=false);
  };

  inline double Random::devNormal(double mean,double stddev,Generator* gen)
  {
    if (stddev<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
    return mean + stddev*devNormal(gen);
  }

  inline double Random::devNormalRatUnif(double mean,double stddev,
					Generator* gen)
  {
    if (stddev<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
    return mean + stddev*devNormalRatUnif(gen);
  }

  inline double Random::devExp(Generator* gen)
  {
    // The negative log of a U[0,1] variable has exponential distribution
    if (gen==0) throw NoGeneratorException();
    return -log(gen->get());
  }

  inline double Random::devExp(double mean,Generator* gen)
  {
    if (mean<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
    if (gen==0) throw NoGeneratorException();
    return -mean*log(gen->get());
  }

  inline double Random::devGamma(double a,double b,Generator* gen)
  {
    if (b<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
    return b*devGamma(a,gen);
  }

  inline double Random::devT(double a,double mean,double stddev,
			    Generator* gen)
  {
    // Note that the mean does not exist for a<=1, the variance does not
    // exist for a<=2. We return a linear combination of a standard deviate
    // even in those cases!
    if (stddev<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
    return mean + stddev*devT(a,gen);
  }

  inline int Random::devBernoulli(double p,Generator* gen)
  {
    if (gen==0) throw NoGeneratorException();
    if (p<0.0 || p>1.0) throw InvalidParameterException(EXCEPT_MSG(""));
    return (gen->get() <= p)?1:0;
  }

  inline int Random::devGeometric(double p,Generator* gen)
  {
    // See Devroye, p.500
    if (gen==0) throw NoGeneratorException();
    if (p<=0.0 || p>=1.0) throw InvalidParameterException(EXCEPT_MSG(""));
    return (int) ceil(log(gen->get())/log(1-p));
  }

  inline int Random::devUniInt(int a,int b,Generator* gen)
  {
    int n=b-a+1;

    if (n<1) throw InvalidParameterException(EXCEPT_MSG("Invalid range"));

    return std::min(a+(int) floor(gen->get()*((double) n)),b);
  }

  inline unsigned Random::devDiscrete(const BaseVector<double>& distr,
				      Generator* gen)
  {
    register int i=0,n=distr.size();
    register double sum=0.0,unif=gen->get();

    for (; i<n; i++)
      if (unif<(sum+=exp(distr[i]))) break;
    if (i==n) i--; // happens with prob. 1 (if 'distr' normalized)

    return i;
  }

  // R code from STATMOD package (Smyth):
  //
  // rinvgauss <- function(n, mu, lambda = 1)
  // {
  // 	if(any(mu<=0)) stop("mu must be positive")
  // 	if(any(lambda<=0)) stop("lambda must be positive")
  // 	if(length(n)>1) n <- length(n)
  // 	if(length(mu)>1 && length(mu)!=n) mu <- rep(mu,length=n)
  // 	if(length(lambda)>1 && length(lambda)!=n) lambda <- rep(lambda,length=n)
  // 	y2 <- rchisq(n,1)
  // 	u <- runif(n)
  // 	r2 <- mu/(2*lambda)*(2*lambda+mu*y2+sqrt(4*lambda*mu*y2+mu^2*y2^2))
  // 	r1 <- mu^2/r2
  // 	ifelse(u < mu/(mu+r1), r1, r2)
  // }
  inline double Random::devInvGauss(double mu,double lam,Generator* gen)
  {
    register double y2tm,r2dm,r1dm;

    if (mu<=0.0 || lam<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
    y2tm=mu*0.5*devGamma(0.5,gen); // Chi^2, dof 1, Gamma(1/2,1/2)
    r2dm=1.0+(y2tm+sqrt(y2tm*(4.0*lam+y2tm)))/(2.0*lam);
    r1dm=1.0/r2dm;

    return mu*((gen->get()<1.0/(1.0+r1dm))?r1dm:r2dm);
  }

  inline double Random::devLogConcave(double (*lfunc) (double),double m,
				      double lfm,double logc,int dtype,
				      Generator* gen)
  {
    register double x,lz,fact;

    if (dtype<0 || dtype>2) throw InvalidParameterException(EXCEPT_MSG(""));
    do {
      x=2.0*gen->get(); lz=log(gen->get());
      if (x>1.0) {
	lz+=(x=log(x-1.0)); x=1.0-x;
      }
      switch (dtype) {
      case 0:
	fact=exp(-lfm-logc);
	break;
      case 1:
	fact=(gen->get()>0.5)?exp(-lfm-logc-M_LN2):(-exp(-lfm-logc-M_LN2));
	break;
      case 2:
	fact=(gen->get()>0.5)?exp(-lfm-logc):(-exp(-lfm-logc));
	break;
      }
      x=x*fact+m;
    } while (lz>lfunc(x)-lfm);

    return x;
  }
//ENDNS

#endif
