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
 * Desc.:  Header class LaplaceSingleLikehood
 * ------------------------------------------------------------------- */

#ifndef LAPLACESINGLELIKEHOOD_H
#define LAPLACESINGLELIKEHOOD_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/quad/default.h"
#include "lhotse/quad/SingleLikehoodFactor.h"
#include "lhotse/specfun/Specfun.h"

//BEGINNS(quad)

#define TAIL_APPROX_THRES 5.0 // See 'compLogPart'

  /**
   * Represents the Laplace factor
   *   t_i(u) = exp(-tau_i |u|).
   * Sites are not assoc. with points of a dataset. The tau_i parameters can
   * be different or all shared. In the latter case, 'ltau' contains the
   * single shared entry only.
   * <p>
   * NOTE: The normalization constant tau_i/2 of the Laplace distrib. is
   * ommitted here in general!
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class LaplaceSingleLikehood : public SingleLikehoodFactor
  {
  protected:
    // Members

    int n;         // Number of sites
    StVector ltau; // Parameters log(tau_i)

  public:
    int debugNo;

  public:
    // Public methods

    /**
     * Default constructor. All tau_i=1 initially.
     *
     * @param pn Size of "dataset" (number of sites)
     */
    LaplaceSingleLikehood(int pn) : n(pn),debugNo(-1) {
      if (pn<1) throw InvalidParameterException(EXCEPT_MSG(""));
      ltau.zeros(1);
    }

    int size() const {
      return n;
    }

    /**
     * Returns false always.
     */
    bool checkDataset(const DataIndex& data) const {
      return false;
    }

    bool suppFracSites() const {
      return true;
    }

    /**
     * The 'ind' argument is used only to select the tau_i entry.
     */
    double compLogPart(int ind,double mean,double var,double* alpha=0,
		       double* nu=0,double frac=1.0) const {
      if (ind<0 || ind>=n) throw OutOfRangeException(EXCEPT_MSG("ind"));
      return compLogPartStatic(mean,var,exp(ltau[(ltau.size()>1)?ind:0]),
			       alpha,nu,frac,0,debugNo==ind);
    }

    /**
     * Static variant, where tau is passed explicitly. If 'exabs' is used,
     * E[-|u|] is returned in '*exabs', E[.] over the tilted distribution.
     * Can be used only if 'alpha', 'nu' are used as well.
     * NOTE: '*exabs' correct only if 'frac'==1!
     */
    static double compLogPartStatic(double mean,double var,double tau=1.0,
				    double* alpha=0,double* nu=0,
				    double frac=1.0,double* exabs=0,
				    bool doDebug=false);

    /**
     * Variant of 'compLogPart' which returns E[-|u|] in 'exAbs', where
     * E[.] is over the tilted distribution \hat{P}_i.
     */
    double compLogPartVar(int ind,double mean,double var,double& exAbs,
			  double frac=1.0) const;

    /**
     * Sets the log(tau_i) entries. If 'pltau' is a vector, it must be of
     * size 'n'. If it is a scalar, all tau_i are shared.
     *
     * @param pltau S.a.
     */
    void setLogTau(const StVector& pltau) {
      if (pltau.size()!=n && pltau.size()!=1)
	throw WrongDimensionException(EXCEPT_MSG(""));
      ltau=pltau;
    }

    void setLogTau(double pltau) {
      ltau.fill(1,pltau);
    }

    const StVector& getLogTau() const {
      return ltau;
    }
  };

  /**
   * For DEBUG!!! Gaussian with variance 1/(2 tau).
   */
  class DebugLaplaceLikehood : public SingleLikehoodFactor
  {
  protected:
    // Members

    int n;         // Number of sites
    StVector ltau; // Parameters log(tau_i)

  public:
    // Public methods

    /**
     * Default constructor. All tau_i=1 initially.
     *
     * @param pn Size of "dataset" (number of sites)
     */
    DebugLaplaceLikehood(int pn) : n(pn) {
      if (pn<1) throw InvalidParameterException(EXCEPT_MSG(""));
      ltau.zeros(1);
    }

    int size() const {
      return n;
    }

    /**
     * Returns false always.
     */
    bool checkDataset(const DataIndex& data) const {
      return false;
    }

    bool suppFracSites() const {
      return true;
    }

    /**
     * The 'ind' argument is used only to select the tau_i entry.
     */
    double compLogPart(int ind,double mean,double var,double* alpha=0,
		       double* nu=0,double frac=1.0) const;

    /**
     * Sets the log(tau_i) entries. If 'pltau' is a vector, it must be of
     * size 'n'. If it is a scalar, all tau_i are shared.
     *
     * @param pltau S.a.
     */
    void setLogTau(const StVector& pltau) {
      if (pltau.size()!=n && pltau.size()!=1)
	throw WrongDimensionException(EXCEPT_MSG(""));
      ltau=pltau;
    }

    void setLogTau(double pltau) {
      ltau.fill(1,pltau);
    }

    const StVector& getLogTau() const {
      return ltau;
    }
  };

  // Inline methods

  /*
   * We use the tail approximation
   *   1 - Phi(x) \approx N(x)/x (1 - u (1 - 3 u (1 - 5 u (1 - 7 u)))),
   *   u =1/(x^2),
   * where Phi, N are cdf., pdf. of standard Gaussian, whenever the
   * argument x >= TAIL_APPROX_THRES. Apart from that, we use 'log1p'
   * to compute log(1+x) for potentially small x.
   */
  inline double
  LaplaceSingleLikehood::compLogPartStatic(double mean,double var,double tau,
					   double* alpha,double* nu,
					   double frac,double* exabs,
					   bool doDebug)
  {
    double temp,a,sqa,h,absh,argm,argp,logi0,ratio,part1,temp2,tailm,tailp;

    if (frac<=0.0) throw InvalidParameterException(EXCEPT_MSG("frac"));
    if (tau<=0.0) throw InvalidParameterException(EXCEPT_MSG("tau"));
    tau*=frac;
    h=mean*tau; a=var*tau*tau; sqa=tau*sqrt(var);
    absh=fabs(h); argm=sqa-absh/sqa; argp=sqa+absh/sqa; // x_+, x_-
    if (doDebug)
      cout << "-----" << endl << "h=" << h << ",a=" << a << ",|h|=" << absh
	   << ",argm=" << argm << ",argp=" << argp << endl;

    // 'part1' = log(1-Phi(x_-)) = F(x_-)
    if (argm>=TAIL_APPROX_THRES) {
      // Use tail approximation
      temp=1.0/argm/argm;
      // g(x_-) = 1 + 'tailm'
      tailm=-temp*(1.0-3.0*temp*(1.0-5.0*temp*(1.0-7.0*temp)));
      part1=-0.5*argm*argm+(temp=-log(argm)+log1p(tailm)-0.5*(M_LN2+M_LNPI));
      logi0=temp-0.5*mean*mean/var; // Dominant part of log(I_0)
      if (doDebug)
	cout << "part1=" << part1 << "[TA], tailm=" << tailm << endl;
    } else {
      part1=Specfun::logCdfNormal(-argm);
      logi0=0.5*a-absh+part1; // Dominant part of log(I_0)
      if (doDebug)
	cout << "part1=" << part1 << "[EX]" << endl;
    }
    if (argp>=TAIL_APPROX_THRES) {
      temp=1.0/argp/argp;
      // g(x_+) = 1 + 'tailp'
      tailp=-temp*(1.0-3.0*temp*(1.0-5.0*temp*(1.0-7.0*temp)));
      if (doDebug)
	cout << "tailp=" << tailp << "[TA]" << endl;
    }
    if (argm>=TAIL_APPROX_THRES)
      // Use tail approximation for both parts of ratio
      ratio=(a-absh)*(1.0+tailp)/(a+absh)/(1.0+tailm);
    else if (argp>=TAIL_APPROX_THRES)
      // Use tail approximation for 'argp' part of ratio
      ratio=exp(-0.5*mean*mean/var+absh-0.5*a-log(argp)+log1p(tailp)-
		0.5*(M_LN2+M_LNPI)-part1);
    else
      // No tail approximations
      ratio=exp(2.0*absh+Specfun::logCdfNormal(-argp)-part1);
    if (doDebug)
      cout << "ratio=" << ratio << endl;
    logi0+=log1p(ratio); // log(I_0) complete
    if (alpha!=0) {
      if (nu==0) throw InvalidParameterException(EXCEPT_MSG(""));
      temp=GSL_SIGN(h)*(1.0-2.0/(1.0+ratio));
      *alpha=temp*tau;
      temp2=exp(-0.5*mean*mean/var-logi0)*M_SQRT2/M_SQRTPI;
      *nu=(temp*temp-1.0+temp2/sqa)*tau*tau;
      if (exabs!=0)
	*exabs=tau*var+temp*mean-temp2*sqrt(var);
    }

    return logi0;
  }

  inline double
  LaplaceSingleLikehood::compLogPartVar(int ind,double mean,double var,
					double &exAbs,double frac)
    const
  {
    double temp,tau,a,sqa,h,absh,argm,argp,logi0,ratio,part1,temp2,tailm,
      tailp,alpha;

    if (ind<0 || ind>=n) throw OutOfRangeException(EXCEPT_MSG(""));
    if (frac<=0.0) throw InvalidParameterException(EXCEPT_MSG("frac"));
    tau=exp(ltau[(ltau.size()>1)?ind:0])*frac;
    h=mean*tau; a=var*tau*tau; sqa=tau*sqrt(var);
    absh=fabs(h); argm=sqa-absh/sqa; argp=sqa+absh/sqa; // x_+, x_-

    // 'part1' = log(1-Phi(x_-)) = F(x_-)
    if (argm>=TAIL_APPROX_THRES) {
      // Use tail approximation
      temp=1.0/argm/argm;
      // g(x_-) = 1 + 'tailm'
      tailm=-temp*(1.0-3.0*temp*(1.0-5.0*temp*(1.0-7.0*temp)));
      part1=-0.5*argm*argm-log(argm)+log1p(tailm)-0.5*(M_LN2+M_LNPI);
    } else {
      part1=log1p(-Specfun::cdfNormal(argm));
    }
    logi0=0.5*a-absh+part1; // Dominant part of log(I_0)
    if (argp>=TAIL_APPROX_THRES) {
      temp=1.0/argp/argp;
      // g(x_+) = 1 + 'tailp'
      tailp=-temp*(1.0-3.0*temp*(1.0-5.0*temp*(1.0-7.0*temp)));
    }
    if (argm>=TAIL_APPROX_THRES) {
      // Use tail approximation for both parts of ratio
      ratio=((a-absh)*(1.0+tailp))/((a+absh)*(1.0+tailm));
    } else if (argp>=TAIL_APPROX_THRES) {
      // Use tail approximation for 'argp' part of ratio
      temp=absh/sqa;
      ratio=exp(-0.5*temp*temp+absh-0.5*a-log(argp)+log1p(tailp)-
		0.5*(M_LN2+M_LNPI)-part1);
    } else
      // No tail approximations
      ratio=exp(2.0*absh+log1p(-Specfun::cdfNormal(argp))-part1);
    logi0+=log1p(ratio); // log(I_0) complete
    alpha=GSL_SIGN(h)*(1.0-2.0/(1.0+ratio));
    temp=absh/sqa;
    exAbs=tau*var+alpha*mean-
      exp(-0.5*temp*temp-logi0)*sqrt(var)*M_SQRT2/M_SQRTPI;

    return logi0;
  }

  inline double
  DebugLaplaceLikehood::compLogPart(int ind,double mean,double var,
				    double* alpha,double* nu,double frac)
    const
  {
    double temp,tau,a,sqa,h,absh,argm,argp,logi0,ratio,part1,temp2,tailm,
      tailp;

    if (ind<0 || ind>=n) throw OutOfRangeException(EXCEPT_MSG(""));
    if (frac<=0.0) throw InvalidParameterException(EXCEPT_MSG("frac"));
    tau=exp(ltau[(ltau.size()>1)?ind:0])*frac;

    temp=2.0*tau*var;
    logi0=-0.5*log1p(temp)-tau*mean*mean/(1.0+temp);
    if (alpha!=0) {
      *alpha=-2.0*tau*mean/(1.0+temp);
      *nu=2.0*tau/(1.0+temp);
    }

    return logi0;
  }
//ENDNS

#endif
