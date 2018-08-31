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
 * Desc.:  Header class Specfun
 * ------------------------------------------------------------------- */

#ifndef SPECFUN_H
#define SPECFUN_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

/*
 * TODO:
 * - Check HAVE_LIBGSL. But this whole class does not compile without GSL!
 * - Kicked out NR code. Some methods not req. in the moment are not
 *   implemented. Use GSL functions to implement them!
 * - Some methods not properly debugged or robust to underflow, etc.!
 */

#include "lhotse/specfun/default.h"
#include <gsl/gsl_sf.h>     // GSL (special functions)
#include <gsl/gsl_sf_psi.h> // GSL (special functions)
#include <gsl/gsl_math.h>   // GSL (general)
#include <gsl/gsl_errno.h>  // GSL (errors)

//BEGINNS(specfun)

#define GSL_EXCEPT(stat) if (stat!=0) { ArrayHandle<char> msg(strlen(gsl_strerror(stat))+6); sprintf(msg.p(),"GSL: %s",gsl_strerror(stat)); throw NumericalException(EXCEPT_MSG(msg.p())); }
#define GSL_EXCEPT_1ARG(stat,x) if (stat!=0) { ArrayHandle<char> msg(strlen(gsl_strerror(stat))+40); sprintf(msg.p(),"GSL: %s. Argument=%f",gsl_strerror(stat),x); throw NumericalException(EXCEPT_MSG(msg.p())); }

#define GSL_INIT if (!gslInit) { gsl_set_error_handler_off(); gslInit=true; }

#define CDFNORMAL_THRES -6.0 // See 'cdfNormal', ...

  /**
   * Provides static methods for special univariate functions, such as
   * pdf's, cdf's of standard distributions.
   * <p>
   * NOTE: This class requires GSL, it will not compile if GSL is not
   * present!
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class Specfun
  {
  protected:
    // Static members

    static bool gslInit; // GSL already initialized?

  public:
    // Public static methods

    /**
     * Computes natural log of Gamma(z) for z>0. Note that if z is a
     * natural number, then z! = Gamma(z+1).
     * <p>
     * Uses GSL.
     */
    static double logGamma(double z) {
      if (z<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
      //if (z<1e-6) cout << "z=" << z << endl; // DEBUG!!
      GSL_INIT;
      gsl_sf_result res;
      int stat=gsl_sf_lngamma_e(z,&res);
      GSL_EXCEPT_1ARG(stat,z);
      return res.val;
    }

    /**
     * Computes Beta function Beta(z,w) for z,w>0.
     */
    static double beta(double z,double w) {
      return exp(logGamma(z)+logGamma(w)-logGamma(z+w));
    }

    /**
     * Computes incomplete Gamma function gamma(a,x) (cdf of Gamma(a,1)
     * variate), x>=0,a>0.
     * ATTENTION: Not implemented!
     */
    static double incGamma(double a,double x) {
      throw NotImplemException(EXCEPT_MSG("IMPLEMENT USING GSL!"));
    }

    /**
     * Computes the digamma function Psi(z), which is the derivative of
     * the log Gamma function, for z>0.
     */
    static double digamma(double z) {
      if (z<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
      GSL_INIT;
      gsl_sf_result res;
      int stat=gsl_sf_psi_e(z,&res);
      GSL_EXCEPT_1ARG(stat,z);
      return res.val;
    }

    /**
     * Computes the trigamma function Psi'(z), the derivative of the
     * digamma function, so the second deriv. of log Gamma, for 'z'>0.
     */
    static double trigamma(double z) {
      if (z<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
      GSL_INIT;
      gsl_sf_result res;
      int stat=gsl_sf_psi_1_e(z,&res);
      GSL_EXCEPT_1ARG(stat,z);
      return res.val;
    }

    /**
     * DOWNW. COMPAT.
     */
    static double derivDigamma(double z) {
      return trigamma(z);
    }

    static double cdfGamma(double a,double b,double x) {
      if (b<=0) throw InvalidParameterException(EXCEPT_MSG(""));
      return (x>0.0) ? incGamma(a,x/b) : 0.0;
    }

    static double pdfGamma(double a,double b,double x) {
      static double olda=-1.0,oldb=-1.0;
      static double abterm,am1;

      if (b<=0 || a<=0) throw InvalidParameterException(EXCEPT_MSG(""));
      if (x<=0.0) return 0.0;
      if (a!=olda || b!=oldb) {
	olda=a; oldb=b;
	abterm=logGamma(a)+a*log(b); am1=a-1.0;
      }
      return exp(am1*log(x) - abterm - x/b);
    }

    static double cdfChiSq(int degFree,double x) {
      if (degFree<=0) throw InvalidParameterException(EXCEPT_MSG(""));
      return (x>0.0) ? incGamma(0.5*(double) degFree,0.5*x) : 0.0;
    }

    /*
     * Methods related to cdf. of Gaussian ('cdfNormal', 'logCdfNormal',
     * 'derivLogCdfNormal'):
     * We make use of the corr. GSL function if x >= CDFNORMAL_THRES.
     * Otherwise, we use the approximation (coming from the asymptotic
     * expansion):
     *   Phi(x) \approx N(-x)/(-x) (1 - y (1 - 3 y (1 - 5 y (1 - 7 y)))),
     *   y = x^{-2},
     * where N(x) is 'pdfNormal'.
     */

    /**
     * Uses GSL.
     *
     * @param x Arg.
     * @return CDF of Gaussian N(0,1) at x
     */
    static double cdfNormal(double x) {
      if (x>=CDFNORMAL_THRES) {
	GSL_INIT;
	gsl_sf_result res;
	int stat=gsl_sf_erf_Q_e(x,&res);
	GSL_EXCEPT_1ARG(stat,x);
	return 1.0-res.val;
      } else {
	double y=1.0/x/x;
	return pdfNormal(x)/(-x)*
	  (1.0-y*(1.0-3.0*y*(1.0-5.0*y*(1.0-7.0*y))));
      }
    }

    static double cdfNormal(double mean,double stddev,double x) {
      if (stddev<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
      return cdfNormal((x-mean)/stddev);
    }

    /**
     * Uses GSL.
     * NOTE: Throws exception for underflow!
     *
     * @param x Arg.
     * @return PDF of Gaussian N(0,1) at x
     */
    static double pdfNormal_GSL(double x) {
      GSL_INIT;
      gsl_sf_result res;
      int stat=gsl_sf_erf_Z_e(x,&res);
      GSL_EXCEPT_1ARG(stat,x);
      return res.val;
    }

    static double pdfNormal(double x) {
      return exp(-0.5*(M_LN2+M_LNPI+x*x));
    }

    static double pdfNormal(double mean,double stddev,double x) {
      if (stddev<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
      return pdfNormal((x-mean)/stddev)/stddev;
    }

    static double logPdfNormal(double x) {
      return -0.5*(M_LN2+M_LNPI+x*x);
    }

    static double cdfExp(double x) {
      return (x>0.0) ? 1.0-exp(-x) : 0.0;
    }

    static double cdfExp(double mean,double x) {
      if (mean<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
      return (x>0.0) ? 1.0-exp(-x/mean) : 0.0;
    }

    static double pdfExp(double x) {
      return (x>0.0) ? exp(-x) : 0.0;
    }

    static double pdfExp(double mean,double x) {
      if (mean<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
      return (x>0.0) ? exp(-x/mean)/mean : 0.0;
    }

    /**
     * Computes cdf of Poisson variate with rate gamma
     */
    static double cdfPoisson(double gamma,int x) {
      if (x>0) return 1.0-incGamma((double) x,gamma);
      else if (x==0) return exp(-gamma);
      else return 0.0;
    }

    static double pdfPoisson(double gamma,int x) {
      if (gamma<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
      double dx=(double) x;
      return (x>=0) ? exp(-gamma + dx*log(gamma) - logGamma(dx+1.0)) : 0.0;
    }

    /**
     * Computes incomplete Beta function (cdf of Beta(a,b) variate) for
     * a,b>0,0<=x<=1.
     * <p>
     * ATTENTION: Not implemented!
     */
    static double cdfBeta(double a,double b,double x) {
      throw NotImplemException(EXCEPT_MSG("IMPLEMENT USING GSL!"));
    }

    static double pdfBeta(double a,double b,double x) {
      static double olda,oldb,atemp;

      if (a<=0.0 || b<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
      if (x<=0.0 || x>=1.0) return 0.0;
      if (olda!=a || oldb!=b) {
	// Setup
	atemp=logGamma(a)+logGamma(b)-logGamma(a+b);
	olda=a; oldb=b;
      }
      return ((a-1.0)*log(x)+(b-1.0)*log(1.0-x)-atemp);
    }

    /**
     * Computes cdf of t variate with alpha degrees of freedom, alpha>0
     */
    static double cdfT(double alpha,double x) {
      if (alpha<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
      if (x>=0.0)
	return 1.0 - 0.5*cdfBeta(alpha*0.5,0.5,alpha/(alpha+x*x));
      else
	return 0.5*cdfBeta(alpha*0.5,0.5,alpha/(alpha+x*x));
    }

    static double pdfT(double alpha,double x) {
      static double oldalpha,atemp1,atemp2;

      if (alpha<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
      if (alpha!=oldalpha) {
	// Setup
	oldalpha=alpha;
	atemp2=0.5*(alpha+1.0);
	atemp1=logGamma(atemp2)-logGamma(0.5*alpha)-
	  0.5*(log(alpha)+M_LNPI);
      }
      return exp(atemp1-atemp2*log(1.0+x*x/alpha));
    }

    /**
     * Computes p value of double-sided t test with alpha degrees of
     * freedom, i.e. the probability of |t| > x, t t-variate with alpha
     * degrees of freedom
     */
    static double pDoubleSidedT(double alpha,double x) {
      if (alpha<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
      return cdfBeta(alpha*0.5,0.5,alpha/(alpha+x*x));
    }

    /**
     * Computes cdf of Binomial variate with parameters n,p
     */
    static double cdfBinomial(int n,double p,int x) {
      if (n<=0 || p<0.0 || p>1.0) throw InvalidParameterException(EXCEPT_MSG(""));
      if (x>=n) return 1.0;
      else if (x>0) return cdfBeta((double) x,(double) (n-x+1),p);
      else if (x==0) return (p!=1.0) ? exp(((double) n)*log(1-p)) : 0.0;
      else return 0.0;
    }

    /**
     * If Phi(z) denotes the c.d.f. of N(0,1), this method computes
     * log Phi(z).
     * Uses GSL.
     *
     * @param z Argument
     * @return  log Phi(z)
     */
    static double logCdfNormal(double z) {
      if (z>=CDFNORMAL_THRES) {
	// Uses GSL to compute log(erfc(x)), then
	// log(Phi(z)) = log(erfc(-z/sqrt(2))) - log(2)
	GSL_INIT;
	gsl_sf_result res;
	double x=-z/M_SQRT2;
	int stat=gsl_sf_log_erfc_e(x,&res);
	GSL_EXCEPT_1ARG(stat,x);
	return res.val-M_LN2;
      } else {
	// Uses tail approximation of Phi(z)
	double y=1.0/z/z;
	double tail=-y*(1.0-3.0*y*(1.0-5.0*y*(1.0-7.0*y)));
	return logPdfNormal(z)-log(-z)+log1p(tail);
      }
    }

    /**
     * If Phi(z) denotes the c.d.f. of N(0,1), this method computes
     *   (d/dz) log Phi(z) = N(z)/Phi(z).
     * If this is called f(z), the function f(-z) is called hazard function
     * in statistics.
     * Uses GSL.
     *
     * @param z Argument
     * @return  (d/dz) log Phi(z)
     */
    static double derivLogCdfNormal(double z) {
      if (z>=CDFNORMAL_THRES) {
	// f(z) = h(-z), h hazard function (GSL)
	GSL_INIT;
	gsl_sf_result res;
	int stat=gsl_sf_hazard_e(-z,&res);
	GSL_EXCEPT_1ARG(stat,-z);
	return res.val;
      } else {
	// Uses tail approx. of Phi(z)
	double y=1.0/z/z;
	return -z/(1.0-y*(1.0-3.0*y*(1.0-5.0*y*(1.0-7.0*y))));
      }
    }
  };
#undef GSL_EXCEPT
#undef GSL_EXCEPT_1ARG
#undef GSL_INIT
#undef CDFNORMAL_THRES
//ENDNS

#endif
