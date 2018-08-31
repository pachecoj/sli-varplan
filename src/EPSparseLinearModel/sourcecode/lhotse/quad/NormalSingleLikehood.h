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
 * Desc.:  Header abstract class NormalSingleLikehood
 * ------------------------------------------------------------------- */

#ifndef NORMALSINGLELIKEHOOD_H
#define NORMALSINGLELIKEHOOD_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/quad/default.h"
#include "lhotse/quad/RegressSingleLikehood.h"

//BEGINNS(quad)
  /**
   * Normal likelihood
   *   t(u) = N(y | u,sigma^2)
   * Has single hyperparameter log(sigma^2). Dataset must have
   * regression targets (if there are >1 response interfaces, the first
   * one is used).
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class NormalSingleLikehood : public RegressSingleLikehood
  {
  public:
    // Public methods

    /**
     * Default constructor.
     *
     * @param dsamp  Data sample
     * @param lsinit Init. value for log(sigma^2) (def.: 0)
     */
    NormalSingleLikehood(const Handle<DataIndex>& dsamp,double lsinit=0.0) :
      RegressSingleLikehood(dsamp,lsinit) {}

    double compExpectLog(int ind,double mean,double var,double* lmean=0,
			 double* lvar=0) const;

    double compExpDeriv(int ind,int comp,double mean,double var) const;

  protected:
    double compLogPartInt(double targ,double mean,double var,
			  double* alpha=0,double* nu=0,double frac=1.0) const;
  };

  // Inline methods

  inline double
  NormalSingleLikehood::compLogPartInt(double targ,double mean,double var,
				       double* alpha,double* nu,double frac)
    const
  {
    double sigsq=getVar();
    double s=1.0/(var+sigsq),logpart;

    if (frac!=1.0) throw InvalidParameterException(EXCEPT_MSG("Fractional sites not supported"));
    targ-=mean;
    logpart=-0.5*(targ*targ*s+M_LN2+M_LNPI+log(var+sigsq));
    if (alpha!=0) {
      if (nu==0) throw InvalidParameterException(EXCEPT_MSG("alpha,nu"));
      *alpha=targ*s; *nu=s;
    }

    return logpart;
  }

  inline double
  NormalSingleLikehood::compExpectLog(int ind,double mean,double var,
				      double* lmean,double* lvar) const
  {
    double temp=getTarget(ind)-mean,elog,isigsq=exp(-getLogVar());

    elog=-0.5*(M_LN2+M_LNPI+getLogVar()+(var+temp*temp)*isigsq);
    if (lmean!=0) {
      if (lvar==0) throw InvalidParameterException(EXCEPT_MSG(""));
      *lmean=sqrt(var)*temp*isigsq;
      *lvar=elog-var*isigsq;
    }

    return elog;
  }

  inline double
  NormalSingleLikehood::compExpDeriv(int ind,int comp,double mean,
				     double var) const
  {
    double alpha,nu,temp;

    if (comp!=0) throw OutOfRangeException(EXCEPT_MSG(""));
    compLogPartInt(getTarget(ind),mean,var,&alpha,&nu);
    temp=getTarget(ind)-mean-var*alpha;
    return 0.5*((temp*temp+var*(1.0-nu*var))/getVar()-1.0);
  }
//ENDNS

#endif
