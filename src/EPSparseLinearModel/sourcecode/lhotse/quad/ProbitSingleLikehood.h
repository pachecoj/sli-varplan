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
 * Desc.:  Header abstract class ProbitSingleLikehood
 * ------------------------------------------------------------------- */

#ifndef PROBITSINGLELIKEHOOD_H
#define PROBITSINGLELIKEHOOD_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/quad/default.h"
#include "lhotse/quad/SingleGaussHermQuadrature.h"
#include "lhotse/quad/BinClassSingleLikehood.h"
#include "lhotse/specfun/Specfun.h"

//BEGINNS(quad)
  /**
   * Helper for 'ProbitSingleLikehood::compExpectLog'.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class ProbitSingleQuad : public SingleGaussHermQuadrature
  {
    const BinClassSingleLikehood* likeh;
    bool targ;
    double bias;

  public:
    ProbitSingleQuad(const BinClassSingleLikehood* lh) : likeh(lh),
    targ(true),bias(0.0) {}

    void setInd(int i) {
      if (i<0 || i>=likeh->size())
	throw OutOfRangeException(EXCEPT_MSG(""));
      targ=likeh->getTarget(i); bias=likeh->getBias();
    }

  protected:
    double compLogT(double u) const {
      return Specfun::logCdfNormal(targ?(u+bias):(-u-bias));
    }
  };

  /**
   * Quadrature tools for binary probit noise model
   *   t(u) = P(y | u) = \Phi(y (u+b)).
   * Here, \Phi is the c.d.f. of N(0,1), b is an intercept parameter and
   * y=+1,-1. b is the single hyperparameter.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class ProbitSingleLikehood : public BinClassSingleLikehood
  {
  protected:
    // Members

    Handle<ProbitSingleQuad> quad; // for 'compExpectLog'

  public:
    // Public methods

    /**
     * Default constructor.
     *
     * @param dsamp    Data sample
     * @param binit    Init. value for b (def.: 0)
     * @param isVariat Is 'compExpectLog' required? Def.: false
     */
    ProbitSingleLikehood(const Handle<DataIndex>& dsamp,double binit=0.0,
			 bool isVariat=false) :
      BinClassSingleLikehood(dsamp,binit) {
      if (isVariat)
	quad.changeRep(new ProbitSingleQuad(this));
    }

    double compExpDeriv(int ind,int comp,double mean,double var) const;

    double compExpectLog(int ind,double mean,double var,double* lmean=0,
			 double* lvar=0) const {
      quad->setInd(ind);
      return quad->compExpectLog(mean,var,lmean,lvar);
    }

    void compLogPredProb(const StVector& means,const StVector& vars,
			 StVector& logprob);

  protected:
    // Internal methods

    double compLogPartInt(bool targ,double mean,double var,double* alpha=0,
			  double* nu=0,double frac=1.0) const;
  };

  // Inline methods

  inline double
  ProbitSingleLikehood::compLogPartInt(bool targ,double mean,double var,
				       double* alpha,double* nu,double frac)
    const
  {
    double s=1.0/sqrt(1.0+var),z,logpart,temp;
    double y=targ?1.0:-1.0;

    if (frac!=1.0) throw InvalidParameterException(EXCEPT_MSG("Fractional sites not supported"));
    mean+=getBias();
    z=y*mean*s;
    logpart=Specfun::logCdfNormal(z);
    if (alpha!=0) {
      if (nu==0) throw InvalidParameterException(EXCEPT_MSG("alpha,nu"));
      *alpha=temp=y*s*exp(Specfun::logPdfNormal(z)-logpart);
      *nu=temp*(temp+mean/(1.0+var));
    }

    return logpart;
  }

  /*
   * b is location parameter, can use 'compExpDerivLocation'
   */
  inline double
  ProbitSingleLikehood::compExpDeriv(int ind,int comp,double mean,double var)
    const
  {
    if (comp!=0) throw OutOfRangeException(EXCEPT_MSG(""));

    return compExpDerivLocation(ind,mean,var);
  }
//ENDNS

#endif
