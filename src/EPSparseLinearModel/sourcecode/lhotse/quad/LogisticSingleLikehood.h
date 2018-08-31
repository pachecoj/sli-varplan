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
 * Desc.:  Header abstract class LogisticSingleLikehood
 * ------------------------------------------------------------------- */

#ifndef LOGISTICSINGLELIKEHOOD_H
#define LOGISTICSINGLELIKEHOOD_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/quad/default.h"
#include "lhotse/quad/SingleGaussHermQuadrature.h"
#include "lhotse/quad/BinClassSingleLikehood.h"

//BEGINNS(quad)
  /**
   * Helper for 'LogisticSingleLikehood', to make use of Gaussian
   * quadrature.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class LogisticSingleQuad : public SingleGaussHermQuadrature
  {
    const BinClassSingleLikehood* likeh;
    bool targ;
    double bias;

  public:
    LogisticSingleQuad(const BinClassSingleLikehood* lh) : likeh(lh),
    targ(true),bias(0.0) {}

    void setTarg(bool t) {
      targ=t; bias=likeh->getBias();
    }

    void setInd(int i) {
      if (i<0 || i>=likeh->size())
	throw OutOfRangeException(EXCEPT_MSG(""));
      setTarg(likeh->getTarget(i));
    }

  protected:
    double compLogT(double u) const {
      double arg=targ?(u+bias):(-u-bias);
      if (arg>=0.0)
	return -log(1.0+exp(-arg));
      else
	return arg-log(1.0+exp(arg));
    }
  };

  /**
   * Quadrature tools for binary logistic noise model
   *   t(u) = P(y | u) = \sigma(y (u+b)), \sigma(a) = [1 + exp(-a)]^(-1)
   * Here, b is an intercept parameter and y=+1,-1. b is the single
   * hyperparameter.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class LogisticSingleLikehood : public BinClassSingleLikehood
  {
  protected:
    // Members

    Handle<LogisticSingleQuad> quad; // for quadrature

  public:
    // Public methods

    /**
     * Default constructor.
     *
     * @param dsamp    Data sample
     * @param binit    Init. value for b (def.: 0)
     */
    LogisticSingleLikehood(const Handle<DataIndex>& dsamp,double binit=0.0) :
      BinClassSingleLikehood(dsamp,binit) {
      quad.changeRep(new LogisticSingleQuad(this));
    }

    /*
     * b is location parameter, can use 'compExpDerivLocation'
     */
    double compExpDeriv(int ind,int comp,double mean,double var) const {
      if (comp!=0) throw OutOfRangeException(EXCEPT_MSG(""));
      return compExpDerivLocation(ind,mean,var);
    }

    double compExpectLog(int ind,double mean,double var,double* lmean=0,
			 double* lvar=0) const {
      quad->setInd(ind);
      return quad->compExpectLog(mean,var,lmean,lvar);
    }

  protected:
    // Internal methods

    double compLogPartInt(bool targ,double mean,double var,double* alpha=0,
			  double* nu=0,double frac=1.0) const {
      if (frac!=1.0) throw InvalidParameterException(EXCEPT_MSG("Fractional sites not supported"));
      quad->setTarg(targ);
      return quad->compLogPart(mean,var,alpha,nu);
    }
  };
//ENDNS

#endif
