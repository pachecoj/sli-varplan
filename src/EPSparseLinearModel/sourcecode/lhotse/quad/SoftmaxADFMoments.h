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
 * Desc.:  Header class SoftmaxADFMoments
 * ------------------------------------------------------------------- */

#ifndef SOFTMAXADFMOMENTS_H
#define SOFTMAXADFMOMENTS_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/quad/default.h"
#include "lhotse/quad/MultiClassQuadrature.h"
#include "lhotse/matrix/StVector.h"

//BEGINNS(quad)
  /**
   * ADF moment projections for softmax multi-class classification
   * noise model, inherits from 'MultiClassQuadrature'.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class SoftmaxADFMoments : public virtual MultiClassQuadrature
  {
  protected:
    // Members

    int cls;       // Current target (0,...,'getDim()'-1)
    StVector bias; // Intercept parameters
    double gamma;  // current re-weight. factor (1.0 -> none)
    mutable StVector vvec;
    //StMatrix rmat; // DEBUG!!
    //int debugInd; // DEBUG!

  public:
    // Public methods

    /**
     * Constructor. 'cls' is set to 0, 'bias' init. with 0.
     *
     * @param dim Number of classes
     */
    SoftmaxADFMoments(int dim) : cls(0),gamma(1.0) {
      if (dim<1) throw InvalidParameterException(EXCEPT_MSG(""));
      bias.zeros(dim);
    }

    void setCls(int cl) {
      if (cl<0 || cl>=getDim())
	throw InvalidParameterException("Invalid class target 'cl'");
      if (cl!=cls) resetTransfer(); // reset
      cls=cl;
    }

    int getCls() const {
      return cls;
    }

    void setBias(const StVector& bvec) {
      if (bvec.size()!=getDim())
	throw WrongDimensionException("'bvec' has wrong size");
      bias=bvec;
      resetTransfer(); // reset
    }

    const StVector& getBias() const {
      return bias;
    }

    void resetBias() {
      bias.zeros(getDim());
      resetTransfer(); // reset
    }

    /**
     * @param val New re-weighting factor in (0,1]
     */
    virtual void setGamma(double val=1.0) {
      if (val<=0.0 || val>1.0)
	throw InvalidParameterException(EXCEPT_MSG(""));
      if (val!=gamma) {
	gamma=val;
	resetTransfer();
      }
    }

    /* DEBUG:
    void setRMat(const StMatrix& r) {
      rmat.reassign(r);
      }*/

    /* DEBUG:
    void setInd(int i) {
      if (i<0 || i>=rmat.cols())
	throw InvalidParameterException(EXCEPT_MSG(""));
      debugInd=i;
      }*/

  protected:
    // Internal methods

    /**
     * The function is
     *   v_{cls} - log 1^T exp(v), v = gamma*u + bias
     * (bias==0 if 'bias' not given). Log-sum-exp computed in stable
     * way.
     *
     * @param u Eval. point
     * @return  Function value f(u)
     */
    double compF(const StVector& u) const {
      const StVector* vP=&u;
      if (gamma!=1.0) {
	vvec.prod(u,gamma);
	vP=&vvec;
      }
      vvec.addprod(*vP,1.0,bias);
      vP=&vvec;
      return (*vP)[cls]-vP->stableLogSum();
    }

    /* DEBUG:
    double compF(const StVector& u) const {
      double temp=rmat(RangeFull::get(),debugInd)->inner(u);
      return -temp*temp;
      }*/

    /**
     * The function is
     *   v - (log 1^T exp(v))*1, v = gamma*u + bias
     * (bias==0 if 'bias' not given).
     *
     * @param u Eval. point
     * @return  Function value f(u)
     */
    void compFMulti(const StVector& u,StVector& val) const {
      if (gamma==1.0)
	val=u;
      else
	val.prod(u,gamma);
      val.addprod(1.0,bias);
      val.addscal(-val.stableLogSum());
    }

    int getFNum() const {
      return getDim();
    }
  };
//ENDNS

#endif
