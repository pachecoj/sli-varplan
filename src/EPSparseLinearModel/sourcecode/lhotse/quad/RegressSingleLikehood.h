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
 * Desc.:  Header abstract class RegressSingleLikehood
 * ------------------------------------------------------------------- */

#ifndef REGRESSSINGLELIKEHOOD_H
#define REGRESSSINGLELIKEHOOD_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/quad/default.h"
#include "lhotse/quad/SingleLikehoodFactor.h"
#include "lhotse/data/VLInterDoubleStatic.h"

//BEGINNS(quad)
  /**
   * Abstract subclass of 'SingleLikehoodFactor' for univariate regression
   * noise models t(u). These have a hyperparameter log(sigma^2), where
   * sigma^2 is the noise variance. Dataset must have regression targets (if
   & there are >1 response interfaces, the first one is used)
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class RegressSingleLikehood : public SingleLikehoodFactor
  {
  protected:
    // Members

    StVector hyperPars;     // Contains log(sigma^2) only
    Handle<DataIndex> data; // Data sample

  public:
    // Public methods

    /**
     * Default constructor.
     *
     * @param dsamp  Data sample
     * @param lsinit Init. value for log(sigma^2) (def.: 0)
     */
    RegressSingleLikehood(const Handle<DataIndex>& dsamp,double lsinit=0.0) :
      data(dsamp) {
      // Check dataset
      if (!checkDataset(*dsamp))
	throw WrongDataTypeException(EXCEPT_MSG("Dataset must have regression targets"));
      hyperPars.fill(1,lsinit);
    }

    int size() const {
      return data->size();
    }

    bool checkDataset(const DataIndex& data) const {
      return (DYNCAST(const VLInterDoubleStatic,data.output().p())!=0);
    }

    /**
     * @param i Datapoint index
     * @return  Value of target y
     */
    double getTarget(int i) const {
      if (i<0 || i>=size())
	throw OutOfRangeException(EXCEPT_MSG(""));
      return DYNCAST(const VLInterDoubleStatic,data->output().p())->access(data->map(i));
    }

    /**
     * Computes log of pred. probability for targets 'targs', pairs of mean,
     * variance from vectors 'means', 'vars'.
     *
     * @param targs   S.a.
     * @param means   S.a.
     * @param vars    S.a.
     * @param logprob Log probs. ret. here
     */
    void compLogPredProb(const StVector& targs,const StVector& means,
			 const StVector& vars,StVector& logprob);

    int getNumHyperpars() const {
      return 1;
    }

    const StVector& getHyperpars() const {
      return hyperPars;
    }

    void setHyperpars(const StVector& vec) {
      if (vec.size()!=1)
	throw WrongDimensionException(EXCEPT_MSG(""));
      hyperPars=vec;
    }

    /**
     * @return Value of sigma^2
     */
    double getVar() const {
      return exp(hyperPars[0]);
    }

    /**
     * @return Value of log(sigma^2)
     */
    double getLogVar() const {
      return hyperPars[0];
    }

    double compLogPart(int ind,double mean,double var,double* alpha=0,
		       double* nu=0,double frac=1.0) const
    {
      return compLogPartInt(getTarget(ind),mean,var,alpha,nu,frac);
    }

    bool canPrint() const {
      return true;
    }

    void print(ostream& os) const {
      os << "RegressSingleLikehood: sigma^2=" << getVar() << endl;
    }

  protected:
    // Internal methods

    /**
     * Does job of 'compLogPart', but for fixed target 'targ'.
     */
    virtual double compLogPartInt(double targ,double mean,double var,
				  double* alpha=0,double* nu=0,
				  double frac=1.0) const = 0;
  };

  inline void
  RegressSingleLikehood::compLogPredProb(const StVector& targs,
					 const StVector& means,
					 const StVector& vars,
					 StVector& logprob)
  {
    int i,n=targs.size();

    if (means.size()!=n || vars.size()!=n)
      throw WrongDimensionException(EXCEPT_MSG(""));
    logprob.zeros(n);
    for (i=0; i<n; i++)
      logprob[i]=compLogPartInt(targs[i],means[i],vars[i]);
  }
//ENDNS

#endif
