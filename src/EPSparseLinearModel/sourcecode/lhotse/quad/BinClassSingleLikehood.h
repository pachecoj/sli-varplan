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
 * Desc.:  Header abstract class BinClassSingleLikehood
 * ------------------------------------------------------------------- */

#ifndef BINCLASSSINGLELIKEHOOD_H
#define BINCLASSSINGLELIKEHOOD_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/quad/default.h"
#include "lhotse/quad/SingleLikehoodFactor.h"
#include "lhotse/data/VLInterCategoryStatic.h"

//BEGINNS(quad)
  /**
   * Abstract subclass of 'SingleLikehoodFactor' for binary classification
   * noise models t(u). These models have an intercept parameter b which is
   * added to u. The hyperpar. vector consists of b. The dataset must have
   * binary targets, which are translated 0 -> -1, 1 -> +1.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class BinClassSingleLikehood : public virtual SingleLikehoodFactor
  {
  protected:
    // Members

    StVector hyperPars;     // Contains b only
    Handle<DataIndex> data; // Data sample

  public:
    /**
     * Default constructor
     *
     * @param dsamp Data sample
     * @param binit Init. value for b (def.: 0)
     */
    BinClassSingleLikehood(const Handle<DataIndex>& dsamp,double binit=0.0) :
      data(dsamp) {
      // Check dataset
      if (!checkDataset(*dsamp))
	throw WrongDataTypeException(EXCEPT_MSG("Dataset must have binary classification targets"));
      hyperPars.fill(1,binit);
    }

    int size() const {
      return data->size();
    }

    bool checkDataset(const DataIndex& data) const {
      const VLInterCategoryStatic* cif=
	DYNCAST(const VLInterCategoryStatic,data.output().p());
      return (cif!=0 &&
	      DYNCAST(const CategoryType,cif->getType())->getNumC()==2);
    }

    /**
     * @param i Datapoint index
     * @return  Is corresponding target == +1?
     */
    bool getTarget(int i) const {
      if (i<0 || i>=size())
	throw OutOfRangeException(EXCEPT_MSG(""));
      return (DYNCAST(const VLInterCategoryStatic,data->output().p())->access(data->map(i))==1);
    }

    const Handle<DataIndex>& getData() const {
      return data;
    }

    /**
     * Computes log of pred. probability for class +1 for pairs of mean,
     * variance from vectors 'means', 'vars'.
     *
     * @param means   S.a.
     * @param vars    S.a.
     * @param logprob Log probs. ret. here
     */
    void compLogPredProb(const StVector& means,const StVector& vars,
			 StVector& logprob);

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

    double getBias() const {
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
      os << "BinClassSingleLikehood: bias=" << getBias() << endl;
    }

  protected:
    // Internal methods

    /**
     * Does job of 'compLogPart', but for fixed target 'targ' (+1 is true,
     * -1 is false).
     */
    virtual double compLogPartInt(bool targ,double mean,double var,
				  double* alpha=0,double* nu=0,
				  double frac=1.0) const = 0;
  };

  inline void
  BinClassSingleLikehood::compLogPredProb(const StVector& means,
					  const StVector& vars,
					  StVector& logprob)
  {
    int i,n=means.size();

    if (vars.size()!=n)
      throw WrongDimensionException(EXCEPT_MSG(""));
    logprob.zeros(n);
    for (i=0; i<n; i++)
      logprob[i]=compLogPartInt(true,means[i],vars[i]);
  }
//ENDNS

#endif
