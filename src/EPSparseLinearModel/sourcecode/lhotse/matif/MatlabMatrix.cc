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
 * Module: matif
 * Desc.:  Definition class MatlabMatrix
 * ------------------------------------------------------------------- */

#include "lhotse/matif/MatlabMatrix.h"
#include "lhotse/matrix/BaseVector.h"
#include "lhotse/matrix/BaseMatrix.h"
#include "lhotse/matrix/BaseLinVec.h"
#include "lhotse/matrix/BaseLinMat.h"
#include "lhotse/matrix/ArrayUtils.h"
#include "lhotse/matrix/StVector.h"

//BEGINNS(matif)
  // Static members

  bool MatlabMatrix::persFlag(false);

  // Public methods

  MatlabMatrix::MatlabMatrix(const BaseVector<int>& vec,bool oown,bool incr,
			     bool col) : isOwn(oown)
  {
    int n=vec.size();

    mat=col?mxCreateDoubleMatrix(n,1,mxREAL) :
      mxCreateDoubleMatrix(1,n,mxREAL);
    StVector msk;
    msk.reassign(mxGetPr(mat),n);
    msk.convert(vec);
    if (incr) msk.addscal(1.0);
    if (persFlag) mexMakeArrayPersistent(mat); // persistent
    persFlag=false; // reset
  }

  MatlabMatrix::MatlabMatrix(const BaseMatrix<double>& src,bool oown) :
    isOwn(oown)
  {
    mat=mxCreateDoubleMatrix(src.rows(),src.cols(),mxREAL);
    if (persFlag) mexMakeArrayPersistent(mat); // persistent
    BaseLinMat<double> linm;
    src.getLinMat(linm,false);
    copyContent(linm);
    persFlag=false; // reset
  }

  MatlabMatrix::MatlabMatrix(const BaseLinMat<double>& src,bool oown) :
    isOwn(oown)
  {
    //if (WriteBackMat<double>::debugFlag) // DEBUG!
    //  mexPrintf("MatlabMatrix: m=%d,n=%d\n",src.m,src.n);
    mat=mxCreateDoubleMatrix(src.m,src.n,mxREAL);
    if (persFlag) mexMakeArrayPersistent(mat); // persistent
    copyContent(src);
    persFlag=false; // reset
  }

  MatlabMatrix::MatlabMatrix(const BaseLinVec<double>& src,bool col,
			     bool oown) : isOwn(oown)
  {
    mat=mxCreateDoubleMatrix(col?src.n:1,col?1:src.n,mxREAL);
    if (persFlag) mexMakeArrayPersistent(mat); // persistent
    src.copy(mxGetPr(mat));
    persFlag=false; // reset
  }

  // Internal methods

  void MatlabMatrix::copyContent(const BaseLinMat<double>& linm)
  {
    if (linm.stride==linm.m) {
      // Flat buffer
      memmove((void*) mxGetPr(mat),(void*) linm.buff,linm.m*linm.n*
	      sizeof(double));
    } else {
      // Strided buffer
      double* trgP=mxGetPr(mat);
      const double* srcP=linm.buff;
      for (int i=0; i<linm.n; i++,srcP+=linm.stride,trgP+=linm.m)
	memmove((void*) trgP,(void*) srcP,linm.m*sizeof(double));
    }
  }
//ENDNS
