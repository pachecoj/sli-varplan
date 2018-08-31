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
 * Module: matrix
 * Desc.:  Definition of class StVector
 * ------------------------------------------------------------------- */

#include "lhotse/matrix/StVector.h"
//#include "lhotse/rando/DistribSamp.h"
//#include "lhotse/rando/Generator.h"
#ifdef MATLAB_MEX
#include "lhotse/matif/MatlabMatrix.h"
#include "lhotse/matif/mex_for_cpp.h"
#endif

// Local macro for args to 'checkTS'
#define TSARG __FILE__,__LINE__

//BEGINNS(matrix)
  // Public methods

#ifdef MATLAB_MEX
  void StVector::maskMatlab(const MatlabMatrix& mat,const Range& rng)
  {
    reassign(mat.buff(),mat.m()*mat.n(),rng);
  }
#endif

  double StVector::maxRelDiff(const BaseVector<double>& b) const
  {
    double mxel=0.0,aelem,belem,temp;
    int i;

    if (b.size()!=n) throw WrongDimensionException(EXCEPT_MSG(""));
    for (i=0; i<n; i++) {
      belem=b[i]; aelem=operator[](i);
      if (aelem!=0.0)
	temp=fabs((belem-aelem)/aelem);
      else
	temp=0.0;
      if (belem!=0.0)
	temp=std::max(temp,fabs((aelem-belem)/belem));
      if (temp>mxel) mxel=temp;
    }

    return mxel;
  }

  /**
   * The current file format is simply the one of 'BaseVector<double>' (note
   * that this includes type information).
   * This method can also load old FF version with the "@StVector" tag.
   * <p>
   * Old formats:
   * - FF version 0:
   *   - tag: "@StVector"
   *   - FF version number [int]
   *   - ['BaseVector<double>::loadInt']
   * - Old format without FFV number:
   *   NOTE: Such files may not be loaded correctly on architectures with
   *   little endian byte order!
   *   - tag: "StVector"
   *   - ['BaseVector<double>::loadInt']
   *
   * @param is    Binary input file stream
   * @param noTag If true, the tag is not loaded. In this case, we assume
   *              that the format is not an old one
   */
  ifstream& StVector::load(ifstream& is,bool noTag)
  {
    int i=0;

    checkTS(TSARG);
    if (!noTag) {
      ArrayHandle<string> tags(3);
      tags[0]="@BaseVector";
      tags[1]="@StVector";
      tags[2]="StVector";
      FileUtils::loadHeaderMulti(is,tags,i,true); // do not load FFV number
    }
    if (i==0) {
      // Current format
      return BaseVector<double>::load(is,true);
    } else {
      // Old format
      if (i==1) {
	NumberFormats<int>::load(is,&i,1,1,0,4); // FF version number
	if (i!=0) throw FileFormatException("Unknown FF version number");
      }
      int newn;
      NumberFormats<int>::load(is,&newn,1,1,0,4);
      if (newn<0) throw FileFormatException("Negative vector length");
      ensureCapacity(newn);
      NumberFormats<double>::load(is,buff,n);
    }

    return is;
  }
//ENDNS

#undef TSARG
