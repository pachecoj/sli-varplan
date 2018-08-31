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
 * Desc.:  Definition class MatlabTools
 * ------------------------------------------------------------------- */

#include "lhotse/matif/MatlabTools.h"
#include "lhotse/matif/MatIF_baseclass.h"

//BEGINNS(matif)
  // Extern global variables

  extern MAP_TYPE(int,MatIF_baseclass*) objList;

  // Public static methods

  MatlabMatrix* MatlabTools::readMatrix(const char* name,const mxArray* arg,
					int rows,int cols)
  {
    char msg[200];
    if (!mxIsDouble(arg)) {
      sprintf(msg,"Expect double matrix for %s",name);
      throw MatIFException(msg);
    }
    if (rows!=0 && mxGetM(arg)!=rows) {
      sprintf(msg,"Expect %d rows for %s",rows,name);
      throw MatIFException(msg);
    }
    if (cols!=0 && mxGetN(arg)!=cols) {
      sprintf(msg,"Expect %d columns for %s",cols,name);
      throw MatIFException(msg);
    }
    return new MatlabMatrix(arg);
  }

  MatlabMatrix* MatlabTools::readVector(const char* name,const mxArray* arg,
					int n,int stat)
  {
    char msg[200];
    if (!mxIsDouble(arg)) {
      sprintf(msg,"Expect double vector for %s",name);
      throw MatIFException(msg);
    }
    int len;
    if (mxGetN(arg)==0 || mxGetM(arg)==0) len=0;
    else if (stat==1) {
      if (mxGetN(arg)!=1) {
	sprintf(msg,"Expect column vector for %s",name);
	throw MatIFException(msg);
      }
      len=mxGetM(arg);
    } else if (stat==2) {
      if (mxGetM(arg)!=1) {
	sprintf(msg,"Expect row vector for %s",name);
	throw MatIFException(msg);
      }
      len=mxGetN(arg);
    } else {
      if (mxGetM(arg)==1) len=mxGetN(arg);
      else if (mxGetN(arg)==1) len=mxGetM(arg);
      else {
	sprintf(msg,"Expect vector for %s",name);
	throw MatIFException(msg);
      }
    }
    if (n!=0 && len!=n) {
      sprintf(msg,"Expect vector of size %d for %s",n,name);
      throw MatIFException(msg);
    }
    return new MatlabMatrix(arg);
  }

  MatlabMatrix* MatlabTools::readIntVector(const char* name,
					   const mxArray* arg,int n,
					   int stat,
					   const Interval<double>* iv)
  {
    char msg[200];
    if (!mxIsDouble(arg)) {
      sprintf(msg,"Expect double vector for %s",name);
      throw MatIFException(msg);
    }
    int len;
    if (mxGetN(arg)==0 || mxGetM(arg)==0) len=0;
    else if (stat==1) {
      if (mxGetN(arg)!=1) {
	sprintf(msg,"Expect column vector for %s",name);
	throw MatIFException(msg);
      }
      len=mxGetM(arg);
    } else if (stat==2) {
      if (mxGetM(arg)!=1) {
	sprintf(msg,"Expect row vector for %s",name);
	throw MatIFException(msg);
      }
      len=mxGetN(arg);
    } else {
      if (mxGetM(arg)==1) len=mxGetN(arg);
      else if (mxGetN(arg)==1) len=mxGetM(arg);
      else {
	sprintf(msg,"Expect vector for %s",name);
	throw MatIFException(msg);
      }
    }
    if (n!=0 && len!=n) {
      sprintf(msg,"Expect vector of size %d for %s",n,name);
      throw MatIFException(msg);
    }
    int i;
    double temp;
    const double* buff=mxGetPr(arg);
    for (i=0; i<len; i++) {
      temp=*(buff++);
      if (temp!=floor(temp)) {
	sprintf(msg,"Components of %s must be integer",name);
	throw MatIFException(msg);
      }
    }
    if (iv!=0 &&
	iv->check(mxGetPr(arg),len)!=0) {
      sprintf(msg,"Components of %s violate range",name);
      throw MatIFException(msg);
    }
    return new MatlabMatrix(arg);
  }

#ifdef LHOTSE_MATLAB_IF
  MatIF_baseclass* MatlabTools::procHandle(int h,const char* name,
					   const char* clname)
  {
    char msg[100];

    MAP_CONSTITER(int,MatIF_baseclass*) it=objList.find(h);
    if (it==objList.end()) {
      sprintf(msg,"Object with handle %s==%d does not exist",name,h);
      throw MatIFException(msg);
    }
    if (!it->second->isInstanceOf(clname)) {
      sprintf(msg,"Object %s is no instance of type %s",name,clname);
      throw MatIFException(msg);
    }
    return it->second;
  }
#endif
//ENDNS
