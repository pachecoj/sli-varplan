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
 * Desc.:  Header class MatlabTools
 * ------------------------------------------------------------------- */

#ifndef MATLABTOOLS_H
#define MATLABTOOLS_H

#include "lhotse/matif/default.h"
#include "lhotse/matif/MatlabMatrix.h"
#include "lhotse/matif/MatlabString.h"

//BEGINNS(matif)
  /**
   * Contains static helper methods, e.g. to process argument lists.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class MatlabTools
  {
  public:
    // Static methods

    /**
     * Reads argument with name 'name' from 'arg'. Type must be double
     * real matrix.
     * If 'rows'!=0, number of rows must be =='rows'. Same for 'cols'.
     * In case of error, a 'MatIFException' with error message is thrown.
     * NOTE: If 'rows' or 'cols' !=0, empty matrix not accepted.
     * NOTE: Returned matrix is NOT a copy, alloc. non-persistent, not
     * owned by matrix object.
     *
     * @param name Arg. name
     * @param arg  Argument
     * @param rows S.a. Optional
     * @param cols "
     * @return     Matrix
     */
    static MatlabMatrix* readMatrix(const char* name,const mxArray* arg,
				    int rows=0,int cols=0);

    /**
     * Reads argument with name 'name' from 'arg'. Type must be double
     * real vector.
     * If 'n'!=0, the length must be =='n'. 'stat' is the vector status:
     * 0: either (def.), 1: column, 2: row.
     * NOTE: Returned matrix is NOT a copy, alloc. non-persistent, not
     * owned by matrix object.
     *
     * @param name Arg. name
     * @param arg  Argument
     * @param n    S.a. Def.: 0
     * @param stat "
     * @return     Vector
     */
    static MatlabMatrix* readVector(const char* name,const mxArray* arg,
				    int n=0,int stat=0);

    /**
     * Reads argument with name 'name' from 'arg'. Type must be double
     * real vector. Elements must be integers.
     * If 'n'!=0, the length must be =='n'. 'stat' is the vector status:
     * 0: either (def.), 1: column, 2: row.
     * If the interval 'iv' is given, each element must lie in 'iv'.
     * NOTE: Returned matrix is NOT a copy, alloc. non-persistent, not
     * owned by matrix object.
     *
     * @param name Arg. name
     * @param arg  Argument
     * @param n    S.a. Def.: 0
     * @param stat "
     * @param iv   S.a. Def: 0
     * @return     Vector
     */
    static MatlabMatrix* readIntVector(const char* name,const mxArray* arg,
				       int n=0,int stat=0,
				       const Interval<double>* iv=0);

    /**
     * Reads argument with name 'name' from 'arg'. Type must be double
     * real scalar. If 'iv' given, ref. to interval which value must fall
     * in.
     *
     * @param name Arg. name
     * @param arg  Argument
     * @param iv   Optional. Def.: 0
     * @return     Scalar
     */
    static double readScalar(const char* name,const mxArray* arg,
			     const Interval<double>* iv=0) {
      char msg[200];
      if (!mxIsDouble(arg) || mxGetM(arg)!=1 || mxGetN(arg)!=1) {
	sprintf(msg,"Expect double scalar for %s",name);
	throw MatIFException(msg);
      }
      double val=*mxGetPr(arg);
      if (iv!=0 && !(*iv)(val)) {
	sprintf(msg,"Value of %s out of range",name);
	throw MatIFException(msg);
      }
      return val;
    }

    /**
     * Version of 'readScalar' which accepts integers only.
     *
     * @param name Arg. name
     * @param arg  Argument
     * @param iv   Optional. Def.: 0
     * @return     Scalar (int)
     */
    static int readScalInt(const char* name,const mxArray* arg,
			   const Interval<int>* iv=0) {
      double val=readScalar(name,arg),temp;
      char msg[200];
      if ((temp=floor(val))!=val) {
	sprintf(msg,"Expect integer scalar for %s",name);
	throw MatIFException(msg);
      }
      int i=(int) temp;
      if (iv!=0 && !(*iv)(i)) {
	sprintf(msg,"Value of %s out of range",name);
	throw MatIFException(msg);
      }
      return i;
    }

    /**
     * Version of 'readScalar' which accepts bools only. A bool must be
     * 0 (false) or 1 (true).
     *
     * @param name Arg. name
     * @param arg  Argument
     * @return     Scalar (bool)
     */
    static bool readScalBool(const char* name,const mxArray* arg) {
      double val=readScalar(name,arg),temp;
      int i;
      char msg[200];
      bool err=((temp=floor(val))!=val);
      i=(int) temp;
      err=err||(i!=0 && i!=1);
      if (err) {
	sprintf(msg,"Expect boolean scalar for %s",name);
	throw MatIFException(msg);
      }
      return (i==1);
    }

    /**
     * Checks whether 'nrhs'<'minr'. If so, throws MatIFException.
     */
    static void checkNrhs(int nrhs,int minr) {
      if (nrhs<minr) {
	char msg[200];
	sprintf(msg,"Expect at least %d input arguments",minr);
	throw MatIFException(msg);
      }
    }

#ifdef LHOTSE_MATLAB_IF
    /**
     * Helper for processing object handles. Looks up 'h' in global
     * 'objList', throws exception using var name 'name' for handle
     * if not found. Next, checks whether object is instance of class
     * 'clname', throws exc. if not (all 'MatIFException'). Finally,
     * returns 'MatIF_baseclass*' for object.
     *
     * @param h      Handle
     * @param name   S.a.
     * @param clname S.a.
     */
    static MatIF_baseclass* procHandle(int h,const char* name,
				       const char* clname);
#endif
  };
//ENDNS

#endif
