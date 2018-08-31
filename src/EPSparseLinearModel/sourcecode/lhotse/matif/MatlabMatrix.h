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
 * Desc.:  Header class MatlabMatrix
 * ------------------------------------------------------------------- */

#ifndef MATLABMATRIX_H
#define MATLABMATRIX_H

#include "lhotse/matif/default.h"
#include "lhotse/matrix/predecl.h"

//BEGINNS(matif)
  /**
   * Wraps a 'mxArray' object (double real matrix), allows easy access
   * (e.g. to dimensions).
   * We own the wrapped object iff 'isOwn'==true. In this case, the object is
   * dealloc. in the destructor. Optionally, the constructor can make the
   * object persistent. If a matrix is used as return argument of a MEX
   * function, we use 'isOwn'==false.
   * <p>
   * Objects are typically used to access in arguments and to create
   * return arguments.
   * To wrap an input argument, create a non-persistent object which is not an
   * owner:
   *   MatlabMatrix inmat(prhs[XXX]);
   * Use 'maskMatlab' to mask 'MatlabMatrix' objects as 'BaseMatrix<double>'
   * or 'BaseVector<double>'.
   * <p>
   * Return arguments must be non-persistent and not owners:
   *   MatlabMatrix retmat(ROWS,COLS,false);
   *   // ...
   *   plhs[XXX]=retmat.matobj();
   * To return a matrix stored in BaseMatrix<double> 'mat':
   *   MatlabMatrix retmat(mat,false);
   *   plhs[XXX]=retmat.matobj();
   * This creates a physical copy of 'mat'.
   * <p>
   * If a BaseMatrix<double> is to be returned as argument and its size
   * need not change during lifetime, it is more efficient to create it
   * as 'MatlabMatrix' and mask it as 'BaseMatrix<double>' or 'StMatrix':
   *   MatlabMatrix retmat(ROWS,COLS,false);
   *   BaseMatrix<double> mat;
   *   mat.maskMatlab(retmat);
   *   // ... (size of 'mat' must not change, it is mask)
   *   plhs[XXX]=retmat.matobj();
   * <p>
   * Return arguments of Matlab functions called from C++ (using
   * 'mexCallMATLAB') should be wrapped immediately in a 'MatlabMatrix'
   * which owns the object:
   *   mexCallMATLAB(1,lhs,...);
   *   MatlabMatrix ret1(lhs[0],true);
   *   // ... 'ret1' disposes off 'lhs[0]' at end of scope
   * NOTE: If the ret. argument is to be ret. itself by the MEX function,
   * the wrapping 'MatlabMatrix' need not own the object.
   * <p>
   * Apart from this, objects can be used without wrapping input/return
   * arguments, e.g. to call 'mexXXX' functions. Constructors from STATSIM
   * matrix types are given.
   * <p>
   * Masking flat buffer:
   * A simple way to "sell" a STATSIM object as 'mxArray' without copying
   * (do NOT do this for return arguments, they have to be copies!) is
   * given by the constructor from flat buffer. The resulting object is
   * always persistent and no owner (to make sure neither the destructor
   * nor Matlab's cleanup touches the buffer).
   * NOTE: Matlab 'mxArray' does not allow for strided buffers, they have
   * to be contiguous.
   * <p>
   * Persistent objects:
   * 'mat' and this enclosing object have to be made persistent at the same
   * time, which is why the static member 'persFlag' is needed. The object is
   * made persistent iff the new operator (which is overloaded) is called with
   * true:
   *   persobj=new(true) MatlabMatrix(10,10);
   * 'new' sets 'persFlag' to this bool argument. Upon return the constructor
   * reads the static flag and gives 'mat' the same status.
   * ==> Not thread-safe! A better way would be to make the flag a non-static
   *     member, but constructor will first call 'new', then init. all
   *     non-static members to their def. values, so cannot pass information
   *     from 'new'!
   * ==> Option to use independent flags in 'new', constructor, but can lead
   *     to nasty memory inconsistencies!
   * NOTE: 'persFlag' is always reset to false and is init. with false.
   * Also, def. value for 'new' flag is false. Means that non-dynamic
   * construction or 'new' without flag results in non-persistent object.
   * <p>
   * Use of persistent objects:
   * It is safer to use persistent 'MatlabMatrix' objects only if abs.
   * necessary. Better to use persistent 'BaseMatrix', 'StMatrix', ...
   * objects and mask them as 'MatlabMatrix' when communicating with
   * Matlab.
   * <p>
   * Use of non-persistent objects:
   * For non-persistent objects, all memory regions and matrices alloc.
   * through Matlab's MM will be dealloc. if the MEX function is interrupted.
   * However, the destructors are NOT called!
   * For persistent objects, an interrupt does not have any effect. If they
   * are not linked into the global static structure, their reference is
   * lost (memory leak).
   * 1. If an assertion fails, do NOT call 'mexErrMsgTxt' directly, but throw
   *    an exception. This is caught in the MEX function, but before that
   *    destructors are called.
   * 2. Objects which should remain persistent must be linked into the global
   *    static structure as soon as they are created.
   * 3. Temporary objects: use non-persistent object whenever the simple
   *    dealloc. of all memory regions/matrices is equivalent to the
   *    destructors job. Can work even for objects containing subobjects if all
   *    of them are non-persistent and the order of dealloc. does not matter.
   *    For other temp. objects: make persistent, which means memory leak in
   *    the worst case.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class MatlabMatrix
  {
  protected:
    // Static members

    static bool persFlag;

    // Members

    mxArray* mat;  // wrapped matrix object
    bool isOwn;    // 'mat' our object?

  public:
    // Constructors

    /**
     * Empty matrix (init. of 'mxCreateDoubleMatrix') of size 'm'-by-'n'.
     * This is ours iff 'oown'==true (def.: true).
     * NOTE: If the matrix is to be used as return value, use 'oown'=false.
     *
     * @param m    Number rows
     * @param n    Number cols
     * @param oown S.a. Def.: true
     */
    MatlabMatrix(int m,int n,bool oown=true) : isOwn(oown) {
      mat=mxCreateDoubleMatrix(m,n,mxREAL);
      if (persFlag) mexMakeArrayPersistent(mat); // persistent
      persFlag=false; // reset
    }

    /**
     * 1-by-1 matrix representing the scalar 'val'.
     * NOTE: If the matrix is to be used as return value, use 'oown'=false.
     * NOTE: Cast 'val' to 'double' to avoid confusion with other
     * constructors.
     *
     * @param val Scalar value
     * @param oown S.a. Def.: true
     */
    MatlabMatrix(double val,bool oown=true) : isOwn(oown) {
      mat=mxCreateDoubleMatrix(1,1,mxREAL);
      if (persFlag) mexMakeArrayPersistent(mat); // persistent
      persFlag=false; // reset
      double* buff=mxGetPr(mat);
      *buff=val;
    }

    /**
     * Wraps 'mxArray' object 'mxm'. This will be our own iff 'oown'==true
     * (def.: false).
     * Use to wrap object which is owned by someone else, f.ex. an input
     * argument.
     *
     * @param mxm  S.a.
     * @param oown S.a. Def.: false
     */
    MatlabMatrix(const mxArray* mxm,bool oown=false) : isOwn(oown) {
      mat=(mxArray*) mxm;
      if (persFlag) mexMakeArrayPersistent(mat); // persistent
      //mexPrintf("mat-cr1: %d, pers=%d\n",mat,persFlag); // DEBUG!
      persFlag=false; // reset
    }

    /**
     * Copy constructor.
     * Draws physical copy of the matrix wrapped in 'mobj', which is our
     * own iff 'oown'==true (def.: true). To draw a copy of an existing
     * 'mxArray' object 'mxm', use:
     *   MatlabMatrix temp(mxm); // temp does not own mxm
     *   mcopy=new(true) MatlabMatrix(temp); // persistent, owns copy
     *
     * @param mobj Source object
     */
    MatlabMatrix(const MatlabMatrix& mobj,bool oown=true) : isOwn(oown) {
      int m=mobj.m(),n=mobj.n();
      mat=mxCreateDoubleMatrix(m,n,mxREAL);
      memmove((char*) mxGetPr(mat),(const char*) mobj.buff(),
	      m*n*sizeof(double));
      if (persFlag) mexMakeArrayPersistent(mat); // persistent
      //mexPrintf("mat-cr2: %d, pers=%d\n",mat,persFlag); // DEBUG!
      persFlag=false; // reset
    }

    /**
     * Constructor for mask on memory region (size m*n).
     * The res. 'mat' does not belong to us and is always persistent,
     * indep. of 'persFlag'.
     *
     * @param mbuff Memory region (contiguous)
     * @param m     Number rows
     * @param n     Number cols
     */
    MatlabMatrix(double* mbuff,int m,int n) : isOwn(false) {
      mat=mxCreateDoubleMatrix(0,0,mxREAL); // empty
      mxSetPr(mat,mbuff);
      mxSetM(mat,m); mxSetN(mat,n);
      mexMakeArrayPersistent(mat); // persistent
      persFlag=false; // reset
    }

    /**
     * Draws physical copy of the matrix in 'mat', which is our own iff
     * 'oown'==true (def).
     * NOTE: To do the same with a vector, use 'mask':
     *   MatlabMatrix(StMatrix::mask(vec));
     *
     * @param src  Source matrix
     * @param oown S.a. Def.: true
     */
    MatlabMatrix(const BaseMatrix<double>& src,bool oown=true);

    /**
     * Same as constructor from 'BaseMatrix<double>', but the argument is
     * of type 'BaseLinMat<double>'.
     *
     * @param src  Source matrix
     * @param oown S.a. Def.: true
     */
    MatlabMatrix(const BaseLinMat<double>& src,bool oown=true);

    /**
     * Same as constructor from 'BaseMatrix<double>', but the argument is
     * of type 'BaseLinVec<double>'. It is treated as column vector iff
     * 'col'==true.
     *
     * @param src  Source vector
     * @param col  S.a. Def.: true
     * @param oown S.a. Def.: true
     */
    MatlabMatrix(const BaseLinVec<double>& src,bool col=true,bool oown=true);

    /**
     * Draws physical copy of index vector in 'vec', which is our own iff
     * 'oown'==true (def). If 'incr' is true (def), 1 is added to each
     * element (conversion C++ base 0 to Matlab base 1).
     * The vector is column iff 'col'==true (def).
     *
     * @param vec  Source vector
     * @param oown S.a. Def.: true
     * @param incr S.a. Def.: true
     * @param col  S.a. Def.: true
     */
    MatlabMatrix(const BaseVector<int>& vec,bool oown=true,bool incr=true,
		 bool col=true);

    /**
     * Destructor. 'mat' is dealloc.iff 'isOwn'==true.
     */
    virtual ~MatlabMatrix() {
      if (isOwn) {
	//mexPrintf("mat-dest: %d...",mat); // DEBUG!
	mxDestroyArray(mat);
	//mexPrintf("OK\n");
      }
    }

    // Overload new, delete

    void* operator new(size_t sz,bool pflag=false) {
      void* ptr=mxMalloc(sz);
      persFlag=pflag;
      if (pflag) mexMakeMemoryPersistent(ptr);
      //mexPrintf("mat-new: %d\n",ptr); // DEBUG!
      return ptr;
    }

    void operator delete(void* ptr) {
      //mexPrintf("mat-delete: %d...",ptr); // DEBUG!
      mxFree(ptr);
      //mexPrintf("OK\n");
    }

    // Public methods

    /**
     * @return Number of rows
     */
    int m() const {
      return mxGetM(mat);
    }

    /**
     * @return Number of cols
     */
    int n() const {
      return mxGetN(mat);
    }

    /**
     * @return Pointer to flat buffer (col-major)
     */
    double* buff() const {
      return mxGetPr(mat);
    }

    /**
     * @return Pointer to 'mxArray' object
     */
    mxArray* matobj() const {
      return mat;
    }

    /**
     * @return Is 'mat' our own?
     */
    bool ourOwn() const {
      return isOwn;
    }

  protected:
    // Internal methods

    /**
     * Helper for constructors. Copies content of linear matrix object
     * 'linm' into Matlab matrix object 'mat'.
     *
     * @param linm S.a.
     */
    void copyContent(const BaseLinMat<double>& linm);
  };
//ENDNS

#endif
