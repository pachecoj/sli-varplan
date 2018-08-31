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
 * Module: GLOBAL
 * Desc.:  Header class MatlabDebug
 * ------------------------------------------------------------------- */

#ifndef MATLABDEBUG_H
#define MATLABDEBUG_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/global.h" // global header
#include "lhotse/matrix/predecl.h"
#if defined(MATLABDEBUG_USEMEX) && defined(MATLAB_MEX)
#include "lhotse/matif/mex_for_cpp.h"
#endif

USING(matrix);

// Macros: Shortcuts

/*
 * - SETMAT:    Double matrix (StMatrix)
 * - SETVEC:    Double vector (StVector)
 * - SETIND:    Index (BaseVector<int>). Matlab code to add 1 to each element
 *              is automatically inserted!
 * - SETIVEC:   Like 'SETIND', but no addition of 1
 * - SETDARR:   Double array -> Matlab vector
 * - SETSYMM:   Symmetric matrix (looks at str. patt.)
 * - SETSCAL:   Scalar (arg. is casted to double)
 * - SETCHOL:   Cholesky factor (lower triangle)
 *
 * - COMP:      Compare two Matlab variables component by component
 * - OUTP:      Write string into Matlab script -> Matlab commands
 *              NOTE: See 'output'!
 * - PTEXT:     Have Matlab print the string
 * - PVAR:      Have Matlab output the variable
 * - DOCOMP:    Compares C++ against Matlab variable (meta macro)
 * - DOCOMP_WB: As DOCOMP, but C++ object copied to Matlab variable after
 *              successful comparison
 * - PUTMAT:    Write Matlab variable into C++ matrix
 */
#define SETMAT(var,mat) MatlabDebug::setMatrix(var,mat,0,__FILE__,__LINE__)
#define STOREMAT(mat,fname) MatlabDebug::setMatrix("",mat,fname)
#define SETVEC(var,mat) MatlabDebug::setVector(var,mat,0,__FILE__,__LINE__)
#define SETIVEC(var,mat) MatlabDebug::setIntVector(var,mat,0,false,__FILE__,__LINE__)
#define SETIND(var,mat) MatlabDebug::setIntVector(var,mat,0,true,__FILE__,__LINE__)
#define STOREVEC(mat,fname) MatlabDebug::setVector("",mat,fname)
#define SETDARR(var,arr,size) MatlabDebug::setDoubleArray(var,arr,size,0,__FILE__,__LINE__)
#define STOREDARR(arr,size,fname) MatlabDebug::setDoubleArray("",arr,size,fname)
#define SETCHOL(var,fact) MatlabDebug::setCholFact(var,fact,__FILE__,__LINE__)
#define SETSYMM(var,mat) MatlabDebug::setSymm(var,mat,__FILE__,__LINE__)
#define SETSCAL(var,val) MatlabDebug::setScalar(var,(double) (val),__FILE__,__LINE__)
#define COMP(var1,var2) MatlabDebug::compare(var1,var2,__FILE__,__LINE__)
#define PTEXT(txt) MatlabDebug::print(txt)
#define PVAR(var) MatlabDebug::printVar(var)
#define OUTP(txt) MatlabDebug::output(txt,__FILE__,__LINE__)
#define DOCOMP(obj,expr,var) MatlabDebug::doCompare(obj,expr,var,false,0,__FILE__,__LINE__)
#define DOCOMP2(obj,expr,var,tname) MatlabDebug::doCompare(obj,expr,var,false,tname,__FILE__,__LINE__)
#define DOCOMP_WB(obj,expr,var) MatlabDebug::doCompare(obj,expr,var,true,0,__FILE__,__LINE__)
#define DOCOMP2_WB(obj,expr,var,tname) MatlabDebug::doCompare(obj,expr,var,true,tname,__FILE__,__LINE__)
#define PUTMAT(var,mat) MatlabDebug::putMatrix(var,mat,__FILE__,__LINE__)

/**
 * Provides static methods for Matlab debugging of LHOTSE code, a very
 * powerful way of finding bugs.
 * The idea is to replicate all code by equivalent Matlab commands, passing
 * all intermediate LHOTSE variables to Matlab and comparing the two
 * versions. The debug script execution stops once a significant difference
 * is detected. A message indicates where in the LHOTSE code the difference
 * occurs. Intermediate messages can be printed as well.
 * When a significant difference is detected, its value is printed together
 * with the code file and line number. Furthermore, the global Matlab
 * variables DIFFVAR1, DIFFVAR2 are init. with the differing arguments
 * (they are not printed).
 * <p>
 * There are two modes:
 * - file-based
 * - MEX based
 * The file-based mode is much less efficient, but more robust and can be
 * used for code which is not a MEX function and which is run independently
 * of Matlab. It is limited to rather short stretches of code.
 * The MEX based mode is much more efficient and can be used for arbitrary
 * amounts of code. It can only be run for a MEX function, and is used iff
 * both MATLAB_MEX and MATLABDEBUG_USEMEX are defined.
 * <p>
 * Running status:
 * Maintained in 'status':
 * 0: Not active
 * 1: Active, doing comparisons
 * 2: Active, but skipping:
 *    - comparisons in 'compare', 'doCompare'
 *    - 'print', 'printVar'
 *    - assignments in 'doCompare'
 *    Not skipped are:
 *    - explicit assignments via 'setXXX'
 *    - Matlab commands in 'output', 'doCompare'
 *    - 'putMatrix'
 *    Useful to do Matlab debug computations without the comparisons.
 * 3: Checking for NaN:
 *    Special mode. Ignores everything, except for commands which import
 *    LHOTSE variables into Matlab ('setXXX', 'doCompare', ...). For
 *    these, the import is done and the Matlab variable is checked for NaN
 *    entries. Nothing else is done.
 *    ATTENTION: Implemented for MEX-based mode only!
 * <p>
 * 'isActive' returns true iff active ('status'!=0). 'activate' sets status
 * to 1, 'deactivate' status to 0 (DOWNW. COMPAT.). Use 'setStatus' instead.
 * In file-based mode, base filename must be passed with first activation.
 * NOTE: Most methods here check status and do nothing if deactivated. The
 * wrapping in
 *   if (MatlabDebug::isActive()) {
 *     // ...
 *   }
 * is not necessary. Code not needed can be wrapped in
 * #ifdef MATLAB_DEBUG_OLD ... #endif. This flag is never set, so such code is
 * not even compiled.
 * <p>
 * File-based mode:
 * A self-contained debug script is produced in 'debugscript.m'. Small
 * matrices are written as text, larger as binary objects:
 * 'debugscript#.mat'
 * where # is a running number. These are not Matlab MAT files, but native
 * LHOTSE objects. The script reads these objects using 'fread'.
 * Matrices and vectors are stored as binary objects iff they contain more than
 * VEC_TEXT_LIMIT elements, otherwise they are stored in text format in
 * 'debugscript.m'. Symmetric matrices are always stored as binary files. A
 * filename can be passed to
 * these methods in which case this name overrides the automatically generated
 * one. In this case, NO Matlab code to read this object is written into the
 * debugscript!
 * <p>
 * Be careful: The script file can grow very large. The script output can
 * be switched on and off at any time. It is on by default.
 * <p>
 * NOTE: We require several Matlab functions which are part of LHOTSE:
 * - LOADSTMATRIX,LOADBASEVECTOR: load LHOTSE objects
 * - RELDIFF: computes symmetric relative difference
 * Make sure these are in the Matlab path when the debug script is run.
 * <p>
 * MEX based mode:
 * Used iff MATLAB_MEX and MATLABDEBUG_USEMEX are both defined. In this
 * case, we use 'mexEvalString' to evaluate the commands in the Matlab
 * workspace (the one calling the MEX function), 'mexPutVariable' to
 * put LHOTSE variables into the caller workspace. In this mode, a
 * comparison which does not result in a signif. difference does not
 * cause an output. The MEX function is aborted once a significant
 * difference is measured. File name and line number are printed in the
 * error message.
 * NOTE: In the MEX based mode, features here can be controlled from Matlab
 * (see MATLABDEBUG_xxx methods in 'MatIF_baseclass').
 * <p>
 * Common variable prefix:
 * In MEX based mode, all variables are kept in the caller workspace, in
 * file-based mode in the global workspace. To avoid interference with the
 * variables there, all debug variables should have an unusual common prefix
 * ('commonPrefix').
 * ATTENTION: This prefix is added automatically to names passed to
 * 'setXXX', but this cannot be done for expressions passed to 'output'.
 * However, all occurences of "%" in expressions passed to 'output' are
 * replaced by 'commonPrefix'.
 * In file based mode, the def. value for 'commonPrefix' is "".
 * <p>
 * Variables "local" to a part of code:
 * All Matlab debug variables are kept in the same workspace (in either
 * mode), so there can be conflicts between several pieces of debug code.
 * One can make variables "local" by extending the common prefix:
 *   string oldPref=MatlabDebug::getCommonPrefix();
 *   MatlabDebug:setCommonPrefix(oldPref+"XXX_");
 *   // local code ...
 *   MatlabDebug::setCommonPrefix(oldPref);
 * This is highly recommended!
 * NOTE: Changing the common prefix makes all variables local.
 * ==> TODO: Allow some variables to be global, but still have an automatic
 *     prefix!
 * <p>
 * Rules for writing Matlab commands in OUTP, 'output':
 * In order to use Matlab debugging in either mode, the following rules must
 * be followed:
 * - all variables names which are used in 'setXXX', 'compare', 'doCompare'
 *   have to be prefixed by "%" in the expression passed to OUTP. The "%" is
 *   automatically replaced by 'commonPrefix'
 *   NOTE: In file-based mode, 'commonPrefix' is "".
 * - only use pure names as variable names in 'setXXX', 'compare', 'doCompare'.
 *   For example, something like
 *     DOCOMP(obj,"%mat(%b,%a)=...","mat(b,a)");
 *   does not work, because the variables "b", "a" do not have the common
 *   prefix. Use:
 *     DOCOMP(obj,"%tmat2=...","tmat2");
 *     OUTP("%mat(%b,%a)=%tmat2;");
 * - do not use the default comp. variable names used by 'doCompare' ("tscal",
 *   "tmat","tvec")
 * - NOTE: Variables are kept in the caller workspace. They are lost if the
 *   MEX function is called several times. To keep variables around, declare
 *   them global (using 'output') before using them!
 * <p>
 * Threshold for comparisons:
 * In both modes, the Matlab function RELDIFF is used to compute the
 * difference which is significant iff larger than 'diffThres'. This is
 * a Matlab variable in file-based mode, but a member here in MEX based
 * mode (see 'setDiffThres').
 * <p>
 * ATTENTION:
 * Some sources of errors:
 * - Matlab debug variables do not behave like local variables corr. to the
 *   C++ method they are used within. They are variables in the Matlab
 *   caller workspace
 *   ==> Use Matlab debug variables locally only, or use names not used
 *       anywhere else
 * - Matlab debug variables belong to the Matlab caller workspace. To share
 *   a variable between several callers, make it global
 * - 'setXXX' cannot be used to assign values to fields in a structure.
 *   For example,
 *     SETMAT("a.b",mat);
 *   is illegal, use instead:
 *     SETMAT("tmat",mat); OUTP("%a.b=%tmat;");
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
class MatlabDebug
{
protected:
  // Members
  static int status;               // Running status
#if !defined(MATLABDEBUG_USEMEX) || !defined(MATLAB_MEX)
  static Handle<ofstream> os;
  static int runNo;                // Running number for binary objects
  static string baseName;          // Base filename of the task
  static char nameBuff[200];       // Buffer for running filenames
#else
  static double diffThres;         // Threshold for comparisons
#endif
  static string commonPrefix;
  static char cmdBuff[2000];

public:
  // Public static methods

  /**
   * Activates debug facility. The base filename must be given only the
   * first time the method is called. Upon later calls, the argument is
   * ignored (and can be left out).
   * NOTE: 'baseP' is not required in MEX based mode and is ignored if
   * given.
   * NOTE: Status set to 1, use 'setStatus' to change that.
   *
   * @param baseP Base filename (optional, see above)
   */
  static void activate(const string& baseP) {
    if (status==0) {
#if !defined(MATLABDEBUG_USEMEX) || !defined(MATLAB_MEX)
      if (runNo==-1) {
	baseName=baseP;
	sprintf(nameBuff,"%sdebugscript.m",baseName.c_str());
	os.changeRep(new ofstream(nameBuff));
	if (!(*os)) throw FileUtilsException("Cannot open Matlab debug file!");
	os->precision(20); os->setf(std::ios::scientific);
	runNo=0;
      } else {
	cout << "MatlabDebug: WARNING: Base filename already set. 'baseP' argument is ignored." << endl;
      }
#endif
    }
    status=1;
  }

  static void activate() {
    setStatus(1);
  }

  static void deactivate() {
    setStatus(0);
  }

  static bool isActive() {
    return (status!=0);
  }

  /**
   * Set running status (see header comment).
   *
   * @param stat New status
   */
  static void setStatus(int stat) {
    if (status==stat) return;
    if (stat<0 || stat>3)
      throw InvalidParameterException(EXCEPT_MSG(""));
#if !defined(MATLABDEBUG_USEMEX) || !defined(MATLAB_MEX)
    if (status==0) {
      if (runNo==-1)
	throw WrongStatusException("Matlab debug feature not initialised. Pass a base filename!");
      // Re-open Matlab debug script file
      sprintf(nameBuff,"%sdebugscript.m",baseName.c_str());
      os.changeRep(new ofstream(nameBuff,ios::app));
      if (!(*os)) throw FileUtilsException("Cannot open Matlab debug file!");
      os->precision(20); os->setf(std::ios::scientific);
    } else if (stat==0) {
      // Close Matlab debug script file
      os.changeRep(0);
    }
#endif
    status=stat;
  }

  /**
   * @return Running status
   */
  static int getStatus() {
    return status;
  }

  /**
   * @return Common prefix
   */
  static const string& getCommonPrefix() {
    return commonPrefix;
  }

  /**
   * @param New common prefix
   */
  static void setCommonPrefix(const string& cp) {
    commonPrefix=cp;
  }

#if defined(MATLABDEBUG_USEMEX) && defined(MATLAB_MEX)
  /**
   * Set comparison threshold 'diffThres'.
   *
   * @param val New value (positive)
   */
  static void setDiffThres(double val) {
    if (val<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
    diffThres=val;
  }
#endif

  /**
   * We compare matrices by computing the corr. matrix of relative differences
   * (using the LHOTSE Matlab function RELDIFF), find the element of max.
   * abs. value and compare against the threshold 'diffThres'.
   * NOTE: 'file', 'line' used in MEX based mode only!
   * <p>
   * NOTE: The comparison fails if any of the arguments has any NaN elements,
   * even if NaN elements are at the same positions resp.
   */
  static void compare(const string& var1,const string& var2,
		      const char* file=0,int line=0);

  /**
   * Assigns standard matrix to matlab variable. If 'fname' is passed,
   * we ignore the variable name, store the object into a binary file readable
   * by Matlab and insert no code into the Matlab debug script (only in
   * file based mode!).
   * ATTENTION: 'var' must ref. to direct variable, NOT structure field.
   */
  static void setMatrix(const string& var,const StMatrix& mat,
			const char* fname=0,const char* file=0,int line=0);

  /**
   * Same as 'setMatrix', but uses linear helper object of type
   * 'BaseLinMat<double>' (used to interface with BLAS/LAPACK code).
   * ATTENTION: 'var' must ref. to direct variable, NOT structure field.
   * <p>
   * ATTENTION: Not implemented for file-based mode!
   */
  static void setMatrix(const string& var,const BaseLinMat<double>& mat,
			const char* fname=0,const char* file=0,int line=0);

  /**
   * Assigns standard vector to matlab variable. If a filename is passed,
   * we ignore the variable name, store the object into a binary file readable
   * by Matlab and insert no code into the Matlab debug script.
   * ATTENTION: 'var' must ref. to direct variable, NOT structure field.
   */
  static void setVector(const string& var,const StVector& mat,
			const char* fname=0,const char* file=0,int line=0);

  /**
   * Assigns 'vec' to Matlab variable 'var'. If 'fname' is passed, we ignore
   * 'var' and store the object into a binary file readable by Matlab.
   * If 'addOne' is true, Matlab code to add 1 to each element is produced
   * (useful for indexes!).
   * ATTENTION: 'var' must ref. to direct variable, NOT structure field.
   */
  static void setIntVector(const string& var,const BaseVector<int>& vec,
			   const char* fname=0,bool addOne=false,const char* file=0,int line=0);

  /**
   * Assigns standard vector to matlab variable. If 'fname' is passed,
   * we ignore the variable name, store the object into a binary file readable
   * by Matlab and insert no code into the Matlab debug script (only in file
   * based mode!).
   * ATTENTION: 'var' must ref. to direct variable, NOT structure field.
   */
  static void setDoubleArray(const string& var,const double* arr,int size,
			     const char* fname=0,const char* file=0,int line=0);

  /**
   * Assigns symm. matrix to matlab variable. Structure pattern must not
   * be 'strctNormal'.
   * NOT DOWNW. COMPAT.!! Earlier version used 'SymmMatrix<double>'.
   * ATTENTION: 'var' must ref. to direct variable, NOT structure field.
   */
  static void setSymm(const string& var,const StMatrix& mat,
		      const char* file=0,int line=0);

  /**
   * Writes Cholesky factor, given by 'fact' into Matlab var. 'var'.
   * The structure pattern is used. For the 'xxxNDg' types, the diagonal
   * is set to 1. For the 'normal' type, the matrix is assumed lower
   * triangular.
   * NOT DOWNW. COMPAT.!! Earlier version always assumed lower triangular!
   * ATTENTION: 'var' must ref. to direct variable, NOT structure field.
   *
   * @param var  Vatiable name
   * @param fact Lower-tr. factor
   */
  static void setCholFact(const string& var,const StMatrix& fact,
			  const char* file=0,int line=0);

  /**`
   * Assigns scalar to Matlab variable
   * ATTENTION: 'var' must ref. to direct variable, NOT structure field.
   */
  static void setScalar(const string& var,double val,const char* file=0,int line=0);

  /**
   * Copies content of Matlab variable 'var' back into 'mat'. Does nothing
   * in status 3.
   * ATTENTION: Works in MEX-based mode only.
   */
  static void putMatrix(const string& var,StMatrix& mat,const char* file=0,
			int line=0);

  /**
   * Prints text
   */
  static void print(const string& text);

  /**
   * Prints variable
   */
  static void printVar(const string& var);

  /**
   * Writes string directly into script. Use this to run Matlab commands.
   * ATTENTION: Multi-line commands cannot be passed using multiple calls
   * of 'output' in the MEX based mode, use one line only!
   * <p>
   * NOTE: Occurences of "%" are replaced by the common prefix
   * 'commonPrefix' which automatically prefixes all variables names used
   * in 'setXXX' commands (MEX based mode only).
   */
  static void output(const string& txt,const char* file=0,int line=0);

  /**
   * Reads matrix from binary file produced by 'setMatrix'. This method is not
   * needed during the normal Matlab debug actions. The name of the binary file
   * must be passed.
   * This method works even if the Matlab debug feature is deactivated, but
   * only in file based mode.
   */
  static void readMatrix(StMatrix& mat,const char* fname);

  /**
   * Higher level method to set matrix / column vector / scalar as Matlab
   * expression and compare against corr. LHOTSE variable 'obj':
   * - call OUTP with 'expr'
   * - call SETxxx passing 'obj' (appr. form of xxx is chosen). The var.
   *   name can be passed in 'tempName'. otherwise is TMAT for matrix, TVEC
   *   for vector, TSCAL for scalar
   * - call COMP. Name of Matlab object is 'varName'
   * NOTE: 'expr' can be empty, in which case OUTP is not called.
   * <p>
   * If 'copy'==true (def.: false), 'obj' (as Matlab variable) is copied
   * to 'varName' after a successful comparison.
   *
   * @param obj      LHOTSE object (matrix / col. vector / scalar)
   * @param expr     S.a.
   * @param varName  S.a.
   * @param copy     S.a. Def.: false
   * @param tempName Optional, s.a.
   */
  static void doCompare(const StMatrix& obj,const string& expr,
			const string& varName,bool copy=false,
			const char* tempName=0,const char* file=0,int line=0);
  static void doCompare(const StVector& obj,const string& expr,
			const string& varName,bool copy=false,
			const char* tempName=0,const char* file=0,int line=0);
  static void doCompare(const BaseVector<int>& obj,const string& expr,
			const string& varName,bool copy=false,
			const char* tempName=0,const char* file=0,int line=0);
  static void doCompare(double obj,const string& expr,
			const string& varName,bool copy=false,
			const char* tempName=0,const char* file=0,int line=0);

protected:
  // Internal methods

  /**
   * All occurences of "%" in 'str' are replaced by 'commonPrefix'.
   *
   * @param str Input
   * @return    Output
   */
  static string insertCommonPrefix(const string& str) {
    string::size_type pos;
    string ret(str);
    string pp("%");
    while ((pos=ret.find(pp))!=string::npos)
      ret=ret.replace(pos,1,commonPrefix);

    return ret;
  }

#if !defined(MATLAB_MEX) || !defined(MATLABDEBUG_USEMEX)
  /**
   * File-based mode only. Draws new file name and writes into 'nameBuff'.
   */
  static void drawName() {
    sprintf(nameBuff,"%sdebugscript%d.mat",baseName.c_str(),runNo++);
  }
#endif

#if defined(MATLAB_MEX) && defined(MATLABDEBUG_USEMEX)
  /**
   * MEX based mode only.
   * Puts file name 'file' and line number 'line' into Matlab caller
   * workspace as variables 'matlabdebug_file', 'matlabdebug_line'.
   *
   * @param file File name
   * @param line Line number
   */
  static void putFileAndLine(const char* file,int line);

  /**
   * MEX based mode only.
   * Puts 'diffThres' member into caller workspace as variable of same
   * name.
   */
  static void putDiffThres();

  /**
   * MEX based mode only.
   * Wrapper for MEX 'mexPutVariable'. Required because the MEX interface
   * did only offer 'mexPutVariable' in Matlab 6.5.
   * NOTE: Uses MATLAB_VER65 (see global.h).
   * NOTE: 'var' is prefixed by 'commonPrefix' here, this is not done if
   * 'noprefix'==true.
   *
   * @param works    Workspace
   * @param var      Variable name
   * @param arr      Array
   * @param noprefix S.a. Def.: false
   * @return         Ret. value of 'mexPutVariable'
   */
  static int putVariable(const char* works,const string& var,mxArray* arr,
			 bool noprefix=false);

  /**
   * MEX based mode only.
   * Wrapper for MEX 'mexGetVariable'. Required because the MEX interface
   * did only offer 'mexGetVariable' in Matlab 6.5.
   * NOTE: Uses MATLAB_VER65 (see global.h).
   * NOTE: 'var' is prefixed by 'commonPrefix' here, this is not done if
   * 'noprefix'==true.
   *
   * @param works    Workspace
   * @param var      Variable name
   * @param noprefix S.a. Def.: false
   * @return         Array, or 0 if failed
   */
  static mxArray* getVariable(const char* works,const string& var,
			      bool noprefix=false);
#endif
};

#endif
