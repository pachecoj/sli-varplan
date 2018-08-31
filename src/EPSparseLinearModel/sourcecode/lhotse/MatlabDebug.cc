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
 * Desc.:  Definition class MatlabDebug
 * ------------------------------------------------------------------- */

#include "lhotse/MatlabDebug.h"
#include "lhotse/matrix/StMatrix.h"
#include "lhotse/matrix/StVector.h"
#include "lhotse/matrix/BaseLinMat.h"
#include "lhotse/NumberFormats.h"
#if defined(MATLABDEBUG_USEMEX) && defined(MATLAB_MEX)
#include "lhotse/matif/MatlabMatrix.h"
#endif

// Constants

#define VEC_TEXT_LIMIT 20

// Static members

int MatlabDebug::status(0);
#if !defined(MATLABDEBUG_USEMEX) || !defined(MATLAB_MEX)
Handle<ofstream> MatlabDebug::os(0);
int MatlabDebug::runNo(-1); // not initialised
char MatlabDebug::nameBuff[200];
string MatlabDebug::baseName("");
string MatlabDebug::commonPrefix("");
#else
double MatlabDebug::diffThres(1e-5);
string MatlabDebug::commonPrefix("MATLABDEBUG_");
#endif
char MatlabDebug::cmdBuff[2000];

// Public static methods

void MatlabDebug::compare(const string& var1,const string& var2,
			  const char* file,int line)
{
  if (status==1) {
    string vname1(commonPrefix+var1),vname2(commonPrefix+var2);
#if !defined(MATLABDEBUG_USEMEX) || !defined(MATLAB_MEX)
    // File-based
    if (file!=0) {
      (*os) << "if ~isempty(find(isnan(" << vname1 << "))), error('** " << vname1 << " CONTAINS NaN!! " << file << "(" << line << ")'); end" << endl;
      (*os) << "if ~isempty(find(isnan(" << vname2 << "))), error('** " << vname2 << " CONTAINS NaN!! " << file << "(" << line << ")'); end" << endl;
    } else {
      (*os) << "if ~isempty(find(isnan(" << vname1 << "))), error('** " << vname1 << " CONTAINS NaN!!'); end" << endl;
      (*os) << "if ~isempty(find(isnan(" << vname2 << "))), error('** " << vname2 << " CONTAINS NaN!!'); end" << endl;
    }
    (*os) << "diff=max(max(reldiff(" << vname1 << "," << vname2 << ",1)));" << endl;
    (*os) << "fprintf('** delta(" << vname1 << ',' << vname2 << ")=%f\\n',diff);" <<
      endl;
    if (file!=0)
      (*os) << "if diff>diffThres, error('** SIGNIFICANT! " << file << "(" << line << ")'); end" << endl;
    else
      (*os) << "if diff>diffThres, error('** SIGNIFICANT!'); end" << endl;
#else
    // MEX based
    putFileAndLine(file,line);
    putDiffThres();
    if (file!=0)
      sprintf(cmdBuff,"if ~isempty(find(isnan(%s))), error(sprintf('** NaN DETECTED\\n   Var=%s\\n   File: %%s (line %%d)',MATDEBUGXYZ_file,MATDEBUGXYZ_line)); end; if ~isempty(find(isnan(%s))), error(sprintf('** NaN DETECTED\\n   Var=%s\\n   File: %%s (line %%d)',MATDEBUGXYZ_file,MATDEBUGXYZ_line)); end",vname1.c_str(),vname1.c_str(),vname2.c_str(),vname2.c_str());
    else
      sprintf(cmdBuff,"if ~isempty(find(isnan(%s))), error('** NaN DETECTED\\n   Var=%s'); end; if ~isempty(find(isnan(%s))), error('** NaN DETECTED\\n   Var=%s'); end",vname1.c_str(),vname1.c_str(),vname2.c_str(),vname2.c_str());
    if (mexEvalString(cmdBuff)!=0)
      throw MatIFException("MatlabDebug::compare",file,line);
    sprintf(cmdBuff,"MATDEBUGXYZ_diff=max(max(reldiff(%s,%s,1)));",
	    vname1.c_str(),vname2.c_str());
    if (mexEvalString(cmdBuff)!=0)
      throw MatIFException("MatlabDebug::compare",file,line);
    if (file!=0)
      sprintf(cmdBuff,"if MATDEBUGXYZ_diff>MATDEBUGXYZ_diffThres, diffvar1=%s; diffvar2=%s; error(sprintf('** SIGNIFICANT (diff=%%f)\\n   Vars=%s,%s\\n   File: %%s (line %%d)',MATDEBUGXYZ_diff,MATDEBUGXYZ_file,MATDEBUGXYZ_line)); end",vname1.c_str(),vname2.c_str(),vname1.c_str(),vname2.c_str());
    else
      sprintf(cmdBuff,"if MATDEBUGXYZ_diff>MATDEBUGXYZ_diffThres, diffvar1=%s; diffvar2=%s; error(sprintf('** SIGNIFICANT (diff=%%f)\\n   Vars=%s,%s',MATDEBUGXYZ_diff)); end",vname1.c_str(),vname2.c_str(),vname1.c_str(),vname2.c_str());
    if (mexEvalString(cmdBuff)!=0)
      throw MatIFException("MatlabDebug::compare",file,line);
#endif
  }
}

void MatlabDebug::setMatrix(const string& var,const StMatrix& mat,
			    const char* fname,const char* file,int line)
{
  int i,j;
  bool storeOnly=(fname!=0);

  if (status!=0) {
    string vname(commonPrefix+var);
#if !defined(MATLABDEBUG_USEMEX) || !defined(MATLAB_MEX)
    if (status==3) throw NotImplemException(EXCEPT_MSG(""));
    // File-based
    if (!storeOnly && (mat.rows()*mat.cols()<=VEC_TEXT_LIMIT)) {
      os->precision(20); os->setf(ios::scientific);
      (*os) << vname << "=[";
      for (i=0; i<mat.rows(); i++) {
	for (j=0; j<mat.cols(); j++)
	  (*os) << mat.get(i,j) << " ";
	(*os) << "; " << endl;
      }
      (*os) << "];" << endl;
    } else {
      if (!storeOnly) drawName();
      ofstream bin(storeOnly?fname:nameBuff);
      if (!bin) throw FileUtilsException("Cannot open MatlabDebug binary object file");
      mat.save(bin);
      CLOSESTR(bin);
      if (!storeOnly)
	(*os) << vname << "=loadbasematrix('" << nameBuff << "');" << endl;
    }
#else
    // MEX based mode
    MatlabMatrix mmat(mat);
    if (putVariable("caller",var,mmat.matobj())!=0)
      throw MatIFException(EXCEPT_MSG(""));
    if (status==3) {
      putFileAndLine(file,line);
      if (file!=0)
	sprintf(cmdBuff,"if ~isempty(find(isnan(%s))), error(sprintf('** NaN DETECTED\\n   Var=%s\\n   File: %%s (line %%d)',MATDEBUGXYZ_file,MATDEBUGXYZ_line)); end",vname.c_str(),vname.c_str());
      else
	sprintf(cmdBuff,"if ~isempty(find(isnan(%s))), error('** NaN DETECTED\\n   Var=%s'); end",vname.c_str(),vname.c_str());
      if (mexEvalString(cmdBuff)!=0)
	throw MatIFException("MatlabDebug::setMatrix",file,line);
    }
#endif
  }
}

void MatlabDebug::setMatrix(const string& var,const BaseLinMat<double>& mat,
			    const char* fname,const char* file,int line)
{
  int i,j;
  bool storeOnly=(fname!=0);

  if (status!=0) {
#if !defined(MATLABDEBUG_USEMEX) || !defined(MATLAB_MEX)
    throw NotImplemException(EXCEPT_MSG(""));
#else
    // MEX based mode
    MatlabMatrix mmat(mat);
    if (putVariable("caller",var,mmat.matobj())!=0)
      throw MatIFException(EXCEPT_MSG(""));
    if (status==3) {
      string vname(commonPrefix+var);
      putFileAndLine(file,line);
      if (file!=0)
	sprintf(cmdBuff,"if ~isempty(find(isnan(%s))), error(sprintf('** NaN DETECTED\\n   Var=%s\\n   File: %%s (line %%d)',MATDEBUGXYZ_file,MATDEBUGXYZ_line)); end",vname.c_str(),vname.c_str());
      else
	sprintf(cmdBuff,"if ~isempty(find(isnan(%s))), error('** NaN DETECTED\\n   Var=%s'); end",vname.c_str(),vname.c_str());
      if (mexEvalString(cmdBuff)!=0)
	throw MatIFException("MatlabDebug::setMatrix",file,line);
    }
#endif
  }
}

void MatlabDebug::setVector(const string& var,const StVector& vec,
			    const char* fname,const char* file,int line)
{
  int i;
  bool storeOnly=(fname!=0);

  if (status!=0) {
    string vname(commonPrefix+var);
#if !defined(MATLABDEBUG_USEMEX) || !defined(MATLAB_MEX)
    if (status==3) throw NotImplemException(EXCEPT_MSG(""));
    if (!storeOnly && vec.size()<=VEC_TEXT_LIMIT) {
      os->precision(20); os->setf(ios::scientific);
      (*os) << vname << "=[";
      for (i=0; i<vec.size(); i++)
	(*os) << vec[i] << ';' << endl;
      (*os) << "];" << endl;
    } else {
      if (!storeOnly) drawName();
      ofstream bin(storeOnly?fname:nameBuff);
      if (!bin)
	throw FileUtilsException("Cannot open MatlabDebug binary object file");
      vec.save(bin);
      CLOSESTR(bin);
      if (!storeOnly)
	(*os) << vname << "=loadbasevector('" << nameBuff << "',6);" << endl;
    }
#else
    // MEX based mode
    MatlabMatrix mvec(StMatrix::mask(vec));
    if (putVariable("caller",var,mvec.matobj())!=0)
      throw MatIFException(EXCEPT_MSG(""));
    if (status==3) {
      putFileAndLine(file,line);
      if (file!=0)
	sprintf(cmdBuff,"if ~isempty(find(isnan(%s))), error(sprintf('** NaN DETECTED\\n   Var=%s\\n   File: %%s (line %%d)',MATDEBUGXYZ_file,MATDEBUGXYZ_line)); end",vname.c_str(),vname.c_str());
      else
	sprintf(cmdBuff,"if ~isempty(find(isnan(%s))), error('** NaN DETECTED\\n   Var=%s'); end",vname.c_str(),vname.c_str());
      if (mexEvalString(cmdBuff)!=0)
	throw MatIFException("MatlabDebug::setVector",file,line);
    }
#endif
  }
}

void MatlabDebug::setIntVector(const string& var,const BaseVector<int>& vec,
			       const char* fname,bool addOne,const char* file,int line)
{
  int i;
  bool storeOnly=(fname!=0);

  if (status!=0) {
    string vname(commonPrefix+var);
#if !defined(MATLABDEBUG_USEMEX) || !defined(MATLAB_MEX)
    if (status==3) throw NotImplemException(EXCEPT_MSG(""));
    if (!storeOnly && vec.size()<=VEC_TEXT_LIMIT) {
      os->precision(20); os->setf(ios::scientific);
      (*os) << vname << "=[";
      for (i=0; i<vec.size(); i++)
	(*os) << vec[i] << ';' << endl;
      (*os) << "];" << endl;
    } else {
      if (!storeOnly) drawName();
      ofstream bin(storeOnly?fname:nameBuff);
      if (!bin)
	throw FileUtilsException("Cannot open MatlabDebug binary object file");
      vec.save(bin);
      CLOSESTR(bin);
      if (!storeOnly)
	(*os) << vname << "=loadbasevector('" << nameBuff << "',3);" << endl;
    }
    if (addOne)
      (*os) << vname << "=" << vname << "+1;" << endl;
#else
    // MEX based mode
    if (status==3) addOne=false; // adding might change things
    MatlabMatrix mvec(vec,true,addOne);
    if (putVariable("caller",var,mvec.matobj())!=0)
      throw MatIFException(EXCEPT_MSG(""));
    if (status==3) {
      putFileAndLine(file,line);
      if (file!=0)
	sprintf(cmdBuff,"if ~isempty(find(isnan(%s))), error(sprintf('** NaN DETECTED\\n   Var=%s\\n   File: %%s (line %%d)',MATDEBUGXYZ_file,MATDEBUGXYZ_line)); end",vname.c_str(),vname.c_str());
      else
	sprintf(cmdBuff,"if ~isempty(find(isnan(%s))), error('** NaN DETECTED\\n   Var=%s'); end",vname.c_str(),vname.c_str());
      if (mexEvalString(cmdBuff)!=0)
	throw MatIFException("MatlabDebug::setMatrix",file,line);
    }
#endif
  }
}

void MatlabDebug::setDoubleArray(const string& var,const double* arr,int size,
				 const char* fname,const char* file,int line)
{
  int i;
  bool storeOnly=(fname!=0);

  if (status!=0) {
    string vname(commonPrefix+var);
#if !defined(MATLABDEBUG_USEMEX) || !defined(MATLAB_MEX)
    if (status==3) throw NotImplemException(EXCEPT_MSG(""));
    if (!storeOnly && size<=VEC_TEXT_LIMIT) {
      os->precision(20); os->setf(ios::scientific);
      (*os) << vname << "=[";
      for (i=0; i<size; i++)
	(*os) << arr[i] << ';' << endl;
      (*os) << "];" << endl;
    } else {
      if (!storeOnly) drawName();
      ofstream bin(storeOnly?fname:nameBuff);
      if (!bin)
	throw FileUtilsException("Cannot open MatlabDebug binary object file");
      StVector msk; msk.reassign((double*) arr,size);
      msk.save(bin);
      CLOSESTR(bin);
      if (!storeOnly)
	(*os) << vname << "=loadbasevector('" << nameBuff << "',6);" << endl;
    }
#else
    // MEX based mode
    MatlabMatrix mvec(StMatrix::mask((double*) arr,size,1,size));
    if (putVariable("caller",var,mvec.matobj())!=0)
      throw MatIFException(EXCEPT_MSG(""));
    if (status==3) {
      putFileAndLine(file,line);
      if (file!=0)
	sprintf(cmdBuff,"if ~isempty(find(isnan(%s))), error(sprintf('** NaN DETECTED\\n   Var=%s\\n   File: %%s (line %%d)',MATDEBUGXYZ_file,MATDEBUGXYZ_line)); end",vname.c_str(),vname.c_str());
      else
	sprintf(cmdBuff,"if ~isempty(find(isnan(%s))), error('** NaN DETECTED\\n   Var=%s'); end",vname.c_str(),vname.c_str());
      if (mexEvalString(cmdBuff)!=0)
	throw MatIFException("MatlabDebug::setMatrix",file,line);
    }
#endif
  }
}

void MatlabDebug::setCholFact(const string& var,const StMatrix& fact,
			      const char* file,int line)
{
  if (status!=0) {
    string vname(commonPrefix+var);
    int n=fact.rows();
    if (fact.cols()!=n)
      throw WrongDimensionException(EXCEPT_MSG(""));
#if !defined(MATLABDEBUG_USEMEX) || !defined(MATLAB_MEX)
    if (status==3) throw NotImplemException(EXCEPT_MSG(""));
    drawName();
    ofstream bin(nameBuff);
    if (!bin) throw FileUtilsException("Cannot open MatlabDebug binary object file");
    fact.save(bin);
    CLOSESTR(bin);
    (*os) << vname << "=loadbasematrix('" << nameBuff << "');" << endl;
    switch (fact.getStrctPatt()) {
    case MatStrct::lowNDg:
      (*os) << vname << "(1:" << n+1 << ":" << n*n << ")=1; ";
    case MatStrct::lower:
    case MatStrct::normal:
      (*os) << vname << "=tril(" << vname << ",0);" << endl;
      break;
    case MatStrct::uppNDg:
      (*os) << vname << "(1:" << n+1 << ":" << n*n << ")=1; ";
    case MatStrct::upper:
      (*os) << vname << "=triu(" << vname << ",0);" << endl;
      break;
    }
#else
    // MEX based mode
    MatlabMatrix mmat(fact);
    if (putVariable("caller",var,mmat.matobj())!=0)
      throw MatIFException(EXCEPT_MSG(""));
    char ndgpart[100]; ndgpart[0]=0;
    switch (fact.getStrctPatt()) {
    case MatStrct::lowNDg:
      sprintf(ndgpart,"%s(1:%d:%d)=1; ",vname.c_str(),n+1,n*n);
    case MatStrct::lower:
    case MatStrct::normal:
      sprintf(cmdBuff,"%s%s=tril(%s,0)",ndgpart,vname.c_str(),vname.c_str());
      break;
    case MatStrct::uppNDg:
      sprintf(ndgpart,"%s(1:%d:%d)=1; ",vname.c_str(),n+1,n*n);
    case MatStrct::upper:
      sprintf(cmdBuff,"%s%s=triu(%s,0)",ndgpart,vname.c_str(),vname.c_str());
      break;
    }
    if (mexEvalString(cmdBuff)!=0)
      throw MatIFException(EXCEPT_MSG(""));
    if (status==3) {
      putFileAndLine(file,line);
      string vname(commonPrefix+var);
      if (file!=0)
	sprintf(cmdBuff,"if ~isempty(find(isnan(%s))), error(sprintf('** NaN DETECTED\\n   Var=%s\\n   File: %%s (line %%d)',MATDEBUGXYZ_file,MATDEBUGXYZ_line)); end",vname.c_str(),vname.c_str());
      else
	sprintf(cmdBuff,"if ~isempty(find(isnan(%s))), error('** NaN DETECTED\\n   Var=%s'); end",vname.c_str(),vname.c_str());
      if (mexEvalString(cmdBuff)!=0)
	throw MatIFException("MatlabDebug::setMatrix",file,line);
    }
#endif
  }
}

void MatlabDebug::setSymm(const string& var,const StMatrix& mat,
			  const char* file,int line)
{
  if (status!=0) {
    if (mat.getStrctPatt()==MatStrct::normal)
      throw WrongStatusException(EXCEPT_MSG(""));
    int n=mat.rows();
    if (n!=mat.cols())
      throw WrongDimensionException(EXCEPT_MSG(""));
    string vname(commonPrefix+var);
#if !defined(MATLABDEBUG_USEMEX) || !defined(MATLAB_MEX)
    if (status==3) throw NotImplemException(EXCEPT_MSG(""));
    drawName();
    ofstream bin(nameBuff);
    if (!bin) throw FileUtilsException(EXCEPT_MSG(""));
    mat.save(bin);
    CLOSESTR(bin);
    (*os) << vname << "=loadbasematrix('" << nameBuff << "');" << endl;
#else
    // MEX based mode
    MatlabMatrix mmat(mat);
    if (putVariable("caller",var,mmat.matobj())!=0)
      throw MatIFException(EXCEPT_MSG(""));
#endif
    char ndgpart[100]; ndgpart[0]=0;
    // Assemble command string
    switch (mat.getStrctPatt()) {
    case MatStrct::lowNDg:
      sprintf(ndgpart,"%s((1:%d:%d)')=1; ",vname.c_str(),n+1,n*n);
    case MatStrct::lower:
      sprintf(cmdBuff,"%s%s=tril(%s,0)+tril(%s,-1)';",ndgpart,vname.c_str(),
	      vname.c_str(),vname.c_str());
      break;
    case MatStrct::uppNDg:
      sprintf(ndgpart,"%s((1:%d:%d)')=1; ",vname.c_str(),n+1,n*n);
    case MatStrct::upper:
      sprintf(cmdBuff,"%s%s=triu(%s,0)+triu(%s,1)';",ndgpart,vname.c_str(),
	      vname.c_str(),vname.c_str());
      break;
    }
#if !defined(MATLABDEBUG_USEMEX) || !defined(MATLAB_MEX)
    (*os) << cmdBuff; (*os) << endl;
#else
    if (mexEvalString(cmdBuff)!=0)
      throw MatIFException(EXCEPT_MSG(""));
    if (status==3) {
      putFileAndLine(file,line);
      if (file!=0)
	sprintf(cmdBuff,"if ~isempty(find(isnan(%s))), error(sprintf('** NaN DETECTED\\n   Var=%s\\n   File: %%s (line %%d)',MATDEBUGXYZ_file,MATDEBUGXYZ_line)); end",vname.c_str(),vname.c_str());
      else
	sprintf(cmdBuff,"if ~isempty(find(isnan(%s))), error('** NaN DETECTED\\n   Var=%s'); end",vname.c_str(),vname.c_str());
      if (mexEvalString(cmdBuff)!=0)
	throw MatIFException("MatlabDebug::setMatrix",file,line);
    }
#endif
  }
}

void MatlabDebug::setScalar(const string& var,double val,const char* file,int line)
{
  if (status!=0) {
    string vname(commonPrefix+var);
#if !defined(MATLABDEBUG_USEMEX) || !defined(MATLAB_MEX)
    if (status==3) throw NotImplemException(EXCEPT_MSG(""));
    os->precision(20); os->setf(ios::scientific);
    (*os) << vname << '=' << val << ";" << endl;
#else
    // MEX based mode
    MatlabMatrix mmat(1,1); *(mmat.buff())=val;
    if (putVariable("caller",var,mmat.matobj())!=0)
      throw MatIFException(EXCEPT_MSG(""));
    if (status==3) {
      putFileAndLine(file,line);
      if (file!=0)
	sprintf(cmdBuff,"if isnan(%s), error(sprintf('** NaN DETECTED\\n   Var=%s\\n   File: %%s (line %%d)',MATDEBUGXYZ_file,MATDEBUGXYZ_line)); end",vname.c_str(),vname.c_str());
      else
	sprintf(cmdBuff,"if isnan(%s), error('** NaN DETECTED\\n   Var=%s'); end",vname.c_str(),vname.c_str());
      if (mexEvalString(cmdBuff)!=0)
	throw MatIFException("MatlabDebug::setMatrix",file,line);
    }
#endif
  }
}

void MatlabDebug::putMatrix(const string& var,StMatrix& mat,const char* file,
			    int line)
{
  if (status==1 || status==2) {
#if !defined(MATLABDEBUG_USEMEX) || !defined(MATLAB_MEX)
    throw NotImplemException(EXCEPT_MSG(""));
#else
    // MEX based mode
    const mxArray* mmmat=getVariable("caller",var);
    MatlabMatrix mmat(mmmat,true);
    StMatrix msk; msk.maskMatlab(mmat);
    mat.setStrctPatt(MatStrct::normal);
    mat=msk;
#endif
  }
}

void MatlabDebug::print(const string& text)
{
  if (status==1) {
#if !defined(MATLABDEBUG_USEMEX) || !defined(MATLAB_MEX)
    (*os) << "fprintf('" << text << "\\n');" << endl;
#else
    // MEX based mode
    mexPrintf("%s\n",text.c_str());
#endif
  }
}

void MatlabDebug::printVar(const string& var)
{
  if (status==1) {
    string vname(commonPrefix+var);
#if !defined(MATLABDEBUG_USEMEX) || !defined(MATLAB_MEX)
    (*os) << "fprintf('" << vname << " =\\n'); disp(" << vname << ");" << endl;
#else
    mexPrintf("'%s' =\n",vname.c_str());
    sprintf(cmdBuff,"disp(%s)",vname.c_str());
    if (mexEvalString(cmdBuff)!=0)
      throw MatIFException(EXCEPT_MSG(""));
#endif
  }
}

void MatlabDebug::output(const string& txt,const char* file,int line)
{
  if (status!=0 && status!=3) {
    string expr(insertCommonPrefix(txt));
#if !defined(MATLABDEBUG_USEMEX) || !defined(MATLAB_MEX)
    (*os) << expr << endl;
#else
    if (mexEvalString(expr.c_str())!=0) {
      if (file!=0)
	throw MatIFException("MatlabDebug::output error",file,line);
      else
	throw MatIFException(EXCEPT_MSG(""));
    }
#endif
  }
}

void MatlabDebug::readMatrix(StMatrix& mat,const char* fname)
{
  if (status!=0) {
#if !defined(MATLABDEBUG_USEMEX) || !defined(MATLAB_MEX)
    ifstream bin(fname);
    if (!bin) throw FileUtilsException("Cannot open MatlabDebug binary object file");
    mat.load(bin);
    CLOSESTR(bin);
#else
    throw NotImplemException(EXCEPT_MSG(""));
#endif
  }
}

void MatlabDebug::doCompare(const StMatrix& obj,const string& expr,
			    const string& varName,bool copy,
			    const char* tempName,const char* file,int line)
{
  if (status!=0) {
    if (expr.size()>0) output(expr,file,line);
    if (status!=2) {
      string tName((tempName!=0)?tempName:"tmat");
      setMatrix(tName,obj,0,file,line);
      if (status==1) {
	compare(tName,varName,file,line);
	if (copy) {
	  string cmd("%");
	  cmd+=varName; cmd+="=%"; cmd+=tName; cmd+=";";
	  output(cmd,file,line);
	}
      }
    }
  }
}

void MatlabDebug::doCompare(const StVector& obj,const string& expr,
			    const string& varName,bool copy,
			    const char* tempName,const char* file,int line)
{
  if (status!=0) {
    if (expr.size()>0) output(expr,file,line);
    if (status!=2) {
      string tName((tempName!=0)?tempName:"tvec");
      setVector(tName,obj,0,file,line);
      if (status==1) {
	compare(tName,varName,file,line);
	if (copy) {
	  string cmd("%");
	  cmd+=varName; cmd+="=%"; cmd+=tName; cmd+=";";
	  output(cmd,file,line);
	}
      }
    }
  }
}

void MatlabDebug::doCompare(const BaseVector<int>& obj,const string& expr,
			    const string& varName,bool copy,
			    const char* tempName,const char* file,int line)
{
  if (status!=0) {
    if (expr.size()>0) output(expr,file,line);
    if (status!=2) {
      string tName((tempName!=0)?tempName:"tvec");
      setIntVector(tName,obj,0,false,file,line);
      if (status==1) {
	compare(tName,varName,file,line);
	if (copy) {
	  string cmd("%");
	  cmd+=varName; cmd+="=%"; cmd+=tName; cmd+=";";
	  output(cmd,file,line);
	}
      }
    }
  }
}

void MatlabDebug::doCompare(double obj,const string& expr,
			    const string& varName,bool copy,
			    const char* tempName,const char* file,int line)
{
  if (status!=0) {
    if (expr.size()>0) output(expr,file,line);
    if (status!=2) {
      string tName((tempName!=0)?tempName:"tscal");
      setScalar(tName,obj);
      if (status==1) {
	compare(tName,varName,file,line);
	if (copy) {
	  string cmd("%");
	  cmd+=varName; cmd+="=%"; cmd+=tName; cmd+=";";
	  output(cmd,file,line);
	}
      }
    }
  }
}

// Internal methods

#if defined(MATLAB_MEX) && defined(MATLABDEBUG_USEMEX)
void MatlabDebug::putFileAndLine(const char* file,int line)
{
  if (file!=0) {
    mxArray* str=mxCreateString(file);
    if (putVariable("caller","MATDEBUGXYZ_file",str,true)!=0)
      throw MatIFException(EXCEPT_MSG(""));
    mxDestroyArray(str);
#ifdef MATLAB_VER65
    mxArray* val=mxCreateDoubleScalar((double) line);
#else
    mxArray* val=mxCreateScalarDouble((double) line);
#endif
    if (putVariable("caller","MATDEBUGXYZ_line",val,true)!=0)
      throw MatIFException(EXCEPT_MSG(""));
    mxDestroyArray(val);
  }
}

void MatlabDebug::putDiffThres()
{
#ifdef MATLAB_VER65
  mxArray* val=mxCreateDoubleScalar(diffThres);
#else
  mxArray* val=mxCreateScalarDouble(diffThres);
#endif
  if (putVariable("caller","MATDEBUGXYZ_diffThres",val,true)!=0)
    throw MatIFException(EXCEPT_MSG(""));
  mxDestroyArray(val);
}

int MatlabDebug::putVariable(const char* works,const string& var,mxArray* arr,
			     bool noprefix)
{
  string vname(noprefix?var:(commonPrefix+var));
#ifdef MATLAB_VER65
  return mexPutVariable(works,vname.c_str(),arr);
#else
  // Works for Matlab <6.5
  mxSetName(arr,vname.c_str()); // set name of 'arr'
  return mexPutArray(arr,works);
#endif
}

mxArray* MatlabDebug::getVariable(const char* works,const string& var,
				  bool noprefix)
{
  string vname(noprefix?var:(commonPrefix+var));
#ifdef MATLAB_VER65
  return mexGetVariable(works,vname.c_str());
#else
  // Works for Matlab <6.5
  return mexGetArray(vname.c_str(),works);
#endif
}
#endif
#undef VEC_TEXT_LIMIT
