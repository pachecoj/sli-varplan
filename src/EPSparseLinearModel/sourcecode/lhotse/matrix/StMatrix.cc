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
 * Desc.:  Definition of class StMatrix
 * ------------------------------------------------------------------- */

#include "lhotse/matrix/StMatrix.h"
#include "lhotse/matrix/cblas_for_cpp.h"
#include "lhotse/FileUtils.h"
#include "lhotse/NumberFormats.h"
#include "lhotse/MachDep.h"
#ifdef MATLAB_MEX
#include "lhotse/matif/MatlabMatrix.h"
#endif
#ifdef MATLAB_DEBUG_OLD
#include "lhotse/MatlabDebug.h"
#endif

// Local macros
#define TSARG __FILE__,__LINE__

#define MASKVEC(hmsk,rngr,rngc) \
  (hmsk).changeRep(DYNCAST(StVector,subrvalInt((rngr),(rngc))))
#define MASKMVEC(hmsk,mat,rngr,rngc) \
  (hmsk).changeRep(DYNCAST(StVector,(mat).subrvalInt((rngr),(rngc))))
#define MASKMAT(hmsk,rngr,rngc) \
  (hmsk).changeRep(DYNCAST(StMatrix,subrvalInt((rngr),(rngc))))
#define MASKMMAT(hmsk,mat,rngr,rngc) \
  (hmsk).changeRep(DYNCAST(StMatrix,(mat).subrvalInt((rngr),(rngc))))
#define CHECKMAT(mat,r,c) do { (mat).checkTS(TSARG); \
  if (((r)!=-1 && (mat).rows()!=(r)) || ((c)!=-1 && (mat).cols()!=(c))) \
    throw WrongDimensionException(EXCEPT_MSG("")); } while (0)
#define CHECKVEC(vec,s) do { (vec).checkTS(TSARG); \
  if ((s)!=-1 && (vec).size()!=(s)) \
    throw WrongDimensionException(EXCEPT_MSG("")); } while (0)
#define COLPOS(i) ((colInd==0)?(i):colInd[i])
#define ROWPOS(i) ((rowInd==0)?(i):rowInd[i])

// Machine precision for type double: 2^(-52) for 52 bits fraction
#define DOUBLE_MACHPRECISION (2.220446049250313e-16)

/*
 * Declaration of LAPACK routines
 * Just a temporary solution here! If many LAPACK routines are used, design
 * a proper wrapper with support of Fortran calling convention.
 * NOTE: Such a C interface exists for BLAS, and we use it (see 'FastUtils'),
 * but there is no standard C interface to LAPACK, or to LINPACK.
 */

extern "C" int
F77_FUNC(dsyevr,DSYEVR) (const char* jobz,const char* range,const char* uplo,
			 const int* n,double* a,
			 const int* lda,double* vl,double* vu,int* il,int *iu,
			 double* abstol,int* m,double* w,double* z,int* ldz,
			 int* isuppz,double* work,int* lwork,int* iwork,
			 int* liwork,int* info);
extern "C" int
F77_FUNC(dpotrf,DPOTRF) (const char* uplo,const int* n,double* a,
			 const int* lda,int* info);
extern "C" int
F77_FUNC(dpotri,DPOTRI) (const char* uplo,const int* n,double* a,
			 const int* lda,int* info);
extern "C" int
F77_FUNC(dgetrf,DGETRF) (const int* m,const int* n,double* a,const int* lda,
			 int* ipiv,int* info);
extern "C" int
F77_FUNC(dgetrs,DGETRS) (const char* trans,const int* n,const int* nrhs,
			 const double* a,const int* lda,const int* ipiv,
			 double* b,const int* ldb,int* info);

// LINPACK routines (external Fortran code)

#ifdef HAVE_FORTRAN
extern "C" void
F77_FUNC(dchex,DCHEX) (double* r,int* ldr,int* p,int* k,int* l,double* z,
		       int* ldz,int* nz,double* c,double* s,int* job);
#endif

//extern "C" void
//F77_FUNC(dchud,DCHUD) (double* r,int* ldr,int* p,double* x,double* z,int* ldz,
//		       int* nz,double* y,double* rho,double* c,double* s);

//extern "C" void
//F77_FUNC(dchdd,DCHDD) (double* r,int* ldr,int* p,double* x,double* z,int* ldz,
//		       int* nz,double* y,double* rho,double* c,double* s,
//		       int* info);

//BEGINNS(matrix)
  // Public methods

#ifdef MATLAB_MEX
  void StMatrix::maskMatlab(MatlabMatrix& mat,const Range& rngR,
			    const Range& rngC)
  {
    if (rngR.checkRange(mat.m()) || rngC.checkRange(mat.n()))
      throw OutOfRangeException(EXCEPT_MSG(""));
    reassign(mat.buff(),mat.m(),
	     rngR.isOpen()?Range(rngR.getStart(),mat.m()-1):rngR,
	     rngC.isOpen()?Range(rngC.getStart(),mat.n()-1):rngC);
  }
#endif

  ifstream& StMatrix::load(ifstream& is,bool noTag)
  {
    int i=0;

    checkTS(TSARG);
    if (!noTag) {
      ArrayHandle<string> tags(2);
      tags[0]="@BaseMatrix";
      tags[1]="@StMatrix";
      FileUtils::loadHeaderMulti(is,tags,i,true); // do not load FFV number
    }
    if (i==0) {
      // Current format
      return BaseMatrix<double>::load(is,true);
    } else {
      // Old format
      NumberFormats<int>::load(is,&i,1,1,0,4); // FFVN
      if (i<0 || i>1) throw FileFormatException("Unknown FF version number");
      if (i==1) {
	// @StMatrix tag has been used briefly for this type (FFVN 1)
	NumberFormats<int>::load(is,&i,1,1,0,4); // type code
	if (i!=TypeCodeConsts::typeDouble)
	  throw FileFormatException("StMatrix: Type code must be double");
	return loadInt(is);
      } else {
	// Old row-major format
	int newm,newn;
	NumberFormats<int>::load(is,&newm,1,1,0,4);
	if (newm<0)
	  throw FileFormatException("Error in file format");
	NumberFormats<int>::load(is,&newn,1,1,0,4);
	if (newn<0)
	  throw FileFormatException("Error in file format");
	ensureCapacity(newm,newn);
	// Read row by row
	StVector vec(n);
	ArrayHandle<double> vecBuff=vec.getFlatBuff();
	Handle<StVector> msk;
	for (int i=0; i<m; i++) {
	  NumberFormats<double>::load(is,vecBuff,n);
	  MASKVEC(msk,i,Range(0,n-1));
	  *msk=vec;
	}
      }
    }
  }

  void StMatrix::fillStrct(int rows,int cols,double a)
  {
    checkTS(TSARG);
    updateStrctPatt();
    if (rows!=m || cols!=n || strpatt==MatStrct::normal) {
      fill(rows,cols,a);
      updateStrctPatt();
    } else {
      Handle<StVector> msk;
      int i,off;
      if (strpatt==MatStrct::lower || strpatt==MatStrct::lowNDg) {
	off=(strpatt==MatStrct::lower)?0:1;
	for (i=0; i<n-off; i++) {
	  MASKVEC(msk,Range(i+off,n-1),i);
	  msk->fill(n-i-off,a);
	}
      } else {
	off=(strpatt==MatStrct::upper)?0:1;
	for (i=off; i<n; i++) {
	  MASKVEC(msk,Range(0,i-off),i);
	  msk->fill(i+1-off,a);
	}
      }
    }
  }

  // For indexed mask version: see 'sumProd' comments
  void StMatrix::mulVec(StVector& avec,const StVector& bvec,bool trans,
			double alpha,double beta) const
  {
    updateStrctPatt();
    if (strpatt!=MatStrct::normal && strpatt!=MatStrct::upper &&
	strpatt!=MatStrct::lower)
      throw WrongStatusException(EXCEPT_MSG("strpatt"));
    checkTS(TSARG);
    avec.checkTS(TSARG);
    bvec.checkTS(TSARG);
    int opn=trans?m:n,opm=trans?n:m;
    if (opn!=bvec.size()) throw DimMismatchException(EXCEPT_MSG("bvec"));
    if (beta!=0.0) {
      if (opm!=avec.size()) throw DimMismatchException(EXCEPT_MSG("avec"));
    } else {
      avec.fill(opm,0.0);
    }
    BaseLinVec<double> tempA,tempB;
    avec.getLinVec(tempA,true);
    bvec.getLinVec(tempB,false);
    if (!isMaskInd() || strpatt!=MatStrct::normal) {
      BaseLinMat<double> tempM;
      getLinMat(tempM,false,strpatt);
      FastUtils::matVec(tempA,tempM,trans,tempB,alpha,beta);
    } else {
      if (beta!=0.0)
	FastUtils::smul(tempA.getBuff(),beta,opm,tempA.step);
      int i;
      if (!trans) {
	for (i=0; i<opn; i++)
	  FastUtils::addsmul(tempA.getBuff(),buff+stride*COLPOS(i),
			     tempB[i]*alpha,opm,tempA.step,1,0,rowInd.p());
      } else {
	for (i=0; i<opm; i++)
	  tempA[i]+=
	    alpha*FastUtils::inner(buff+stride*COLPOS(i),tempB.getBuff(),opn,
				   1,tempB.step,rowInd.p());
      }
    }
  }

  void StMatrix::triMulVec(StVector& avec,bool trans) const
  {
    updateStrctPatt();
    if (strpatt==MatStrct::normal)
      throw WrongStatusException(EXCEPT_MSG("Need triangular matrix pattern"));
    // Matrix square here
    checkTS(TSARG); avec.checkTS(TSARG);
    if (n!=avec.size()) throw DimMismatchException(EXCEPT_MSG("avec"));
    BaseLinMat<double> tempM;
    getLinMat(tempM,false,strpatt);
    BaseLinVec<double> tempA;
    avec.getLinVec(tempA,true);
    FastUtils::trimatVec(tempA,tempM,trans);
  }

  void StMatrix::backsubst(StVector& avec,bool trans) const
  {
    updateStrctPatt();
    if (strpatt==MatStrct::normal)
      throw WrongStatusException(EXCEPT_MSG("strpatt"));
    // Matrix square here
    checkTS(TSARG); avec.checkTS(TSARG);
    if (n!=avec.size()) throw DimMismatchException(EXCEPT_MSG("avec"));
    BaseLinMat<double> tempM;
    getLinMat(tempM,false,strpatt);
    BaseLinVec<double> tempA;
    avec.getLinVec(tempA,true);
    FastUtils::backsubst(tempA,tempM,trans);
  }

  void StMatrix::backsubst(StVector& avec,const StVector& dvec,bool trans)
    const
  {
    updateStrctPatt();
    if (strpatt==MatStrct::normal)
      throw WrongStatusException(EXCEPT_MSG("strpatt"));
    // Matrix square here
    checkTS(TSARG); avec.checkTS(TSARG);
    if (n!=avec.size() || n!=dvec.size())
      throw DimMismatchException(EXCEPT_MSG("avec,dvec"));
    BaseLinMat<double> tempM;
    getLinMat(tempM,false,strpatt);
    BaseLinVec<double> tempA,tempD;
    avec.getLinVec(tempA,true);
    dvec.getLinVec(tempD,false);
    if (!trans) {
      FastUtils::div(tempA.getBuff(),tempD.getBuff(),n,tempA.step,tempD.step);
      FastUtils::backsubst(tempA,tempM,false);
    } else {
      FastUtils::backsubst(tempA,tempM,true);
      FastUtils::div(tempA.getBuff(),tempD.getBuff(),n,tempA.step,tempD.step);
    }
  }

  void StMatrix::rankOne(const StVector& bvec,double alpha)
  {
    updateStrctPatt();
    if (strpatt==MatStrct::normal) {
      if (m!=n) throw WrongDimensionException(EXCEPT_MSG("Must be square"));
    } else if (strpatt!=MatStrct::upper && strpatt!=MatStrct::lower)
      throw WrongStatusException(EXCEPT_MSG("strpatt"));
    checkTS(TSARG); bvec.checkTS(TSARG);
    if (n!=bvec.size()) throw DimMismatchException(EXCEPT_MSG("bvec"));
    BaseLinMat<double> tempM;
    getLinMat(tempM,true,strpatt);
    BaseLinVec<double> tempB;
    bvec.getLinVec(tempB,false);
    FastUtils::symRankOne(tempM,tempB,alpha);
  }

  void StMatrix::rankOne(const StVector& bvec,double alpha,const StVector& cvec)
  {
    updateStrctPatt();
    if (strpatt!=MatStrct::normal)
      throw WrongStatusException(EXCEPT_MSG("strpatt"));
    checkTS(TSARG); bvec.checkTS(TSARG); cvec.checkTS(TSARG);
    if (m!=bvec.size() || n!=cvec.size())
      throw DimMismatchException(EXCEPT_MSG("bvec,cvec"));
    BaseLinMat<double> tempM;
    getLinMat(tempM,true,strpatt);
    BaseLinVec<double> tempB,tempC;
    bvec.getLinVec(tempB,false);
    cvec.getLinVec(tempC,false);
    FastUtils::rankOne(tempM,tempB,tempC,alpha);
  }

  void StMatrix::addVec(const StVector& avec,double alpha,bool col)
  {
    int i;
    Handle<StVector> msk;
    checkTS(TSARG); avec.checkTS(TSARG);
    updateStrctPatt();
    if (strpatt!=MatStrct::normal)
      throw WrongStatusException(EXCEPT_MSG("strpatt"));
    if ((col && avec.size()!=m) || (!col && avec.size()!=n))
      throw WrongDimensionException(EXCEPT_MSG(""));
    if (col) {
      for (i=0; i<n; i++) {
	msk.changeRep(DYNCAST(StVector,subrvalInt(RangeFull::get(),i)));
	msk->addprod(alpha,avec);
      }
    } else {
      for (i=0; i<m; i++) {
	msk.changeRep(DYNCAST(StVector,subrvalInt(i,RangeFull::get())));
	msk->addprod(alpha,avec);
      }
    }
  }

  void StMatrix::symRankTwo(const StVector& bvec,const StVector& cvec,
			    double alpha)
  {
    updateStrctPatt();
    if (strpatt!=MatStrct::upper && strpatt!=MatStrct::lower)
      throw WrongStatusException(EXCEPT_MSG("strpatt"));
    checkTS(TSARG); bvec.checkTS(TSARG); cvec.checkTS(TSARG);
    if (n!=bvec.size() || n!=cvec.size())
      throw DimMismatchException(EXCEPT_MSG("bvec,cvec"));
    BaseLinMat<double> tempM;
    getLinMat(tempM,true,strpatt);
    BaseLinVec<double> tempB,tempC;
    bvec.getLinVec(tempB,false);
    cvec.getLinVec(tempC,false);
    FastUtils::symRankTwo(tempM,tempB,tempC,alpha);
  }

  void StMatrix::backsubst(StMatrix& amat,bool left,bool trans,double alpha)
    const
  {
    updateStrctPatt(); amat.updateStrctPatt();
    if (strpatt==MatStrct::normal)
      throw WrongStatusException(EXCEPT_MSG("strpatt"));
    if (amat.strpatt!=MatStrct::normal)
      throw WrongStatusException(EXCEPT_MSG("amat.strpatt"));
    // Matrix square here
    checkTS(TSARG); amat.checkTS(TSARG);
    if ((left && n!=amat.m) || (!left && n!=amat.n))
      throw DimMismatchException(EXCEPT_MSG("amat"));
    BaseLinMat<double> tempM,tempA;
    getLinMat(tempM,false,strpatt);
    amat.getLinMat(tempA,true);
    FastUtils::backsubstMat(tempA,tempM,trans,left,alpha);
  }

  void StMatrix::backsubst(StMatrix& amat,const StVector& dvec,bool left,
			   bool trans) const
  {
    updateStrctPatt(); amat.updateStrctPatt();
    if (strpatt==MatStrct::normal)
      throw WrongStatusException(EXCEPT_MSG("strpatt"));
    if (amat.strpatt!=MatStrct::normal)
      throw WrongStatusException(EXCEPT_MSG("amat.strpatt"));
    // Matrix square here
    checkTS(TSARG); amat.checkTS(TSARG); dvec.checkTS(TSARG);
    if ((left && n!=amat.m) || (!left && n!=amat.n))
      throw DimMismatchException("amat");
    if (n!=dvec.size()) throw DimMismatchException("dvec");
    BaseLinMat<double> tempM,tempA;
    BaseLinVec<double> tempD;
    getLinMat(tempM,false,strpatt);
    amat.getLinMat(tempA,true,amat.strpatt);
    dvec.getLinVec(tempD,false);
    if (left^trans) {
      FastUtils::mulDiag(tempA,tempD,tempA,true,left);
      FastUtils::backsubstMat(tempA,tempM,trans,left);
    } else {
      FastUtils::backsubstMat(tempA,tempM,trans,left);
      FastUtils::mulDiag(tempA,tempD,tempA,true,left);
    }
  }

  void StMatrix::mul(const StMatrix& amat,bool atrans,const StMatrix& bmat,
		     bool btrans,double alpha,double beta)
  {
    if (this==&amat || this==&bmat)
      throw EqualArgsException(EXCEPT_MSG(""));
    updateStrctPatt(); amat.updateStrctPatt(); bmat.updateStrctPatt();
    if (beta!=0.0 && strpatt!=MatStrct::normal)
      throw WrongStatusException(EXCEPT_MSG("strpatt"));
    if (amat.strpatt==MatStrct::uppNDg || amat.strpatt==MatStrct::lowNDg)
      throw WrongStatusException(EXCEPT_MSG("amat.strpatt"));
    if (bmat.strpatt==MatStrct::uppNDg || bmat.strpatt==MatStrct::lowNDg)
      throw WrongStatusException(EXCEPT_MSG("bmat.strpatt"));
    checkTS(TSARG); amat.checkTS(TSARG); bmat.checkTS(TSARG);
    // Sizes
    int needm,needn;
    if (atrans) {
      needm=amat.n;
      if (btrans) {
	needn=bmat.m;
	if (bmat.n!=amat.m) throw DimMismatchException(EXCEPT_MSG(""));
      } else {
	needn=bmat.n;
	if (bmat.m!=amat.m) throw DimMismatchException(EXCEPT_MSG(""));
      }
    } else {
      needm=amat.m;
      if (btrans) {
	needn=bmat.m;
	if (bmat.n!=amat.n) throw DimMismatchException(EXCEPT_MSG(""));
      } else {
	needn=bmat.n;
	if (bmat.m!=amat.n) throw DimMismatchException(EXCEPT_MSG(""));
      }
    }
    if (beta!=0.0) {
      if (m!=needm || n!=needn)
	throw DimMismatchException(EXCEPT_MSG(""));
    } else {
      zeros(needm,needn);
      strpatt=MatStrct::normal;
    }
    // Select BLAS-III routine
    BaseLinMat<double> tempM,tempA,tempB;
    getLinMat(tempM,true); // always 'normal', write back
    Handle<StMatrix> tempMat;
    if (amat.strpatt==MatStrct::normal) {
      bmat.getLinMat(tempB,false,bmat.strpatt);
      if (bmat.strpatt==MatStrct::normal) {
	// Both normal
	amat.getLinMat(tempA,false);
	FastUtils::matMat(tempM,tempA,atrans,tempB,btrans,alpha,beta);
      } else {
	// A normal, B symm.
	if (!atrans) {
	  amat.getLinMat(tempA,false);
	} else {
	  // Have to transpose A
	  tempMat.changeRep(new StMatrix());
	  tempMat->trans(amat);
	  tempMat->getLinMat(tempA,false);
	}
	FastUtils::symmatMat(tempM,tempA,tempB,alpha,beta);
      }
    } else {
      amat.getLinMat(tempA,false,amat.strpatt);
      if (bmat.strpatt==MatStrct::normal) {
	// A symm., B normal
	if (!btrans) {
	  bmat.getLinMat(tempB,false);
	} else {
	  // Have to transpose B
	  tempMat.changeRep(new StMatrix());
	  tempMat->trans(bmat);
	  tempMat->getLinMat(tempB,false);
	}
	FastUtils::symmatMat(tempM,tempA,tempB,alpha,beta);
      } else {
	// Both symmetric -> copy B
	tempMat.changeRep(new StMatrix());
	tempMat->assignInt(bmat,bmat.strpatt); // copy
	tempMat->makeSymm(bmat.strpatt==MatStrct::lower); // make symm.
	tempMat->getLinMat(tempB,false); // as normal
	FastUtils::symmatMat(tempM,tempA,tempB,alpha,beta);
      }
    }
  }

  void StMatrix::symMul(const StMatrix& amat,bool atrans,double alpha,
			double beta)
  {
    amat.updateStrctPatt();
    if (beta!=0.0) {
      updateStrctPatt();
      if (strpatt!=MatStrct::upper && strpatt!=MatStrct::lower)
	throw WrongStatusException(EXCEPT_MSG("strpatt"));
    }
    if (amat.strpatt!=MatStrct::normal)
      throw WrongStatusException(EXCEPT_MSG("amat.strpatt"));
    checkTS(TSARG); amat.checkTS(TSARG);
    // Sizes
    int needn=atrans?amat.n:amat.m;
    if (beta!=0.0) {
      if (n!=needn) throw WrongDimensionException(EXCEPT_MSG("amat"));
    } else {
      ensureCapacity(needn,needn);
      if (strpatt!=MatStrct::upper) strpatt=MatStrct::lower;
      fillStrct(needn,needn);
    }
    // BLAS-III routine
    BaseLinMat<double> tempM,tempA;
    getLinMat(tempM,true,strpatt);
    amat.getLinMat(tempA,false,MatStrct::normal);
    FastUtils::symRankK(tempM,tempA,atrans,alpha,beta);
  }

  void StMatrix::symMulDiag(const StMatrix& amat,bool atrans,
			    const StVector& dvec,double alpha,double beta)
  {
    int i,j;
    double temp;
    bool increm=(beta!=0.0);

    amat.updateStrctPatt();
    if (increm) {
      updateStrctPatt();
      if (strpatt!=MatStrct::upper && strpatt!=MatStrct::lower)
	throw WrongStatusException(EXCEPT_MSG("strpatt"));
    }
    if (amat.strpatt!=MatStrct::normal)
      throw WrongStatusException(EXCEPT_MSG("amat.strpatt"));
    checkTS(TSARG); amat.checkTS(TSARG);
    // Sizes
    if (dvec.size()!=(atrans?amat.m:amat.n))
      throw DimMismatchException(EXCEPT_MSG(""));
    int needn=atrans?amat.n:amat.m;
    if (increm) {
      if (n!=needn) throw WrongDimensionException(EXCEPT_MSG(""));
    } else {
      ensureSize(needn,needn);
      if (strpatt!=MatStrct::upper) strpatt=MatStrct::lower;
    }
    // Loop
    BaseLinMat<double> linA;
    amat.getLinMat(linA,false);
    ArrayHandle<double> tempArr;
    ArrayHandle<double> dArr=dvec.getFlatBuff();
    const double* iP,*jP;
    if (atrans) {
      tempArr.changeRep(amat.m);
      for (i=0,iP=linA.buff.p(); i<n; i++,iP+=linA.stride) {
	FastUtils::prod(tempArr.p(),iP,dArr.p(),amat.m);
	FastUtils::smul(tempArr.p(),alpha,amat.m);
	for (j=i,jP=iP; j<n; j++,jP+=linA.stride) {
	  temp=FastUtils::inner(tempArr.p(),jP,amat.m);
	  if (strpatt==MatStrct::upper) {
	    if (increm) temp+=beta*get(i,j);
	    set(i,j,temp);
	  } else {
	    if (increm) temp+=beta*get(j,i);
	    set(j,i,temp);
	  }
	}
      }
    } else {
      tempArr.changeRep(amat.n);
      for (i=0,iP=linA.buff.p(); i<n; i++,iP++) {
	FastUtils::prod(tempArr.p(),iP,dArr.p(),amat.n,1,linA.stride);
	FastUtils::smul(tempArr.p(),alpha,amat.n);
	for (j=i,jP=iP; j<n; j++,jP++) {
	  temp=FastUtils::inner(tempArr.p(),jP,amat.n,1,linA.stride);
	  if (strpatt==MatStrct::upper) {
	    if (increm) temp+=beta*get(i,j);
	    set(i,j,temp);
	  } else {
	    if (increm) temp+=beta*get(j,i);
	    set(j,i,temp);
	  }
	}
      }
    }
  }

  void StMatrix::symMul(const StMatrix& amat,bool atrans,
			const StMatrix& bmat,
			double alpha,double beta)
  {
    updateStrctPatt(); amat.updateStrctPatt(); bmat.updateStrctPatt();
    if (beta!=0.0 && (strpatt!=MatStrct::upper && strpatt!=MatStrct::lower))
      throw WrongStatusException(EXCEPT_MSG("strpatt"));
    if (amat.strpatt!=MatStrct::normal || bmat.strpatt!=MatStrct::normal)
      throw WrongStatusException(EXCEPT_MSG("amat/bmat.strpatt"));
    checkTS(TSARG); amat.checkTS(TSARG); bmat.checkTS(TSARG);
    // Sizes
    if (amat.m!=bmat.m || amat.n!=bmat.n)
      throw DimMismatchException(EXCEPT_MSG(""));
    int needn=atrans?amat.n:amat.m;
    if (beta!=0.0) {
      if (n!=needn) throw WrongDimensionException(EXCEPT_MSG("amat"));
    } else {
      ensureCapacity(needn,needn);
      if (strpatt!=MatStrct::upper) strpatt=MatStrct::lower;
      fillStrct(needn,needn);
    }
    // BLAS-III routine
    BaseLinMat<double> tempM,tempA,tempB;
    getLinMat(tempM,true,strpatt);
    amat.getLinMat(tempA,false,amat.strpatt);
    bmat.getLinMat(tempB,false,bmat.strpatt);
    FastUtils::symRank2K(tempM,tempA,atrans,tempB,alpha,beta);
  }

  void StMatrix::triMul(StMatrix& amat,bool mtrans,bool left,double alpha)
    const
  {
    updateStrctPatt(); amat.updateStrctPatt();
    if (strpatt==MatStrct::normal)
      throw WrongStatusException(EXCEPT_MSG("strpatt"));
    if (amat.strpatt!=MatStrct::normal)
      throw WrongStatusException(EXCEPT_MSG("amat.strpatt"));
    checkTS(TSARG); amat.checkTS(TSARG);
    // Sizes
    if ((left && n!=amat.m) || (!left && n!=amat.n))
      throw DimMismatchException(EXCEPT_MSG(""));
    // BLAS-III routine
    BaseLinMat<double> tempM,tempA;
    getLinMat(tempM,false,strpatt);
    amat.getLinMat(tempA,true,amat.strpatt);
    FastUtils::trimatMat(tempA,tempM,mtrans,left,alpha);
  }

  void StMatrix::mulDiag(const StVector& avec,const StMatrix& bmat,
			 bool inv,bool incr)
  {
    bmat.updateStrctPatt(); updateStrctPatt();
    if (bmat.strpatt==MatStrct::uppNDg || bmat.strpatt==MatStrct::lowNDg)
      throw WrongStatusException(EXCEPT_MSG("bmat.strpatt"));
    if (avec.size()!=bmat.m)
      throw DimMismatchException(EXCEPT_MSG(""));
    checkTS(TSARG); avec.checkTS(TSARG); bmat.checkTS(TSARG);
    if (!incr) {
      ensureCapacity(bmat.m,bmat.n);
      strpatt=bmat.strpatt;
    } else if (m!=bmat.m || n!=bmat.n || strpatt!=bmat.strpatt)
      throw InvalidParameterException(EXCEPT_MSG(""));
    if (!isMaskInd() && !bmat.isMaskInd()) {
      BaseLinMat<double> tempM,tempB;
      getLinMat(tempM,true,strpatt);
      bmat.getLinMat(tempB,false,bmat.strpatt);
      BaseLinVec<double> tempA;
      avec.getLinVec(tempA,false);
      FastUtils::mulDiag(tempM,tempA,tempB,inv,true,incr);
    } else {
      // One of this, B is indexed mask
      int i;
      Handle<StVector> tMsk,bMsk;
      if (!incr) {
	if (!inv) {
	  for (i=0; i<n; i++) {
	    tMsk.changeRep(DYNCAST(StVector,sublvalInt(RangeFull::get(),i)));
	    bMsk.changeRep(DYNCAST(StVector,
				   bmat.subrvalInt(RangeFull::get(),i)));
	    tMsk->prod(*bMsk,avec);
	  }
	} else {
	  for (i=0; i<n; i++) {
	    tMsk.changeRep(DYNCAST(StVector,sublvalInt(RangeFull::get(),i)));
	    bMsk.changeRep(DYNCAST(StVector,
				   bmat.subrvalInt(RangeFull::get(),i)));
	    tMsk->div(*bMsk,avec);
	  }
	}
      } else {
	if (!inv) {
	  for (i=0; i<n; i++) {
	    tMsk.changeRep(DYNCAST(StVector,sublvalInt(RangeFull::get(),i)));
	    bMsk.changeRep(DYNCAST(StVector,
				   bmat.subrvalInt(RangeFull::get(),i)));
	    tMsk->addprod(*bMsk,avec);
	  }
	} else {
	  for (i=0; i<n; i++) {
	    tMsk.changeRep(DYNCAST(StVector,sublvalInt(RangeFull::get(),i)));
	    bMsk.changeRep(DYNCAST(StVector,
				   bmat.subrvalInt(RangeFull::get(),i)));
	    tMsk->adddiv(*bMsk,avec);
	  }
	}
      }
    }
  }

  void StMatrix::mulDiag(const StMatrix& bmat,const StVector& avec,
			 bool inv,bool incr)
  {
    bmat.updateStrctPatt();
    if (bmat.strpatt==MatStrct::uppNDg || bmat.strpatt==MatStrct::lowNDg)
      throw WrongStatusException(EXCEPT_MSG("bmat.strpatt"));
    if (avec.size()!=bmat.n)
      throw DimMismatchException(EXCEPT_MSG(""));
    checkTS(TSARG); avec.checkTS(TSARG); bmat.checkTS(TSARG);
    if (!incr) {
      ensureCapacity(bmat.m,bmat.n);
      strpatt=bmat.strpatt;
    } else if (m!=bmat.m || n!=bmat.n || strpatt!=bmat.strpatt)
      throw InvalidParameterException(EXCEPT_MSG(""));
    if (!isMaskInd() && !bmat.isMaskInd()) {
      BaseLinMat<double> tempM,tempB;
      getLinMat(tempM,true,strpatt);
      bmat.getLinMat(tempB,false,bmat.strpatt);
      BaseLinVec<double> tempA;
      avec.getLinVec(tempA,false);
      FastUtils::mulDiag(tempM,tempA,tempB,inv,false,incr);
    } else {
      // One of this, B is indexed mask
      int i;
      Handle<StVector> tMsk,bMsk;
      double ascal;
      if (!incr) {
	for (i=0; i<n; i++) {
	  ascal=(!inv)?avec[i]:1.0/avec[i];
	  tMsk.changeRep(DYNCAST(StVector,sublvalInt(RangeFull::get(),i)));
	  bMsk.changeRep(DYNCAST(StVector,
				 bmat.subrvalInt(RangeFull::get(),i)));
	  tMsk->prod(*bMsk,ascal);
	}
      } else {
	for (i=0; i<n; i++) {
	  ascal=(!inv)?avec[i]:1.0/avec[i];
	  tMsk.changeRep(DYNCAST(StVector,sublvalInt(RangeFull::get(),i)));
	  bMsk.changeRep(DYNCAST(StVector,
				 bmat.subrvalInt(RangeFull::get(),i)));
	  tMsk->addprod(ascal,*bMsk);
	}
      }
    }
  }

  void StMatrix::mulDiagSym(const StMatrix& amat,const StVector& dvec,
			    double alpha,double beta)
  {
    int i;
    Handle <StVector> srcMsk,trgMsk,dMsk;

    amat.updateStrctPatt(); updateStrctPatt();
    checkTS(TSARG); dvec.checkTS(TSARG); amat.checkTS(TSARG);
    if (amat.strpatt!=MatStrct::upper && amat.strpatt!=MatStrct::lower)
      throw WrongStatusException(EXCEPT_MSG("amat.strpatt"));
    if (beta==0.0) {
      ensureCapacity(amat.n,amat.n);
      setStrctPatt(amat.strpatt);
    } else {
      if (this==&amat)
	throw EqualArgsException(EXCEPT_MSG(""));
      if (strpatt!=MatStrct::upper && strpatt!=MatStrct::lower)
	throw WrongStatusException(EXCEPT_MSG("strpatt"));
      if (n!=amat.n)
	throw WrongDimensionException(EXCEPT_MSG(""));
    }
    for (i=0; i<n; i++) {
      Range rng(i);
      if (strpatt==MatStrct::lower)
	trgMsk.changeRep(DYNCAST(StVector,sublvalInt(rng,i)));
      else
	trgMsk.changeRep(DYNCAST(StVector,sublvalInt(i,rng)));
      if (amat.strpatt==MatStrct::lower)
	srcMsk.changeRep(DYNCAST(StVector,amat.subrvalInt(rng,i)));
      else
	srcMsk.changeRep(DYNCAST(StVector,amat.subrvalInt(i,rng)));
      dMsk.changeRep(DYNCAST(StVector,dvec.subrvalInt(rng)));
      if (beta==0) {
	trgMsk->prod(*srcMsk,*dMsk);
	trgMsk->prod(alpha*dvec[i]);
      } else {
	trgMsk->prod(beta);
	trgMsk->addsprod(*trgMsk,*srcMsk,*dMsk,alpha*dvec[i]);
      }
    }
  }

  void StMatrix::addMulDiagSym(const StMatrix& amat,const StVector& dvec,
			       double alpha,double beta)
  {
    int i;
    Handle <StVector> srcMsk,trgMsk,dMsk;

    amat.updateStrctPatt(); updateStrctPatt();
    checkTS(TSARG); dvec.checkTS(TSARG); amat.checkTS(TSARG);
    if (amat.strpatt!=MatStrct::upper && amat.strpatt!=MatStrct::lower)
      throw WrongStatusException(EXCEPT_MSG("amat.strpatt"));
    if (beta==0.0) {
      ensureCapacity(amat.n,amat.n);
      setStrctPatt(amat.strpatt);
    } else {
      if (this==&amat)
	throw EqualArgsException(EXCEPT_MSG(""));
      if (strpatt!=MatStrct::upper && strpatt!=MatStrct::lower)
	throw WrongStatusException(EXCEPT_MSG("strpatt"));
      if (n!=amat.n)
	throw WrongDimensionException(EXCEPT_MSG(""));
    }
    StVector tempVec(n);
    for (i=0; i<n; i++) {
      Range rng(i);
      if (strpatt==MatStrct::lower)
	trgMsk.changeRep(DYNCAST(StVector,sublvalInt(rng,i)));
      else
	trgMsk.changeRep(DYNCAST(StVector,sublvalInt(i,rng)));
      if (amat.strpatt==MatStrct::lower)
	srcMsk.changeRep(DYNCAST(StVector,amat.subrvalInt(rng,i)));
      else
	srcMsk.changeRep(DYNCAST(StVector,amat.subrvalInt(i,rng)));
      dMsk.changeRep(DYNCAST(StVector,dvec.subrvalInt(rng)));
      tempVec.addscal(*dMsk,dvec[i]);
      tempVec.prod(alpha);
      if (beta==0)
	trgMsk->prod(*srcMsk,tempVec);
      else {
	trgMsk->prod(beta);
	trgMsk->addprod(*srcMsk,tempVec);
      }
    }
  }

  void StMatrix::addsmul(const StMatrix& bmat,double alpha)
  {
    int i;
    Handle<StVector> trgM,srcM;

    updateStrctPatt(); bmat.updateStrctPatt();
    if (strpatt!=bmat.strpatt)
      throw WrongStatusException(EXCEPT_MSG("strpatt"));
    if (m!=bmat.m || n!=bmat.n) throw DimMismatchException(EXCEPT_MSG(""));
    checkTS(TSARG); bmat.checkTS(TSARG);
    // Column by column
    if (strpatt==MatStrct::normal) {
      if (!isMaskInd() && !bmat.isMaskInd()) {
	double* trgP=buff;
	const double* srcP=bmat.buff;
	if (stride==m && bmat.stride==m)
	  FastUtils::addsmul(trgP,srcP,alpha,m*n); // both contiguous
	else
	  for (i=0; i<n; i++,trgP+=stride,srcP+=bmat.stride)
	    FastUtils::addsmul(trgP,srcP,alpha,m);
      } else {
	Range rng(0,m-1);
	for (i=0; i<n; i++) {
	  MASKVEC(trgM,rng,i);
	  MASKMVEC(srcM,bmat,rng,i);
	  trgM->addprod(alpha,*srcM);
	}
      }
    } else if (strpatt==MatStrct::upper) {
      // Both upper triangular
      if (!isMaskInd() && !bmat.isMaskInd()) {
	double* trgP=buff;
	const double* srcP=bmat.buff;
	(*trgP)+=alpha*(*srcP);
	for (i=1; i<n; i++) {
	  trgP+=stride; srcP+=bmat.stride;
	  FastUtils::addsmul(trgP,srcP,alpha,i+1);
	}
      } else {
	set(0,0,get(0,0)+alpha*bmat.get(0,0));
	for (i=1; i<n; i++) {
	  Range rng(0,i);
	  MASKVEC(trgM,rng,i);
	  MASKMVEC(srcM,bmat,rng,i);
	  trgM->addprod(alpha,*srcM);
	}
      }
    } else if (strpatt==MatStrct::uppNDg) {
      // Both upper triangular, no diag.
      if (!isMaskInd() && !bmat.isMaskInd()) {
	double* trgP=buff;
	const double* srcP=bmat.buff;
	for (i=1; i<n; i++) {
	  trgP+=stride; srcP+=bmat.stride;
	  FastUtils::addsmul(trgP,srcP,alpha,i);
	}
      } else {
	for (i=1; i<n; i++) {
	  Range rng(0,i-1);
	  MASKVEC(trgM,rng,i);
	  MASKMVEC(srcM,bmat,rng,i);
	  trgM->addprod(alpha,*srcM);
	}
      }
    } else if (strpatt==MatStrct::lower) {
      // Both lower triangular
      if (!isMaskInd() && !bmat.isMaskInd()) {
	double* trgP=buff;
	const double* srcP=bmat.buff;
	for (i=0; i<n-1; i++,trgP+=(stride+1),srcP+=(bmat.stride+1)) {
	  FastUtils::addsmul(trgP,srcP,alpha,m-i);
	}
	(*trgP)+=alpha*(*srcP);
      } else {
	for (i=0; i<n-1; i++) {
	  Range rng(i,m-1);
	  MASKVEC(trgM,rng,i);
	  MASKMVEC(srcM,bmat,rng,i);
	  trgM->addprod(alpha,*srcM);
	}
	set(n-1,n-1,get(n-1,n-1)+alpha*bmat.get(n-1,n-1));
      }
    } else if (strpatt==MatStrct::lowNDg) {
      // Both lower triangular, no diag.
      if (!isMaskInd() && !bmat.isMaskInd()) {
	double* trgP=buff+1;
	const double* srcP=bmat.buff+1;
	for (i=0; i<n-1; i++,trgP+=(stride+1),srcP+=(bmat.stride+1)) {
	  FastUtils::addsmul(trgP,srcP,alpha,m-i-1);
	}
      } else {
	for (i=0; i<n-1; i++) {
	  Range rng(i+1,m-1);
	  MASKVEC(trgM,rng,i);
	  MASKMVEC(srcM,bmat,rng,i);
	  trgM->addprod(alpha,*srcM);
	}
      }
    }
  }

  void StMatrix::addsmul(double bmat,double alpha)
  {
    int i;
    double fact=alpha*bmat;
    Handle<StVector> trgM;

    updateStrctPatt();
    checkTS(TSARG);
    // Column by column
    if (strpatt==MatStrct::normal) {
      if (!isMaskInd()) {
	double* trgP=buff;
	if (stride==m)
	  FastUtils::addscal(trgP,fact,m*n); // contiguous
	else
	  for (i=0; i<n; i++,trgP+=stride)
	    FastUtils::addscal(trgP,fact,m);
      } else {
	Range rng(0,m-1);
	for (i=0; i<n; i++) {
	  MASKVEC(trgM,rng,i);
	  trgM->addscal(fact);
	}
      }
    } else if (strpatt==MatStrct::upper) {
      // Both upper triang.
      if (!isMaskInd()) {
	double* trgP=buff;
	(*trgP)+=fact;
	for (i=1; i<n; i++) {
	  trgP+=stride;
	  FastUtils::addscal(trgP,fact,i+1);
	}
      } else {
	set(0,0,get(0,0)+fact);
	for (i=1; i<n; i++) {
	  MASKVEC(trgM,Range(0,i),i);
	  trgM->addscal(fact);
	}
      }
    } else if (strpatt==MatStrct::uppNDg) {
      // Both upper triang., no diag.
      if (!isMaskInd()) {
	double* trgP=buff;
	for (i=1; i<n; i++) {
	  trgP+=stride;
	  FastUtils::addscal(trgP,fact,i);
	}
      } else {
	for (i=1; i<n; i++) {
	  MASKVEC(trgM,Range(0,i-1),i);
	  trgM->addscal(fact);
	}
      }
    } else if (strpatt==MatStrct::lower) {
      // Both lower triang.
      if (!isMaskInd()) {
	double* trgP=buff;
	for (i=0; i<n-1; i++,trgP+=(stride+1)) {
	  FastUtils::addscal(trgP,fact,m-i);
	}
	(*trgP)+=fact;
      } else {
	for (i=0; i<n-1; i++) {
	  MASKVEC(trgM,Range(i,m-1),i);
	  trgM->addscal(fact);
	}
	set(n-1,n-1,get(n-1,n-1)+fact);
      }
    } else if (strpatt==MatStrct::lowNDg) {
      // Both lower triang., no diag.
      if (!isMaskInd()) {
	double* trgP=buff+1;
	for (i=0; i<n-1; i++,trgP+=(stride+1)) {
	  FastUtils::addscal(trgP,fact,m-i-1);
	}
      } else {
	for (i=0; i<n-1; i++) {
	  MASKVEC(trgM,Range(i+1,m-1),i);
	  trgM->addscal(fact);
	}
      }
    }
  }

  double StMatrix::traceProd(const StMatrix& amat) const
  {
    int i;
    double sum;
    Handle<StVector> mskM,mskA;

    amat.updateStrctPatt(); updateStrctPatt();
    if (m!=amat.rows() || n!=amat.cols())
      throw WrongDimensionException(EXCEPT_MSG(""));
    if (strpatt==MatStrct::uppNDg || strpatt==MatStrct::lowNDg)
      throw WrongStatusException(EXCEPT_MSG("strpatt"));
    if (amat.strpatt==MatStrct::uppNDg || amat.strpatt==MatStrct::lowNDg)
      throw WrongStatusException(EXCEPT_MSG("amat.strpatt"));
    if (strpatt!=amat.strpatt)
      throw WrongStatusException("Structure patterns must be the same");
    checkTS(TSARG); amat.checkTS(TSARG);
    if (strpatt==MatStrct::normal) {
      // Both normal
      if (!isMaskInd() && !amat.isMaskInd()) {
	// Both flat
	if (m==stride && m==amat.stride) {
	  // Both contiguous
	  return FastUtils::inner(buff,amat.buff,m*n);
	} else {
	  const double* mP=buff,*aP=amat.buff;
	  for (i=0,sum=0.0; i<n; i++,mP+=stride,aP+=amat.stride)
	    sum+=FastUtils::inner(mP,aP,m);
	  return sum;
	}
      } else {
	Range rng(0,m-1);
	for (i=0,sum=0.0; i<n; i++) {
	  MASKVEC(mskM,rng,i);
	  MASKMVEC(mskA,amat,rng,i);
	  sum+=mskM->inner(*mskA);
	}
	return sum;
      }
    } else if (strpatt==MatStrct::lower) {
      // Both symmetric (lower triangle)
      if (!isMaskInd() && !amat.isMaskInd()) {
	// Both flat
	const double* mP=buff,*aP=amat.buff;
	for (i=0,sum=0.0; i<n-1; i++,mP+=stride,aP+=amat.stride) {
	  sum+=((*(mP++))*(*(aP++)));
	  sum+=2.0*FastUtils::inner(mP,aP,m-i-1);
	}
	sum+=((*mP)*(*aP));
	return sum;
      } else {
	for (i=0,sum=0.0; i<n-1; i++) {
	  sum+=(get(i,i)*amat.get(i,i));
	  Range rng(i+1,m-1);
	  MASKVEC(mskM,rng,i);
	  MASKMVEC(mskA,amat,rng,i);
	  sum+=2.0*mskM->inner(*mskA);
	}
	sum+=(get(n-1,n-1)*amat.get(n-1,n-1));
	return sum;
      }
    } else {
      // Both symmetric (upper triangle)
      if (!isMaskInd() && !amat.isMaskInd()) {
	// Both flat
	const double* mP=buff,*aP=amat.buff;
	sum=(*mP)*(*aP);
	for (i=1; i<n; i++) {
	  mP+=stride; aP+=amat.stride;
	  sum+=2.0*FastUtils::inner(mP,aP,i);
	  sum+=((*(mP+i))*(*(aP+i)));
	}
	return sum;
      } else {
	sum=get(0,0)*amat.get(0,0);
	for (i=1; i<n; i++) {
	  sum+=(get(i,i)*amat.get(i,i));
	  Range rng(0,i-1);
	  MASKVEC(mskM,rng,i);
	  MASKMVEC(mskA,amat,rng,i);
	  sum+=2.0*mskM->inner(*mskA);
	}
	return sum;
      }
    }
  }

  void StMatrix::sum(StVector& vec,bool rowMaj,bool increm) const
  {
    int i;
    Handle<StVector> msk;

    checkTS(TSARG); vec.checkTS(TSARG);
    if (!rowMaj) {
      if (!increm || m!=vec.size()) vec.zeros(m);
      Range rng(0,m-1);
      for (i=0; i<n; i++) {
	MASKVEC(msk,rng,i);
	vec.addprod(1.0,*msk);
      }
    } else {
      if (!increm || n!=vec.size()) vec.zeros(n);
      Range rng(0,n-1);
      for (i=0; i<m; i++) {
	MASKVEC(msk,i,rng);
	vec.addprod(1.0,*msk);
      }
    }
  }

  void StMatrix::smul(const StMatrix& amat,double alpha)
  {
    int i;
    Handle<StVector> mskM,mskA;

    checkTS(TSARG); amat.checkTS(TSARG);
    amat.updateStrctPatt();
    ensureCapacity(amat.rows(),amat.cols());
    strpatt=amat.strpatt;
    if (strpatt==MatStrct::normal) {
      // Normal
      if (!isMaskInd() && !amat.isMaskInd()) {
	// Both flat
	if (m==stride && m==amat.stride) {
	  // Both contiguous
	  FastUtils::smul(buff,amat.buff,alpha,m*n);
	} else {
	  double* mP=buff;
	  const double* aP=amat.buff;
	  for (i=0; i<n; i++,mP+=stride,aP+=amat.stride)
	    FastUtils::smul(mP,aP,alpha,m);
	}
      } else {
	Range rng(0,m-1);
	for (i=0; i<n; i++) {
	  MASKVEC(mskM,rng,i);
	  MASKMVEC(mskA,amat,rng,i);
	  mskM->prod(*mskA,alpha);
	}
      }
    } else if (strpatt==MatStrct::lower) {
      // Lower triangle
      if (!isMaskInd() && !amat.isMaskInd()) {
	// Both flat
	double* mP=buff;
	const double* aP=amat.buff;
	for (i=0; i<n; i++,mP+=stride+1,aP+=amat.stride+1)
	  FastUtils::smul(mP,aP,alpha,m-i);
      } else {
	for (i=0; i<n; i++) {
	  Range rng(i,m-1);
	  MASKVEC(mskM,rng,i);
	  MASKMVEC(mskA,amat,rng,i);
	  mskM->prod(*mskA,alpha);
	}
      }
    } else if (strpatt==MatStrct::upper) {
      // Upper triangle
      if (!isMaskInd() && !amat.isMaskInd()) {
	// Both flat
	double* mP=buff;
	const double* aP=amat.buff;
	for (i=0; i<n; i++,mP+=stride,aP+=amat.stride)
	  FastUtils::smul(mP,aP,alpha,i+1);
      } else {
	for (i=0; i<n; i++) {
	  Range rng(0,i);
	  MASKVEC(mskM,rng,i);
	  MASKMVEC(mskA,amat,rng,i);
	  mskM->prod(*mskA,alpha);
	}
      }
    } else if (strpatt==MatStrct::lowNDg) {
      // Lower triangle unit diag.
      if (!isMaskInd() && !amat.isMaskInd()) {
	// Both flat
	double* mP=buff+1;
	const double* aP=amat.buff+1;
	for (i=0; i<n-1; i++,mP+=stride+1,aP+=amat.stride+1)
	  FastUtils::smul(mP,aP,alpha,m-i-1);
      } else {
	for (i=0; i<n-1; i++) {
	  Range rng(i+1,m-1);
	  MASKVEC(mskM,rng,i);
	  MASKMVEC(mskA,amat,rng,i);
	  mskM->prod(*mskA,alpha);
	}
      }
    } else {
      // Upper triangle unit diag.
      if (!isMaskInd() && !amat.isMaskInd()) {
	// Both flat
	double* mP=buff+stride;
	const double* aP=amat.buff+amat.stride;
	for (i=1; i<n; i++,mP+=stride,aP+=amat.stride)
	  FastUtils::smul(mP,aP,alpha,i);
      } else {
	for (i=1; i<n; i++) {
	  Range rng(0,i-1);
	  MASKVEC(mskM,rng,i);
	  MASKMVEC(mskA,amat,rng,i);
	  mskM->prod(*mskA,alpha);
	}
      }
    }
  }

  void StMatrix::smul(double alpha)
  {
    int i;
    Handle<StVector> mskM;

    checkTS(TSARG);
    updateStrctPatt();
    if (strpatt==MatStrct::normal) {
      // Normal
      if (!isMaskInd()) {
	// Both flat
	if (m==stride) {
	  // Both contiguous
	  FastUtils::smul(buff,alpha,m*n);
	} else {
	  double* mP=buff;
	  for (i=0; i<n; i++,mP+=stride)
	    FastUtils::smul(mP,alpha,m);
	}
      } else {
	Range rng(0,m-1);
	for (i=0; i<n; i++) {
	  MASKVEC(mskM,rng,i);
	  mskM->prod(alpha);
	}
      }
    } else if (strpatt==MatStrct::lower) {
      // Lower triangle
      if (!isMaskInd()) {
	// Both flat
	double* mP=buff;
	for (i=0; i<n; i++,mP+=stride+1)
	  FastUtils::smul(mP,alpha,m-i);
      } else {
	for (i=0; i<n; i++) {
	  MASKVEC(mskM,Range(i,m-1),i);
	  mskM->prod(alpha);
	}
      }
    } else if (strpatt==MatStrct::upper) {
      // Upper triangle
      if (!isMaskInd()) {
	// Both flat
	double* mP=buff;
	for (i=0; i<n; i++,mP+=stride)
	  FastUtils::smul(mP,alpha,i+1);
      } else {
	for (i=0; i<n; i++) {
	  MASKVEC(mskM,Range(0,i),i);
	  mskM->prod(alpha);
	}
      }
    } else if (strpatt==MatStrct::lowNDg) {
      // Lower triangle unit diag.
      if (!isMaskInd()) {
	// Both flat
	double* mP=buff+1;
	for (i=0; i<n-1; i++,mP+=stride+1)
	  FastUtils::smul(mP,alpha,m-i-1);
      } else {
	Handle<StVector> mskM;
	for (i=0; i<n-1; i++) {
	  MASKVEC(mskM,Range(i+1,m-1),i);
	  mskM->prod(alpha);
	}
      }
    } else {
      // Upper triangle unit diag.
      if (!isMaskInd()) {
	// Both flat
	double* mP=buff+stride;
	for (i=1; i<n; i++,mP+=stride)
	  FastUtils::smul(mP,alpha,i);
      } else {
	Handle<StVector> mskM;
	for (i=1; i<n; i++) {
	  MASKVEC(mskM,Range(0,i-1),i);
	  mskM->prod(alpha);
	}
      }
    }
  }

  void StMatrix::prod(const StMatrix& amat,const StMatrix& bmat)
  {
    int i;
    Handle<StVector> mskM,mskA,mskB;

    amat.updateStrctPatt(); bmat.updateStrctPatt();
    if (amat.strpatt!=bmat.strpatt)
      throw InvalidParameterException(EXCEPT_MSG("amat, bmat structure pattern must be the same"));
    if (amat.rows()!=bmat.rows() || amat.cols()!=bmat.cols())
      throw DimMismatchException(EXCEPT_MSG("amat, bmat must have same size"));
    checkTS(TSARG); amat.checkTS(TSARG); bmat.checkTS(TSARG);
    ensureCapacity(amat.rows(),amat.cols());
    strpatt=amat.strpatt;
    if (strpatt==MatStrct::normal) {
      // Normal
      if (!isMaskInd() && !amat.isMaskInd() && !bmat.isMaskInd()) {
	// All flat
	if (m==stride && m==amat.stride && m==bmat.stride) {
	  // All contiguous
	  FastUtils::prod(buff,amat.buff,bmat.buff,m*n);
	} else {
	  double* mP=buff;
	  const double* aP=amat.buff,
	    *bP=bmat.buff;
	  for (i=0; i<n; i++,mP+=stride,aP+=amat.stride,bP+=bmat.stride)
	    FastUtils::prod(mP,aP,bP,m);
	}
      } else {
	Range rng(0,m-1);
	for (i=0; i<n; i++) {
	  MASKVEC(mskM,rng,i);
	  MASKMVEC(mskA,amat,rng,i);
	  MASKMVEC(mskB,bmat,rng,i);
	  mskM->prod(*mskA,*mskB);
	}
      }
    } else if (strpatt==MatStrct::lower) {
      // Lower triangle
      if (!isMaskInd() && !amat.isMaskInd() && !bmat.isMaskInd()) {
	// Both flat
	double* mP=buff;
	const double* aP=amat.buff,
	  *bP=bmat.buff;
	for (i=0; i<n; i++,mP+=stride+1,aP+=amat.stride+1,
	       bP+=bmat.stride+1)
	  FastUtils::prod(mP,aP,bP,m-i);
      } else {
	for (i=0; i<n; i++) {
	  Range rng(i,m-1);
	  MASKVEC(mskM,rng,i);
	  MASKMVEC(mskA,amat,rng,i);
	  MASKMVEC(mskB,bmat,rng,i);
	  mskM->prod(*mskA,*mskB);
	}
      }
    } else if (strpatt==MatStrct::upper) {
      // Upper triangle
      if (!isMaskInd() && !amat.isMaskInd() && !bmat.isMaskInd()) {
	// Both flat
	double* mP=buff;
	const double* aP=amat.buff,
	  *bP=bmat.buff;
	for (i=0; i<n; i++,mP+=stride,aP+=amat.stride,bP+=bmat.stride)
	  FastUtils::prod(mP,aP,bP,i+1);
      } else {
	Handle<StVector> mskM,mskA,mskB;
	for (i=0; i<n; i++) {
	  Range rng(0,i);
	  MASKVEC(mskM,rng,i);
	  MASKMVEC(mskA,amat,rng,i);
	  MASKMVEC(mskB,bmat,rng,i);
	  mskM->prod(*mskA,*mskB);
	}
      }
    } else if (strpatt==MatStrct::lowNDg) {
      // Lower triangle unit diag.
      if (!isMaskInd() && !amat.isMaskInd() && !bmat.isMaskInd()) {
	// Both flat
	double* mP=buff+1;
	const double* aP=amat.buff+1,
	  *bP=bmat.buff+1;
	for (i=0; i<n-1; i++,mP+=stride+1,aP+=amat.stride+1,bP+=bmat.stride+1)
	  FastUtils::prod(mP,aP,bP,m-i-1);
      } else {
	Handle<StVector> mskM,mskA,mskB;
	for (i=0; i<n-1; i++) {
	  Range rng(i+1,m-1);
	  MASKVEC(mskM,rng,i);
	  MASKMVEC(mskA,amat,rng,i);
	  MASKMVEC(mskB,bmat,rng,i);
	  mskM->prod(*mskA,*mskB);
	}
      }
    } else {
      // Upper triangle unit diag.
      if (!isMaskInd() && !amat.isMaskInd() && !bmat.isMaskInd()) {
	// Both flat
	double* mP=buff+stride;
	const double* aP=amat.buff+amat.stride,
	  *bP=bmat.buff+bmat.stride;
	for (i=1; i<n; i++,mP+=stride,aP+=amat.stride,bP+=bmat.stride)
	  FastUtils::prod(mP,aP,bP,i);
      } else {
	Handle<StVector> mskM,mskA,mskB;
	for (i=1; i<n; i++) {
	  Range rng(0,i-1);
	  MASKVEC(mskM,rng,i);
	  MASKMVEC(mskA,amat,rng,i);
	  MASKMVEC(mskB,bmat,rng,i);
	  mskM->prod(*mskA,*mskB);
	}
      }
    }
  }

  void StMatrix::prod(const StMatrix& bmat)
  {
    int i;
    Handle<StVector> mskM,mskB;

    bmat.updateStrctPatt(); updateStrctPatt();
    if (strpatt!=bmat.strpatt)
      throw InvalidParameterException(EXCEPT_MSG("Structure pattern must be the same"));
    if (m!=bmat.rows() || n!=bmat.cols())
      throw DimMismatchException(EXCEPT_MSG("Matrices must have same size"));
    checkTS(TSARG); bmat.checkTS(TSARG);
    if (strpatt==MatStrct::normal) {
      // Normal
      if (!isMaskInd() && !bmat.isMaskInd()) {
	// Both flat
	if (m==stride && m==bmat.stride) {
	  // Both contiguous
	  FastUtils::prod(buff,bmat.buff,m*n);
	} else {
	  double* mP=buff;
	  const double* bP=bmat.buff;
	  for (i=0; i<n; i++,mP+=stride,bP+=bmat.stride)
	    FastUtils::prod(mP,bP,m);
	}
      } else {
	Range rng(0,m-1);
	for (i=0; i<n; i++) {
	  MASKVEC(mskM,rng,i);
	  MASKMVEC(mskB,bmat,rng,i);
	  mskM->prod(*mskB);
	}
      }
    } else if (strpatt==MatStrct::lower) {
      // Lower triangle
      if (!isMaskInd() && !bmat.isMaskInd()) {
	// Both flat
	double* mP=buff;
	const double* bP=bmat.buff;
	for (i=0; i<n; i++,mP+=stride+1,bP+=bmat.stride+1)
	  FastUtils::prod(mP,bP,m-i);
      } else {
	for (i=0; i<n; i++) {
	  Range rng(i,m-1);
	  MASKVEC(mskM,rng,i);
	  MASKMVEC(mskB,bmat,rng,i);
	  mskM->prod(*mskB);
	}
      }
    } else if (strpatt==MatStrct::upper) {
      // Upper triangle
      if (!isMaskInd() && !bmat.isMaskInd()) {
	// Both flat
	double* mP=buff;
	const double* bP=bmat.buff;
	for (i=0; i<n; i++,mP+=stride,bP+=bmat.stride)
	  FastUtils::prod(mP,bP,i+1);
      } else {
	Handle<StVector> mskM,mskB;
	for (i=0; i<n; i++) {
	  Range rng(0,i);
	  MASKVEC(mskM,rng,i);
	  MASKMVEC(mskB,bmat,rng,i);
	  mskM->prod(*mskB);
	}
      }
    } else if (strpatt==MatStrct::lowNDg) {
      // Lower triangle unit diag.
      if (!isMaskInd() && !bmat.isMaskInd()) {
	// Both flat
	double* mP=buff+1;
	const double* bP=bmat.buff+1;
	for (i=0; i<n-1; i++,mP+=stride+1,bP+=bmat.stride+1)
	  FastUtils::prod(mP,bP,m-i-1);
      } else {
	Handle<StVector> mskM,mskB;
	for (i=0; i<n-1; i++) {
	  Range rng(i+1,m-1);
	  MASKVEC(mskM,rng,i);
	  MASKMVEC(mskB,bmat,rng,i);
	  mskM->prod(*mskB);
	}
      }
    } else {
      // Upper triangle unit diag.
      if (!isMaskInd() && !bmat.isMaskInd()) {
	// Both flat
	double* mP=buff+stride;
	const double* bP=bmat.buff+bmat.stride;
	for (i=1; i<n; i++,mP+=stride,bP+=bmat.stride)
	  FastUtils::prod(mP,bP,i);
      } else {
	Handle<StVector> mskM,mskB;
	for (i=1; i<n; i++) {
	  Range rng(0,i-1);
	  MASKVEC(mskM,rng,i);
	  MASKMVEC(mskB,bmat,rng,i);
	  mskM->prod(*mskB);
	}
      }
    }
  }

  // See 'sumProd' comments
  void StMatrix::sumSquares(StVector& vec,bool rowMaj,bool increm) const
  {
    int i;

    checkTS(TSARG); vec.checkTS(TSARG);
    if (n==0) return;
    int sz,cols;
    if (!rowMaj) {
      sz=m; cols=n;
    } else {
      sz=n; cols=m;
    }
    if (!increm)
      vec.zeros(sz);
    else if (vec.size()!=sz)
      throw WrongDimensionException(EXCEPT_MSG(""));
    BaseLinVec<double> linv;
    vec.getLinVec(linv,true);
    const double* m1P=buff;
    const int* mcind=0,*mrind=0;
    ArrayHandle<int> hmrind;
    int mrstep,mcstep;
    if (rowMaj) {
      mrstep=stride; mcstep=1;
      if (isMaskInd()) {
	mcind=rowInd.p(); mrind=colInd.p();
	if (mrind!=0) {
	  // Multiply index by stride
	  hmrind.changeRep(sz);
	  for (i=0; i<sz; i++) hmrind[i]=mrind[i]*mrstep;
	  mrind=hmrind;
	  mrstep=1;
	}
      }
    } else {
      mrstep=1; mcstep=stride;
      if (isMaskInd()) {
	mrind=rowInd.p(); mcind=colInd.p();
      }
    }
    if (mcind==0) {
      for (i=0; i<cols; i++,m1P+=mcstep)
	FastUtils::addprod(linv.getBuff(),m1P,sz,1.0,linv.step,mrstep,0,mrind);
    } else {
      for (i=0; i<cols; i++) {
	m1P=buff+mcstep*(*(mcind++));
	FastUtils::addprod(linv.getBuff(),m1P,sz,1.0,linv.step,mrstep,0,mrind);
      }
    }
  }

  void StMatrix::sumSquares(StVector& vec,const StVector& wvec,bool rowMaj,
			    bool increm) const
  {
    int i;

    checkTS(TSARG); vec.checkTS(TSARG);
    if (n==0) return;
    int sz,cols;
    if (!rowMaj) {
      sz=m; cols=n;
    } else {
      sz=n; cols=m;
    }
    if (wvec.size()!=cols)
      throw WrongDimensionException(EXCEPT_MSG("wvec"));
    if (!increm)
      vec.zeros(sz);
    else if (vec.size()!=sz)
      throw WrongDimensionException(EXCEPT_MSG(""));
    BaseLinVec<double> linv;
    vec.getLinVec(linv,true);
    const double* m1P=buff;
    const int* mcind=0,*mrind=0;
    ArrayHandle<int> hmrind;
    int mrstep,mcstep;
    if (rowMaj) {
      mrstep=stride; mcstep=1;
      if (isMaskInd()) {
	mcind=rowInd.p(); mrind=colInd.p();
	if (mrind!=0) {
	  // Multiply index by stride
	  hmrind.changeRep(sz);
	  for (i=0; i<sz; i++) hmrind[i]=mrind[i]*mrstep;
	  mrind=hmrind;
	  mrstep=1;
	}
      }
    } else {
      mrstep=1; mcstep=stride;
      if (isMaskInd()) {
	mrind=rowInd.p(); mcind=colInd.p();
      }
    }
    if (mcind==0) {
      for (i=0; i<cols; i++,m1P+=mcstep)
	FastUtils::addprod(linv.getBuff(),m1P,sz,wvec[i],linv.step,mrstep,0,
			   mrind);
    } else {
      for (i=0; i<cols; i++) {
	m1P=buff+mcstep*(*(mcind++));
	FastUtils::addprod(linv.getBuff(),m1P,sz,wvec[i],linv.step,mrstep,0,
			   mrind);
      }
    }
  }

  /*
   * Idea is consider "virtual" matrices M_v, A_v. A row of M_v is
   * traversed either using step size 'mrstep' or index 'mrind' (not both;
   * if 'mrind'!=0, then 'mrstep'==1). If 'mcind'==0, 'mcstep' is the stride
   * to jump to the next column, otherwise the offset for column i is
   * 'mcind[i]*mcstep'. Same for A_v.
   * M_v is M if 'rowMaj'==false, M^T otherwise.
   * A_v is A^T if ex. one of 'rowMaj', 'atrans' true, A otherwise.
   * We then compute diag( M_v^T A_v ), runing column by column over M_v,
   * A_v.
   */
  void StMatrix::sumProd(const StMatrix& a,StVector& vec,bool atrans,
			 bool rowMaj,bool increm) const
  {
    int i;

    checkTS(TSARG); vec.checkTS(TSARG); a.checkTS(TSARG);
    if ((!atrans && (a.rows()!=m || a.cols()!=n)) ||
	(atrans && (a.cols()!=m || a.rows()!=n)))
      throw DimMismatchException(EXCEPT_MSG(""));
    if (n==0) return;
    int sz,cols;
    if (!rowMaj) {
      sz=m; cols=n;
    } else {
      sz=n; cols=m;
    }
    if (!increm)
      vec.zeros(sz);
    else if (vec.size()!=sz)
      throw WrongDimensionException(EXCEPT_MSG(""));
    BaseLinVec<double> linv;
    vec.getLinVec(linv,true);
    const double* abuff=a.buff;
    const double* m1P=buff,*m2P=abuff;
    const int* mcind=0,*mrind=0,*acind=0,*arind=0;
    ArrayHandle<int> hmrind;
    int mrstep,mcstep;
    if (rowMaj) {
      mrstep=stride; mcstep=1;
      if (isMaskInd()) {
	mcind=rowInd.p(); mrind=colInd.p();
	if (mrind!=0) {
	  // Multiply index by stride
	  hmrind.changeRep(sz);
	  for (i=0; i<sz; i++) hmrind[i]=mrind[i]*mrstep;
	  mrind=hmrind;
	  mrstep=1;
	}
      }
    } else {
      mrstep=1; mcstep=stride;
      if (isMaskInd()) {
	mrind=rowInd.p(); mcind=colInd.p();
      }
    }
    ArrayHandle<int> harind;
    int arstep,acstep;
    if (atrans^rowMaj) {
      arstep=a.getStride_INT(); acstep=1;
      if (a.isMaskInd()) {
	acind=a.rowInd.p(); arind=a.colInd.p();
	if (arind!=0) {
	  // Multiply index by stride
	  harind.changeRep(sz);
	  for (i=0; i<sz; i++) harind[i]=arind[i]*arstep;
	  arind=harind;
	  arstep=1;
	}
      }
    } else {
      arstep=1; acstep=a.getStride_INT();
      if (a.isMaskInd()) {
	arind=a.rowInd.p(); acind=a.colInd.p();
      }
    }
    if (mcind==0) {
      if (acind==0) {
	for (i=0; i<cols; i++,m1P+=mcstep,m2P+=acstep)
	  FastUtils::addprod(linv.getBuff(),m1P,m2P,sz,1.0,linv.step,mrstep,
			     arstep,0,mrind,arind);
      } else {
	for (i=0; i<cols; i++,m1P+=mcstep) {
	  m2P=abuff+acstep*(*(acind++));
	  FastUtils::addprod(linv.getBuff(),m1P,m2P,sz,1.0,linv.step,mrstep,
			     arstep,0,mrind,arind);
	}
      }
    } else {
      if (acind==0) {
	for (i=0; i<cols; i++,m2P+=acstep) {
	  m1P=buff+mcstep*(*(mcind++));
	  FastUtils::addprod(linv.getBuff(),m1P,m2P,sz,1.0,linv.step,mrstep,
			     arstep,0,mrind,arind);
	}
      } else {
	for (i=0; i<cols; i++) {
	  m1P=buff+mcstep*(*(mcind++));
	  m2P=abuff+acstep*(*(acind++));
	  FastUtils::addprod(linv.getBuff(),m1P,m2P,sz,1.0,linv.step,mrstep,
			     arstep,0,mrind,arind);
	}
      }
    }
  }

  void StMatrix::makeSymmAuto()
  {
    updateStrctPatt();
    if (strpatt==MatStrct::normal)
      throw WrongStatusException(EXCEPT_MSG("strpatt"));
    checkTS(TSARG);
    switch (strpatt) {
    case MatStrct::uppNDg:
      // unit diag.
      ((StVector&) diag())=1.0;
      // no break intended!
    case MatStrct::upper:
      makeSymm(false);
      break;
    case MatStrct::lowNDg:
      // unit diag.
      ((StVector&) diag())=1.0;
      // no break intended!
    case MatStrct::lower:
      makeSymm(true);
      break;
    }
    strpatt=MatStrct::normal;
  }

  void StMatrix::sym()
  {
    int i;

    checkTS(TSARG);
    if (m==0) return;
    if (m!=n) throw WrongDimensionException(EXCEPT_MSG("Must be square"));
    Handle<StVector> mskS,mskT;
    for (i=0; i<n-1; i++) {
      Range rng(i+1);
      MASKVEC(mskS,rng,i);
      MASKVEC(mskT,i,rng);
      mskT->addprod(1.0,*mskS);
      mskT->prod(0.5);
      *mskS=*mskT;
    }
  }

  int StMatrix::cholDecomp(bool clear)
  {
    if (isMaskInd())
      throw MaskNotImplException(EXCEPT_MSG(""));
    if (n!=m || n==0) throw WrongDimensionException(EXCEPT_MSG(""));
    checkTS(TSARG); updateStrctPatt();
    bool upper;
    if (strpatt==MatStrct::upper) upper=true;
    else if (strpatt==MatStrct::lower) upper=false;
    else throw WrongStatusException(EXCEPT_MSG(""));
    if (clear) {
      // Clear opposite triangle
      strpatt=upper?MatStrct::lowNDg:MatStrct::uppNDg;
      fillStrct(n,n,0.0);
      strpatt=upper?MatStrct::upper:MatStrct::lower; // reset
    }
    // LAPACK DPOTRF
    const char* uplo=upper?"U":"L";
    int info;
    F77_FUNC(dpotrf,DPOTRF) (uplo,&n,buff,&stride,&info);

    return info;
  }

  int StMatrix::invCholInPlace()
  {
    if (n!=m || n==0) throw WrongDimensionException(EXCEPT_MSG(""));
    if (isMaskInd())
      throw MaskNotImplException(EXCEPT_MSG(""));
    checkTS(TSARG); updateStrctPatt();
    bool upper=(strpatt==MatStrct::upper);
    if (!upper && strpatt!=MatStrct::lower)
      throw WrongStatusException(EXCEPT_MSG("strpatt"));
    // LAPACK DPOTRI
    const char* uplo=upper?"U":"L";
    int info;
    F77_FUNC(dpotri,DPOTRI) (uplo,&n,buff,&stride,&info);

    return info;
  }

  /*
   * The i-th comp. of diag A^{-1} is squared norm of i-th col. of L^{-1}.
   * We use 'invCholInPlace' (first part) to compute L^{-1} from L.
   * In 'invCholInPlace' L is overwritten col. by col. by L^{-1}. The
   * computation of [L^{-1}]_{j,i}, j>i, requires
   * elements of cols k>i of L and elements of [L^{-1}]_{k,i}, k<j. To
   * compute diag A^{-1} WITHOUT overwriting L, we only need one aux.
   * vector 'diagV', storing the current column i of L^{-1} in pos.
   * i,...,n-1, while the result builds up in 0,...,i-1.
   */
  void StMatrix::diagInvChol(const StVector& p,StVector& diagV) const
  {
    int i,j;
    const double* rowB;
    double sum,temp,actEl;

    if (n!=m || n==0) throw WrongDimensionException(EXCEPT_MSG(""));
    if (p.size()!=n) throw DimMismatchException(EXCEPT_MSG(""));
    if (isMaskInd())
      throw MaskNotImplException(EXCEPT_MSG(""));
    checkTS(TSARG); p.checkTS(TSARG); diagV.checkTS(TSARG);
    BaseLinVec<double> lindiag;
    diagV.zeros(n);
    diagV.getLinVec(lindiag,true);
    for (i=0; i<n; i++) {
      temp=lindiag[i]=1.0/p[i];
      actEl=(temp*temp); // accum. for diag. elem.
      rowB=buff+(i*stride+i+1);
      for (j=i+1; j<n; j++,rowB++) {
	// NOTE: 'lindiag' uses Fortran convention for negative
	// lindiag.step. Access lindiag[i..j-1] here
	sum=-FastUtils::inner(rowB,lindiag.getBuff()+
			      ((lindiag.step>0)?i:(j-n))*lindiag.step,j-i,
			      stride,lindiag.step);
	temp=diagV[j]=sum/p[j];
	actEl+=(temp*temp);
      }
      diagV[i]=actEl;
    }
  }

  int StMatrix::luDecomp(BaseVector<int>& pind)
  {
    if (isMaskInd())
      throw MaskNotImplException(EXCEPT_MSG(""));
    if (n!=m || n==0) throw WrongDimensionException(EXCEPT_MSG(""));
    checkTS(TSARG); pind.checkTS(TSARG);
    // !!! Def. fill value broken, have to spec fill value (0) !!
    pind.fill(n,0); // ensure that buffer large enough
    BaseLinVec<int> pilin;
    pind.getLinVec(pilin,true,true,true); // has to be flat
    // LAPACK DGETRF
    int info;
    F77_FUNC(dgetrf,DGETRF) (&n,&n,buff,&stride,pilin.getBuff(),&info);

    return info;
  }

  void StMatrix::luSolve(const BaseVector<int>& pind,const StMatrix& b,
			 StMatrix& x,bool trans) const
  {
    if (isMaskInd())
      throw MaskNotImplException(EXCEPT_MSG(""));
    if (n!=m || n==0 || n!=b.rows() || n!=pind.size())
      throw WrongDimensionException(EXCEPT_MSG(""));
    checkTS(TSARG); pind.checkTS(TSARG); b.checkTS(TSARG); x.checkTS(TSARG);
    x=b; // copy
    BaseLinMat<double> xlin;
    x.getLinMat(xlin,true);
    BaseLinVec<int> pilin;
    pind.getLinVec(pilin,false,true,true); // has to be flat
    // LAPACK DGETRS
    int info,nrhs=b.cols();
    const char* strans=trans?"T":"N";
    F77_FUNC(dgetrs,DGETRS) (strans,&n,&nrhs,buff,&stride,pilin.getBuff(),
			     xlin.buff,&xlin.stride,&info);
  }

  double StMatrix::logDetLU(const BaseVector<int>& pind,double& sign) const
  {
    int i;

    if (n!=m || n==0 || n!=pind.size())
      throw WrongDimensionException(EXCEPT_MSG(""));
    checkTS(TSARG); pind.checkTS(TSARG);
    // Only need parity of P
    bool ispeven=true;
    for (i=0; i<n; i++)
      if (pind[i]!=i+1) ispeven=!ispeven;
    // Diagonal of U
    Handle<StVector> dgMsk(DYNCAST(StVector,diagInt(0,RangeFull::get())));
    double ret=0.0,temp;
    for (i=0; i<n; i++) {
      if ((temp=(*dgMsk)[i])<=0.0) {
	ispeven=!ispeven;
	ret+=log(-temp);
      } else
	ret+=log(temp);
    }

    sign=ispeven?1.0:(-1.0);
    return ret;
  }

  /*
   * e_j^T U^-1 L^-1 P^T e_j = (U^-T e_j)^T (L^-1 P^T e_j) = x^T y,
   * U^T x = e_j, L y = P^T e_j.
   * Here, P^T e_j means to apply permutation given by 'pind' on e_j.
   * NOTE: Present implementation does not use the fact that r.h.s.
   * are unit vectors (could save above half the time)!
   */
  void StMatrix::diagInvLU(const BaseVector<int>& pind,StVector& diagV)
  {
    int j,i,pi;
    double temp;

    if (n!=m || n==0)
      throw WrongDimensionException(EXCEPT_MSG(""));
    checkTS(TSARG); pind.checkTS(TSARG);
    if (pind.size()!=n)
      throw WrongDimensionException(EXCEPT_MSG(""));
    uchar spatt=strpatt;
    StVector xvec(n),yvec(n);
    diagV.ensureCapacity(n);
    for (j=0; j<n; j++) {
      xvec.eye(j);
      setStrctPatt(MatStrct::upper);
      backsubst(xvec,true);
      yvec.eye(j);
      for (i=0; i<n; i++) {
	pi=pind[i]-1;
	temp=yvec[i]; yvec[i]=yvec[pi]; yvec[pi]=temp;
      }
      setStrctPatt(MatStrct::lowNDg);
      backsubst(yvec);
      diagV[j]=xvec.inner(yvec);
    }
    setStrctPatt(spatt);
  }

  int StMatrix::eigenVals(StVector& eigs)
  {
    if (isMaskInd())
      throw MaskNotImplException(EXCEPT_MSG(""));
    if (n!=m || n==0) throw WrongDimensionException(EXCEPT_MSG(""));
    checkTS(TSARG); eigs.checkTS();
    updateStrctPatt();
    bool upper=(strpatt==MatStrct::upper);
    if (!upper && strpatt!=MatStrct::lower)
      throw WrongStatusException(EXCEPT_MSG("strpatt"));
    eigs.ensureCapacity(n);
    BaseLinVec<double> lineigs;
    eigs.getLinVec(lineigs,true);
    // LAPACK DSYEVR
    // Query call to find out sizes for working arrays
    int i1,i2,i3,info;
    double d1,d2,d3;
    const char* jobz="N",*range="A",*uplo=upper?"U":"L";
    i1=-1;
    d1=0.0; i3=1;
    F77_FUNC(dsyevr,DSYEVR) (jobz,range,uplo,&n,buff,&stride,&d3,&d3,&i3,&i3,
			     &d1,&i3,lineigs.getBuff(),&d3,&i3,&i3,&d2,&i1,
			     &i2,&i1,&info);
    if (info!=0)
      throw InternalException(EXCEPT_MSG("WS query failed"));
    ArrayHandle<double> work((int) d2); // workspace
    ArrayHandle<int> iwork(i2);
    // Prepare for call
    // We use ABSTOL=0 (appropriate??)
    d1=0.0; i1=work.size(); i2=iwork.size(); i3=1;
    F77_FUNC(dsyevr,DSYEVR) (jobz,range,uplo,&n,buff,&stride,&d3,&d3,&i3,&i3,
			     &d1,&i3,lineigs.getBuff(),&d3,&i3,&i3,work,&i1,
			     iwork,&i2,&info);

    return info;
  }

  /*
   * The argument ISUPPZ is ignored. It seems to supply additional information,
   * NOT support a sparse storage of the Z array.
   */
  int StMatrix::eigenDecomp(StVector& dvec,StMatrix& umat)
  {
    if (isMaskInd())
      throw MaskNotImplException(EXCEPT_MSG(""));
    if (n!=m || n==0) throw WrongDimensionException(EXCEPT_MSG(""));
    checkTS(TSARG); dvec.checkTS(); umat.checkTS();
    updateStrctPatt();
    bool upper=(strpatt==MatStrct::upper);
    if (!upper && strpatt!=MatStrct::lower)
      throw WrongStatusException(EXCEPT_MSG("strpatt"));
    dvec.ensureCapacity(n);
    BaseLinVec<double> lind;
    dvec.getLinVec(lind,true);
    umat.ensureSize(n,n); umat.setStrctPatt(MatStrct::normal);
    // Ensure write-back of 'linu'
    int info;
    {
      BaseLinMat<double> linu;
      umat.getLinMat(linu,true);
      ArrayHandle<int> isuppz(2*n);
      // Query call to find out sizes for working arrays
      int i1,i2,i3;
      double d1,d2,d3;
      const char* jobz="V",*range="A",*uplo=upper?"U":"L";
      i1=-1;
      d1=0.0; i3=1;
      F77_FUNC(dsyevr,DSYEVR) (jobz,range,uplo,&n,buff,&stride,&d3,&d3,&i3,&i3,
			       &d1,&i3,lind.getBuff(),linu.buff,&(linu.stride),
			       isuppz,&d2,&i1,&i2,&i1,&info);
      if (info!=0)
	throw InternalException(EXCEPT_MSG("WS query failed"));
      ArrayHandle<double> work((int) d2); // workspace
      ArrayHandle<int> iwork(i2);
      // Prepare for call
      // We use ABSTOL=0 (appropriate??)
      d1=0.0; i1=work.size(); i2=iwork.size(); i3=1;
      F77_FUNC(dsyevr,DSYEVR) (jobz,range,uplo,&n,buff,&stride,&d3,&d3,&i3,&i3,
			       &d1,&i3,lind.getBuff(),linu.buff,&(linu.stride),
			       isuppz,work,&i1,iwork,&i2,&info);
    }
    // Make sure U is unique (if all eigenvalues distinct): For each col.,
    // the first elem. !=0 must be positive
    StVector colMsk;
    int i,pos;
    for (i=0; i<n; i++) {
      colMsk.reassign(umat(RangeFull::get(),i));
      if ((pos=colMsk.findfst(bind2nd(std::not_equal_to<double>(),0.0)))!=-1)
	if (colMsk[pos]<0.0) colMsk.prod(-1.0); // flip
    }

    return info;
  }

  /*
   * The algorithm iterates a scheme for
   *   D^(i) + a a^T = L D^(i+1) L^T,
   * where L is defined by 2 (m-1)-vectors and a backsubstitution with L
   * is O(n). Initially, D^(0) = D. If a is the first col. of A and B
   * collects the rest (m-by-(n-1)), then
   *   D + A A^T = L^(1) D^(1) L^{(1)T} + B B^T.
   * Let C = L^{(1) -1} B, and continue on D^(1) + C C^T in the same
   * way. The implementation is in place, does not need additional memory.
   *
   * 'lFact' is used as working memory as follows. For initialization, we
   * copy this matrix into the lower part of 'lFact', i.e. col j of this
   * matrix A is copied into row j of 'lFact', at positions m-2,...,2m-3.
   * Also, set D^(0) = D (stored in 'diag').
   * The outer loop is i=1,...,n. In iteration i:
   * - Compute \beta^(i) and D^(i) from p^(i), D^(i-1). D^(i) overwrites
   *   D^(i-1) in 'diag', p^(i) is read from col i, pos. m-2,...,2m-3. For
   *   example, p^(1) is just the first col of A. \beta^(i) is stored in
   *   col i of 'lFact', at positions 0,...,m-2. This defines the new L^(i).
   *   NOTE: The last element of \beta^(i) overwrites the first element of
   *   p^(i) in col i of 'lFact'. That is OK, because this first element is
   *   not needed anymore at that time.
   * - For j=i+1,...,n-1: Solve L^(i) x = p^(j) by backsubstitution, the
   *   result x overwrites p^(j) in col j of 'lFact'.
   * Note that at the end of iteration i, we have:
   * - \beta^(j),p^(j) in col j of 'lFact' for all j<=i
   *   NOTE: In the documentation to this routine, \beta and p are written
   *   as m-vectors, while here they are stored as (m-1)-vectors. This is OK,
   *   because the last element of the \beta and the first element of the p
   *   vectors are not needed anymore.
   * - Sol. of L^(1) ... L^(i) x = a^(j) in col j of 'lFact' for all j>i,
   *   where a^(j) is the j-th col of A.
   */
  void StMatrix::cholDecLowRank(StVector& diag,StMatrix& lFact) const
  {
    int i,j,k;
    double temp,sum,t,tlam,lam;
    double* colB,*pB,*betaB;

    if (diag.size()!=m) throw WrongDimensionException(EXCEPT_MSG(""));
    checkTS(TSARG); diag.checkTS(TSARG); lFact.checkTS(TSARG);
    if (isMaskInd() || lFact.isMaskInd())
      throw MaskObjectException("StMatrix::cholDecLowRank not implemented for indexed masks!");
    lFact.ensureCapacity(2*(m-1),n);
    // Copy this matrix into lower part of 'lFact' -> a^(j)
    ((StMatrix&) lFact(Range(m-2,2*m-3),RangeFull::get()))=*this;

    // Main loop
    for (i=0; i<n; i++) {
      // Compute L^(i), D^(i): Using the t-recurrence
      betaB=lFact.buff+(i*lFact.stride); // \beta^(i)
      pB=betaB+(m-2); // p^(i) (m-vector, will be m-1 when finished)
      t=1.0;
      for (j=0; j<m-1; j++) {
	temp=*(pB++); // p_j
	lam=diag[j]; // \lambda_j
	tlam=t*lam+temp*temp; // t_j \lambda_j
	diag[j]=tlam/t; // \tilde{\lambda}_j
	t=tlam/lam; // t_j
	*(betaB++)=temp/tlam; // \beta_j
      }
      temp=*pB;
      tlam=t*diag[m-1]+temp*temp;
      diag[m-1]=tlam/t;
      // Backsubstitutions (columns >i)
      for (j=i+1; j<n; j++) {
	betaB=lFact.buff+(i*lFact.stride);
	pB=betaB+(m-1); // \beta^(i),p^(i) ((m-1)-vectors)
	colB=lFact.buff+(j*lFact.stride+m-2);
	for (k=1,sum=0.0,temp=*(colB++); k<m; k++) {
	  sum+=(temp*(*(betaB++)));
	  temp=((*(colB++))-=(sum*(*(pB++))));
	}
      }
    }
    diag.apply1(ptr_fun<double, double>(sqrt));
  }

  void StMatrix::cholDecLowRankMinus(StVector& diag,StMatrix& lFact) const
  {
    int i,j,k;
    double temp,sum,s,b,tlam;
    double* colB,*pB,*betaB;

    if (diag.size()!=m) throw WrongDimensionException(EXCEPT_MSG(""));
    checkTS(TSARG); diag.checkTS(TSARG); lFact.checkTS(TSARG);
    if (isMaskInd() || lFact.isMaskInd())
      throw MaskObjectException("StMatrix::cholDecLowRankMinus not implemented for indexed masks!");
    lFact.ensureCapacity(2*(m-1),n);
    // Copy this matrix into lower part of 'lFact'
    ((StMatrix&) lFact(Range(m-2,2*m-3),RangeFull::get()))=*this;

    // Main loop
    for (i=0; i<n; i++) {
      // Compute L^(i), D^(i): Using the s-recurrence
      betaB=lFact.buff+(i*lFact.stride); // \beta^(i)
      pB=betaB+(m-2); // p^(i) (m-vector, will be m-1 when finished)
      s=-1.0;
      for (j=0; j<n-1; j++) {
	temp=*(pB++); // p_j
	b=temp*s; // p_j s_{j-1}
	tlam=(diag[j]+=(temp*b)); // \tilde{\lambda}_j
	if (tlam<=0)
	  throw NumericalException("StMatrix::cholDecLowRankMinus: Matrix not pos. definite!");
	temp=*(betaB++)=b/tlam; // \beta_j
	s-=(b*temp); // s_j
      }
      temp=*pB;
      diag[n-1]+=(s*temp*temp);
      if (diag[n-1]<=0)
	throw NumericalException("StMatrix::cholDecLowRankMinus: Matrix not pos. definite!");
      // Backsubstitutions (cols j>i)
      for (j=i+1; j<n; j++) {
	betaB=lFact.buff+(i*lFact.stride);
	pB=betaB+(m-1); // \beta^(i), p^(i)
	colB=lFact.buff+(i*lFact.stride+m-2);
	for (k=1,sum=0.0,temp=*(colB++); k<n; k++) {
	  sum+=(temp*(*(betaB++)));
	  temp=((*(colB++))-=(sum*(*(pB++))));
	}
      }
    }
    diag.apply1(ptr_fun<double, double>(sqrt));
  }

  /*
   * We have L = L^(1) ... L^(n) (D^(n))^{1/2}. L x = b is done by backsubst.
   * with L^(i), i=1,...,n, then mult. with (D^(n))^{-1/2}.
   */
  void StMatrix::cholLRBacksubst(const StVector& diag,const StVector& b,
				 StVector& x) const
  {
    int i,j,realn;
    double temp,sum;
    double* xB;
    const double* pB,*betaB;

    realn=diag.size();
    if (b.size()!=realn || m!=2*(realn-1))
      throw WrongDimensionException(EXCEPT_MSG(""));
    checkTS(TSARG); diag.checkTS(TSARG); b.checkTS(TSARG); x.checkTS(TSARG);
    x=b;
    {
      BaseLinVec<double> xlin;
      x.getLinVec(xlin,true);
      for (i=0; i<n; i++) {
	// Backsubstitution with L^(i)
	xB=(double*) x.buff;
	betaB=buff+(i*stride); pB=betaB+(realn-1);
	for (j=1,sum=0.0,temp=*(xB++); j<realn; j++) {
	  sum+=(temp*(*(betaB++)));
	  temp=((*(xB++))-=(sum*(*(pB++))));
	}
      }
    }
    // xlin -> x written back here
    x.div(diag);
  }

  void StMatrix::cholLRBacksubst(const StVector& diag,const StMatrix& b,
				 StMatrix& x,bool trans) const
  {
    int i,j,realn,bsz;
    double* tvB,*colB;
    const double* pB,*betaB;

    realn=diag.size();
    if ((!trans && b.rows()!=realn) || (trans && b.cols()!=realn) ||
	m!=2*(realn-1))
      throw WrongDimensionException(EXCEPT_MSG(""));
    checkTS(TSARG); diag.checkTS(TSARG); b.checkTS(TSARG); x.checkTS(TSARG);
    if (isMaskInd() || b.isMaskInd() || x.isMaskInd())
      throw MaskObjectException("StMatrix::cholLRBacksubst not implemented for indexed masks!");
    x=b;
    x.setStrctPatt(MatStrct::normal);
    bsz=trans?b.rows():b.cols();
    StVector tempVec(bsz); tvB=(double*) tempVec.buff;

    if (!trans) {
      for (i=0; i<n; i++) {
	// Backsubstitution with L^(i)
	betaB=buff+(i*stride); pB=betaB+(realn-1);
	tempVec.zeros(); // "sum" variables
	colB=(double*) x.buff;
	for (j=1; j<realn; j++) {
	  FastUtils::addsmul(tvB,colB,*(betaB++),bsz,1,x.stride);
	  colB++;
	  FastUtils::addsmul(colB,tvB,-*(pB++),bsz,x.stride);
	}
      }
      x.mulDiag(diag,x,true);
    } else {
      for (i=0; i<n; i++) {
	// Backsubstitution with L^(i)
	betaB=buff+(i*stride); pB=betaB+(realn-1);
	tempVec.zeros(); // "sum" variables
	colB=(double*) x.buff;
	for (j=1; j<realn; j++) {
	  FastUtils::addsmul(tvB,colB,*(betaB++),bsz);
	  colB+=x.stride;
	  FastUtils::addsmul(colB,tvB,-*(pB++),bsz);
	}
      }
      x.mulDiag(x,diag,true);
    }
  }

  /*
   * We have L^T = (D^(n))^{1/2} L^(n)^T ... L^(1)^T . L^T x = b is done by
   * mult. with (D^(n))^{-1/2}, then backsubst. with L^(i)^T, i=n,...,1.
   */
  void StMatrix::cholLRBacksubsT(const StVector& diag,const StVector& b,
				 StVector& x) const
  {
    int i,j,realn;
    double temp,sum;
    double* xB;
    const double* pB,*betaB;

    realn=diag.size();
    if (b.size()!=realn || m!=2*(realn-1))
      throw WrongDimensionException(EXCEPT_MSG(""));
    checkTS(TSARG); diag.checkTS(TSARG); b.checkTS(TSARG); x.checkTS(TSARG);
    if (isMaskInd())
      throw MaskObjectException("StMatrix::cholLRBacksubsT not implemented for indexed masks!");

    x.div(b,diag);
    BaseLinVec<double> xlin;
    x.getLinVec(xlin,true);
    for (i=n-1; i>=0; i--) {
      // Backsubstitution with L^(i)^T
      xB=(double*) x.buff+(realn-1);
      betaB=buff+(i*stride+realn-2); pB=betaB+(realn-1);
      for (j=realn-2,sum=0.0,temp=*(xB--); j>=0; j--) {
	sum+=(temp*(*(pB--)));
	temp=((*(xB--))-=(sum*(*(betaB--))));
      }
    }
  }

  void StMatrix::cholLRBacksubsT(const StVector& diag,const StMatrix& b,
				 StMatrix& x,bool trans) const
  {
    int i,j,realn,bsz;
    double* tvB,*colB;
    const double* pB,*betaB;

    realn=diag.size();
    if ((!trans && b.rows()!=realn) || (trans && b.cols()!=realn) ||
	m!=2*(realn-1))
      throw WrongDimensionException(EXCEPT_MSG(""));
    checkTS(TSARG); diag.checkTS(TSARG); b.checkTS(TSARG); x.checkTS(TSARG);
    if (isMaskInd() || b.isMaskInd() || x.isMaskInd())
      throw MaskObjectException("StMatrix::cholLRBacksubst not implemented for indexed masks!");
    x=b;
    x.setStrctPatt(MatStrct::normal);
    bsz=trans?b.rows():b.cols();
    StVector tempVec(bsz); tvB=(double*) tempVec.buff;

    if (!trans) {
      x.mulDiag(diag,x,true);
      for (i=m-1; i>=0; i--) {
	// Backsubstitution with L^(i)^T
	betaB=buff+(i*stride+realn-2); pB=betaB+(realn-1);
	colB=(double*) x.buff+(realn-1);
	for (j=realn-2; j>=0; j--) {
	  FastUtils::addsmul(tvB,colB,*(pB--),bsz,1,x.stride);
	  colB--;
	  FastUtils::addsmul(colB,tvB,-*(betaB--),bsz,x.stride);
	}
      }
    } else {
      x.mulDiag(x,diag,true);
      for (i=m-1; i>=0; i--) {
	// Backsubstitution with L^(i)^T
	betaB=buff+(i*stride+realn-2); pB=betaB+(realn-1);
	colB=(double*) x.buff+((realn-1)*x.stride);
	for (j=realn-2; j>=0; j--) {
	  FastUtils::addsmul(tvB,colB,*(pB--),bsz);
	  colB-=x.stride;
	  FastUtils::addsmul(colB,tvB,-*(betaB--),bsz);
	}
      }
    }
  }

  /*
   * We have \tilde{L} = L^(1) ... L^(n) (D^(n))^{1/2}. We multiply the
   * conv. factor L by L^(i), i=1,...,n, finally by the diagonal matrix
   * (D^(n))^{1/2}. Each multiplication with a L^(i) can be done in
   * place.
   */
  void StMatrix::cholLRUpdateFact(const StVector& diag,StMatrix& fact,
				  StVector& dfact) const
  {
    int i,j,k,realn;
    double rho,sum,temp;
    double* rowB;
    const double* pB,*betaB;

    realn=diag.size();
    if (dfact.size()!=realn || m!=2*(realn-1))
      throw WrongDimensionException(EXCEPT_MSG(""));
    if (fact.rows()!=realn || fact.cols()!=realn)
      throw WrongDimensionException(EXCEPT_MSG(""));
    checkTS(TSARG); diag.checkTS(TSARG); fact.checkTS(TSARG); dfact.checkTS(TSARG);
    if (isMaskInd() || fact.isMaskInd())
      throw MaskObjectException("StMatrix::cholLRUpdateFact not implemented for indexed masks!");

    for (i=0; i<n; i++) {
      // Multiplication with L^(i)
      for (j=1; j<realn; j++) {
	betaB=buff+(i*stride+j-1); pB=betaB+(n-1);
	rho=dfact[j]*(*(pB--));
	rowB=fact.buff+(j*fact.stride+j-1);
	for (k=i-1,sum=0.0; k>=0; k--) {
	  sum+=rho;
	  temp=*rowB; rho=temp*(*(pB--));
	  *rowB=temp+sum*(*(betaB--));
	  rowB-=fact.stride;
	}
      }
    }
    // Multiplication with D^{(n) 1/2}
    StVector tempVec;
    Handle<StVector> mskDg(DYNCAST(StVector,fact.diagInt(0,RangeFull::get())));
    tempVec=*mskDg; // copy diagonal
    (*mskDg)=dfact;
    uchar oldstrpatt=fact.strpatt;
    fact.strpatt=MatStrct::lower;
    fact.mulDiag(fact,diag);
    dfact=*mskDg;
    (*mskDg)=tempVec; // copy back
    fact.strpatt=oldstrpatt;
  }

  /*
   * NOTE: We do not exactly follow the TR
   *   M. Seeger
   *   Low Rank Updates for the Cholesky Decomposition
   * but rather the convention used by LINPACK dchud. The plane rotations
   * are always done on (n+1)-vectors with one auxiliary dimension. In dchud,
   * this is the last one, in the TR it is the first one. The LINPACK
   * convention has the advantage that BLAS drotg can be used more easily.
   *
   * NOTE: dchud can end up with negative entries on diag(L). We correct this
   * problem here, by flipping signs of c_i, s_i whenever it happens.
   */
  void StMatrix::cholUpdRk1(const StVector& vvec,StVector& cvec,StVector& svec,
			    StVector& wkvec,StMatrix* dragZ,
			    const StVector* dragY)
  {
    int i,r=0,stp,sz;
    double temp;
    bool doDrag=(dragZ!=0),islower=(strpatt==MatStrct::lower);
    double* tbuff;

    if (n!=m || n==0) throw WrongDimensionException(EXCEPT_MSG(""));
    checkTS(TSARG); vvec.checkTS(TSARG); wkvec.checkTS(TSARG);
    cvec.checkTS(TSARG); svec.checkTS(TSARG);
    if (!islower && strpatt!=MatStrct::upper)
      throw WrongStatusException(EXCEPT_MSG("Need triangular matrix"));
    if (vvec.size()!=n)
      throw WrongDimensionException(EXCEPT_MSG(""));
    if (isMaskInd())
      throw MaskObjectException("Not implemented for indexed masks!");
    if (doDrag) {
      r=dragZ->rows();
      if (dragZ->cols()!=n || dragY==0 || dragY->size()!=r)
	throw WrongDimensionException(EXCEPT_MSG(""));
      dragZ->checkTS(TSARG); dragY->checkTS(TSARG);
    }
    cvec.zeros(n); svec.zeros(n);

    // Generate Givens rotations, update L
    wkvec=vvec;
    ArrayHandle<double> wkarr(wkvec.getFlatBuff());
    stp=islower?1:stride;
    for (i=0,sz=n,tbuff=buff; i<n-1; i++) {
      // drotg(a,b,c,s): J = [c s; -s c], s.t. J [a; b] = [r; 0]
      // a overwritten by r, b by some other information (NOT 0!)
      cblas_drotg(tbuff,wkarr.p()+i,&cvec[i],&svec[i]);
      // Do not want negative elements on factor diagonal
      if ((temp=*tbuff)<0.0) {
	*tbuff=-temp; cvec[i]=-cvec[i]; svec[i]=-svec[i];
      }
      // drot(x,y,c,s): J = [c s; -s c]. [x_i; y_i] overwritten by
      // J [x_i; y_i], for all i
      // BAD: Slower for upper triangular!
      cblas_drot(--sz,tbuff+stp,stp,wkarr.p()+(i+1),1,cvec[i],svec[i]);
      tbuff+=(stride+1);
    }
    cblas_drotg(tbuff,wkarr.p()+(n-1),&cvec[i],&svec[i]);
    if ((temp=*tbuff)<0.0) {
      *tbuff=-temp; cvec[i]=-cvec[i]; svec[i]=-svec[i];
    }

    // Dragging along
    if (doDrag) {
      wkvec=*dragY; wkarr=wkvec.getFlatBuff();
      BaseLinMat<double> zlin;
      dragZ->getLinMat(zlin,true);
      for (i=0; i<n; i++)
	cblas_drot(r,zlin.buff.p()+(i*zlin.stride),1,wkarr.p(),1,cvec[i],
		   svec[i]);
    }
  }

  int StMatrix::cholDndRk1(const StVector& vvec,StVector& cvec,StVector& svec,
			   BaseVector<int>& flind,StVector& wkvec,bool isp,
			   StMatrix* dragZ,const StVector* dragY)
  {
    int i,r=0,stp,sz;
    double qs;
    bool doDrag=(dragZ!=0),islower=(strpatt==MatStrct::lower);
    double* tbuff;

    if (n!=m || n==0) throw WrongDimensionException(EXCEPT_MSG(""));
    checkTS(TSARG); vvec.checkTS(TSARG); wkvec.checkTS(TSARG);
    cvec.checkTS(TSARG); svec.checkTS(TSARG);
    if (!islower && strpatt!=MatStrct::upper)
      throw WrongStatusException(EXCEPT_MSG("Need triangular matrix"));
    if (vvec.size()!=n)
      throw WrongDimensionException(EXCEPT_MSG(""));
    if (isMaskInd())
      throw MaskObjectException("Not implemented for indexed masks!");
    if (doDrag) {
      r=dragZ->rows();
      if (dragZ->cols()!=n || dragY==0 || dragY->size()!=r)
	throw WrongDimensionException(EXCEPT_MSG(""));
      dragZ->checkTS(TSARG); dragY->checkTS(TSARG);
    }
    cvec.zeros(n); svec.zeros(n);

    // Compute p (if not given)
    wkvec=vvec;
    if (!isp)
      backsubst(wkvec,!islower); // p = L^-1 v
    // Generate Givens rotations
    if (((qs=1.0-wkvec.inner(wkvec))<=0.0) ||
	MachDep::isUndef(qs=sqrt(qs)))
      return 1; // Matrix A' not pos. def.
    ArrayHandle<double> wkarr(wkvec.getFlatBuff());
    for (i=n-1; i>=0; i--) {
      cblas_drotg(&qs,wkarr.p()+i,&cvec[i],&svec[i]);
      // 'qs' must remain positive
      if (qs<0.0) {
	qs=-qs; cvec[i]=-cvec[i]; svec[i]=-svec[i];
      }
    }
    // NOTE: 'qs' should be 1 now
    //cout << "qs=" << qs << endl; // DEBUG!

    // Update L
    wkvec.zeros(n);
    flind.fill(0,0);
    stp=islower?1:stride;
    for (i=n-1,sz=0,tbuff=buff+((n-1)*(stride+1)); i>=0; i--) {
      // BAD: Slower for upper triangular!
      cblas_drot(++sz,wkarr.p()+i,1,tbuff,stp,cvec[i],svec[i]);
      // Do not want negative elements on diagonal
      if (*tbuff<0.0) {
	flind.insert(i,1,0);
	FastUtils::smul(tbuff,-1.0,sz,stp);
      }
      tbuff-=(stride+1);
    }
    // NOTE: Should have v in 'wkvec' now
    //if (flind.size()>0) cout << "flind_sz=" << flind.size() << endl; // DEBUG

    // Dragging along
    if (doDrag) {
      double* zcol;
      double cval,sval,c1,c2;
      int nxi=-1,j;
      wkvec=*dragY; wkarr=wkvec.getFlatBuff();
      BaseLinMat<double> zlin;
      dragZ->getLinMat(zlin,true);
      if (flind.size()>0) {
	nxi=flind[0]; j=0;
      }
      for (i=0; i<n; i++) {
	zcol=zlin.buff.p()+(i*zlin.stride);
	cval=cvec[i]; sval=svec[i];
	FastUtils::addsmul(zcol,wkarr.p(),-sval,r);
	if (nxi==i) {
	  if (++j<flind.size()) nxi=flind[j];
	  c1=-1.0/cval; c2=sval;
	} else {
	  c1=1.0/cval; c2=-sval;
	}
	FastUtils::smul(zcol,c1,r);
	if (i<n-1) {
	  FastUtils::smul(wkarr.p(),cval,r);
	  FastUtils::addsmul(wkarr.p(),zcol,c2,r);
	}
      }
    }

    return 0;
  }

  /*
  void StMatrix::dchud(const StVector& xvec,StMatrix* zmat,
		       const StVector* yvec,ArrayHandle<double>& carr,
		       ArrayHandle<double>& sarr)
  {
    int nz=0;

    // Check arguments
    if (n!=m || n==0) throw WrongDimensionException(EXCEPT_MSG(""));
    checkTS(TSARG); xvec.checkTS(TSARG);
    if (strpatt!=MatStrct::upper)
      throw WrongStatusException(EXCEPT_MSG("Need upper triangular matrix"));
    if (isMaskInd())
      throw MaskObjectException("Not implemented for indexed masks!");
    if (xvec.size()!=n) throw WrongDimensionException(EXCEPT_MSG(""));
    if (zmat!=0) {
      nz=zmat->cols();
      if (zmat->rows()!=n || yvec==0 || yvec->size()!=nz)
	throw InvalidParameterException(EXCEPT_MSG(""));
      zmat->checkTS(TSARG); yvec->checkTS(TSARG);
    }
    ArrayHandle<double> cwarr,swarr;
    if (carr==0) cwarr.changeRep(n);
    else {
      if (carr.size()<n) throw WrongDimensionException(EXCEPT_MSG(""));
      cwarr=carr;
    }
    if (sarr==0) swarr.changeRep(n);
    else {
      if (sarr.size()<n) throw WrongDimensionException(EXCEPT_MSG(""));
      swarr=sarr;
    }

    // Call LINPACK DCHUD
    BaseLinVec<double> xlin,ylin;
    BaseLinMat<double> zlin;
    xvec.getLinVec(xlin,false);
    ArrayHandle<double> rho;
    if (nz>0) {
      zmat->getLinMat(zlin,true);
      yvec->getLinVec(ylin,false);
      rho.changeRep(nz); ArrayUtils<double>::fill(rho.p(),-1.0,nz);
    }
    F77_FUNC(dchud,DCHUD) (buff,&stride,&n,xlin.buff.p(),zlin.buff.p(),
			   &(zlin.stride),&nz,ylin.buff.p(),rho.p(),cwarr.p(),
			   swarr.p());
  }
  */

  /*
  int StMatrix::dchdd(const StVector& xvec,StMatrix* zmat,const StVector* yvec,
		      ArrayHandle<double>& carr,ArrayHandle<double>& sarr)
  {
    int nz=0,info;

    // Check arguments
    if (n!=m || n==0) throw WrongDimensionException(EXCEPT_MSG(""));
    checkTS(TSARG); xvec.checkTS(TSARG);
    if (strpatt!=MatStrct::upper)
      throw WrongStatusException(EXCEPT_MSG("Need upper triangular matrix"));
    if (isMaskInd())
      throw MaskObjectException("Not implemented for indexed masks!");
    if (xvec.size()!=n) throw WrongDimensionException(EXCEPT_MSG(""));
    if (zmat!=0) {
      nz=zmat->cols();
      if (zmat->rows()!=n || yvec==0 || yvec->size()!=nz)
	throw InvalidParameterException(EXCEPT_MSG(""));
      zmat->checkTS(TSARG); yvec->checkTS(TSARG);
    }
    ArrayHandle<double> cwarr,swarr;
    if (carr==0) cwarr.changeRep(n);
    else {
      if (carr.size()<n) throw WrongDimensionException(EXCEPT_MSG(""));
      cwarr=carr;
    }
    if (sarr==0) swarr.changeRep(n);
    else {
      if (sarr.size()<n) throw WrongDimensionException(EXCEPT_MSG(""));
      swarr=sarr;
    }

    // Call LINPACK DCHDD
    BaseLinVec<double> xlin,ylin;
    BaseLinMat<double> zlin;
    xvec.getLinVec(xlin,false);
    ArrayHandle<double> rho;
    if (nz>0) {
      zmat->getLinMat(zlin,true);
      yvec->getLinVec(ylin,false);
      rho.changeRep(nz); ArrayUtils<double>::fill(rho.p(),-1.0,nz);
    }
    F77_FUNC(dchdd,DCHDD) (buff,&stride,&n,xlin.buff.p(),zlin.buff.p(),
			   &(zlin.stride),&nz,ylin.buff.p(),rho.p(),cwarr.p(),
			   swarr.p(),&info);

    return info;
  }
  */

  int StMatrix::cholext(const StVector& avec,double ascal,bool isl,
			ArrayHandle<StMatrix*>& dragX1,
			ArrayHandle<StVector*>& dragB1,
			ArrayHandle<StMatrix*>& dragX2,
			ArrayHandle<StVector*>& dragB2)
  {
    int i,numm,oldn;
    double lscal;
    bool isr;

    if (n!=m || n==0) throw WrongDimensionException(EXCEPT_MSG(""));
    if (isMask()) throw MaskObjectException(EXCEPT_MSG(""));
    if (strpatt==MatStrct::lower) isr=false;
    else if (strpatt==MatStrct::upper) isr=true;
    else
      throw WrongStatusException(EXCEPT_MSG("Need triangular matrix"));
    CHECKVEC(avec,n);
    numm=dragX1.size();
    if (numm>0 && dragB1.size()!=numm)
      throw WrongDimensionException(EXCEPT_MSG(""));
    for (i=0; i<numm; i++) {
      CHECKMAT(*dragX1[i],dragB1[i]->size(),n);
      dragB1[i]->checkTS(TSARG);
    }
    numm=dragX2.size();
    if (numm>0 && dragB2.size()!=numm)
      throw WrongDimensionException(EXCEPT_MSG(""));
    for (i=0; i<numm; i++) {
      CHECKMAT(*dragX2[i],dragB2[i]->size(),n);
      dragB2[i]->checkTS(TSARG);
    }

    Handle<StVector> lvec;
    if (!isl) {
      lvec.changeRep(new StVector(avec));
      backsubst(*lvec,isr); // L^-1 a (R^-T a)
    } else
      lvec.changeRep((StVector*) &avec,false);
    lscal=ascal-lvec->inner(*lvec);
    if (lscal<=0.0) return 1; // new A not pos. def.
    lscal=sqrt(lscal);
    oldn=n;
    expand(n+1,n+1);
    Handle<StVector> lmsk;
    MASKVEC(lmsk,oldn,Range(0,oldn-1));
    if (!isr)
      *lmsk=*lvec;
    else
      lmsk->zeros();
    MASKVEC(lmsk,Range(0,oldn-1),oldn);
    if (!isr)
      lmsk->zeros();
    else
      *lmsk=*lvec;
    set(oldn,oldn,lscal); // L (R) updated
    Handle<StMatrix> mmsk;
    // Dragging along, case 1: X = B L^-T
    for (i=0; i<dragX1.size(); i++) {
      dragX1[i]->expand(dragX1[i]->rows(),oldn+1);
      MASKMMAT(mmsk,*dragX1[i],RangeFull::get(),Range(0,oldn-1)); // old M
      MASKMVEC(lmsk,*dragX1[i],RangeFull::get(),oldn); // new M col.
      *lmsk=*dragB1[i];
      mmsk->mulVec(*lmsk,*lvec,false,-1.0,1.0);
      lmsk->prod(1.0/lscal);
    }
    // Dragging along, case 2: X = B L
    // NOT TESTED!!
    for (i=0; i<dragX2.size(); i++) {
      dragX2[i]->expand(dragX2[i]->rows(),oldn+1);
      MASKMMAT(mmsk,*dragX2[i],RangeFull::get(),Range(0,oldn-1)); // old M
      MASKMVEC(lmsk,*dragX2[i],RangeFull::get(),oldn); // new M col.
      lmsk->prod(*dragB2[i],lscal);
      mmsk->rankOne(*dragB2[i],1.0,*lvec);
    }

    return 0;
  }

  int StMatrix::cholext(const StMatrix& avecs,const StMatrix& ascals,bool isl,
			ArrayHandle<StMatrix*>& dragX1,
			ArrayHandle<StMatrix*>& dragB1,
			ArrayHandle<StMatrix*>& dragX2,
			ArrayHandle<StMatrix*>& dragB2)
  {
    int i,numm,oldn,d;
    bool isr;

    if (n!=m || n==0) throw WrongDimensionException(EXCEPT_MSG(""));
    if (isMask()) throw MaskObjectException(EXCEPT_MSG(""));
    if (strpatt==MatStrct::lower) isr=false;
    else if (strpatt==MatStrct::upper) isr=true;
    else
      throw WrongStatusException(EXCEPT_MSG("Need triangular matrix"));
    CHECKMAT(avecs,n,-1);
    d=avecs.cols();
    CHECKMAT(ascals,d,d);
    ascals.updateStrctPatt();
    if (ascals.getStrctPatt()!=MatStrct::lower &&
	ascals.getStrctPatt()!=MatStrct::upper)
      throw WrongStatusException(EXCEPT_MSG(""));
    numm=dragX1.size();
    if (numm>0 && dragB1.size()!=numm)
      throw WrongDimensionException(EXCEPT_MSG(""));
    for (i=0; i<numm; i++) {
      CHECKMAT(*dragX1[i],dragB1[i]->rows(),n);
      CHECKMAT(*dragB1[i],-1,d);
    }
    numm=dragX2.size();
    if (numm>0 && dragB2.size()!=numm)
      throw WrongDimensionException(EXCEPT_MSG(""));
    for (i=0; i<numm; i++) {
      CHECKMAT(*dragX2[i],dragB2[i]->rows(),n);
      CHECKMAT(*dragB2[i],-1,d);
    }

    Handle<StMatrix> lvecs;
    if (!isl) {
      lvecs.changeRep(new StMatrix(avecs));
      backsubst(*lvecs,true,isr);
    } else
      lvecs.changeRep((StMatrix*) &avecs,false);
    StMatrix lscals(ascals); lscals.setStrctPatt(ascals.getStrctPatt());
    lscals.makeSymmAuto(); lscals.setStrctPatt(MatStrct::lower);
    lscals.symMul(*lvecs,true,-1.0,1.0);
    if (lscals.cholDecomp(true)!=0) return 1; // new A not pos. def.
    oldn=n;
    expand(n+d,n+d);
    Handle<StMatrix> lmsk;
    Range rng1(oldn),rng2(0,oldn-1);
    MASKMAT(lmsk,rng1,rng2);
    if (!isr)
      lmsk->trans(*lvecs);
    else
      lmsk->zeros();
    MASKMAT(lmsk,rng2,rng1);
    if (!isr)
      lmsk->zeros();
    else
      *lmsk=*lvecs;
    MASKMAT(lmsk,rng1,rng1);
    if (!isr)
      *lmsk=lscals;
    else
      lmsk->trans(lscals);
    // Dragging along, case 1: X = B L^-T
    Handle<StMatrix> mmsk;
    for  (i=0; i<dragX1.size(); i++) {
      dragX1[i]->expand(dragX1[i]->rows(),oldn+d);
      MASKMMAT(mmsk,*dragX1[i],RangeFull::get(),rng2); // old M
      MASKMMAT(lmsk,*dragX1[i],RangeFull::get(),rng1); // new cols of M
      *lmsk=*dragB1[i];
      lmsk->mul(*mmsk,*lvecs,-1.0,1.0);
      lscals.backsubst(*lmsk,false,true);
    }
    // Dragging along, case 2: X = B L
    // NOT TESTED!!
    for  (i=0; i<dragX2.size(); i++) {
      dragX2[i]->expand(dragX2[i]->rows(),oldn+d);
      MASKMMAT(mmsk,*dragX2[i],RangeFull::get(),rng2); // old M
      MASKMMAT(lmsk,*dragX2[i],RangeFull::get(),rng1); // new cols of M
      *lmsk=*dragB2[i];
      lscals.triMul(*lmsk,false,false);
      mmsk->mul(*dragB2[i],*lvecs,1.0,1.0);
    }

    return 0;
  }

#define CHOLUPD_CHECKDRAGARRAY(arr,n) do { for (i=0; i<(arr).size(); i++) { \
  (arr)[i]->checkTS(TSARG); \
  if ((arr)[i]->cols()!=n) throw WrongDimensionException(EXCEPT_MSG("")); \
  } } while (0)

  int StMatrix::cholUpdRk1_old(const StVector& vec,double alpha,bool isp,
			   ArrayHandle<StMatrix*>& dragX1,
			   ArrayHandle<StMatrix*>& dragX2,StVector* wkvec,
			   ArrayHandle<double>& wkarrn1,
			   ArrayHandle<double>& wkarrn2)
  {
    int i,j,numb,sz;
    double temp,temp2,t,dlam,b,s;
    double* colB,*lamB,*betaB;
    const double* pB;

    if (alpha==0.0) return 0; // nothing todo
    if (n!=m || n==0) throw WrongDimensionException(EXCEPT_MSG(""));
    checkTS(TSARG);
    if (strpatt!=MatStrct::lower)
      throw WrongStatusException(EXCEPT_MSG("Need lower triangular matrix"));
    vec.checkTS(TSARG);
    if (vec.size()!=n)
      throw WrongDimensionException(EXCEPT_MSG(""));
    if (isMaskInd())
      throw MaskObjectException("Not implemented for indexed masks!");
    if (alpha==0.0) return 0;
    CHOLUPD_CHECKDRAGARRAY(dragX1,n);
    CHOLUPD_CHECKDRAGARRAY(dragX2,n);
    // Working arrays
    Handle<StVector> hwvec;
    if (wkvec==0) {
      hwvec.changeRep(new StVector(n));
      wkvec=hwvec;
    } else
      wkvec->zeros(n);
    ArrayHandle<double> warr,betaArr;
    if (wkarrn1==0) warr.changeRep(n);
    else {
      if (wkarrn1.size()<n) throw WrongDimensionException(EXCEPT_MSG(""));
      warr=wkarrn1;
    }
    if (wkarrn2==0) betaArr.changeRep(n);
    else {
      if (wkarrn2.size()<n) throw WrongDimensionException(EXCEPT_MSG(""));
      betaArr=wkarrn2;
    }
    if (alpha<0.0)
      return cholUpdRk1Negative_old(vec,alpha,isp,dragX1,dragX2,wkvec,warr,
				    betaArr);

    // Backsubstitution (if 'pVec' not given)
    Handle<StVector> pvHand;
    const StVector* pVec;
    if (isp) pVec=&vec;
    else {
      pvHand.changeRep(new StVector(vec)); // copy of 'vec'
      backsubst(*pvHand);
      pVec=pvHand;
    }
    ArrayHandle<double> plin(pVec->getFlatBuff()); // flat version

    // I + alpha p p^T = \tilde{L} \tilde{\Lambda} \tilde{L}^T. Here,
    // \tilde{L} is represented by p and \beta
    ArrayHandle<double> lamArr(wkvec->getFlatBuff());
    t=1.0/alpha;
    pB=plin.p(); betaB=(double*) betaArr.p();
    lamB=(double*) lamArr.p();
    for (j=0; j<n-1; j++) {
      temp=*(pB++);
      temp2=t+temp*temp;
      if ((*(lamB++)=temp2/t)<1e-30) return 1; // Numerical breakdown
      t=temp2;
      *(betaB++)=temp/temp2;
    }
    temp=*pB;
    if ((*(lamB)=1.0+temp*temp/t)<1e-30) return 1; // Numerical breakdown

    // Update of L. Note that \tilde{l}_{i,j} = p_i \beta_j for i>j and
    // \tilde{l}_{i,i}=1
    wkvec->apply1(ptr_fun<double,double>(sqrt)); // sqrt
    ArrayUtils<double>::fill(warr.p(),0.0,n);
    pB=plin.p()+(n-1); betaB=(double*) betaArr.p()+(n-2);
    lamB=(double*) lamArr.p()+(n-1);
    warr[n-1]=(temp=this->get(n-1,n-1))*(*(pB--));
    double* vfrom=warr.p()+(n-2);
    this->set(n-1,n-1,temp*(*(lamB--)));
    for (j=n-2,sz=2; j>=0; j--,sz++) {
      colB=buff+(j*(stride+1)); // (j,j)
      /*
      if (j>0) FastUtils::addsmul(vto,vfrom,colB,*(pB--),sz);
      FastUtils::addsmul(colB+1,vfrom+1,*(betaB--),sz-1);
      FastUtils::smul(colB,*(lamB--),sz);
      vexch=vto; vto=vfrom-1; vfrom=vexch-1; // Exchange and decr.
      */
      // Avoid slow 'addsmul':
      temp=*(betaB--); temp2=*(pB--);
      FastUtils::addsmul(colB+1,vfrom+1,temp,sz-1);
      if (j>0) {
	FastUtils::smul(vfrom+1,1.0-temp*temp2,sz-1);
	FastUtils::addsmul(vfrom--,colB,temp2,sz);
      }
      FastUtils::smul(colB,*(lamB--),sz);
    }

    // Dragging along, case 1: X = B L^-T
    if ((numb=dragX1.size())>0) {
      // Need buffer vector
      int maxm=0;
      for (i=0; i<numb; i++)
	if ((j=dragX1[i]->rows())>maxm) maxm=j;
      if (wkarrn1.size()<maxm) wkarrn1.changeRep(maxm);
      for (i=0; i<numb; i++) {
	StMatrix* bAct=dragX1[i];
	bAct->setStrctPatt(MatStrct::normal);
	{
	  BaseLinMat<double> blin;
	  bAct->getLinMat(blin,true);
	  // Backsubstitution with \tilde{L}
	  betaB=(double*) betaArr.p();
	  pB=plin.p();
	  int bm=bAct->rows();
	  double* sumB=(double*) wkarrn1.p(),*colB;
	  ArrayUtils<double>::fill(sumB,0.0,bm);
	  colB=blin.buff;
	  for (j=1; j<n; j++) {
	    FastUtils::addsmul(sumB,colB,*(betaB++),bm);
	    colB+=blin.stride; pB++;
	    FastUtils::addsmul(colB,sumB,-(*pB),bm);
	  }
	}
	// blin -> bAct written back here
	bAct->mulDiag(*bAct,*wkvec,true); // right-mult. with D^-1/2
      }
    }

    // NOT TESTED!!
    // Dragging along, case 2: X = B L
    if ((numb=dragX2.size())>0) {
      // Need buffer vectors
      int maxm=0;
      for (i=0; i<numb; i++)
	if ((j=dragX2[i]->rows())>maxm) maxm=j;
      if (wkarrn1.size()<maxm) wkarrn1.changeRep(maxm);
      if (wkarrn2.size()<maxm) wkarrn2.changeRep(maxm);
      for (i=0; i<numb; i++) {
	StMatrix* bAct=dragX2[i];
	bAct->setStrctPatt(MatStrct::normal);
	{
	  BaseLinMat<double> blin;
	  bAct->getLinMat(blin,true);
	  // Leftmult. with \tilde{L}^T
	  betaB=((double*) betaArr.p())+(n-1);
	  pB=plin.p()+(n-1);
	  int bm=bAct->rows();
	  double* sumB=(double*) wkarrn1.p(),
	    *rhoB=(double*) wkarrn2.p(),*colB;
	  ArrayUtils<double>::fill(sumB,0.0,bm);
	  colB=blin.buff+((n-1)*blin.stride);
	  FastUtils::smul(rhoB,colB,*(pB--),bm);
	  for (j=n-2; j>=0; j--) {
	    FastUtils::addsmul(sumB,rhoB,1.0,bm);
	    colB-=blin.stride;
	    FastUtils::smul(rhoB,colB,*(pB--),bm);
	    betaB--;
	    FastUtils::addsmul(colB,sumB,*betaB,bm);
	  }
	}
	// blin -> bAct written back here
	bAct->mulDiag(*bAct,*wkvec); // right-mult. with D^1/2
      }
    }

    return 0; // OK
  }

  /*
   * Method of Stewart for Cholesky downdating (C3 in Gill et al), i.e.
   * alpha < 0.
   *
   * U^T U = R^T R - v v^T for upper triangular.
   * Find orthogonal Q s.t. Q [ 0^T; R ] = [v^T; U].
   * Here, Q = J_1 ... J_n, J_k Givens rotation on comp. 0,k. The J_k are
   * determined by the relation
   *   J_1 ... J_n q = delta_1, q = [beta b],
   * where R^T p = v, b = (-alpha)^(1/2) p, beta = sqrt(1 - b^T b).
   * J_k is repres. by scalars c,s with c^2 + s^2 = 1.
   */
  int StMatrix::cholUpdRk1Negative_old(const StVector& vec,double alpha,
				       bool isp,
				   ArrayHandle<StMatrix*>& dragX1,
				   ArrayHandle<StMatrix*>& dragX2,
				   StVector* wkvec,
				   ArrayHandle<double>& wkarrn1,
				   ArrayHandle<double>& wkarrn2)
  {
    int i,j,k,numb,bm;
    double temp,q,qs,gamma,salph,cscal,sscal,xelem;
    double* xP,*colP,*col2P;

    // No checks done here!
    ArrayHandle<double> carr(wkvec->getFlatBuff()),sarr(wkarrn2);
    ArrayHandle<double> workVec(wkarrn1);

    // Backsubstitution (if 'pVec' not given): L p = v
    Handle<StVector> pvHand;
    const StVector* pVec;
    if (isp) pVec=&vec;
    else {
      pvHand.changeRep(new StVector(vec)); // copy of 'vec'
      backsubst(*pvHand);
      pVec=pvHand;
    }
    ArrayHandle<double> plin(pVec->getFlatBuff()); // flat version

    // Compute Givens rotations c_i, s_i
    // TODO: Use BLAS drotg !!!
    qs=1.0+alpha*pVec->inner(*pVec);
    if (qs<=0.0) return 1; // res. matrix not positive definite
    q=sqrt(qs);
    salph=sqrt(-alpha);
    for (i=n-1; i>=0; i--) {
      temp=plin[i]*salph;
      qs+=(temp*temp);
      gamma=sqrt(qs);
      carr[i]=q/gamma; sarr[i]=temp/gamma;
      q=gamma;
    }

    // Update of L
    ArrayUtils<double>::fill(workVec,0.0,n);
    for (i=n-1; i>=0; i--)
      cblas_drot(n-i,workVec.p()+i,1,buff+(i*(stride+1)),1,carr[i],sarr[i]);

    // Dragging along, case 1: X = B L^-T
    j=0;
    for (i=0; i<dragX1.size(); i++)
      j=std::max(j,dragX1[i]->rows());
    for (i=0; i<dragX2.size(); i++)
      j=std::max(j,dragX2[i]->rows());
    if (j>wkarrn1.size()) wkarrn1.changeRep(j);
    if ((numb=dragX1.size())>0) {
      for (i=0; i<numb; i++) {
	StMatrix* bAct=dragX1[i];
	bAct->setStrctPatt(MatStrct::normal);
	BaseLinMat<double> blin;
	bAct->getLinMat(blin,true);
	bm=bAct->rows();
	ArrayUtils<double>::fill(wkarrn1,0.0,bm); // beta values
	col2P=blin.buff;
	for (j=0; j<n; j++) {
	  cscal=carr[j]; sscal=sarr[j];
	  colP=col2P;
	  xP=wkarrn1.p();
	  for (k=0; k<bm; k++) {
	    *(colP++)=temp=((*colP)-sscal*(xelem=*xP))/cscal;
	    *(xP++)=cscal*xelem-sscal*temp;
	  }
	  col2P+=blin.stride;
	}
      }
    }

    // NOT TESTED!!
    // Dragging along, case 2: X = B L
    if ((numb=dragX2.size())>0) {
      for (i=0; i<numb; i++) {
	StMatrix* bAct=dragX2[i];
	bAct->setStrctPatt(MatStrct::normal);
	BaseLinMat<double> blin;
	bAct->getLinMat(blin,true);
	bm=bAct->rows();
	ArrayUtils<double>::fill(wkarrn1,0.0,bm); // beta values
	col2P=blin.buff+((n-1)*blin.stride);
	for (j=n-1; j>=0; j--) {
	  cscal=carr[j]; sscal=sarr[j];
	  colP=col2P;
	  xP=wkarrn1.p();
	  for (k=0; k<bm; k++) {
	    *(colP++)=cscal*(temp=*colP)-sscal*(xelem=*xP);
	    *(xP++)=cscal*xelem+sscal*temp;
	  }
	  col2P-=blin.stride;
	}
      }
    }

    return 0; // OK
  }

  /*
   * Method of Goldfarb for indef. rank 2 update
   *
   * If w = L^T u, z = L^-1 v, then
   *  A' = L (I + z w^T) (I + w z^T) L^T.
   * Now, I + z w^T = \tilde{L} Q^T, Q orthonormal, \tilde{L} lower
   * triangular, therefore L' = L \tilde{L}. Q is the product of 2(n-1)
   * Givens rotations. Only \tilde{L} is required.
   * \tilde{L} has a simple form depending on O(n) parameters, so updating
   * L is O(n^2), dragging along is O(n) per vector.
   */
  /*
  int StMatrix::cholUpdRk2Goldfarb(const StVector& uvec,const StVector& vvec,
				   bool isw,bool isz,
				   ArrayHandle<StMatrix*>& dragX1,
				   ArrayHandle<StMatrix*>& dragX2)
  {
    int i,j,k,sz,numb,bm;
    double temp,cscal,sscal,belem,bargam,barlam,lamelem,rho,welem;
    double* colP,*betaP,*sP,*gammaP,*lamP,*sbP,*sgP,*wvP;
    const double* wP,*zP;

    if (n!=m || n==0) throw WrongDimensionException(EXCEPT_MSG(""));
    checkTS(TSARG);
    if (strpatt!=MatStrct::lower)
      throw WrongStatusException(EXCEPT_MSG("Need lower triangular matrix"));
    uvec.checkTS(TSARG); vvec.checkTS(TSARG);
    if (uvec.size()!=n || vvec.size()!=n)
      throw WrongDimensionException(EXCEPT_MSG(""));
    if (isMaskInd())
      throw MaskObjectException("Not implemented for indexed masks!");
    // Allocate workspace
    ArrayHandle<double> beta(n),gamma(n),lambda(n),workVec(n);
    j=n;
    CHOLUPD_CHECKDRAGARRAY(dragX1,n);
    for (i=0; i<dragX1.size(); i++)
      if ((k=dragX1[i]->rows())>j) j=k;
    CHOLUPD_CHECKDRAGARRAY(dragX2,n);
    for (i=0; i<dragX2.size(); i++)
      if ((k=dragX2[i]->rows())>j) j=k;
    ArrayHandle<double> sbeta(j),sgamma(j);

    // Compute w,z (if not given)
    Handle<StVector> wvHand,zvHand;
    const StVector* wVec,*zVec;
    if (isw) wVec=&uvec;
    else {
      wvHand.changeRep(new StVector(uvec)); // copy of 'uvec'
      wVec=wvHand;
      triMulVec((StVector&) *wVec,true); // w = L^T u
    }
    if (isz) zVec=&vvec;
    else {
      zvHand.changeRep(new StVector(vvec)); // copy of 'vvec'
      zVec=zvHand;
      backsubst((StVector&) *zVec); // z = L^-1 v
    }
    ArrayHandle<double> wlin(wVec->getFlatBuff()); // flat versions
    ArrayHandle<double> zlin(zVec->getFlatBuff());

    // Skip trailing elem. of w close to 0 --> use 'sz' inst. of 'n'
    temp=DOUBLE_MACHPRECISION*wVec->norm2();
    for (sz=n-1; sz>=0 && fabs(wlin[sz])<=temp; sz--);
    if (sz<0)
      return 0; // w==0: nothing to do
    else
      sz++; // size of non-zero w prefix
    if (sz<n) {
      ArrayUtils<double>::fill(((double*) beta.p())+sz,0.0,n-sz);
      ArrayUtils<double>::fill(((double*) gamma.p())+sz,0.0,n-sz);
      ArrayUtils<double>::fill(((double*) lambda.p())+sz,1.0,n-sz);
    }

    // Recurrence 1 (first sz-1 Givens rotations)
    // \bar{beta} in 'beta', s in 'gamma'. c not required
    betaP=((double*) beta.p())+(sz-1);
    sP=((double*) gamma.p())+(sz-1);
    wP=wlin.p()+(sz-1);
    belem=1.0/(*(wP--));
    for (i=sz-2; i>=0; i--) {
      temp=(*(wP--))*belem; // r_i
      *(--sP)=sscal=1.0/sqrt(1.0+temp*temp);
      belem*=sscal;
      *(betaP--)=temp*belem;
    }
    *betaP=belem;
    // Recurrence 2 (2nd sz-1 Givens rotations)
    betaP=(double*) beta.p();
    gammaP=(double*) gamma.p();
    lamP=(double*) lambda.p();
    wP=wlin.p(); zP=zlin.p();
    // 'belem'==\bar{beta}, 'bargam'==\bar{gamma}, 'barlam'==\bar{lambda}
    if ((belem=*betaP)==0.0) return 1; // breakdown
    if (MachDep::isUndef(bargam=1.0/belem)) return 1; // breakdown
    for (i=0; i<sz-1; i++) {
      barlam=belem*(*(wP++))+bargam*(*(zP++));
      sscal=*gammaP; // s
      if (fabs(sscal)>fabs(barlam)) {
	temp=barlam/fabs(sscal); rho=sqrt(1.0+temp*temp);
	lamelem=fabs(sscal)*rho;
	sscal=(sscal<0.0)?(-1.0/rho):(1.0/rho);
	cscal=temp/rho;
      } else {
	temp=sscal/fabs(barlam); rho=sqrt(1.0+temp*temp);
	lamelem=fabs(barlam)*rho;
	cscal=(barlam<0.0)?(-1.0/rho):(1.0/rho);
	sscal=temp/rho;
      }
      if (MachDep::isUndef(*(lamP++)=lamelem)) return 1; // breakdown
      temp=*(betaP+1); // next \bar{beta}
      if (MachDep::isUndef(*(betaP++)=cscal*belem-sscal*temp))
	return 1; // breakdown
      if (MachDep::isUndef(*(gammaP++)=cscal*bargam)) return 1; // breakdown
      bargam*=sscal;
      belem=sscal*belem+cscal*temp; // new value for 'belem'
    }
    *lamP=belem*(*wP)+bargam*(*zP);

    // Update of L
    betaP=((double*) beta.p())+(n-1);
    gammaP=((double*) gamma.p())+(n-1);
    lamP=((double*) lambda.p())+(n-1);
    wP=wlin.p()+(n-1); zP=zlin.p()+(n-1);
    ArrayUtils<double>::fill((double*) sbeta.p(),0.0,n); // reset accus
    ArrayUtils<double>::fill((double*) sgamma.p(),0.0,n);
    sbP=((double*) sbeta.p())+(n-1);
    sgP=((double*) sgamma.p())+(n-1);
    wvP=((double*) workVec.p())+(n-1); // column copies
    colP=buff+((n-1)*(stride+1)); // (n-1,n-1)
    *wvP=*colP; (*colP)*=(*(lamP--));
    for (j=n-2; j>=0; j--) {
      FastUtils::addsmul(sbP,wvP,*(wP--),n-j-1);
      FastUtils::addsmul(sgP,wvP,*(zP--),n-j-1);
      wvP--;
      colP-=(stride+1); // (j,j)
      ArrayUtils<double>::copy(wvP,colP,n-j); // copy column
      FastUtils::smul(colP,*(lamP--),n-j);
      FastUtils::addsmul(colP+1,sbP,*(--betaP),n-j-1);
      FastUtils::addsmul(colP+1,sgP,*(--gammaP),n-j-1);
      sbP--; sgP--;
    }

    // Dragging along, case 1: X = B L^-T
    if ((numb=dragX1.size())>0) {
      for (i=0; i<numb; i++) {
	StMatrix* bAct=dragX1[i];
	bAct->setStrctPatt(MatStrct::normal);
	BaseLinMat<double> blin;
	bAct->getLinMat(blin,true);
	bm=bAct->rows();
	sbP=(double*) sbeta.p(); sgP=(double*) sgamma.p();
	ArrayUtils<double>::fill(sbP,0.0,bm); // reset accus
	ArrayUtils<double>::fill(sgP,0.0,bm);
	colP=blin.buff;
	FastUtils::smul(colP,1.0/lambda[0],bm);
	for (j=1; j<n; j++) {
	  FastUtils::addsmul(sbP,colP,beta[j-1],bm);
	  FastUtils::addsmul(sgP,colP,gamma[j-1],bm);
	  colP+=blin.stride;
	  FastUtils::addsmul(colP,sbP,-wlin[j],bm);
	  FastUtils::addsmul(colP,sgP,-zlin[j],bm);
	  FastUtils::smul(colP,1.0/lambda[j],bm);
	}
      }
    }

    // Dragging along, case 2: X = B L
    if ((numb=dragX2.size())>0) {
      for (i=0; i<numb; i++) {
	StMatrix* bAct=dragX2[i];
	bAct->setStrctPatt(MatStrct::normal);
	BaseLinMat<double> blin;
	bAct->getLinMat(blin,true);
	bm=bAct->rows();
	sbP=(double*) sbeta.p(); sgP=(double*) sgamma.p();
	ArrayUtils<double>::fill(sbP,0.0,bm); // reset accus
	ArrayUtils<double>::fill(sgP,0.0,bm);
	colP=blin.buff+((n-1)*blin.stride);
	FastUtils::smul(colP,lambda[n-1],bm);
	for (j=n-2; j>=0; j--) {
	  FastUtils::addsmul(sbP,colP,wlin[j+1],bm);
	  FastUtils::addsmul(sgP,colP,zlin[j+1],bm);
	  colP-=blin.stride;
	  FastUtils::smul(colP,lambda[j],bm);
	  FastUtils::addsmul(colP,sbP,beta[j],bm);
	  FastUtils::addsmul(colP,sgP,gamma[j],bm);
	}
      }
    }

    return 0; // OK
  }
  */
#undef CHOLUPD_CHECKDRAGARRAY

  /*
   * There is a bug in LINPACK DCHEX. In some cases, the returned factor
   * R' has R'(j,j) < 0 for some k<=j<=l. In this case, row j of R' (and of
   * X' if given) has to be multiplied by -1.
   */
  void StMatrix::cholUpdExch(int k,int l,int job,StMatrix* dragX)
  {
#ifndef HAVE_FORTRAN
    throw NotImplemException(EXCEPT_MSG("NO EXTERNAL FORTRAN CODE"));
#else
    int farg1,farg2,farg3,j;

    // Check arguments
    if (n!=m || n==0) throw WrongDimensionException(EXCEPT_MSG(""));
    checkTS(TSARG);
    if (strpatt!=MatStrct::upper)
      throw WrongStatusException(EXCEPT_MSG("Need upper triangular matrix"));
    if (isMaskInd())
      throw MaskObjectException("Not implemented for indexed masks!");
    if (k<0 || k>=l || l>=n)
      throw InvalidParameterException(EXCEPT_MSG(""));
    if (job<1 || job>2)
      throw InvalidParameterException(EXCEPT_MSG(""));
    if (dragX!=0) {
      dragX->checkTS(TSARG);
      if (dragX->rows()!=n)
	throw WrongDimensionException(EXCEPT_MSG(""));
    }
    // Allocate workspace
    ArrayHandle<double> cvec(n),svec(n);
    // Call LINPACK DCHEX. Work around bug
    farg1=k+1; farg2=l+1;
    if (dragX!=0) {
      BaseLinMat<double> xlin;
      dragX->getLinMat(xlin,true);
      F77_FUNC(dchex,DCHEX) (buff,&stride,&n,&farg1,&farg2,xlin.buff.p(),
			     &(xlin.stride),&(xlin.n),cvec.p(),svec.p(),&job);
    } else {
      farg3=0;
      F77_FUNC(dchex,DCHEX) (buff,&stride,&n,&farg1,&farg2,0,&farg3,&farg3,
			     cvec.p(),svec.p(),&job);
    }
    Handle<StVector> msk;
    for (j=k; j<=l; j++)
      if (get(j,j)<0) {
	// Work around bug
	Handle<StVector> msk;
	MASKVEC(msk,j,Range(j,n-1));
	msk->prod(-1.0);
	if (dragX!=0) {
	  MASKMVEC(msk,*dragX,j,RangeFull::get());
	  msk->prod(-1.0);
	}
      }
#endif
  }

  /*
   * OLD CODE!! NOT STABLE!
   * We update M explicitly alongside its CF (in 'cholM') and refresh the CF
   * from M after each 'd' iterations.
   */
  /*
  int StMatrix::cholUpdRkd(const StMatrix& vmat,double alpha,bool isp,
			   const StMatrix* cholM,
			   ArrayHandle<StMatrix*>& bMat,StMatrix* buffMat)
  {
    int i,j,numb,d,count;
    double temp,dlam;
    double* rowB,*lamB;

    if (n!=m || n==0) throw WrongDimensionException(EXCEPT_MSG(""));
    checkTS(TSARG);
    if (strpatt!=MatStrct::lower)
      throw WrongStatusException(EXCEPT_MSG("Need lower triangular matrix"));
    vmat.checkTS(TSARG);
    if (vmat.rows()!=n)
      throw WrongDimensionException(EXCEPT_MSG(""));
    d=vmat.cols();
    if (cholM!=0 && (cholM->rows()!=d || cholM->cols()!=d))
      throw WrongDimensionException(EXCEPT_MSG(""));
    if (isMaskInd())
      throw NotImplemException(EXCEPT_MSG("Not for indexed mask"));
    if (buffMat!=0 && buffMat->isMask())
      throw MaskObjectException(EXCEPT_MSG("buffMat"));
    if (alpha==0.0) return 0;
    double sign=(alpha>0.0)?1.0:-1.0; // positive or negative update?
    alpha=sqrt(fabs(alpha));
    numb=bMat.size();
    for (i=0; i<numb; i++) {
      bMat[i]->checkTS(TSARG);
      if (bMat[i]->cols()!=n)
	throw WrongDimensionException(EXCEPT_MSG("bMat"));
    }

    // Initialization
    StMatrix qfact(d),mmat(d);
    if (cholM!=0) {
      qfact.smul(*cholM,alpha);
      mmat.setStrctPatt(MatStrct::lower);
      mmat.symMul(qfact,false); // |alpha| M
    } else {
      qfact.zeros(); ((StVector&) qfact.diag())=alpha;
      mmat.zeros(); ((StVector&) mmat.diag())=alpha*alpha;
      mmat.setStrctPatt(MatStrct::lower);
    }
    qfact.setStrctPatt(MatStrct::lower);
    // Backsubstitution (if 'vmat' does not contain P)
    Handle<StMatrix> pvHand;
    const StMatrix* pMat;
    if (isp)
      pMat=&vmat;
    else {
      pvHand.changeRep(new StMatrix(vmat)); // copy of 'vmat'
      backsubst(*pvHand,true);
      pMat=pvHand;
    }

    // Build \tilde{L}, D' (note: D==I)
    StVector zvec(d),tempvec(d);
    ArrayHandle<Handle<StVector> > pMsk(n);
    ArrayHandle<StVector> beta(n);
    ArrayHandle<BaseLinVec<double> > pArr(n);
    ArrayHandle<BaseLinVec<double> > betaArr(n);
    for (i=0; i<n; i++) {
      beta[i].zeros(d);
      beta[i].getLinVec(betaArr[i],false);
      pMsk[i].changeRep(DYNCAST(StVector,
				pMat->subrvalInt(i,RangeFull::get())));
      pMsk[i]->getLinVec(pArr[i],false);
    }
    StVector lamDiag(n);
    ArrayHandle<StMatrix*> dummyArr; // for call to 'cholUpdRk1'
    for (j=0,count=0; j<n-1; j++) {
      zvec=*pMsk[j]; qfact.triMulVec(zvec,true); // z
      if ((lamDiag[j]=dlam=1.0+sign*zvec.inner(zvec))<1e-30)
	return 1; // Numerical break-down
      beta[j].prod(zvec,sign/dlam); qfact.triMulVec(beta[j],false);
      // Update of M
      mmat.mulVec(tempvec,*pMsk[j]); // M p_i
      mmat.rankOne(tempvec,-sign/dlam);
      // Update of Q
      if (++count>d) {
	// Refresh
	qfact=mmat;
	if (qfact.cholDecomp(true)!=0) return 1; // break-down
	count=0;
      } else {
	if (qfact.cholUpdRk1(zvec,-sign/dlam,true,dummyArr)!=0)
	  return 1; // break-down
      }
    }
    zvec=*pMsk[n-1]; qfact.triMulVec(zvec,true); // z
    if ((lamDiag[n-1]=dlam=1.0+sign*zvec.inner(zvec))<1e-30)
      return 1; // Numerical break-down

    // Update of L
    lamDiag.apply1(ptr_fun(sqrt)); // sqrt
    ArrayHandle<double> lamDA(lamDiag.getFlatBuff());
    ArrayHandle<double> sigma(tempvec.getFlatBuff()),rho(zvec.getFlatBuff());
    (*buff)*=(*lamDA);
    for (i=1; i<n; i++) {
      lamB=(double*) lamDA.p()+i;
      rowB=buff+(i*(stride+1)); // (i,i)
      FastUtils::smul(rho.p(),pArr[i].getBuff(),*rowB,d,1,pArr[i].step); // rho
      ArrayUtils<double>::fill(sigma.p(),0.0,d);
      (*rowB)*=(*lamB);
      for (j=i-1; j>=0; j--) {
	FastUtils::addsmul(sigma.p(),rho.p(),1.0,d);
	rowB-=stride; // (i,j)
	temp=*rowB;
	FastUtils::smul(rho.p(),pArr[j].getBuff(),temp,d,1,pArr[j].step);
	*rowB=(temp+FastUtils::inner(sigma.p(),betaArr[j].getBuff(),d))*
	  (*(--lamB));
      }
    }
    // Dragging vectors along (if 'bMat' given)
    if (numb>0) {
      ArrayHandle<double> colMsk;
      BaseLinVec<double> collin;
      BaseLinMat<double> siglin;
      // Allocate buffer (if not given via 'buffMat')
      ArrayHandle<double> sigBuff;
      int maxm=0;
      for (i=0; i<numb; i++)
	if (bMat[i]->rows()>maxm)
	  maxm=bMat[i]->rows();
      if (buffMat!=0)
	buffMat->zeros(maxm,d);
      else
	sigBuff.changeRep(maxm*d);
      for (i=0; i<numb; i++) {
	StMatrix* bAct=bMat[i];
	bAct->setStrctPatt(MatStrct::normal);
	{
	  BaseLinMat<double> blin;
	  bAct->getLinMat(blin,true); // with write-back
	  // Backsubstitution with \tilde{L}
	  int bm=bAct->rows();
	  // Prepare sigma matrix
	  if (buffMat!=0) {
	    buffMat->zeros(bm,d);
	    buffMat->getLinMat(siglin,false);
	  } else {
	    ArrayUtils<double>::fill(sigBuff.p(),0.0,bm*d);
	    siglin.init(bm,d,sigBuff,bm);
	  }
	  double* colB=blin.buff;
	  for (j=1; j<n; j++) {
	    colMsk.changeRep(colB,bm,false);
	    collin.init(bm,colMsk);
	    FastUtils::rankOne(siglin,collin,betaArr[j-1],1.0);
	    colB+=blin.stride;
	    colMsk.changeRep(colB,bm,false);
	    collin.init(bm,colMsk);
	    FastUtils::matVec(collin,siglin,false,pArr[j],-1.0,1.0);
	  }
	}
	// blin -> bAct written back here
	bAct->mulDiag(*bAct,lamDiag,true); // right-mult. with D^-1/2
      }
    }

    return 0; // OK
  }
  */

  /*
  int StMatrix::cholDecUpdateOld(const StVector& vec,bool add,bool isp,
				 ArrayHandle<StMatrix*>& bMat)
  {
    int i,j,numb;
    double temp,s,t,tlam,b;
    double* rowB,*lamB,*betaB;
    const double* pB;

    if (n!=m || n==0) throw WrongDimensionException(EXCEPT_MSG(""));
    checkTS(TSARG);
    if (strpatt!=MatStrct::lower)
      throw WrongStatusException(EXCEPT_MSG("Need lower triangular matrix"));
    vec.checkTS(TSARG);
    if (vec.size()!=n)
      throw WrongDimensionException(EXCEPT_MSG(""));
    if (isMaskInd())
      throw MaskObjectException("StMatrix::cholDecUpdate not implemented for indexed masks!");
    numb=bMat.size();
    for (i=0; i<numb; i++) {
      bMat[i]->checkTS(TSARG);
      if (bMat[i]->cols()!=n)
	throw WrongDimensionException(EXCEPT_MSG("bMat"));
    }

    // Backsubstitution (if 'pVec' not given)
    Handle<StVector> pvHand;
    const StVector* pVec;
    if (isp)
      pVec=&vec;
    else {
      pvHand.changeRep(new StVector(vec)); // copy of 'vec'
      backsubst(*pvHand);
      pVec=pvHand;
    }
    ArrayHandle<double> plin(pVec->getFlatBuff()); // flat version

    // I +- p p^T = \tilde{L} \tilde{\Lambda} \tilde{L}^T. Here, \tilde{L}
    //   is represented by p and \beta.
    StVector lamDiag(n);
    ArrayHandle<double> betaArr(n-1),lamArr(lamDiag.getFlatBuff());
    if (add) {
      // I + p p^T: t recurrence
      t=1.0;
      pB=plin.p(); betaB=(double*) betaArr.p();
      lamB=(double*) lamArr.p();
      for (j=0; j<n-1; j++) {
	temp=*(pB++);
	tlam=t+temp*temp; // t_j
	*(lamB++)=tlam/t; // \tilde{\lambda}_j
	t=tlam;
	*(betaB++)=temp/tlam; // \beta_j
      }
      temp=*pB;
      *lamB=(t+temp*temp)/t;
    } else {
      // I - p p^T: s recurrence
      s=-1.0;
      pB=plin.p(); betaB=(double*) betaArr.p();
      lamB=(double*) lamArr.p();
      for (j=0; j<n-1; j++) {
	temp=*(pB++); // p_j
	b=temp*s; // p_j s_{j-1}
	tlam=(*(lamB++)=1.0+(temp*b)); // \tilde{\lambda}_j
	if (tlam<=0.0)
	  return 1; // Numerical breakdown
	temp=*(betaB++)=b/tlam; // \beta_j
	s-=(b*temp); // s_j
      }
      temp=*pB;
      *lamB=tlam=1.0+s*temp*temp;
      if (tlam<=0.0)
	return 1; // Numerical breakdown
    }
    // Update of L. Note that \tilde{l}_{i,j} = p_i \beta_j for i>j and
    //  \tilde{l}_{i,i}=1.
    lamDiag.apply1(ptr_fun(sqrt)); // sqrt
    (*buff)*=lamArr[0];
    for (i=1; i<n; i++) {
      pB=plin.p()+i;
      lamB=(double*) lamArr.p()+i;
      rowB=buff+(i*(stride+1)); // (i,i)
      b=(*rowB)*(*pB);
      (*rowB)*=(*lamB);
      betaB=(double*) betaArr.p()+i;
      for (j=i-1,s=0.0; j>=0; j--) {
	s+=b;
	rowB-=stride; // (i,j)
	temp=*rowB;
	b=temp*(*(--pB));
	*rowB=(temp+s*(*(--betaB)))*(*(--lamB));
      }
    }
    // Dragging vectors along (if 'bMat' given). See 'cholLRBacksubst'
    if (numb>0) {
      // Need buffer vector
      int maxm=0;
      for (i=0; i<numb; i++)
	if (bMat[i]->rows()>maxm) maxm=bMat[i]->rows();
      ArrayHandle<double> sumArr(maxm);
      for (i=0; i<numb; i++) {
	StMatrix* bAct=bMat[i];
	bAct->setStrctPatt(MatStrct::normal);
	{
	  BaseLinMat<double> blin;
	  bAct->getLinMat(blin,true);
	  // Backsubstitution with \tilde{L}
	  betaB=(double*) betaArr.p();
	  pB=plin.p();
	  int bm=bAct->rows();
	  double* sumB=(double*) sumArr.p(),*colB;
	  ArrayUtils<double>::fill(sumB,0.0,bm);
	  colB=blin.buff;
	  for (j=1; j<n; j++) {
	    FastUtils::addsmul(sumB,colB,*(betaB++),bm);
	    colB+=blin.stride; pB++;
	    FastUtils::addsmul(colB,sumB,-(*pB),bm);
	  }
	}
	// blin -> bAct written back here
	bAct->mulDiag(*bAct,lamDiag,true); // right-mult. with D^-1/2
      }
    }

    return 0; // OK
  }
  */

#define EIG_TOL 1e-5
  /*
   * We use the same update algorithm than 'cholDecLowRank'/
   * 'cholDecLowRankMinus'.
   * The safeguard is used by Gill and Murray in their BFGS implementation.
   * We have
   *   A' = L (D + s p p^T) L^T, L p = v,
   * which is pos. def. iff
   *   1 + s p^T D^{-1} p > 0.
   * Thus, if s < (eps-1)/(p^T D^{-1} p), we use
   *   s' = (eps-1)/(p^T D^{-1} p)
   * instead.
   */
  /*
  void StMatrix::cholDecUpdateSpecOld(StVector& diag,const StVector& vec,
				      double scal,bool safe,
				      const StVector* pVec)
  {
    int i,j;
    double temp,s,t,tlam,b,lam;
    double* rowB,*betaB;
    const double* pB;

    if (n!=m || n==0) throw WrongDimensionException(EXCEPT_MSG(""));
    if (diag.size()!=n || vec.size()!=n) throw WrongDimensionException(EXCEPT_MSG(""));
    checkTS(TSARG); diag.checkTS(TSARG); vec.checkTS(TSARG);
    if (pVec!=0) {
      pVec->checkTS(TSARG);
      if (pVec->size()!=n) throw WrongDimensionException(EXCEPT_MSG(""));
    }
    if (isMaskInd())
      throw MaskObjectException("StMatrix::cholDecUpdateSpec not implemented for indexed masks!");
    if (scal==0.0) return;
    uchar spatt=getStrctPatt();
    setStrctPatt(MatStrct::lowNDg);

    // Backsubstitution
    Handle<StVector> pvHand;
    if (pVec==0) {
      pvHand.changeRep(new StVector());
      *pvHand=vec;
      backsubst(*pvHand); // L p = v
      pVec=pvHand;
    }
    ArrayHandle<double> plin(pVec->getFlatBuff());

    // Safeguard
    StVector betaVec(n);
    if (safe && scal<0.0) {
      betaVec.div(*pVec,diag);
      temp=pVec->inner(betaVec);
      if (1.0+scal*temp<EIG_TOL) {
	cout << "cholDecUpdateSpec: Need to raise scal from " << scal; // DEBUG!!
	scal=(EIG_TOL-1.0)/temp;
	cout << " to " << scal << endl; // DEBUG!!
      }
    }

    // D + s p p^T = \tilde{L} \tilde{\Lambda} \tilde{L}^T. Here, \tilde{L}
    // is represented by p and \beta.
    ArrayHandle<double> betaArr=betaVec.getFlatBuff();
    t=1.0/scal;
    pB=plin.p(); betaB=(double*) betaArr.p();
    for (j=0; j<n-1; j++) {
      temp=*(pB++); // p_j
      lam=diag[j]; // \lambda_j
      tlam=t*lam+temp*temp; // t_j \lambda_j
      diag[j]=tlam/t; // \tilde{\lambda}_j
      t=tlam/lam;
      *(betaB++)=temp/tlam; // \beta_j
    }
    temp=*pB;
    tlam=t*diag[n-1]+temp*temp;
    diag[n-1]=tlam/t;

    // Update of L. Note that \tilde{l}_{i,j} = p_i \beta_j for i>j and
    // \tilde{l}_{i,i}=1.
    for (i=1; i<n; i++) {
      pB=plin.p()+i;
      b=*pB;
      rowB=buff+(i*(stride+1)); // pos. (i,i)
      betaB=(double*) betaArr.p()+i;
      for (j=i-1,s=0.0; j>=0; j--) {
	s+=b;
	rowB-=stride; // (i,j)
	temp=*rowB;
	b=temp*(*(--pB)); // l_{i,j} p_j
	*rowB=temp+s*(*(--betaB));
      }
    }
    setStrctPatt(spatt); // restore
  }
#undef EIG_TOL
  */

#define EIG_TOL 1e-5
  /*
   * Same update algorithm as in 'cholDecUpdate' (rank 1 version).
   * The safeguard is used by Gill and Murray in their BFGS implementation.
   * We have
   *   A' = L (D + s p p^T) L^T, L p = v,
   * which is pos. def. iff
   *   1 + s p^T D^{-1} p > 0.
   * If 'safe'==true and s < (eps-1)/(p^T D^{-1} p), we abort the method
   * without changing any variables (except 'pVec'). Here, eps=='EIG_TOL'.
   */
  bool StMatrix::cholDecUpdateSpec(StVector& diag,const StVector& vec,
				   double scal,bool isp,bool safe,
				   StVector* pVec,StVector* betaVec,
				   StMatrix* mMat)
  {
    int j;
    double temp,t,tlam,lam;
    double* betaB;
    const double* pB;
    bool updateL=(betaVec==0),doDrag=(mMat!=0);

    checkTS(TSARG); diag.checkTS(TSARG); vec.checkTS(TSARG);
    if (n!=m || n==0) throw WrongDimensionException(EXCEPT_MSG(""));
    if (diag.size()!=n || vec.size()!=n) throw WrongDimensionException(EXCEPT_MSG(""));
    if (isMaskInd())
      throw MaskObjectException("StMatrix::cholDecUpdateSpec not implemented for indexed masks!");
    if (!updateL && !isp && pVec==0)
      throw InvalidParameterException(EXCEPT_MSG("Need 'pVec' for return"));
    if (doDrag) {
      mMat->checkTS(TSARG);
      if (mMat->cols()!=n) throw WrongDimensionException(EXCEPT_MSG(""));
    }
    if (scal==0.0) throw InvalidParameterException(EXCEPT_MSG(""));
    uchar spatt=getStrctPatt();
    setStrctPatt(MatStrct::lowNDg); // L is unit lower triangular

    // Compute p, setup beta
    Handle<StVector> pHand,betaHand;
    if (isp) {
      if (pVec!=0) *pVec=vec;
      else pVec=(StVector*) &vec;
    } else {
      if (pVec!=0) *pVec=vec;
      else {
	pHand.changeRep(new StVector(vec));
	pVec=pHand;
      }
      backsubst(*pVec); // p = L^-1 v
    }
    if (betaVec==0) {
      betaHand.changeRep(new StVector());
      betaVec=betaHand;
    }
    betaVec->zeros(n);
    ArrayHandle<double> pArr(pVec->getFlatBuff());
    ArrayHandle<double> betaArr(betaVec->getFlatBuff());

    // Safeguard
    if (safe && scal<0.0) {
      betaVec->div(*pVec,diag);
      temp=pVec->inner(*betaVec);
      if (1.0+scal*temp<EIG_TOL) {
	// Abort method
	setStrctPatt(spatt);
	return true;
      }
    }

    // D + s p p^T = \tilde{L} \tilde{\Lambda} \tilde{L}^T. Here, \tilde{L}
    // is represented by p and \beta.
    t=1.0/scal;
    pB=pArr.p(); betaB=(double*) betaArr.p();
    for (j=0; j<n-1; j++) {
      temp=*(pB++); // p_j
      lam=diag[j]; // \lambda_j
      tlam=t*lam+temp*temp; // t_j \lambda_j
      diag[j]=tlam/t; // \tilde{\lambda}_j
      t=tlam/lam;
      *(betaB++)=temp/tlam; // \beta_j
    }
    temp=*pB;
    tlam=t*diag[n-1]+temp*temp;
    diag[n-1]=tlam/t; // D' complete

    // Update of L
    if (updateL)
      cholDecUpdateSpecComplete(pArr,betaArr);

    // Dragging vectors along
    if (doDrag)
      cholDecUpdateSpecDrag(*mMat,pArr,betaArr);

    setStrctPatt(spatt); // restore
    return false;
  }
#undef EIG_TOL

  void StMatrix::cholDecUpdateSpecComplete(const ArrayHandle<double>& pArr,
					   const ArrayHandle<double>& betaArr)
  {
    int i,j;
    double temp,b,s;
    double* rowB,*betaB;
    const double* pB;

    checkTS(TSARG);
    if (n!=m || n==0) throw WrongDimensionException(EXCEPT_MSG(""));
    if (isMaskInd())
      throw MaskObjectException("StMatrix::cholDecUpdateSpecComplete not implemented for indexed mask");
    if (pArr.size()<n || betaArr.size()<n-1)
      throw WrongDimensionException(EXCEPT_MSG(""));

    // Update of L. Note that \tilde{l}_{i,j} = p_i \beta_j for i>j and
    // \tilde{l}_{i,i}=1.
    for (i=1; i<n; i++) {
      pB=pArr.p()+i;
      b=*pB;
      rowB=buff+(i*(stride+1)); // pos. (i,i)
      betaB=(double*) betaArr.p()+i;
      for (j=i-1,s=0.0; j>=0; j--) {
	s+=b;
	rowB-=stride; // (i,j)
	temp=*rowB;
	b=temp*(*(--pB)); // l_{i,j} p_j
	*rowB=temp+s*(*(--betaB));
      }
    }
  }

  void StMatrix::cholDecUpdateSpecDrag(StMatrix& mMat,
				       const ArrayHandle<double>& pArr,
				       const ArrayHandle<double>& betaArr,
				       bool trans)
  {
    int j;

    mMat.checkTS(TSARG);
    int n=mMat.cols(),m=mMat.rows();
    if (pArr.size()<n || betaArr.size()<n-1)
      throw WrongDimensionException(EXCEPT_MSG(""));

    // Need buffer vector
    ArrayHandle<double> sumArr(m);
    mMat.setStrctPatt(MatStrct::normal);
    {
      BaseLinMat<double> mlin;
      mMat.getLinMat(mlin,true); // with write-back
      if (!trans) {
	// Backsubstitution with L''
	const double* betaB=(double*) betaArr.p(),*pB=pArr.p();
	double* sumB=(double*) sumArr.p(),*colB;
	ArrayUtils<double>::fill(sumB,0.0,m);
	colB=mlin.buff;
	for (j=1; j<n; j++) {
	  FastUtils::addsmul(sumB,colB,*(betaB++),m);
	  colB+=mlin.stride; pB++;
	  FastUtils::addsmul(colB,sumB,-(*pB),m);
	}
      } else {
	// Backsubstitution with (L'')^T
	const double* betaB=(double*) betaArr.p()+(n-1),*pB=pArr.p()+(n-1);
	double* sumB=(double*) sumArr.p(),*colB;
	ArrayUtils<double>::fill(sumB,0.0,m);
	colB=mlin.buff+((n-1)*mlin.stride);
	for (j=n-2; j>=0; j--) {
	  FastUtils::addsmul(sumB,colB,*(pB--),m);
	  colB-=mlin.stride; betaB--;
	  FastUtils::addsmul(colB,sumB,-(*betaB),m);
	}
      }
    }
    // mlin -> mMat written back here (at end of scope)
  }
//ENDNS
