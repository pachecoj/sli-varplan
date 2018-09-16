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
 * Project source file
 * Module: eplin
 * Desc.:  Definition of class ExpectPropLinear
 * ------------------------------------------------------------------- */

#include "src/eplin/ExpectPropLinear.h"
#include "src/eplin/EPSingleLaplaceManager.h"
#include "lhotse/rando/Random.h"
#include "lhotse/specfun/Specfun.h"
#include "lhotse/PrintToStdout.h"
//#include "lhotse/DebugMacros.h" // DEBUG
#include "lhotse/quad/LaplaceSingleLikehood.h" // DEBUG!!
//#include "lhotse/FileUtils.h" // DEBUG!
//#include "lhotse/NumberFormats.h" // DEBUG!

//BEGINNS(eplin)

#define CHECKCONSISTENT do { if (prior->cols()!=sites->size()) \
  throw InvalidParameterException(EXCEPT_MSG("")); } while (0)

#define UPDATE_REPRES do { if (!up2date) refresh(); } while (0)

#define PI_LOWER_BOUND 1e-8

#define DISCRELDIFF(a,b) (fabs((a)-(b))/std::max(fabs(a),std::max(fabs(b),1e-3)))

/* DEBUG:
#define DEBUGCOMP(a,txt) do { \
  NumberFormats<double>::load(*deb_is,&debScal); \
  if (fabs(debScal-(a))>1e-10) \
    cout << "DIFF[" << (txt) << "]: val=" << (a) << ",file=" << debScal \
         << endl; \
  } while(0)
#define DEBUGCOMPVEC(a,txt) do { \
  debVec.load(*deb_is); \
  debVec.addprod(-1.0,(a)); debVec.apply1(ptr_fun(fabs)); \
  if ((debScal=debVec.max())>5e-12) \
    cout << "DIFF[" << (txt) << "]: maxabsdiff=" << debScal << endl; \
  } while(0)
#define DEBUGCOMPMAT(a,txt) do { \
  debMat.load(*deb_is); \
  debMat.addsmul((a),-1.0); debMat.apply1(ptr_fun(fabs)); \
  if ((debScal=debMat.max().max())>5e-12) \
    cout << "DIFF[" << (txt) << "]: maxabsdiff=" << debScal << endl; \
  } while(0)
*/

  ExpectPropLinear::ExpectPropLinear(double peps,const Handle<EPPrior>& pprior,
				     const Handle<EPSingleSiteManager>& psites,
				     int nonDegM,int istat,
				     const Handle<StMatrix>& buffmat1,
				     const Handle<StMatrix>& buffmat2,
				     const Handle<StMatrix>& buffmat3,
				     const Handle<StVector>& buffvec) :
    eps(peps),prior(pprior),sites(psites),tempmat(buffmat1),up2date(false),
    tempmat2(buffmat2),tempmat3(buffmat3),tempvec(buffvec),verbose(0)
    //deb_cnt(0) // DEBUG
    //deb_is(0) // DEBUG!
  {
    int n=pprior->cols(),m=pprior->rows();

    CHECKCONSISTENT;
    if (nonDegM==-1) nonDegM=n/2;
    else if (nonDegM<0 || nonDegM>n)
      throw InvalidParameterException(EXCEPT_MSG("nonDegM"));
    nonDegenM=nonDegM;
    if (peps<0.0 || (m<n && peps<PI_LOWER_BOUND))
      throw InvalidParameterException(EXCEPT_MSG(""));
    if (buffmat1==0) tempmat.changeRep(new StMatrix());
    if (buffmat2==0) tempmat2.changeRep(new StMatrix());
    if (buffmat3==0) tempmat3.changeRep(new StMatrix());
    if (buffvec==0) tempvec.changeRep(new StVector());
    if (istat==0) {
      sites->setSiteB(0.0); sites->setSitePi(peps);
    }
  }

  void ExpectPropLinear::refresh()
  {
    CHECKCONSISTENT;
    if (isDegenerate()) {
      // Degenerate repres.
      tempvec->div(1.0,sites->getSitePi()); // Pi^-1
      prior->outerProd(*tempvec,factL);
      factL.addseye(1.0); // I + X Pi^-1 X^T
      factL.setStrctPatt(MatStrct::lower);
      if (factL.cholDecomp(true)!=0) { // L
	throw NumericalException(EXCEPT_MSG(""));
      }
      tempvec->addprod(sites->getSiteB(),1.0,prior->getBVec());
      tempvec->div(sites->getSitePi());
      prior->mult(*tempvec,gamma);
      factL.backsubst(gamma); // gamma
    } else {
      // Non-degenerate repres.
      prior->outerTProd(factL); // X^T X
      factL.diag()->addprod(1.0,sites->getSitePi()); // A^-1
      factL.setStrctPatt(MatStrct::lower);
      if (factL.cholDecomp(true)!=0) // L
	throw NumericalException(EXCEPT_MSG(""));
      gamma.addprod(sites->getSiteB(),1.0,prior->getBVec());
      factL.backsubst(gamma); // gamma
    }
    up2date=true;
  }

  bool ExpectPropLinear::updateSite(int i,double thres,double frac,
				    double* delta,int utype)
  {
    double ch,crho,beta,nu,pi,b,pinew,bnew,delpi,temp,b0,temp2,temp3;
    double h,rho,hpr,rhopr;
    StMatrix matmsk;
    StVector vmsk;
    bool doFrac=(frac!=1.0);
    int n=prior->cols(),m=prior->rows();
    double sigsq=prior->getSigSq();
    bool modifPi=false;
    double tau;

    CHECKCONSISTENT;
    if (i<0 || i>=n) throw OutOfRangeException(EXCEPT_MSG(""));
    if (thres<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
    if (frac<=0.0 || frac>1.0)
      throw InvalidParameterException(EXCEPT_MSG("frac"));
    if (doFrac && !sites->suppFractional(i))
      throw InvalidParameterException(EXCEPT_MSG("LH factor does not support fractional updates"));
    if (utype<0 || utype>3)
      throw InvalidParameterException(EXCEPT_MSG(""));
    if (doFrac && utype>1)
      throw NotImplemException(EXCEPT_MSG(""));
    if (utype>1) {
      const EPSingleLaplaceManager* lsites=
	DYNCAST(const EPSingleLaplaceManager,sites.p());
      if (lsites==0)
	throw WrongStatusException(EXCEPT_MSG("Need Laplace sites!"));
      tau=lsites->getTau();
    }
    UPDATE_REPRES;
    // DEBUG:
    //DCOMP_STCMSCAL((double) i,"i"); DCOMP_STCMSCAL(thres,"thres");
    //DCOMP_STCMVEC(gamma,"gamma(INIT)");
    //DCOMP_STCMVEC((StVector&) factL.diag(),"dgL(INIT)");
    // END DEBUG

    // Compute marginal (cavity) moments -> 'ch', 'crho'
    pi=sites->getSitePi()[i]; b=sites->getSiteB()[i];
    //DCOMP_STCMSCAL(pi,"pi"); DCOMP_STCMSCAL(b,"b"); // DEBUG
    if (utype!=0 && b!=0.0)
      throw InternalException(EXCEPT_MSG("")); // Sanity check
    if (verbose>1)
      mycout << "i=" << i << ": pi=" << pi << ",b=" << b << "\n";
    if (isDegenerate()) {
      // Degenerate repres.
      // NOTE: In MEX mode, '*tempvec' is a mask vector of size n
#ifndef MATLAB_MEX
      tempvec->zeros(m);
#endif
      vmsk.reassign(*tempvec,Range(0,m-1));
      prior->getCol(i,vmsk);
      factL.backsubst(vmsk); // v = L^-1 x_i
      //DCOMP_STCMVEC(vmsk,"tvec(CAV-DEG)"); // DEBUG
      temp=vmsk.inner(vmsk); // |v|^2
      b0=(prior->getBVec())[i];
      temp2=vmsk.inner(gamma);
      if (!doFrac) {
	crho=1.0/temp-1.0/pi;
	ch=(b0-temp2)/temp+b/pi;
      } else {
	temp3=1.0-temp/pi;
	crho=temp3/pi;
	crho/=(temp3=1.0-frac*temp3);
	ch=((b0-temp2)/temp3+b)/pi;
      }
      h=(b0+b-temp2)/pi; rho=(1.0-temp/pi)/pi;
      if (verbose>1) {
	mycout << "marg: mean=" << h << ",rho=" << rho << "\n";
	if (utype<2)
	  mycout << "cav:  mean=" << ch << ",rho=" << crho << "\n";
      }
      if (utype<2 && crho<=0.0) {
	mycout << "i=" << i << ": Cavity variance not positive!\n";
	return false; // Break
      }
    } else {
      // Non-degen. repres.
      tempvec->zeros(n); (*tempvec)[i]=1.0;
      factL.backsubst(*tempvec); // v = L^-1 delta_i
      //DCOMP_STCMVEC(*tempvec,"tvec(CAV)"); // DEBUG
      rho=tempvec->inner(*tempvec); // Marg. moments
      h=tempvec->inner(gamma);
      if (verbose>1)
	mycout << "marg: mean=" << h << ",rho=" << rho << "\n";
      if (utype<2) {
	if ((temp=1.0-frac*pi*rho)<=0.0) {
	  mycout << "i=" << i << ": Cavity variance not positive!\n";
	  mycout << "temp=" << temp << ",frac=" << frac << ",pi=" << pi
		 << ",rho=" << rho << "\n";
	  return false; // Break
	}
	ch=(h-frac*rho*b)/temp; crho=rho/temp;
	if (verbose>1)
	  mycout << "cav:  mean=" << ch << ",rho=" << crho << "\n";
      }
    }
    // DEBUG:
    //DCOMP_STCMSCAL(h,"h"); DCOMP_STCMSCAL(rho,"rho");
    //if (utype<2) { DCOMP_STCMSCAL(ch,"ch"); DCOMP_STCMSCAL(crho,"crho"); }

    // Update: pinew, bnew, hpr, rhopr
    bnew=0.0; // for all 'utype' except 0
    switch (utype) {
    case 0:
      // Standard ADF update
      // ADF projection -> beta, nu. Don't need log Z itself
      sites->compMoments(i,sigsq,ch,crho,&beta,&nu,frac);
      //DCOMP_STCMSCAL(beta,"beta"); DCOMP_STCMSCAL(nu,"nu"); // DEBUG
      hpr=ch+beta*crho*sigsq; rhopr=crho*(1.0-crho*sigsq*nu);
      if (verbose>1)
	mycout << "beta=" << beta << ",nu=" << nu << "\n";
      temp=1.0-nu*crho*sigsq;
      pinew=(1.0-frac)*pi+nu*sigsq/temp;
      bnew=(1.0-frac)*b+sigsq*(ch*nu+beta)/temp;
      break;
    case 1:
      // ADF update, but b_i=0
      sites->compMoments(i,sigsq,ch,crho,&beta,&nu,frac);
      //DCOMP_STCMSCAL(beta,"beta"); DCOMP_STCMSCAL(nu,"nu"); // DEBUG
      // Need to additional projection, because b_i cannot be used
      hpr=ch+beta*crho*sigsq; rhopr=crho*(1.0-crho*sigsq*nu);
      temp=rhopr+hpr*hpr/sigsq; // beta (in paper)
      temp=0.5*(crho+sqrt(crho*crho+4.0*ch*ch*temp/sigsq))/temp; // kappa
      pinew=(temp-1.0)/rho/frac;
      hpr=ch/temp; rhopr=crho/temp;
      break;
    case 2:
    case 3:
      // Tipping SBL update / Girolami variational
      beta=rho+h*h/sigsq; // E[a_i^2]/sigma^2
      pinew=(utype==2)?(0.5*(sqrt(9.0+4.0*beta*tau*tau)-3.0)/beta):
	(tau/sqrt(beta));
      break;
    }
    // DEBUG:
    //DCOMP_STCMSCAL(pinew,"pinew"); DCOMP_STCMSCAL(bnew,"bnew");
    //if (utype<2) { DCOMP_STCMSCAL(hpr,"hpr"); DCOMP_STCMSCAL(rhopr,"rhopr"); }
    // END DEBUG
    if (verbose>1) {
      mycout << "pinew=" << pinew << ",bnew=" << bnew << "\n";
      if (utype<2)
	mycout << "tilt: mean=" << hpr << ",rho=" << rhopr << "\n";
    }

    // Update of representation
    if (isDegenerate()) {
      // Degenerate repres.
      if (pinew<PI_LOWER_BOUND) {
	if (verbose>0)
	  mycout << "i=" << i << ": pinew fix to " << PI_LOWER_BOUND << "\n";
	pinew=PI_LOWER_BOUND; modifPi=true;
      }
      if (pinew<pi) {
	delpi=1.0-pinew/pi; temp=pinew;
      } else {
	delpi=pi/pinew-1.0; temp=pi;
      }
      if (fabs(delpi)<thres*temp) {
	// Do not do the update (would be too small)
	delpi=0.0; pinew=pi; bnew=b;
      }
      //DCOMP_STCMSCAL(delpi,"delpi(DEG)"); // DEBUG
      if (delpi!=0.0) {
	matmsk.reassign(gamma,false); // Mask gamma as row vector
	if (delpi>0.0) {
	  // pinew < pi
	  prior->getCol(i,vmsk); // x_i
	  temp=sqrt(pinew); // sqrt(min(pinew,pi))
	  // sqrt(1/pinew-1/pi) = temp2/temp:
	  vmsk.prod((temp2=sqrt(delpi))/temp);
	  temp3=(bnew==b)?
	    (b0+b)*temp2/temp:
	    (b0*temp2+(bnew-b*pinew/pi)/temp2)/temp;
	  //DCOMP_STCMVEC(vmsk,"tvec(POS-DEG)"); // DEBUG
	  //DCOMP_STCMSCAL(temp3,"temp3(POS-DEG)"); // DEBUG
#ifndef MATLAB_MEX
	  factL.cholUpdRk1(vmsk,tempvec2,tempvec3,vmsk,&matmsk,
			   &((StVector&) StVector::mask(&temp3,1)));
#else
	  factL.cholUpdRk1(vmsk,tempvec2(Range(0,m-1)),tempvec3(Range(0,m-1)),
			   vmsk,&matmsk,
			   &((StVector&) StVector::mask(&temp3,1)));
#endif
	} else {
	  // pi < pinew
	  // 'vmsk' contains L^-1 x_i
	  temp=sqrt(pi); // sqrt(min(pinew,pi))
	  // sqrt(-1/pinew+1/pi) = temp2/temp:
	  vmsk.prod((temp2=sqrt(-delpi))/temp);
	  temp3=(bnew==b)?
	    (b0+b)*temp2/temp:
	    (b0*temp2+(b-bnew*pi/pinew)/temp2)/temp;
	  //DCOMP_STCMVEC(vmsk,"tvec(NEG-DEG)"); // DEBUG
	  //DCOMP_STCMSCAL(temp3,"temp3(NEG-DEG)"); // DEBUG
#ifndef MATLAB_MEX
	  if (factL.cholDndRk1(vmsk,tempvec2,tempvec3,tempind,vmsk,true,
			       &matmsk,
			       &((StVector&) StVector::mask(&temp3,1)))!=0)
	    throw NumericalException(EXCEPT_MSG(""));
#else
	  if (factL.cholDndRk1(vmsk,tempvec2(Range(0,m-1)),
			       tempvec3(Range(0,m-1)),tempind,vmsk,true,
			       &matmsk,
			       &((StVector&) StVector::mask(&temp3,1)))!=0)
	    throw NumericalException(EXCEPT_MSG(""));
#endif
	}
      }
    } else {
      // Non-degen. repres.
      temp=(n<m)?0.0:PI_LOWER_BOUND; // Could still be degen. case!
      if (pinew<temp) {
	if (verbose>0)
	  mycout << "i=" << i << ": pinew fix to " << temp << "\n";
	pinew=temp; modifPi=true;
      }
      delpi=pinew-pi;
      if (fabs(delpi)<thres) {
	// Do not do the update (would be too small)
	delpi=0.0; pinew=pi; bnew=b;
      }
      //DCOMP_STCMSCAL(delpi,"delpi"); // DEBUG
      if (delpi!=0.0) {
	// DEBUG: COMPARE AGAINST TRUTH
	//deb_factL=factL; deb_factL.setStrctPatt(MatStrct::normal);
	//deb_mat.setStrctPatt(MatStrct::lower);
	//deb_mat.symMul(deb_factL,false); // L L^T
	// END DEBUG
	matmsk.reassign(gamma,false);
	if (delpi>0.0) {
	  tempvec->zeros(n); (*tempvec)[i]=temp=sqrt(delpi);
	  temp=(bnew-b)/temp;
	  //DCOMP_STCMVEC(*tempvec,"tvec(POS)"); // DEBUG
	  //DCOMP_STCMSCAL(temp,"temp(POS)"); // DEBUG
	  factL.cholUpdRk1(*tempvec,tempvec2,tempvec3,*tempvec,&matmsk,
			   &((StVector&) StVector::mask(&temp,1)));
	} else {
	  // '*tempvec' contains v = L^-1 delta_i
	  tempvec->prod(temp=sqrt(-delpi));
	  temp=(b-bnew)/temp;
	  //DCOMP_STCMVEC(*tempvec,"tvec(NEG)"); // DEBUG
	  //DCOMP_STCMSCAL(temp,"temp(NEG)"); // DEBUG
	  if (factL.cholDndRk1(*tempvec,tempvec2,tempvec3,tempind,*tempvec,
			       true,&matmsk,
			       &((StVector&) StVector::mask(&temp,1)))!=0)
	    throw NumericalException(EXCEPT_MSG(""));
	}
	/* DEBUG: COMPARE AGAINST EXACT
	deb_vvec.zeros(n); deb_vvec[i]=1.0;
	deb_mat.rankOne(deb_vvec,delpi);
	if (deb_mat.cholDecomp(true)!=0)
	throw NumericalException(EXCEPT_MSG(""));
	deb_gamma.addprod(sites->getSiteB(),1.0,prior->getBVec());
	deb_gamma[i]+=(bnew-b);
	deb_mat.backsubst(deb_gamma);
	//if (delpi>0.0) cout << "POS" << endl;
	//else cout << "NEG" << endl;
	double deb_mx;
	if ((deb_mx=temp=deb_mat.maxRelDiff(factL))>1e-7)
	cout << "factL: " << temp << endl;
	if ((temp=deb_gamma.maxRelDiff(gamma))>1e-7)
	cout << "gamma: " << temp << endl;
	if ((deb_mx=std::max(deb_mx,temp))>1e-5 && ++deb_cnt<10) {
	// Store variables
	cout << "STORE DEBUG FILES" << endl;
	ofstream os;
	char sbuff[80];
	sprintf(sbuff,"debug%d-factl.stm",deb_cnt);
	FileUtils::openFileWrite(sbuff,os);
	deb_factL.save(os); CLOSESTR(os);
	sprintf(sbuff,"debug%d-gamma.stv",deb_cnt);
	FileUtils::openFileWrite(sbuff,os);
	deb_gamma.insert(delpi); deb_gamma.insert(bnew-b);
	deb_gamma.save(os); CLOSESTR(os);
	}
	// END DEBUG
	*/
      }
    }
    if (delpi!=0.0) {
      // DEBUG
      //DCOMP_STCMVEC(gamma,"gamma(UPD)");
      //DCOMP_STCMVEC((StVector&) factL.diag(),"dgL(UPD)");
      // END DEBUG
      if (delta!=0) {
	if (modifPi || utype>1) {
	  temp=1.0+rho*(pinew-pi);
	  hpr=(h+rho*(bnew-b))/temp; rhopr=rho/temp;
	}
	temp=sqrt(sigsq*rhopr); temp2=sqrt(sigsq*rho);
	*delta=std::max(DISCRELDIFF(hpr,h),DISCRELDIFF(temp,temp2));
	//DCOMP_STCMSCAL(*delta,"delta");
      }
      sites->setSitePi(i,pinew);
      if (utype==0) sites->setSiteB(i,bnew);
    } else if (delta!=0) *delta=0.0;

    return (delpi!=0.0);
  }

  /*
   * ATTENTION: Must not use 'tempvecXXX'. Otherwise, call from
   * 'compCritGrad' goes wrong!
   */
  void ExpectPropLinear::compLogZ(const StVector& means,const StVector& vars,
				  StVector& logz,StVector* exabs)
  {
    int i,n=prior->cols();
    double h,rho,mrho,temp,sigsq=prior->getSigSq(),beta,nu;
    StVector vecmsk;

    if (means.size()!=n || vars.size()!=n)
      throw WrongDimensionException(EXCEPT_MSG(""));
    if (exabs!=0 &&
	DYNCAST(const EPSingleLaplaceManager,sites.p())==0)
      throw InvalidParameterException(EXCEPT_MSG("Need Laplace sites for 'exabs'"));
    UPDATE_REPRES;
    logz.zeros(n);
    if (exabs!=0) exabs->zeros(n);
    // DEBUG:
    //double tau=DYNCAST(const EPSingleLaplaceManager,sites.p())->getTau(),
    //  temp2;
    //cout << "tau=" << tau << ",sigsq=" << sigsq << endl;
    // END DEBUG
    for (i=0; i<n; i++) {
      mrho=vars[i]/sigsq;
      if ((temp=1.0-mrho*sites->getSitePi()[i])<=0.0)
	throw NumericalException(EXCEPT_MSG(""));
      rho=mrho/temp;
      h=(means[i]-mrho*sites->getSiteB()[i])/temp;
      if (exabs==0)
	logz[i]=sites->compMoments(i,sigsq,h,rho);
      else {
	vecmsk.reassign(*exabs,Range(i,i));
	logz[i]=sites->compMoments(i,sigsq,h,rho,&beta,&nu,1.0,&vecmsk);
      }
      // DEBUG:
      //temp=tau/sqrt(sigsq);
      //if (exabs==0)
      //logz[i]=LaplaceSingleLikehood::compLogPartStatic(h*temp,rho*tau*tau,
      //						 1.0)+log(0.5*temp);
      //else {
      //logz[i]=LaplaceSingleLikehood::compLogPartStatic(h*temp,rho*tau*tau,
      //						 1.0,&beta,&nu,1.0,
      //						 &temp2)+log(0.5*temp);
      //(*exabs)[i]=temp2/temp;
      //}
      // END DEBUG
    }
  }

  double
  ExpectPropLinear::compCritGrad(StMatrix& gradX,double& gradSig,
				 double& gradTau,int meth,StVector* wkvecm,
				 StVector* wkvecn,int gradStat)
  {
    int i,j,m=prior->rows(),n=prior->cols();
    double temp,crit,h,rho,pi,b,sigsq=prior->getSigSq(),tau;
    StVector* means,*vars,*logz,*exabs=0,*uvec;
    Handle<StVector> hwkvecn,hwkvecm;

    // Initialization
    if (meth<0 || meth>3) throw InvalidParameterException(EXCEPT_MSG("meth"));
    const EPSingleLaplaceManager* lsites=
      DYNCAST(const EPSingleLaplaceManager,sites.p());
    if (lsites==0)
      throw WrongStatusException(EXCEPT_MSG("Need Laplace sites!"));
    tau=lsites->getTau();
    tempvec->zeros(n); logz=tempvec.p();
    tempvec2.zeros(n); means=&tempvec2;
    tempvec3.zeros(n); vars=&tempvec3;
    if (wkvecm==0) {
      hwkvecm.changeRep(new StVector(m)); wkvecm=hwkvecm.p();
    } else
      wkvecm->zeros(m);
    uvec=wkvecm;
    if (gradStat>0) {
      if (meth<2) {
	if (wkvecn==0) {
	  hwkvecn.changeRep(new StVector(n)); wkvecn=hwkvecn.p();
	} else
	  wkvecn->zeros(n);
	exabs=wkvecn;
      }
    }
    crit=0.0;
    if (gradStat>0) {
      gradSig=0.0; gradTau=0.0;
      if (gradStat>1) gradX.zeros(m,n); 
    }

    // Accumulate criterion and gradient
    if (isDegenerate() && meth>1)
      throw NotImplemException(EXCEPT_MSG("NOT IMPLEMENTED FOR DEGEN."));
    // NOTE: 'means', 'vars' ref. to 'tempvec2', 'tempvec3', which are not
    // used in 'getMargMoments':
    // DEBUG:
    //if (!doDebug)
    getMargMoments(*means,*vars); // h_i, sigma^2 rho_i
    //else
    //epDebug->getMargMoments(*means,*vars);
#ifdef MATLAB_DEBUG_OLD
    DOCOMP(*means,"%means=%factL'\\%gamma;","means");
    DOCOMP(*vars,"%vars=diag(%factL'\\(%factL\\eye(%n)))*%sigsq;","vars");
#endif // MATLAB_D
    // NOTE: 'compLogZ' does not use any 'tempvecXXX':
    if (gradStat==0) {
      if (meth<2) compLogZ(*means,*vars,*logz); // log Z_i values
    } else {
      // Compute gradient
      if (meth<2) {
	compLogZ(*means,*vars,*logz,exabs); // log Z_i values
	gradTau-=(exabs->sum()+((double) n)*sqrt(sigsq)/tau);
      } else {
	// SBL / Girolami
	logz->div(1.0,sites->getSitePi());
	gradTau+=(tau*logz->sum()-(temp=((double) n)/tau));
	if (meth==2) gradTau-=temp;
      }
      // DEBUG:
      //if (!doDebug) {
      prior->mult(*means,*uvec);
      uvec->prod(-1.0);
      prior->getVVec(*uvec,true); // f = u - X h
      //} else {
      //*uvec=epDebug->uvec;
      //epDebug->xmat.mulVec(*uvec,*means,false,-1.0,1.0);
      //}
#ifdef MATLAB_DEBUG_OLD
      DOCOMP(*uvec,"%fvec=%uvec-%xmat*%means;","fvec");
#endif // MATLAB_D
      if (gradStat>1) {
	gradX.rankOne(*uvec,-1.0/sigsq,*means);
#ifdef MATLAB_DEBUG_OLD
	DOCOMP(gradX,"%gradX=%gradX-%fvec*%means'/%sigsq;","gradX");
#endif // MATLAB_D
	prior->getXMat(*tempmat); // X
	if (isDegenerate()) {
	  factL.backsubst(*tempmat,true,false);
	  factL.backsubst(*tempmat,true,true);
	  tempmat->mulDiag(*tempmat,sites->getSitePi(),true);
	} else {
	  // DEBUG:
	  //if (!doDebug) {
	  factL.backsubst(*tempmat,false,true);
	  factL.backsubst(*tempmat,false,false);
	  //} else {
	  //epDebug->factL.backsubst(*tempmat,false,true);
	  //epDebug->factL.backsubst(*tempmat,false,false);
	  //}
	}
	gradX.addsmul(*tempmat,1.0); // += X Sigma
#ifdef MATLAB_DEBUG_OLD
	DOCOMP(gradX,"%gradX=%gradX+(%xmat/%factL')/%factL;","gradX");
#endif // MATLAB_D
      }
      if (meth<2) {
	// DEBUG:
	//if (!doDebug)
	gradSig+=0.5*(uvec->inner(*uvec)+sigsq*((double) (n-m))
		      -vars->inner(sites->getSitePi()));
	//else
	//  gradSig+=0.5*(uvec->inner(*uvec)+epDebug->sigsq*((double) (n-m))
	//		-vars->inner(epDebug->sitePi));
#ifdef MATLAB_DEBUG_OLD
	DOCOMP(gradSig,"%gradSig=%gradSig+0.5*(%fvec'*%fvec+%sigsq*(%n-%m)-%vars'*%sitePi);","gradSig");
#endif // MATLAB_D
      } else {
	// SBL / Girolami: only for non-degenerate!
	prior->getVVec(*uvec); // u
	gradSig+=
	  0.5*(uvec->inner(*uvec)-gamma.inner(gamma)-sigsq*((double) m));
      }
    }
    // Compute criterion
    prior->getVVec(*uvec); // u
    if (meth<2) {
      // DEBUG:
      //if (!doDebug)
      crit+=(-logz->sum()+0.5*uvec->inner(*uvec)/sigsq-
	     0.5*((double) (n-m))*(M_LN2+M_LNPI+log(sigsq)));
      //else
      //crit+=(-logz->sum()+0.5*epDebug->uvec.inner(epDebug->uvec)/sigsq-
      //       0.5*((double) (n-m))*(M_LN2+M_LNPI+log(sigsq)));
#ifdef MATLAB_DEBUG_OLD
      SETSCAL("tscal",logz.sum());
      DOCOMP(crit,"%crit=%crit+(-%tscal+0.5*%uvec'*%uvec/%sigsq-0.5*(%n-%m)*log(2*pi*%sigsq));","crit");
#endif // MATLAB_D
      for (i=0; i<n; i++) {
	h=(*means)[i]; rho=(*vars)[i]/sigsq;
	// DEBUG:
	//if (!doDebug) {
	pi=sites->getSitePi()[i]; b=sites->getSiteB()[i];
	//} else {
	//pi=epDebug->sitePi[i]; b=epDebug->siteB[i];
	//}  
	crit+=(temp=0.5*(log1p(-pi*rho)-
			 (pi*h*h-2.0*h*b+rho*b*b)/(sigsq-pi*(*vars)[i])));
#ifdef MATLAB_DEBUG_OLD
	SETSCAL("i",i+1);
	DOCOMP(crit,"%h=%means(%i); %rho=%vars(%i)/%sigsq; %p=%sitePi(%i); %b=%siteB(%i); %crit=%crit+0.5*(log(1-%p*%rho)-(%p*%h*%h-2*%h*%b+%rho*%b*%b)/(%sigsq-%p*%vars(%i)));","crit");
#endif // MATLAB_D
      }
    } else {
      // SBL / Girolami
      logz->div(1.0,sites->getSitePi());
      crit+=0.5*(tau*tau*logz->sum()+
		 ((double) m)*(M_LN2+M_LNPI+log(sigsq))+
		 uvec->inner(*uvec)/sigsq);
      if (meth==2)
	crit+=(1.5*sites->getSitePi().logDet()-
	       ((double) n)*(2.0*log(tau)-M_LN2));
      else
	crit-=0.5*((double) n)*(M_LNPI-M_LN2+2.0*log(tau));
    }
    if (isDegenerate()) {
      crit+=0.5*(sites->getSitePi().logDet()+factL.logDetChol());
      logz->addprod(prior->getBVec(),1.0,sites->getSiteB());
      crit-=0.5*logz->inner(*means)/sigsq;
    } else {
      // DEBUG:
      //if (!doDebug) {
      crit-=0.5*gamma.inner(gamma)/sigsq;
      crit+=0.5*factL.logDetChol();
      /*
	} else {
	crit-=0.5*epDebug->gamma.inner(epDebug->gamma)/epDebug->sigsq;
	#ifdef MATLAB_DEBUG_OLD
	DOCOMP(crit,"%crit=%crit-0.5*%gamma'*%gamma/%sigsq;","crit");
	#endif // MATLAB_D
	crit+=0.5*epDebug->factL.logDetChol();
	#ifdef MATLAB_DEBUG_OLD
	DOCOMP(crit,"%crit=%crit+sum(log(diag(%factL)));","crit");
	#endif // MATLAB_D
	}
      */
    }
#ifdef MATLAB_DEBUG_OLD
    DOCOMP(crit,"","crit");
    if (gradStat>0) {
      DOCOMP(gradSig,"","gradSig");
      if (gradStat>1)
	DOCOMP(gradX,"","gradX");
    }
#endif // MATLAB_D
    // EP: Convert 'gradSig', 'gradTau' from old to new convention
    // Old: tilde{tau} = tau/sigma was tau
    if (gradStat>0 && meth<2) {
      temp=sqrt(sigsq);
      gradSig+=0.5*gradTau*tau*temp;
      gradTau/=temp;
    }

    return crit;
  }

  void ExpectPropLinear::priorChange()
  {
    int n=prior->cols(),m=prior->rows(); // m already incremented, X'
    StVector xvec;
    double vscal;

    CHECKCONSISTENT;
    UPDATE_REPRES;
    prior->getVVec(*tempvec); vscal=(*tempvec)[m-1]; // v_*
    prior->getRow(m-1,*tempvec);
    xvec.reassign(*tempvec); // x (new row of X')
    if (isDegenerate()) {
      // Degenerate repres. (after update, therefore also before)
      xvec.div(sites->getSitePi()); // Pi^-1 x
      // X' Pi^-1 x, so last elem. is x^T Pi^-1 x:
      prior->mult(xvec,tempvec2);
      double beta=tempvec2[m-1]+1.0; // 1 + x^T Pi^-1 x
      tempvec2.remove(m-1);
      factL.backsubst(tempvec2); // l = L^-1 X Pi^-1 x
      gamma.addprod(vscal,tempvec2); // += v_* l
      *tempmat=StMatrix::mask(gamma,false);
      StVector newb;
      newb.fill(1,xvec.inner(prior->getBVec())+xvec.inner(sites->getSiteB()));
      ArrayHandle<StMatrix*> dragx(1);
      dragx[0]=tempmat.p();
      ArrayHandle<StVector*> dragb(1);
      dragb[0]=&newb;
      if (factL.cholext(tempvec2,beta,true,dragx,dragb)!=0)
	throw NumericalException(EXCEPT_MSG(""));
      gamma=(StVector&) (*tempmat)(0);
    } else if (m==nonDegenM)
      // Switch from degenerate to non-degen. repres.
      refresh();
    else {
      // Non-degen. repres. (before and after update)
      StMatrix matmsk;
      matmsk.reassign(gamma,false); // gamma as row vector
      factL.cholUpdRk1(xvec,tempvec2,tempvec3,xvec,&matmsk,
		       &((StVector&) StVector::mask(&vscal,1)));
    }
  }

  void ExpectPropLinear::reset()
  {
    int n=prior->cols();

    CHECKCONSISTENT;
    sites->setSiteB(0.0); sites->setSitePi(eps);
    refresh();
  }

  void ExpectPropLinear::initRepres(const StVector& newb,const StVector& newpi,
				    const StMatrix* newl,
				    const StVector* newgam)
  {
    int n=prior->cols(),m=prior->rows();
    int sz=isDegenerate()?m:n;

    CHECKCONSISTENT;
    if (newb.size()!=n || newpi.size()!=n)
      throw WrongDimensionException(EXCEPT_MSG(""));
    double lbnd=(m<n)?PI_LOWER_BOUND:0.0;
    if (!newpi.checkBounds(Interval<double>(lbnd,0.0,IntVal::ivClosed,
					    IntVal::ivInf)))
      throw InvalidParameterException(EXCEPT_MSG("newpi"));
    sites->setSiteB(newb); sites->setSitePi(newpi);
    if (newl!=0) {
      if (newgam==0 || newl->rows()!=sz || newl->cols()!=sz ||
	  newgam->size()!=sz)
	throw WrongDimensionException(EXCEPT_MSG(""));
      factL=*newl; factL.setStrctPatt(MatStrct::lower);
      gamma=*newgam;
    } else {
      if (newgam!=0) throw InvalidParameterException(EXCEPT_MSG(""));
      refresh();
    }
    up2date=true;
  }

  void ExpectPropLinear::getMargMoments(StVector& means)
  {
    CHECKCONSISTENT;
    UPDATE_REPRES;
    if (isDegenerate()) {
      // Degenerate repres.
#ifndef MATLAB_MEX
      *tempvec=gamma;
      factL.backsubst(*tempvec,true);
      prior->multT(*tempvec,means); // X^T L^-T gamma
#else
      // '*tempvec' is a mask, must not change size!
      StVector vecmsk;
      vecmsk.reassign(*tempvec,Range(0,prior->rows()-1));
      vecmsk=gamma;
      factL.backsubst(vecmsk,true);
      prior->multT(vecmsk,means); // X^T L^-T gamma
#endif
      means.addprod(prior->getBVec(),-1.0,means);
      means.addprod(1.0,sites->getSiteB());
      means.div(sites->getSitePi());
    } else {
      // Non-degen. case
      means=gamma;
      factL.backsubst(means,true); // L^-T gamma
    }
  }

  /*
   * ATTENTION: Must not use 'tempvec2', 'tempvec3'. Otherwise, call
   * from 'compCritGrad' fails!
   *
   * NOTE: Modified prev. implementation for degenerate repres. Now,
   * 'prior->getXMat' is not called, and no temp. matrix is used.
   * Rather, we do a m-loop, one X^T MVM per iter.
   */
  void ExpectPropLinear::getMargMoments(StVector& means,StVector& vars)
  {
    double sigsq=prior->getSigSq();

    CHECKCONSISTENT;
    UPDATE_REPRES;
    if (isDegenerate()) {
      // Degenerate repres.
      int i,m=prior->rows(),n=prior->cols();
      StVector tvec;
      // Compute e1, e2 (m-loop)
      means.zeros(n); vars.zeros(n);
      for (i=0; i<m; i++) {
	tvec.zeros(m); tvec[i]=1.0;
	factL.backsubst(tvec,true);
	prior->multT(tvec,*tempvec); // v = X^T L^-T delta_i
	means.addprod(gamma[i],*tempvec);
	vars.addprod(*tempvec,*tempvec);
      }
      // Mean, variances from e1, e2
      means.addprod(prior->getBVec(),-1.0,means);
      means.addprod(1.0,sites->getSiteB());
      means.div(sites->getSitePi());
      vars.div(sites->getSitePi());
      vars.addscal(-1.0);
      vars.div(sites->getSitePi());
      vars.prod(-sigsq);
    } else {
      // Non-degen. repres.
      means=gamma;
      factL.backsubst(means,true); // Means
      factL.diagInvChol(vars);
      vars.prod(sigsq); // Variances
    }
  }

  // FLORIAN
  void ExpectPropLinear::addMargVarRank(BaseVector<int>& ranks)
  {
    int n=prior->cols();
    StVector vars;
    BaseVector<int> h(n,0);

    CHECKCONSISTENT;
    if (ranks.size()!=n) throw WrongDimensionException(EXCEPT_MSG(""));
    UPDATE_REPRES;
    if (isDegenerate()) {
      // Degenerate repres.
      throw NotImplemException(EXCEPT_MSG("Not implemented"));
    } else {
      // Non-degen. case
      factL.diagInvChol(vars); // Variances / sigma^2
    }
    vars.sort(false,&h);
    for (int i = 0; i < n; i++) ranks[i] += h[i];
  }
  // END FLORIAN

  void ExpectPropLinear::sample(int num,StMatrix& smat,Handle<Generator>& gen)
  {
    int i,j,n=prior->cols(),m=prior->rows();
    double temp;

    CHECKCONSISTENT;
    if (num<1) throw InvalidParameterException(EXCEPT_MSG(""));
    UPDATE_REPRES;
    // Sample N(0,sigma^2) matrix -> S
    smat.apply0(mybindarg(ptr_fun(Random::devNormalRatUnif),gen.p()),n,num);
    smat.smul(sqrt(prior->getSigSq()));
    if (isDegenerate()) {
      // Degenerate repres.
      // See paper for how this works
      tempvec->div(1.0,sites->getSitePi()); // Pi^-1
      prior->outerProd(*tempvec,*tempmat2);
      tempmat2->addseye(1.0); // I + X Pi^-1 X^T
      tempmat2->setStrctPatt(MatStrct::lower);
      tempvec->apply1(sites->getSitePi(),ptr_fun<double, double>(sqrt));
      tempvec->div(1.0,*tempvec);
      smat.mulDiag(*tempvec,smat); // Pi^-1/2 S
      if (tempmat2->eigenDecomp(*tempvec,*tempmat)!=0) // U D U^T
	throw NumericalException(EXCEPT_MSG(""));
      for (i=0; i<m; i++) { // R from D
	temp=sqrt((*tempvec)[i]);
	(*tempvec)[i]=1.0/(temp*(temp+1.0));
      }
      // U, R in 'tempmat', 'tempvec'
      tempmat2->setStrctPatt(MatStrct::normal);
      prior->mult(smat,*tempmat2); // X ...
      tempmat3->mul(*tempmat,true,*tempmat2,false); // U^T ...
      tempmat3->mulDiag(*tempvec,*tempmat3); // R ...
      tempmat2->mul(*tempmat,false,*tempmat3,false); // U ...
      prior->multT(*tempmat2,*tempmat3); // X^T ...
      tempvec->div(1.0,sites->getSitePi());
      tempmat3->mulDiag(*tempvec,*tempmat3); // Pi^-1 ...
      smat.addsmul(*tempmat3,-1.0); // Has covar. sigma^2 Sigma now
      getMargMoments(*tempvec); // h
      smat.addVec(*tempvec); // Add means h
    } else {
      // Non-degen. repres.
      smat.addVec(gamma);
      factL.backsubst(smat,true,true);
    }
  }

  void ExpectPropLinear::compScores(const StVector& v,double& ent,
				    double& logprob,double& exl1)
  {
    int i,n=prior->cols(),m=prior->rows();
    double logdet,sqa,temp,temp2;
    StVector means,vars;
    double sigsq=prior->getSigSq();

    CHECKCONSISTENT;
    if (v.size()!=n) throw WrongDimensionException(EXCEPT_MSG(""));
    getMargMoments(means,vars); // h, vars
    // Expected L_1
    for (i=0,exl1=0.0; i<n; i++) {
      sqa=sqrt(vars[i]); temp=(v[i]-means[i])/sqa;
      exl1+=sqa*(2.0*Specfun::pdfNormal(temp)+
		 temp*(2.0*Specfun::cdfNormal(temp)-1.0));
    }
    means.addprod(v,-1.0,means); // v-h
    if (isDegenerate()) {
      logdet=-sites->getSitePi().logDet()-factL.logDetChol(); // log|Sigma|
      prior->mult(means,vars); // X (v-h)
      temp=-0.5*(means.inner(means,&(sites->getSitePi()))+
		 vars.inner(vars)/sigsq);
    } else {
      logdet=-factL.logDetChol(); // log|Sigma|
      factL.triMulVec(means,true); // L^T (v-h)
      temp=-0.5*means.inner(means)/sigsq;
    }
    temp2=(double) n;
    logdet+=(temp2*(M_LN2+M_LNPI+log(sigsq)));
    ent=0.5*(logdet+temp2);
    logprob=temp-0.5*logdet;
  }

  // DEBUG!!
  void ExpectPropLinear::debugGetPostCov(StMatrix& amat,bool inv)
  {
    if (inv || isDegenerate()) {
      prior->outerTProd(amat);
      amat.diag()->addprod(1.0,sites->getSitePi()); // X^T X + Pi
    }
    if (!inv) {
      if (isDegenerate()) {
	amat.setStrctPatt(MatStrct::lower);
	if (amat.cholDecomp(true)!=0)
	  throw NumericalException(EXCEPT_MSG(""));
      } else {
	UPDATE_REPRES;
	amat=factL; amat.setStrctPatt(MatStrct::lower);
      }
      amat.invChol(); amat.setStrctPatt(MatStrct::normal);
    }
    amat.smul(prior->getSigSq());
  }

  /*
   * If x is InvGauss(mu,1), then lam*x is InvGauss(mu*lam,lam).
   * For very large mu, InvGauss(mu,1) is essentially the same as
   * InvGamma(1/2,1/2), so x = 2/y with y being Gamma(1/2,1).
   */
  void ExpectPropLinear::mcmcSamplePi(const StVector& curra,
				      Handle<Generator>& gen)
  {
    int i,n=prior->cols();
    double tau,temp,sigsq=prior->getSigSq();
    const EPSingleLaplaceManager* lsites=
      DYNCAST(const EPSingleLaplaceManager,sites.p());

    // Sample inverse Gaussian variates
    if (lsites==0)
      throw WrongStatusException(EXCEPT_MSG("Need Laplace sites!"));
    tau=lsites->getTau();
    tempvec2.apply1(curra,ptr_fun<double, double>(fabs));
    tempvec2.prod(tau/sqrt(sigsq)); // 1/mu (for lam=1)
    for (i=0; i<n; i++)
      if ((temp=tempvec2[i])<1e-30)
	tempvec2[i]=2.0/Random::devGamma(0.5,gen);
      else
	tempvec2[i]=Random::devInvGauss(1.0/temp,gen);
    tempvec2.prod(tau*tau);
    // Update repres.
    mcmcInitRepres(tempvec2);
  }

  double ExpectPropLinear::mcmcEvalFullCond(const StVector& a)
  {
    int n=prior->cols();
    double temp,sigsq=prior->getSigSq();

    CHECKCONSISTENT;
    if (a.size()!=n) throw InvalidParameterException(EXCEPT_MSG(""));
    if (isDegenerate())
      throw NotImplemException(EXCEPT_MSG(""));
    UPDATE_REPRES;

    // Implemented for non-degen. repres. only right now
    temp=0.5*((double) n);
    *tempvec=a; factL.triMulVec(*tempvec,true);
    tempvec->addprod(-1.0,gamma); // L^T (a-h) = L^T a - gamma

    return -temp*(M_LN2+M_LNPI+log(sigsq))+0.5*factL.logDetChol()-
      0.5*tempvec->inner(*tempvec)/sigsq;
  }
//ENDNS
