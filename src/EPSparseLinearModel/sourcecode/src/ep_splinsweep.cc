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
 * Module: GLOBAL
 * Desc.:  MEX function for doing an EP sweep for the sparse linear
 *         model (Laplace prior), fixed noise variance sigma^2.
 *         Does not use the LHOTSE MEX gateway
 * ------------------------------------------------------------------- */

#include "lhotse/matif/MatlabTools.h"
#include "lhotse/matif/MatlabString.h"
#include "lhotse/matif/MatlabMatrix.h"
#include "lhotse/matrix/StVector.h"
#include "lhotse/matrix/StMatrix.h"
#include "lhotse/matrix/IndexVector.h"
#include "lhotse/rando/GenstateKISS.h"
#include "lhotse/rando/Random.h"
#include "src/eplin/MaskEPPrior.h"
#include "src/eplin/EPSingleLaplaceManager.h"
#include "src/eplin/ExpectPropLinear.h"
#include "lhotse/PrintToStdout.h"
//#include "lhotse/DebugMacros.h" // DEBUG
#ifdef MATLAB_DEBUG_OLD
#include "lhotse/MatlabDebug.h"
#endif // MATLAB_D

/**
 * MEX function EP_SPLINSWEEP
 *
 * Does a single EP sweep over selected sites of a sparse linear model with
 * Laplace prior. The likelihood is N(U | X*A, SIGSQ*EYE), and each prior
 * factor has the form
 *   (TT/2)*EXP( -TT*ABS(A(i)) ),
 * where TT = TAU/SQRT(SIGSQ). For all details, see 'ExpectPropLinear' in
 * LHOTSE project module 'eplin'.
 *
 * Input:
 * - X:       Design matrix [m-by-n]
 * - U:       Vector of responses [m]
 * - SIGSQ:   Noise variance
 * - TAU:     Scale parameter of Laplace prior
 * - SITEPI:  Site parameters pi_i [n] (*)
 * - SITEB:   Site parameters b_i [n] (*)
 * - UPDIND:  Index of sites to be updated. Entries in 1:n. Can have any
 *            size.
 *            Can be empty, for example if only LOGZ, LOGTILZ are to
 *            be computed.
 * - L:       Cholesky factor of representation [n-by-n / m-by-m] (*)
 * - GAMMA:   Vector of representation [n / m] (*)
 * - PITHRES: Threshold det. when updates are skipped because too small
 *            (see 'ExpectPropLinear::updateSite')
 * - REINIT:  Boolean flag. If true, the representation is recomputed
 *            from the site pars. before updates are done
 * - WKVEC1:  Scratch vector [n] (*)
 * - WKVEC2:  Scratch vector [n] (*)
 * - WKVEC3:  Scratch vector [n] (*)
 * - METH:    Approximate inference method to be used:
 *            0: Standard EP [default]
 *            2: Tipping SBL
 *            3: Girolami variational, VMFB
 * - REFRESH: If true, the repres. is recomp. from scratch after
 *            updates are done. Optional. Def.: true
 * - FRAC:    Value in (0,1]. If FRAC<1, a fractional EP update is done.
 *            See below. Can only be used with METH==0. Def.: 1
 *
 * Return:
 * - DELTA:   Optional. Maximum rel. marginal change over all sites
 *            (see 'ExpectPropLinear::updateSite')
 * - NUMUPD:  Optional. Number of updates done (updates are not done if
 *            change in pi_i is below PITHRES)
 *
 * We use the undocumented feature that the content of input arguments
 * can be overwritten to get a "call by reference". This is done for the
 * (*) arguments, whose content is overwritten. We cannot change the
 * size that way, so these matrices must be of correct size, even if their
 * input content is not used here.
 *
 * Different representations:
 * If m<n, you can use the non-degenerate representation (L n-by-n, GAMMA n)
 * or the degenerate one (L m-by-m, GAMMA m). If m>=n, only the non-degen.
 * can be used.
 * The kind of representation to be used is determined by the size of
 * L, GAMMA passed (even if REINIT==1).
 * The non-degen. is more numerically stable in general, and is more
 * efficient than the degen. one for about m >= n/2. If m << n, use the
 * degen. one to save time. Otherwise, or if max(n,m) is not large, we
 * recommend using the non-degen. one.
 *
 * REFRESH and REINIT:
 * The low-rank updates of the representation after each site introduce
 * numerical errors. It is wise to "refresh" the representation by recomp.
 * it from scratch after each n updates.
 * If REFRESH==true, a refresh is done at the end. If REINIT==true, the
 * repres. is recomp. at the beginning. The default usage of this function
 * is to use a permutation of 1:n in UPDIND, together with REFRESH==true
 * and REINIT==false, but REINIT==true for the first sweep in order to
 * initialize the representation.
 * NOTE: A refresh after each sweep is fairly conservative, but
 * subdominant anyway. A complete EP sweep takes far longer than a
 * refresh, although they both have the same complexity.
 *
 * Fractional EP updates:
 * If FRAC<1, a fractional EP update is done instead of a standard one.
 * This leads to different fixed points. It can improve robustness of the
 * method. In some cases, fractional EP leads to better solutions than
 * standard EP (if both converge properly).
 */
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  try {
    int i,j,m,n,numupd;
    double temp,sigsq,tau,thres,delta,h,rho,pi,b;
    Handle<MatlabMatrix> mxmat,muvec,msitepi,msiteb,ml,mgamma,mlogz,mlogtz,
      mwkvec1,mwkvec2,mwkvec3,mupdind;
    StMatrix xmat,lfact;
    StVector uvec,sitePi,siteB,gamma,wkvec1,wkvec2,wkvec3,updind;
    bool degenRep=false,reinit,refresh=true;
    int nonDegenM;
    int meth=0;
    double frac=1.0;

    // Read arguments
    MatlabTools::checkNrhs(nrhs,14);
    if (nlhs>2)
      throw MatIFException("Too many return arguments");
    mxmat.changeRep(MatlabTools::readMatrix("X",prhs[0]));
    if ((m=mxmat->m())==0 || (n=mxmat->n())==0)
      throw MatIFException("X must not be empty");
    xmat.maskMatlab(*mxmat);
    muvec.changeRep(MatlabTools::readVector("U",prhs[1],m));
    uvec.maskMatlab(*muvec);
    sigsq=MatlabTools::readScalar("SIGSQ",prhs[2],&DefIVal<double>::posit());
    tau=MatlabTools::readScalar("TAU",prhs[3],&DefIVal<double>::posit());
    msitepi.changeRep(MatlabTools::readVector("SITEPI",prhs[4],n));
    sitePi.maskMatlab(*msitepi);
    msiteb.changeRep(MatlabTools::readVector("SITEB",prhs[5],n));
    siteB.maskMatlab(*msiteb);
    Interval<double> tiv(1.0,(double) n,IntVal::ivClosed,IntVal::ivClosed);
    mupdind.changeRep(MatlabTools::readIntVector("UPDIND",prhs[6],0,0,&tiv));
    updind.maskMatlab(*mupdind);
    ml.changeRep(MatlabTools::readMatrix("L",prhs[7]));
    if ((i=ml->m())!=ml->n()) throw MatIFException("L wrong size");
    if (i!=n) {
      if (degenRep=(i==m)) {
	if (m>n)
	  throw MatIFException("Non-degenerate representation cannot be used");
      } else
	throw MatIFException("L wrong size");
    }
    lfact.maskMatlab(*ml); lfact.setStrctPatt(MatStrct::lower);
    mgamma.changeRep(MatlabTools::readVector("GAMMA",prhs[8],i));
    gamma.maskMatlab(*mgamma);
    thres=MatlabTools::readScalar("PITHRES",prhs[9],&DefIVal<double>::posit());
    reinit=MatlabTools::readScalBool("REINIT",prhs[10]);
    mwkvec1.changeRep(MatlabTools::readVector("WKVEC1",prhs[11],n));
    wkvec1.maskMatlab(*mwkvec1);
    mwkvec2.changeRep(MatlabTools::readVector("WKVEC2",prhs[12],n));
    wkvec2.maskMatlab(*mwkvec2);
    mwkvec3.changeRep(MatlabTools::readVector("WKVEC3",prhs[13],n));
    wkvec3.maskMatlab(*mwkvec3);
    if (nrhs>14) {
      meth=MatlabTools::readScalInt("METH",prhs[14],&DefIVal<int>::nonneg());
      if (meth==1 || meth>3)
	throw MatIFException("METH must be in 0,2,3");
      if (nrhs>15) {
	refresh=MatlabTools::readScalBool("REFRESH",prhs[15]);
	if (nrhs>16) {
	  Interval<double> tiv2(0.0,1.0,IntVal::ivOpen,IntVal::ivClosed);
	  frac=MatlabTools::readScalar("FRAC",prhs[16],&tiv2);
	  if (frac<1.0 && meth!=0)
	    throw MatIFException("Cannot use fractional sites with METH!=0");
	}
      }
    }

    // Initialization
    nonDegenM=degenRep?n:0;
    Handle<EPPrior> prior(new MaskEPPrior(xmat,uvec,sigsq));
    Handle<EPSingleLaplaceManager>
      sites(new EPSingleLaplaceManager(siteB,sitePi,tau));
    ExpectPropLinear epobj(0.1,prior,sites,nonDegenM,1,
			   HandleZero<StMatrix>::get(),
			   HandleZero<StMatrix>::get(),
			   HandleZero<StMatrix>::get(),
			   Handle<StVector>(&wkvec1,false));
    epobj.provideBuffers(lfact,gamma,wkvec2,wkvec3);
    if (!reinit)
      epobj.initRepres(siteB,sitePi,&lfact,&gamma);
    else
      epobj.refresh();

    // Loop over sites in UPDIND
    //epobj.setVerbose(2); // DEBUG!
    for (j=numupd=0,delta=0.0; j<updind.size(); j++) {
      i=((int) updind[j])-1; // Matlab is base 1
      //mycout << "i=" << i << "\n";
      try {
	if (epobj.updateSite(i,thres,frac,&temp,meth)) numupd++;
	delta=std::max(delta,temp);
      } catch (StandardException ex) {
	//mycout << ex.msg() << "\n"; // DEBUG!
	delta=std::max(delta,1.0);
      }
    }
    if (refresh) epobj.refresh();

    // Return arguments
    // NOTE: All (*) input arguments are dealt with, thanks to them being
    // masks
    if (nlhs>0) {
      plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
      *(mxGetPr(plhs[0]))=delta;
      if (nlhs>1) {
	plhs[1]=mxCreateDoubleMatrix(1,1,mxREAL);
	*(mxGetPr(plhs[1]))=(double) numupd;
      }
    }
  } catch (MatIFException ex) {
    mexErrMsgTxt(ex.msg());
  } catch (StandardException ex) {
    string msg="Caught LHOTSE exception:\n";
    msg+=ex.msg();
    mexErrMsgTxt(msg.c_str());
  } catch (...) {
    mexErrMsgTxt("Caught unspecified exception");
  }
}
