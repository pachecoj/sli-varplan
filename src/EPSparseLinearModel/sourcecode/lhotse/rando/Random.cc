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
 * Module: rando
 * Desc.:  Definition of class Random
 * ------------------------------------------------------------------- */

#include "lhotse/rando/Random.h"
#include "lhotse/specfun/Specfun.h"
#include <algorithm>
#include "lhotse/matrix/BaseVector.h"
#include "lhotse/matrix/BaseMatrix.h"
#include "lhotse/matrix/StVector.h"
#include "lhotse/matrix/StMatrix.h"

//BEGINNS(rando)
  // Static members

  Generator* Random::defGen=0;

  // Public static methods

  double Random::devNormal(Generator* gen)
  {
    // Standard Box-Muller
    static bool oneMore=false;
    static double second;

    if (gen==0) throw NoGeneratorException();
    if (!oneMore) {
      register double v1,v2,rsq,fac;
      do {
	v1=2.0*gen->get()-1.0;
	v2=2.0*gen->get()-1.0;
	rsq=v1*v1+v2*v2;
      } while (rsq>=1.0);
      if (rsq==0.0) {
	// Happens with prob 0!
	second=0.0; oneMore=true;
	return 0.0;
      }
      fac=sqrt(-2.0*log(rsq)/rsq);
      second=v2*fac; oneMore=true;
      return v1*fac;
    } else {
      oneMore=false;
      return second;
    }
  }

  /*
   * See Devroye, p. 194ff, for details. The constants are, for
   * f(x) = exp(-0.5*x*x),
   * b = sup sqrt(x) = 1,
   * a_+ = sup y sqrt(x) = sqrt(2/e)
   * a_- = -(a_+)
   * We generate (u,v) uniformly on A, where
   * A = {(u,v) | 0 <= u <= sqrt(f(v/u))},
   * by rejection from the rectangle [0,b] x [a_-,a_+].
   * Then, sqrt(v/u) has standard normal density \propto f(x)
   *
   * The final acc. test is: x^2 <= -4*log(u).
   * We use the bounds
   *   2*((u-2)^2 - 1) <= -4*log(u) <= -2*(u-1/u)
   * to squeeze the logarithm.
   */
  double Random::devNormalRatUnif(Generator* gen)
  {
    if (gen==0) throw NoGeneratorException();
    static const double ap2=1.71552776992141359294; // 2*a_+
    register double u,x,sx,temp;

    do {
      // u uniformly on (0,b), b=1, v uniformly on (a_-,a_+)
      u=gen->get(); // from U[0,b]
      x=(gen->get()-0.5)*ap2/u; // x=v/u, v from U[a_-,a_+]
      sx=-0.5*x*x;
      // Test quick acceptance condition
      temp=u-2.0;
      if (sx >= 1.0-temp*temp) break; // accept
      /* Test reject conditions: The quick test is on the left, therefore
	 is evaluated first! */
    } while (sx<=u-1.0/u || 2.0*log(u)>sx);
    return x;
  }

  /*
   * For a>=1, we use Best's rejection algorithm from a t distribution with
   * 2 dof's. This is rather heavy-tailed and therefore can be used as
   * dominating curve for a properly scaled Gamma density.
   *
   * For a<1, we use the transformation of an exponential power distribution
   * (EPD) variate. The latter is obtained by the rejection method (see
   * Devroye, p. 304).
   */
  double Random::devGamma(double a,Generator* gen)
  {
    if (gen==0) throw NoGeneratorException();
    if (a<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
    if (a>=1.0) {
      /* Best's rejection algorithm (see Devroye, p.410). Note that b used
	 here has nothing to do with the scale parameter b (the latter is
	 1 here). */
      const double b=a-1.0,c=3.0*a-0.75;
      register double u,vw,w,y,x,z;
      do {
	do {
	  u=gen->get();
	  w=u*(1.0-u); y=(u-0.5)*sqrt(c/w); x=y+b;
	} while (x<0);
	vw=gen->get()*w; // v*w
	y=64.0*w*vw*vw; // 4^3*w^3*v^2
	if (z <= 1.0-2.0*y*y/x) break; // quick acceptance test
      } while (log(z) > 2.0*(b*log(x/b)-y)); // rejection test
      return x;
    } else {
      // See Devroye, p.304
      register double u,e,estar,x,lim;
      do {
	u=gen->get(); estar=-log(gen->get()); // estar is exponential var.
	if (u <= 1-a) {
	  // Here, u is uniform on [0,1-a]
	  x=exp(log(u)/a);
	  lim=estar;
	} else {
	  /* Here, u is uniform on [1-a,1], i.e. (1-u)/a is uniform on [0,1].
	     We can therefore generate an exponential variate e is follows: */
	  e=-log((1-u)/a);
	  x=exp(log(1-a+a*e)/a);
	  lim=estar+e;
	}
      } while (x>lim);
      return x;
    }
  }

  /*
   * Cheng's rejection method (Devroye, p. 438) uses rejection from a
   * Burr XII density. If both a,b are smaller than 1 and a+b <= 1.5,
   * we use Johnk's beta generator which employs the property of theorem
   * 3.4 in Devroye, p. 416.
   */
  double Random::devBeta(double a,double b,Generator* gen)
  {
    static double olda=-1.0,oldb=-1.0;
    static double s,gamma,u; // precomputed variables
    static const double log4=1.38629436111989061883;

    if (gen==0) throw NoGeneratorException();
    if (a<=0.0 || b<=0.0) throw InvalidParameterException(EXCEPT_MSG(""));
    if (a+b > 1.5 || std::max(a,b) >= 1.0) {
      // Cheng's rejection method (see Devroye, p.438)
      // Setup
      if (olda!=a || oldb!=b) {
	olda=a; oldb=b;
	s=a+b;
	double temp=std::min(a,b);
	if (temp <= 1.0)
	  gamma=temp;
	else
	  gamma=sqrt((2.0*a*b-s)/(s-2.0));
	u=a+gamma;
      }
      // Generator
      register double u1,v,y;
      do {
	u1=gen->get();
	v=u1/(gamma*(1.0-u1)); y=a*exp(v);
      } while (s*log(s/(b+y))+u*v-log4 < log(u1*u1*gen->get()));
      return y/(b+y);
    } else {
      // Johnk's generator (see Devroye, p.418)
      register double x,y;
      do {
	x=exp(log(gen->get())/a);
	y=exp(log(gen->get())/b);
      } while (x+y > 1.0);
      return x/(x+y);
    }
  }

  double Random::devT(double a,Generator* gen)
  {
    static double olda=-1.0;
    static double ap2,s,fact;
    if (gen==0) throw NoGeneratorException();
    if (a<=0) throw InvalidParameterException(EXCEPT_MSG(""));
    if (a >= 1) {
      // Ratio-of-uniforms method (see Devroye, p.201)
      if (a!=olda) {
	// Setup
	olda=a;
	s=(a+1.0)/4.0;
	ap2=sqrt(8*a)*exp((s-0.5)*log(a-1.0) - s*log(a+1.0)); // 2*a_+
	fact=exp(s*log(1.0+1.0/a));
      }
      register double u,x,sx,temp;
      bool age3=(a>=3.0);
      do {
	/* u uniformly on (0,b), b=1, v uniformly on (a_-,a_+), x=v/u.
	   As the density is symmetric around 0, a_- = -(a_+). */
	u=gen->get();
	x=(gen->get()-0.5)*ap2/u;
	sx=x*x;
	// Test quick acceptance condition
	temp=u*fact;
	if (sx <= 5.0-4.0*temp) break; // accept
	/* Test reject conditions. Note that the quick test is on the left
	   and will be evaluated first. */
      } while ((age3 && sx>=-3.0+4.0/temp) || sx>a*(exp(-log(u)/s)-1.0));
      return x;
    } else {
      /* If N is N(0,1) and G is Gamma(a/2), N,G independent, then
	 sqrt(2*a)*N/sqrt(G) is t_a distributed (Devroye, p. 445). */
      return devNormalRatUnif(gen)*sqrt(2.0*a/devGamma(0.5*a,gen));
    }
  }

  void Random::randomPerm(int n,BaseVector<int>& iv,Generator* gen,bool jumble)
  {
    // Swapping method, see Devroye p.646
    register int i,pos;
    double id;

    if (n<=0) throw WrongDimensionException(EXCEPT_MSG(""));
    if (gen==0) throw NoGeneratorException();
    if (jumble) {
      if (iv.size()!=n)
	throw WrongDimensionException("'iv' must be of size 'n'");
    } else {
      iv.fill(n,0);
      for (i=1; i<n; i++) iv[i]=i;
    }
    for (i=n-1,id=(double) n; i>=1; i--,id-=1.0) {
      pos=(int) floor(((double) gen->get())*id);
      if (pos==(i+1)) pos=i;
      iv.swap(pos,i); // swap pos <-> i
    }
  }

  void Random::randomSet(int n,int d,BaseVector<int>& iv,Generator* gen,
			 bool shift)
  {
    register int i,pos;
    double id;
    BaseVector<int> workInd; // mask, or working copy

    if (n<=0 || d<=0 || d>n) throw InvalidParameterException(EXCEPT_MSG(""));
    if (gen==0) throw NoGeneratorException();
    if (!shift) {
      workInd.fill(n,0);
      for (i=1; i<n; i++) workInd[i]=i; // identity
    } else {
      if (iv.size()!=n) throw WrongDimensionException("iv");
      workInd.reassign(iv); // mask
    }
    for (i=0,id=(double) n; i<d; i++,id-=1.0) {
      pos=i+((int) floor(((double) gen->get())*id));
      if (pos==n) pos--;
      workInd.swap(pos,i); // swap pos <-> i
    }
    if (!shift)
      iv=workInd(Range(0,d-1)); // copy back
  }

  void Random::devMultiNormal(BaseVector<double>& targ,Generator* gen)
  {
    int i,dim=targ.size();

    if (gen==0) throw NoGeneratorException();
    if (dim==0) throw WrongDimensionException(EXCEPT_MSG(""));
    for (i=0; i<dim; i++)
      targ[i]=devNormalRatUnif(gen);
  }

  void Random::devMultiNormal(BaseVector<double>& targ,
			      const BaseVector<double>& mean,Generator* gen)
  {
    int i,dim=targ.size();

    if (gen==0) throw NoGeneratorException();
    if (dim==0 || dim!=mean.size()) throw WrongDimensionException(EXCEPT_MSG(""));
    for (i=0; i<dim; i++)
      targ[i]=devNormalRatUnif(gen)+mean[i];
  }

  /*
   * If u\sim N(0,I), then x = L u + \mu is N(\mu,\Sigma), where
   * \Sigma = L L' (Cholesky decomp.).
   */
  void Random::devMultiNormal(StVector& targ,
			      const StVector& mean,
			      const StMatrix& chCov,
			      Generator* gen)
  {
    devMultiNormal(targ,gen); // u
    uchar spatt=chCov.getStrctPatt(); // copy
    chCov.setStrctPatt(WriteBackMat<double>::strctLower);
    chCov.triMulVec(targ); // L u
    chCov.setStrctPatt(spatt); // reset
    targ.addprod(1.0,mean); // L u + mu
  }

  /*
   * If u\sim N(0,I), then x = (L')^{-1} (u + L^{-1} r) is N(\mu,\Sigma),
   * where \Sigma^{-1} = L L' and \mu = \Sigma r.
   */
  void Random::devMultiNormalInv(StVector& targ,
				 const StVector& r,
				 const StMatrix& chICov,
				 Generator* gen)
  {
    uchar spatt=chICov.getStrctPatt(); // copy
    chICov.setStrctPatt(WriteBackMat<double>::strctLower);
    targ=r;
    chICov.backsubst(targ); // L^-1 r
    for (int i=0; i<targ.size(); i++)
      targ[i]+=devNormalRatUnif(gen); // add u
    chICov.backsubst(targ,true); // L^-T (...)
    chICov.setStrctPatt(spatt); // reset
  }

  /*
   * It is easy to show that if we sample x_i independently from
   * Gamma(\alpha_i,\beta), i=1,\dots,k, set s=\sum_i x_i and
   * \theta_i = x_i/s, then:
   * - (\theta_1,\dots,\theta_k) is Dirichlet(\alpha)
   * - s is Gamma(\sum_i \alpha_i,\beta)
   * - s and the \theta_i are independent
   */
  void Random::devDirichlet(const BaseVector<double>& alpha,
			    BaseVector<double>& rv,Generator* gen)
  {
    register int i,k=alpha.size();
    register double sum;

    if (k<1) throw InvalidParameterException(EXCEPT_MSG(""));
    rv.fill(k,0.0);
    if (k==1) {
      rv[0]=1.0;
      return;
    }
    for (i=0,sum=0.0; i<k; i++)
      sum+=(rv[i]=devGamma(alpha[i],gen));
    for (i=0; i<k; i++)
      rv[i]/=sum;
  }
//ENDNS
