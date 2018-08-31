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
 * Desc.:  Header class TestUtils (DOWNW. COMPAT.)
 * ------------------------------------------------------------------- */

#ifndef TESTUTILS_H
#define TESTUTILS_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/rando/default.h"
#include "lhotse/specfun/Specfun.h"

//BEGINNS(rando)

  /**
   * THIS CLASS EXISTS FOR DOWNWARD COMPATIBILITY ONLY! DO NOT USE IT.
   * USE Specfun (module specfun) INSTEAD!
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class TestUtils
  {
  public:
    // Public static methods

    static double logGamma(double z) {
      return Specfun::logGamma(z);
    }

    static double beta(double z,double w) {
      return Specfun::beta(z,w);
    }

    static double incGamma(double a,double x) {
      return Specfun::incGamma(a,x);
    }

    static double digamma(double z) {
      return Specfun::digamma(z);
    }

    static double derivDigamma(double z) {
      return Specfun::derivDigamma(z);
    }

    static double cdfGamma(double a,double b,double x) {
      return Specfun::cdfGamma(a,b,x);
    }

    static double pdfGamma(double a,double b,double x) {
      return Specfun::pdfGamma(a,b,x);
    }

    static double cdfChiSq(int degFree,double x) {
      return Specfun::cdfChiSq(degFree,x);
    }

    static double cdfNormal(double x) {
      return Specfun::cdfNormal(x);
    }

    static double cdfNormal(double mean,double stddev,double x) {
      return Specfun::cdfNormal(mean,stddev,x);
    }

    static double pdfNormal_GSL(double x) {
      return Specfun::pdfNormal_GSL(x);
    }

    static double pdfNormal(double x) {
      return Specfun::pdfNormal(x);
    }

    static double pdfNormal(double mean,double stddev,double x) {
      return Specfun::pdfNormal(mean,stddev,x);
    }

    static double logPdfNormal(double x) {
      return Specfun::logPdfNormal(x);
    }

    static double cdfExp(double x) {
      return Specfun::cdfExp(x);
    }

    static double cdfExp(double mean,double x) {
      return Specfun::cdfExp(mean,x);
    }

    static double pdfExp(double x) {
      return Specfun::pdfExp(x);
    }

    static double pdfExp(double mean,double x) {
      return Specfun::pdfExp(mean,x);
    }

    static double cdfPoisson(double gamma,int x) {
      return Specfun::cdfPoisson(gamma,x);
    }

    static double pdfPoisson(double gamma,int x) {
      return Specfun::pdfPoisson(gamma,x);
    }

    static double cdfBeta(double a,double b,double x) {
      return Specfun::cdfBeta(a,b,x);
    }

    static double pdfBeta(double a,double b,double x) {
      return Specfun::pdfBeta(a,b,x);
    }

    static double cdfT(double alpha,double x) {
      return Specfun::cdfT(alpha,x);
    }

    static double pdfT(double alpha,double x) {
      return Specfun::pdfT(alpha,x);
    }

    static double pDoubleSidedT(double alpha,double x) {
      return Specfun::pDoubleSidedT(alpha,x);
    }

    static double cdfBinomial(int n,double p,int x) {
      return Specfun::cdfBinomial(n,p,x);
    }

    static double logCdfNormal(double z) {
      return Specfun::logCdfNormal(z);
    }

    static double derivLogCdfNormal(double z) {
      return Specfun::derivLogCdfNormal(z);
    }
  };
//ENDNS

#endif
