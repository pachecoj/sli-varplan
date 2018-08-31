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
 * Desc.:  Header class ExpPropPars
 * ------------------------------------------------------------------- */

#ifndef EPLIN_EXPPROPPARS_H
#define EPLIN_EXPPROPPARS_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eplin/default.h"

//BEGINNS(eplin)
  /**
   * Collects parameters to configure the EP algorithm.
   * If all abs. changes of site pars. are < 'delta' during a sweep, the
   * algorithm is stopped.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class ExpPropPars
  {
  public:

    int maxiter;    // Max. number of sweeps
    double delta;
    double thres;   // See 'ExpectPropLinear::updateSite'
    double frac;    // Fractional updates. See 'ExpectPropLinear::updateSite'
    double eps;     // Init. value for pi_i (>0)
    int nonDegenM;  // See 'ExpectPropLinear'. 0: always use non-degen. repr.
    int verbose;    // Verbosity level (0: no messages)

    ExpPropPars() : maxiter(30),delta(1e-7),thres(1e-12),frac(1.0),eps(1.0),
		    verbose(0),nonDegenM(0) {}
  };
//ENDNS

#endif
