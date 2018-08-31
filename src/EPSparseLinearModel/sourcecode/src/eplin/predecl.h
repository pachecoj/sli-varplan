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
 * Desc.:  Module predeclarations
 * ------------------------------------------------------------------- */

#ifndef DE_EPLIN_H
#define DE_EPLIN_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

//BEGINNS(eplin)
  class EPSingleSite;
  class EPGaussLaguerreSite;
  class EPLaplaceSite;
  class EPExponentialSite;
  class EPDebugSite;
  class EPLaplaceDebugSite;
//class EPLaplaceOldSite;
//class EPAdaptQuadSite;
  class EPSiteManager;
  class EPSingleSiteManager;
  class EPJointSiteManager;
  class EPVectorSiteManager;
  class EPSingleVectorSiteManager;
  class EPJointVectorSiteManager;
  class EPSingleLaplaceManager;
  class EPJointLaplaceManager;
  class EPMaskVectorSiteManager;
  class EPSingleVectorLaplaceManager;
  class ExpPropPars;
  class EPPrior;
  class MaskEPPrior;
  class SpLMPrior;
  class SpLMPriorDebug;
  class ExpectPropLinear;
  class ExpectPropJointLinear;
  class ExpectPropDirect;
  class SparseLinModel;
  class SparseImageModel;
  class SparseJointImageModel;
//ENDNS

#endif
