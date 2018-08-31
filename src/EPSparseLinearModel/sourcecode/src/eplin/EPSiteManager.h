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
 * Desc.:  Header abstract class EPSiteManager
 * ------------------------------------------------------------------- */

#ifndef EPLIN_EPSITEMANAGER_H
#define EPLIN_EPSITEMANAGER_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eplin/default.h"

//BEGINNS(eplin)
  /**
   * Abstract base of all EP site managers.
   * A site manager gives access to a number of sites, and to corr. local
   * services. These can be moment computations requiring quadrature, or
   * computations depending on vectors assoc. with the site.
   * The parameters of the sites are also maintained here.
   * <p>
   * ATTENTION: Some diamond shapes start with this one. Use virtual
   * inheritance. When coming together, the methods here have to be
   * overwritten again, even if parent classes do so, to break ambiguity.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class EPSiteManager
  {
  public:
    // Public methods

    virtual ~EPSiteManager() {}

    /**
     * @return Number of sites
     */
    virtual int size() const = 0;
  };
//ENDNS

#endif
