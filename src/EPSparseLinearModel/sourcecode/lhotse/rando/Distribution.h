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
 * Desc.:  Header abstract class Distribution
 * ------------------------------------------------------------------- */

#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/rando/default.h"

//BEGINNS(rando)
  /**
   * Abstract base class of probability distribution hierarchy. Childs of
   * this class represent univariate probability distributions and must
   * export methods to evaluate the pdf and cdf of the distribution. If the
   * distribution is discrete, the pdf is meant to be a at points of non-zero
   * mass a and 0 elsewhere.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class Distribution
  {
  public:
    // Public methods

    virtual ~Distribution() {}

    /**
     * Returns value of pdf at given point
     */
    virtual double pdf(double x) const = 0;

    /**
     * Returns value of cdf at given point
     */
    virtual double cdf(double x) const = 0;
  };
//ENDNS

#endif
