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
 * Desc.:  Header abstract class DistribSamp
 * ------------------------------------------------------------------- */

#ifndef DISTRIBSAMP_H
#define DISTRIBSAMP_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/rando/default.h"
#include "lhotse/rando/Generator.h"
#include "lhotse/rando/Distribution.h"
#include "lhotse/rando/Random.h"

//BEGINNS(rando)
  /**
   * Abstract class in probability distribution hierarchy. Represents
   * univariate probability distributions. Apart from methods of Distribution,
   * a method to draw independent samples (given a generator for independent
   * uniform variates) must be provided.
   * <p>
   * Note: If no generator is supplied and no default generator has been set
   * for the DistribSamp instance, the default generator of Random is used
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class DistribSamp : public Distribution
  {
  protected:
    Generator* genDef; // default generator

  public:
    // Public methods

    DistribSamp() : genDef(0) {}

    virtual ~DistribSamp() {}

    /**
     * Produce sample from distribution given a generator
     */
    virtual double sample(Generator* gen=0) const = 0;

    Generator* getGenerator() const {
      return genDef;
    }

    void setGenerator(Generator* gen) {
      genDef=gen;
    }

  protected:
    // Internal methods
    Generator* generator() const {
      return (genDef!=0) ? genDef : Random::getDefGen();
    }
  };
//ENDNS

#endif
