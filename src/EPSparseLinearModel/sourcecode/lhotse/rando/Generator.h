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
 * Desc.:  Header abstract class Generator
 * ------------------------------------------------------------------- */

#ifndef GENERATOR_H
#define GENERATOR_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/rando/default.h"
#include "lhotse/rando/GeneratorState.h"

//BEGINNS(rando)
  /**
   * Abstract base class of PRN generator hierarchy.
   * A PRN generator in the sense of this class is a dynamical system
   * which advances a hidden state (or seed) and produces a double value
   * within [0,1) as output with each call. The state of a generator may
   * be of arbitrary type, but it must be s.t. in any case, the future
   * output of the generator is a deterministic function of the current
   * state. States are subclasses of 'GeneratorState'.
   * NOTE: The generator state must be maintained as 'mutable', so that
   * methods like 'get' can be const. This is required s.t. this class can
   * be a 'NullaryFunc'.
   * <p>
   * Many generators actually generate random integer values and map them
   * into [0,1). In this case, the method 'canGenUInt' should return true
   * (def.: false), and the method 'getUInt' should return a random (uniform)
   * unsigned int value.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class Generator : public NullaryFunc<double>
  {
  public:
    // Public methods

    virtual ~Generator() {}

    /**
     * Returns sample of uniform pseudorandom deviate within [0,1).
     *
     * @return See above
     */
    virtual double get() const = 0;

    double operator()() const {
      return get();
    }

    /**
     * Restarts generator with new state (seed)
     *
     * @param seed New state
     */
    virtual void restart(const GeneratorState* seed) const = 0;

    /**
     * Returns current state, i.e. the one which would be used to generate
     * the next output.
     *
     * @param seed Current state returned here
     */
    virtual void getSeed(GeneratorState* seed) const = 0;

    /**
     * The method 'getUInt' works iff this method returns true.
     *
     * @return Can unsigned int's be generated directly?
     */
    virtual bool canGenUInt() const {
      return false; // def. implementation: No
    }

    /**
     * See header comment. Works only if 'canGenUInt' returns true!
     * NOTE: Must be implemented if 'canGenUInt' returns true!
     *
     * @return See above
     */
    virtual unsigned int getUInt() const {
      throw NotImplemException(EXCEPT_MSG(""));
    }
  };
//ENDNS

#endif
