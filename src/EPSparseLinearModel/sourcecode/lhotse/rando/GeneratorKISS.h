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
 * Desc.:  Header class GeneratorKISS
 * ------------------------------------------------------------------- */

#ifndef GENERATORKISS_H
#define GENERATORKISS_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/rando/default.h"
#include "lhotse/rando/Generator.h"
#include "lhotse/rando/GenstateKISS.h"

//BEGINNS(rando)
  /**
   * Implements the KISS PRN generator, using states of type 'GenstateKISS'.
   * This generator is described in:
   *   Marsaglia and Zaman (1993)
   *   The KISS Generator
   *   Technical Report, Dept. of Statistics, Univ. of Florida
   * It runs three unsigned int hidden chains, two of them are advanced using
   * shift-register generators, one using a congruential generator. The output
   * is generated from the sum modulo 2^32, divided by 2^32. The period is of
   * order 2^95. See the report for details, and also:
   *   Robert and Casella
   *   Monte Carlo Statistical Methods, Chap. 2.1.2,
   * which was the source for this implementation.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class GeneratorKISS : public Generator
  {
  protected:
    // Members

    mutable GenstateKISS state; // Current state

  public:
    // Public methods

    /**
     * Constructor.
     *
     * @param initSeed Init. state
     */
    GeneratorKISS(const GenstateKISS& initSeed) : state(initSeed) {}

    /**
     * Returns sample of uniform pseudorandom deviate within [0,1),
     * and advances state.
     *
     * @return See above
     */
    double get() const {
      return ((double) getUInt())/(((double) 0xFFFFFFFF)+1.0);
    }

    /**
     * Restarts generator with new state (seed)
     *
     * @param seed New state
     */
    void restart(const GeneratorState* seed) const {
      state.copy(seed);
    }

    /**
     * Returns current state, i.e. the one which would be used to generate
     * the next output.
     *
     * @param seed Current state returned here
     */
    void getSeed(GeneratorState* seed) const {
      seed->copy(&state);
    }

    bool canGenUInt() const {
      return true;
    }

    /**
     * Returns PRN deviate, uniform over the domain of unsigned int, which
     * is 0,...,2^32-1.
     *
     * @return See above
     */
    unsigned int getUInt() const {
      uint si=state.stateI(),sj=state.stateJ(),sk=state.stateK();

      sj^=(sj<<17); state.stateJ(sj^=(sj>>15));
      sk=(sk^(sk<<18))&0x7FFFFFFF;
      state.stateK(sk^=(sk>>13));
      state.stateI(si=69069*si+23606797);
      uint comb=si+sj+sk;

      return comb;
    }
  };
//ENDNS

#endif
