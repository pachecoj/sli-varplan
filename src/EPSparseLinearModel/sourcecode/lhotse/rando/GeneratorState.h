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
 * Desc.:  Header abstract class GeneratorState
 * ------------------------------------------------------------------- */

#ifndef GENERATORSTATE_H
#define GENERATORSTATE_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/rando/default.h"

//BEGINNS(rando)
  /**
   * Abstract base class for PRN generator state variables. See 'Generator'
   * for details. The future sequence of outputs of a generator must be a
   * deterministic function of its current state.
   * <p>
   * TODO: Add mandatory 'save','load' methods and use a
   * 'GeneratorStateFactory' where all 'GeneratorState' subclasses are
   * registered s.t. they can be loaded generically.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class GeneratorState
  {
  public:
    /**
     * Draws a copy of state 'src' and stores it in this object.
     *
     * @param src Source state
     */
    virtual void copy(const GeneratorState* src) = 0;

    /**
     * Prints current state together with name of concrete
     * subclass. The information must be sufficient to initialize the
     * state. Use 'printMsgStdout' for the output.
     */
    virtual void print() const = 0;

    /**
     * Initializes state from string representation of the seed. The format
     * depends on the subclass.
     *
     * @param str String repres. of seed
     */
    virtual void initFromString(const char* str) = 0;

    /**
     * Writes string representation of current state into 'str'. The format
     * depends on the subclass, but must be the same as accepted by
     * 'initFromString'.
     *
     * @param str String repres. of seed ret. here
     */
    virtual void getString(string& str) const = 0;
  };
//ENDNS

#endif
