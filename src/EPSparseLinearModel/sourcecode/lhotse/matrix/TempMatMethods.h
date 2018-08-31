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
 * Module: matrix
 * Desc.:  Header abstract class TempMatMethods
 * ------------------------------------------------------------------- */

#ifndef TEMPMATMETHODS_H
#define TEMPMATMETHODS_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/matrix/default.h"

//BEGINNS(matrix)
  /**
   * Helper class, used in the 'TempMatMethods::assignVirtual' variant
   * where the argument is a scalar. The problem is that here, we cannot
   * know the element type of the vector/matrix class. Specializations
   * of this class must contain a member of the element type.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class TempMatMethods_AssignElement_Generic
  {
  public:
    virtual ~TempMatMethods_AssignElement_Generic() {}
  };

  /**
   * 'BaseMatrix', 'BaseVector' inherit from this base class. It is required to
   * allow 'TempXXX' classes to properly destroy 'XXX' objects without knowing
   * them directly. The virtual destructor here ensures that the correct
   * subclass destructor is called.
   * <p>
   * We also define some further abstract methods here, which are used
   * in 'TempXXX' to call 'XXX' methods without actually knowing 'XXX':
   * - assignment: 'assignVirtual' (matrix and scalar variant)
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class TempMatMethods
  {
  public:
    //bool debugFlag;

    /**
     * Constructor
     */
    TempMatMethods() {
      //debugFlag=false;
    }

    /**
     * Destructor
     */
    virtual ~TempMatMethods() {}

    /**
     * Has to be implemented in subclasses. The implementation has to
     * check whether 'arg' has the correct type (if not, this is an
     * internal error: throw 'InternalException'). If so, assign this
     * object to 'arg'.
     * If this object is a vector, 'arg' must be a vector as well.
     * If this object is a matrix, then 'arg' can always be a matrix,
     * but for special subclasses, it may also be a vector, which should
     * then be treated as a column vector.
     * <p>
     * NOTE: The implementation must do exactly the same as the corr.
     * 'operator='!
     *
     * @param arg Source argument
     */
    virtual void assignVirtual(const TempMatMethods* arg) = 0;

    /**
     * Has to be implemented in subclasses. The implementation has to
     * check whether 'arg' has the correct type (if not, this is an
     * internal error: throw 'InternalException'). If so, 'arg' contains
     * a member of the matrix/vector element type. The matrix/vector
     * should be filled with this element.
     * <p>
     * NOTE: The implementation must do exactly the same as the corr.
     * 'operator='!
     *
     * @param arg Source argument
     */
    virtual void
    assignVirtual(const TempMatMethods_AssignElement_Generic* arg) = 0;
  };
//ENDNS

#endif
