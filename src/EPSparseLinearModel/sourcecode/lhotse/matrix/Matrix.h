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
 * Desc.:  Header abstract class Matrix
 * ------------------------------------------------------------------- */

#ifndef MATRIX_H
#define MATRIX_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/matrix/default.h"

//BEGINNS(matrix)
  /**
   * Abstract base class of the hierarchy provided by the matrix module.
   * A Matrix is a two-dimensional field whose elements are of a common
   * type.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class Matrix
  {

  protected:
    // Variables

    int m;      // Number of rows
    int n;      // Number of columns

  public:
    // Constructors

    /**
     * Constructs matrix with 0 rows/columns, no entries
     */
    Matrix() : m(0),n(0) {}

    /**
     * Destructor
     */
    virtual ~Matrix() {}

    // Public methods

    /**
     * Returns number of rows.
     *
     * @return Number of rows
     */
    int rows() const {
      return m;
    }

    /**
     * Returns number of columns.
     *
     * @return Number of columns
     */
    int cols() const {
      return n;
    }
  };
//ENDNS

#endif
