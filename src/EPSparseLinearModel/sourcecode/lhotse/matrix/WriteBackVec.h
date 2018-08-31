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
 * Desc.:  Header abstract class WriteBackVec
 * ------------------------------------------------------------------- */

#ifndef WRITEBACKVEC_H
#define WRITEBACKVEC_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/matrix/default.h"

//BEGINNS(matrix)
  /**
   * Helper class to break mutual dependency between 'BaseVector' and
   * 'BaseLinVec'. 'BaseVector' has to inherit from this class.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  template<class T> class WriteBackVec
  {
  public:
    /**
     * Copy back buffer content from 'buff' into this vector whose size is
     * the same as that of 'buff'.
     *
     * @param buff S.a.
     */
    virtual void writeBackBuff(const ArrayHandle<T>& buff) = 0;
  };
//ENDNS

#endif
