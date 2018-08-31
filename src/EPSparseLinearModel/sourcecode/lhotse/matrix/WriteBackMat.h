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
 * Desc.:  Header abstract class WriteBackMat
 * ------------------------------------------------------------------- */

#ifndef WRITEBACKMAT_H
#define WRITEBACKMAT_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/matrix/default.h"
#include "lhotse/matrix/MatStrct.h"

//BEGINNS(matrix)
  /**
   * Helper class to break mutual dependency between 'BaseMatrix' and
   * 'BaseLinMat'. 'BaseMatrix' has to inherit from this class.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  template<class T> class WriteBackMat : public MatStrct
  {
  public:
    static bool debugFlag; // DEBUG: Used in 'FastUtils<T>'

    /**
     * Copy back buffer content from 'buff' into this matrix whose size is
     * the same as that of 'buff'.
     * NOTE: 'buff' is not only flat, but contiguous, i.e. the striding
     * constant for the underlying matrix is the same as the number of rows.
     * 'strpatt' is the structure pattern for the write-back (see
     * 'strctXXX' constants and BaseLinMat header comment).
     *
     * @param buff    S.a.
     * @param strpatt S.a.
     */
    virtual void writeBackBuff(const ArrayHandle<T>& buff,
			       unsigned char strpatt) = 0;
  };

  template<class T> bool WriteBackMat<T>::debugFlag(false);
//ENDNS

#endif
