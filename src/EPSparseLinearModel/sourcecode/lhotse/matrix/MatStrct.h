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
 * Desc.:  Header abstract class MatStrct
 * ------------------------------------------------------------------- */

#ifndef MATSTRCT_H
#define MATSTRCT_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/matrix/default.h"

//BEGINNS(matrix)
  /**
   * Abstract class, contains structure pattern constants. Could put this
   * into 'WriteBackMat', but the present class is not templatized.
   * See 'StMatrix' for application.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class MatStrct
  {
  public:
    /*
     * Structure pattern (within rectangular matrix):
     * normal:  Normal rectangular
     * upper:   Upper triangle and diagonal
     * lower:   Lower triangle and diagonal
     * uppNDg:  Upper triangle without diagonal
     * lowNDg:  Lower triangle without diagonal
     * NOTE: A symmetric matrix can be stored as 'upper' or 'lower'.
     * NOTE: The 'strctXXX' variants are there for downw. compat.
     * Do not use them anymore!
     */
    static const unsigned char normal=0;
    static const unsigned char upper =1;
    static const unsigned char lower =2;
    static const unsigned char uppNDg=3;
    static const unsigned char lowNDg=4;
    static const unsigned char last  =4;

    // DOWNW. COMPATIBILITY (DO NOT USE):
    static const unsigned char strctNormal=normal;
    static const unsigned char strctUpper =upper;
    static const unsigned char strctLower =lower;
    static const unsigned char strctUppNDg=uppNDg;
    static const unsigned char strctLowNDg=lowNDg;
    static const unsigned char strctLast  =last;
  };
//ENDNS

#endif
