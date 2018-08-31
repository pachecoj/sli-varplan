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
 * Desc.:  Definition of class MatDefMembers
 * ------------------------------------------------------------------- */

#include "lhotse/matrix/MatDefMembers.h"

//BEGINNS(matrix)

  template<> void MatDefMembers_setDefFill(char& a)
  {
    a=0;
  }

  template<> void MatDefMembers_setDefFill(uchar& a)
  {
    a=0;
  }

  template<> void MatDefMembers_setDefFill(int& a)
  {
    a=0;
  }

  template<> void MatDefMembers_setDefFill(uint& a)
  {
    a=0;
  }

  template<> void MatDefMembers_setDefFill(float& a)
  {
    a=0.0;
  }

  template<> void MatDefMembers_setDefFill(double& a)
  {
    a=0.0;
  }

  template<> void MatDefMembers_setDefFill(bool& a)
  {
    a=false;
  }

  template<> void MatDefMembers_setDefFill(long& a)
  {
    a=0;
  }

  template<> void MatDefMembers_setDefFill(ulong& a)
  {
    a=0;
  }
//ENDNS
