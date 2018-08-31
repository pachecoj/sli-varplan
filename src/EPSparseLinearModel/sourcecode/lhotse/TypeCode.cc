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
 * Module: GLOBAL
 * Desc.:  Definition specializations of function getTypeCode
 * ------------------------------------------------------------------- */

#include "lhotse/global.h"
#include "lhotse/TypeCode.h"

const int TypeCodeConsts::typeOther ;
const int TypeCodeConsts::typeChar  ;
const int TypeCodeConsts::typeUChar ;
const int TypeCodeConsts::typeInt   ;
const int TypeCodeConsts::typeUInt  ;
const int TypeCodeConsts::typeFloat ;
const int TypeCodeConsts::typeDouble;
const int TypeCodeConsts::typeBool  ;
const int TypeCodeConsts::typeLong  ;
const int TypeCodeConsts::typeULong ;
const int TypeCodeConsts::typeLast  ;

// Specialisations

template<>
int getTypeCode(const char& a)
{
  return TypeCodeConsts::typeChar;
}

template<>
int getTypeCode(const uchar& a)
{
  return TypeCodeConsts::typeUChar;
}

template<>
int getTypeCode(const int& a)
{
  return TypeCodeConsts::typeInt;
}

template<>
int getTypeCode(const uint& a)
{
  return TypeCodeConsts::typeUInt;
}

template<>
int getTypeCode(const float& a)
{
  return TypeCodeConsts::typeFloat;
}

template<>
int getTypeCode(const double& a)
{
  return TypeCodeConsts::typeDouble;
}

template<>
int getTypeCode(const bool& a)
{
  return TypeCodeConsts::typeBool;
}

template<>
int getTypeCode(const long& a)
{
  return TypeCodeConsts::typeLong;
}

template<>
int getTypeCode(const ulong& a)
{
  return TypeCodeConsts::typeULong;
}
