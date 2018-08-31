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
 * Desc.:  Header class TypeCodeConsts, function getTypeCode
 * ------------------------------------------------------------------- */

/*
 * Test whether template specialization works, as used here!
 */

#ifndef TYPECODE_H
#define TYPECODE_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

/**
 * Superclass of 'TypeCode', holds typecode constants.
 * <p>
 * NOTE: These codes are used in files ==> Do NOT change the code for a
 * type! New type codes can be added at any time.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
class TypeCodeConsts
{
public:
  // Constants

  static const int typeOther    =0;
  static const int typeChar     =1;
  static const int typeUChar    =2;
  static const int typeInt      =3;
  static const int typeUInt     =4;
  static const int typeFloat    =5;
  static const int typeDouble   =6;
  static const int typeBool     =7;
  static const int typeLong     =8;
  static const int typeULong    =9;
  static const int typeLast     =9;
};

/**
 * Template function which returns the typecode depending on its argument.
 * A ref. to a variable of type T has to be provided (the variable is not
 * accessed), to force the compiler to use the correct specialization.
 * <p>
 * ATTENTION (unresolved problem):
 * 'getTypeCode' does return 'typeOther' even for elementary types!
 * Somehow, template specialization does not seem to work properly in
 * our context.
 * ==> WEIRD: It works in smaller test programs!
 * As long as this is not resolved properly, specific workarounds have
 * to be used (see 'BaseMatrix', 'MatDefMembers::getElemTypeCode' for an
 * example).
 *
 * @param a Ref. to T variable
 * @return  Type code, TypeCodeConsts::typeOther for non-default type
 */
template<class T>
int getTypeCode(const T& a)
{
  return TypeCodeConsts::typeOther; // default implementation
}

template<> int getTypeCode(const char& a);
template<> int getTypeCode(const uchar& a);
template<> int getTypeCode(const int& a);
template<> int getTypeCode(const uint& a);
template<> int getTypeCode(const float& a);
template<> int getTypeCode(const double& a);
template<> int getTypeCode(const bool& a);
template<> int getTypeCode(const long& a);
template<> int getTypeCode(const ulong& a);

#endif
