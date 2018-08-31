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
 * Desc.:  Header class CommandParser
 * ------------------------------------------------------------------- */

#ifndef COMMANDPARSER_H
#define COMMANDPARSER_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/global.h" // global header

/**
 * Simple parser to read the command file to control the main program. Upon
 * construction of the object, the control file is processed and a table of
 * keyword-value pairs is built. Values are stored as strings. Values can
 * then be requested by giving their keyword and the value type.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
class CommandParser
{
protected:
  // Types

  typedef MAP_CONSTITER(my_string,my_string) LstIter;

  // Members

  MAP_TYPE(my_string,my_string) keyValLst;

public:
  // Constants
  static const int maxLineLength=10000;  // Maximal line length
  static const int maxListEntries=100; // Max. number of entries per list par.

  static const int typeInt=0;         // Type constants
  static const int typeDouble=1;
  static const int typeString=2;
  static const int typeBool=3;
  static const int typeLong=4;
  static const int typeListInt=5;
  static const int typeListDouble=6;
  static const int typeLast=6;        // Must cont. highest type const.

  // Constructors

  /**
   * Default constructor. The structure of the command file is described in
   * the documentation.
   *
   * @param fname Name of command file
   * @returns     True in case of success. False if parse error
   */
  CommandParser(const my_string& fname);

  ~CommandParser() {}

  // Public methods

  /**
   * Returns value for a given keyword.
   * Type code      Type of 'value'
   * typeInt        int
   * typeDouble     double
   * typeString     my_string
   * typeBool       bool
   * typeLong       long
   * typeListInt    BaseVector<int>
   * typeListDouble BaseVector<double>
   * If 'ival' is given, it is an interval against which the corr. value
   * is checked, or each value in case of a vector type. 'ival' must be
   * of correct 'Interval<T>' type, it is ignored for 'typeString',
   * 'typeBool'. For a violation, a 'OutOfRangeException' is thrown.
   *
   * @param key   Keyword
   * @param type  Typecode (see constants)
   * @param value Pointer to variable where to store value
   * @param intv  S.a. Optional. Def.: 0
   * @exception   KeyNotFoundException,ParseException,OutOfRangeException
   */
  void getValue(const my_string& key,int type,void* value,
		const IntVal* ival=0) const;
};

#endif
