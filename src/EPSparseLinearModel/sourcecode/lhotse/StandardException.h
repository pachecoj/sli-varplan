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
 * Desc.:  Header class StandardException
 * ------------------------------------------------------------------- */

#ifndef STANDARDEXCEPTION_H
#define STANDARDEXCEPTION_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

/**
 * Base class of the exceptions defined by LHOTSE. Maintains an error
 * message.
 * <p>
 * Supporting 'DebugVars':
 * The 'DebugVars' class maintains the debug variables for LHOTSE. One of
 * the facilities is to mess up LHOTSE exceptions: some code in the
 * constructor of 'StandardException' is executed leading to a segmentation
 * fault. In order to avoid a cyclic inclusion, we have to define the
 * constructor in the code file.
 * <p>
 * We also want selected subclasses of 'StandardException' NOT to cause a
 * seg. fault (e.g. 'KeyNotFoundException'). To this end, we define a
 * second non-standard constructor which never causes an exception.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
class StandardException
{
protected:
  // Variables

  string message; // error message
  static string name;

public:
  // Constructors

  /**
   * Constructor.
   * Can pass filename and line number in 'file','line'. Use EXCEPT_MSG
   * macro:
   *   throw ExcName(EXCEPT_MSG("exc. message"));
   *
   * @param nam  Class name
   * @param mess Message string
   * @param file File name of throw
   * @param line Line number of throw
   */
  StandardException(const char* nam="StandardException",const char* mess=0,
		    const char* file=0,int line=0);

  ~StandardException() {}

  // Public methods

  /**
   * @return Exc. message
   */
  const char* msg() const {
    return message.c_str();
  }

  /**
   * @return Class name
   */
  const char* getName() const {
    return name.c_str();
  }
};

#endif
