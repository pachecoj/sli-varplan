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
 * Desc.:  Header class PrintToStdout
 * ------------------------------------------------------------------- */

#ifndef PRINTTOSTDOUT_H
#define PRINTTOSTDOUT_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/global.h"
#include "lhotse/ArrayHandle.h"
#include "lhotse/matrix/BaseVector.h"
#include "lhotse/matrix/BaseMatrix.h"
#ifdef MATLAB_MEX
#include "lhotse/matif/mex_for_cpp.h"
#endif

/**
 * Provides static method for printing text to stdout. Furthermore, the
 * << operator is implemented here.
 * <p>
 * Unifies message printing for normal and Matlab MEX code. Replaces the
 * global function 'printMsgStdout' (deprecated, but can still be used).
 * The global object 'mycout' is declared here, replacing a basic version of
 * 'cout' in code which is used for MEX files.
 * Advanced features of streams such as manipulators cannot be used. For
 * normal code, the << operator calls a corr. << operator function which must
 * be defined for the arg. type.
 * In MEX mode, only certain types are supported.
 * NOTE: Use "\n" instead of 'endl' with 'mycout'.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
class PrintToStdout
{
protected:
  // Members

  static char sbuff[50];

public:
  // Public methods

  /**
   * Static method. Prints string to stdout.
   * NOTE: Not terminated by CR/LF.
   *
   * @param str String to print
   */
  static void print(const char* str) {
#ifndef MATLAB_MEX
    cout << str;
#else
    mexPrintf(str);
#endif
  }

  /**
   * << operator. In normal (non-MEX) mode, a << global function must be
   * defined for the type of the argument. In MEX mode, this works for
   * specific argument types only.
   * NOTE: This is not identical to the usual stream operator. Namely,
   * control tags such as maniuplators cannot be used. Instead of 'endl',
   * use the string "\n".
   *
   * @param arg Argument
   */
#ifndef MATLAB_MEX
  template<class T> PrintToStdout& operator<<(const T& arg) {
    cout << arg;
    return *this;
  }
#else
  PrintToStdout& operator<<(int arg) {
    sprintf(sbuff,"%d",arg); print(sbuff);
    return *this;
  }
  
  PrintToStdout& operator<<(uint arg) {
    sprintf(sbuff,"%d",arg); print(sbuff);
    return *this;
  }
  
  PrintToStdout& operator<<(double arg) {
    sprintf(sbuff,"%f",arg); print(sbuff);
    return *this;
  }
  
  PrintToStdout& operator<<(char arg) {
    sprintf(sbuff,"%c",arg); print(sbuff);
    return *this;
  }
  
  PrintToStdout& operator<<(uchar arg) {
    sprintf(sbuff,"%c",arg); print(sbuff);
    return *this;
  }
  
  PrintToStdout& operator<<(long arg) {
    sprintf(sbuff,"%ld",arg); print(sbuff);
    return *this;
  }
  
  PrintToStdout& operator<<(ulong arg) {
    sprintf(sbuff,"%ld",arg); print(sbuff);
    return *this;
  }
  
  PrintToStdout& operator<<(bool arg) {
    sprintf(sbuff,"%d",(int) arg); print(sbuff);
    return *this;
  }
  
  PrintToStdout& operator<<(const char* arg) {
    print(arg);
    return *this;
  }

  template<class T> PrintToStdout& operator<<(const BaseVector<T>& arg) {
    arg.matlabPrint();
    return *this;
  }

  template<class T> PrintToStdout& operator<<(const BaseMatrix<T>& arg) {
    arg.matlabPrint();
    return *this;
  }
#endif
};

// Global object

extern PrintToStdout mycout;

#endif
