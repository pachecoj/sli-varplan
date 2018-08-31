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
 * Desc.:  Header class MachDep
 * ------------------------------------------------------------------- */

/*
 * TODO:
 * This is bad and incomplete!
 * - Fix range of systems on which LHOTSE should compile
 * - Use system detection in the 'configure' script
 * - Do proper implementation here, using variables passed by 'configure'
 *   (through config.h)
 */

#ifndef MACHDEP_H
#define MACHDEP_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/global.h"   // global header

#ifdef HAVE_LIBGSL
#include <gsl/gsl_math.h> // GSL (general)
#endif

extern "C" {
#include <unistd.h>
#include <errno.h>
#ifdef HAVE_OS_SUN_SOLARIS
#include <ieeefp.h>
#endif
#ifdef HAVE_OS_ALPHA_OSF
#include <fp_class.h>
#include <signal.h>
#endif
}

// Functions

#ifdef HAVE_OS_ALPHA_OSF
void signal_handler(int); // req. for 'setupWorkarounds'
#endif

/**
 * Some utilities which may be machine-dependent are collected here.
 * However, there might be more more machine dependencies lurking somewhere.
 * For example: The file utilities are still machine-dependent: Files written
 * under one system cannot be read under another.
 * NOTE: Some things included here may be "standard" in C. But then, you never
 * know with C...
 * <p>
 * For SPARC SUN SOLARIS GCC 2.7/EGCS:
 * - values.h: Contains numerical limits
 * - ieeefp.h: Contains macros to test for NaN and infinity (for floating
 *   point numbers)
 * - unistd.h: BSD stuff(?). Need it for 'sleep'
 * <p>
 * For DEC ALPHA OSF1:
 * - unistd.h: Same as for Sun Solaris
 * - math.h,fp_class.h: The IEEE standard for floating point numbers
 * <p>
 * Conversion from string to numerical types:
 * This should really be done using a stringstream object and the >>
 * operator, but somehow this is broken (GCC 2.96, Linux i486).
 * ==> CHECK OUT!
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
class MachDep
{
public:
  // Public static methods

  /**
   * For some systems, one might have to do certain things to get stuff
   * running. Call this method at the beginning of the program, e.g. in
   * the LHOTSE initialization function.
   * <p>
   * This is annoying, but seems to be necessary to get things done!
   */
  static void setupWorkarounds() {
#ifdef HAVE_OS_ALPHA_OSF
    // Floating exceptions are thrown at RANDOM times, killing your
    // program. We deactivate the exception handler to avoid this.
    // True floating exceptions, however, are also deactivated by this
    // measure
    cout << "DEC OSF workaround: Deactivating floating exceptions." << endl;
    signal(SIGFPE, signal_handler); // install dummy signal handler
#endif
  }

  /**
   * @param val Floating point value
   * @return    Is 'val' == NaN (not a number)?
   */
  static bool isNaN(double val) {
#ifdef HAVE_LIBGSL
    return (bool) gsl_isnan(val);
#elif defined(HAVE_OS_SUN_SOLARIS)
    // Sun Solaris
    return (bool) isnand(val);
#else
    // Alpha OSF, Linux
    return (bool) isnan(val);
#endif
  }

  /**
   * @param val Floating point value
   * @return    Is 'val' == +infty?
   */
  static bool isPlusInf(double val) {
#ifdef HAVE_LIBGSL
    return (gsl_isinf(val)==+1);
#elif defined(HAVE_OS_SUN_SOLARIS)
    return (fpclass(val)==FP_PINF);
#elif defined(HAVE_OS_ALPHA_OSF)
    return (fp_class(val)==FP_POS_INF);
#else
    return (isinf(val)==1);
#endif
  }

  /**
   * @param val Floating point value
   * @return    Is 'val' == -infty?
   */
  static bool isMinusInf(double val) {
#ifdef HAVE_LIBGSL
    return (gsl_isinf(val)==-1);
#elif defined(HAVE_OS_SUN_SOLARIS)
    return (fpclass(val)==FP_NINF);
#elif defined(HAVE_OS_ALPHA_OSF)
    return (fp_class(val)==FP_NEG_INF);
#else
    return (isinf(val)==-1);
#endif
  }

  /**
   * @param val Floating point value
   * @return    Is 'val' NaN or +- infinity?
   */
  static bool isUndef(double val) {
#ifdef HAVE_LIBGSL
    return (gsl_finite(val)==0);
#else
    return isNaN(val) || isPlusInf(val) || isMinusInf(val);
#endif
  }

  /**
   * Sleep for 'sec' seconds.
   * ATTENTION: Does not work for Windows!
   *
   * @param sec
   */
  static void sleepSecs(uint sec) {
#ifndef HAVE_OS_WINDOWS
    // For Sun Solaris, Alpha OSF, Linux
    sleep(sec); // def. in unistd.h
#else
    throw NotImplemException(EXCEPT_MSG(""));
#endif
  }

  /**
   * Sleep for 'msec' microseconds. 'msec' must be smaller than 1e+6
   * <p>
   * ATTENTION: Does not work for ALPHA OSF, Windows!
   *
   * @param msec
   */
  static void sleepMSecs(uint msec) {
#if defined(HAVE_OS_LINUX) || defined(HAVE_OS_SUN_SOLARIS)
    // For Sun Solaris, Linux
    usleep(msec); // def. in unistd.h
#else
    throw NotImplemException(EXCEPT_MSG(""));
#endif
  }

  /**
   * Convert string into double value. In case of a conversion error, a
   * 'ParseException' is thrown.
   *
   * @param str String representation of double number
   * @return    Double value
   */
  static double stringToDouble(const char* str) {
    errno=0; // reset
    double val=atof(str);
    if (errno!=0) throw ParseException(EXCEPT_MSG(""));

    return val;
  }

  /**
   * Convert string into int value. In case of a conversion error, a
   * 'ParseException' is thrown.
   *
   * @param str String representation of int number
   * @return    Int value
   */
  static int stringToInt(const char* str) {
    errno=0; // reset
    int val=atoi(str);
    if (errno!=0) throw ParseException(EXCEPT_MSG(""));

    return val;
  }

  /**
   * Convert string into long (int) value. In case of a conversion error, a
   * 'ParseException' is thrown.
   *
   * @param str String representation of double number
   * @return    Long (int) value
   */
  static long stringToLong(const char* str) {
    errno=0; // reset
    long val=atol(str);
    if (errno!=0) throw ParseException(EXCEPT_MSG(""));

    return val;
  }
};

#endif
