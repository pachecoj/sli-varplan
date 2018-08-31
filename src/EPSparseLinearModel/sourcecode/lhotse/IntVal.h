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
 * Desc.:  Header abstract class IntVal
 * ------------------------------------------------------------------- */

#ifndef INTVAL_H
#define INTVAL_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

/**
 * Abstract superclass of 'Interval'. Required by polymorphic code, in
 * which the type T of 'Interval' is not known.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
class IntVal
{
public:
  /* Constants
     NOTE: Must lie within 0,...,15! */
  static const int ivOpen  =0;
  static const int ivClosed=1;
  static const int ivInf   =2;
  static const int ivLast  =2;

public:
  // Public methods

  virtual ~IntVal() {}
};

#endif
