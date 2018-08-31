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
 * Module: rando
 * Desc.:  Standard exceptions
 * ------------------------------------------------------------------- */

#ifndef EX_RANDO_H
#define EX_RANDO_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/global.h"

//BEGINNS(rando)
  /**
   * No random generator supplied when calling for nonuniform deviate
   */
  class NoGeneratorException : public StandardException {
  public:
    NoGeneratorException(const char* mess=0) :
      StandardException("NoGeneratorException",mess) {}
  };
//ENDNS

#endif
