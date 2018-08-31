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
 * Desc.:  Module predeclarations
 * ------------------------------------------------------------------- */

#ifndef DE_RANDO_H
#define DE_RANDO_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

//BEGINNS(rando)
  class NoGeneratorException;

  class GeneratorState;
  class Generator;
  class GenstateKISS;
  class GeneratorKISS;
//class RandoNumrecCode;
  class Random;
  class TestUtils; // DOWNW. COMPAT.!!
  class Distribution;
  class DistribSamp;
  class DistrNormal;
  class DistrGamma;
  class DistrPoisson;
  class DistrT;
  class DistrBeta;
  class DistrExp;
  class DistrUnif;
//ENDNS

#endif
