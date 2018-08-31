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
 * Project source file
 * Module: GLOBAL
 * Desc.:  Default LHOTSE initialisation function
 * ------------------------------------------------------------------- */

#ifndef DEFINITFUNC_H
#define DEFINITFUNC_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/global.h"
#include "lhotse/FileUtils.h"
#include "lhotse/MachDep.h"
#include "lhotse/CommandParser.h"
#include "lhotse/rando/GenstateKISS.h"
#include "lhotse/rando/GeneratorKISS.h"
#include "lhotse/rando/Random.h"

USING(rando);

/**
 * This default implementation should be called by all LHOTSE main programs.
 * It does:
 * - initialize memory manager 'FixedMemManager' (if flag USE_OWN_MEMMAN set)
 * - initialise debug facilities (if flag HAVE_DEBUG is set)
 * - test number format of architecture (FileUtils)
 * - setup machine-dep. workarounds (MachDep)
 * - creates and initialises the default PRN generator:
 *   This is (in the moment) the KISS generator 'GeneratorKISS'. The seed
 *   is chosen using the system time if not provided, it is printed to
 *   stdout in any case
 *
 * @param initSeed Seed (init. state) for KISS generator. 0 => Seed computed
 *                 from output of 'time'
 * @param argsP    Optional. Command line arguments object
 */
int LHOTSE_init(const GeneratorState* initSeed=0,CommandParser* argsP=0)
{
#ifdef USE_OWN_MEMMAN
  // Initialize own memory manager iff 'own-mem-manager' is true
  bool ownMM=false;
  if (argsP!=0)
    try {
      argsP->getValue("own-mem-manager",CommandParser::typeBool,&ownMM);
    } catch (KeyNotFoundException ex) {}
  if (ownMM) {
    printMsgStdout("Initialize memory manager for fixed-size objects");
    FixedMemManager::init();
  }
#endif
#ifdef HAVE_DEBUG
  // LHOTSE debug facilities: initialize from control file
  if (argsP!=0) {
    printMsgStdout("LHOTSE debug facilities: Initialize from control file...");
    DebugVars::init(*argsP);
  }
#endif
  // Test number formats of architecture
  FileUtils::testFormats();
  // Install workarounds for machine
  MachDep::setupWorkarounds();
  // Create and initialise default PRNG
  const GenstateKISS* state=DYNCAST(const GenstateKISS,initSeed);
  Handle<GenstateKISS> hstate;
  if (initSeed!=0) {
    if (state==0)
      throw InvalidParameterException(EXCEPT_MSG("Need GenstateKISS"));
    Random::setDefGen(new GeneratorKISS(*state));
  } else {
    hstate.changeRep(new GenstateKISS());
    hstate->initFromTime();
    state=hstate;
    Random::setDefGen(new GeneratorKISS(*state));
  }
  state->print();

  return 1;
}

#endif
