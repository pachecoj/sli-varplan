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
 * Desc.:  Definition class DebugVars
 * ------------------------------------------------------------------- */

#include "lhotse/global.h"
#include "lhotse/DebugVars.h"
#include "lhotse/CommandParser.h"
#include "lhotse/MatlabDebug.h"

// Static members

string DebugVars::matDebBaseFName("");
bool DebugVars::doMessUpExc(false);
bool DebugVars::doPrintExc(false);
string DebugVars::messUpName("");

// Public static methods

void DebugVars::init(CommandParser& args)
{
  int dummy;

  matDebBaseFName="";
  try {
    args.getValue("debug-matlab-base-fname",CommandParser::typeString,
		  &matDebBaseFName);
  } catch (KeyNotFoundException ex) {}
  try {
    args.getValue("debug-matlab-activate",CommandParser::typeInt,
		  &dummy);
    if (dummy==0)
      MatlabDebug::deactivate();
    else {
      MatlabDebug::activate(matDebBaseFName);
      cout << "LHOTSE debug: Matlab debug activated." << endl;
    }
  } catch (KeyNotFoundException ex) {}
  doMessUpExc=false; doPrintExc=false; messUpName="";
  try {
    args.getValue("debug-messup-exc-activate",CommandParser::typeBool,
		  &doMessUpExc);
    if (doMessUpExc)
      cout << "LHOTSE debug: Messing up exceptions activated" << endl;
  } catch (KeyNotFoundException ex) {}
  try {
    args.getValue("debug-messup-exc-name",CommandParser::typeString,
		  &messUpName);
  } catch (KeyNotFoundException ex) {}
  try {
    args.getValue("debug-print-exc-early",CommandParser::typeInt,
		  &dummy);
    doPrintExc=(dummy!=0);
    if (doPrintExc)
      cout << "LHOTSE debug: Printing exc. messages early activated" << endl;
  } catch (KeyNotFoundException ex) {}
}

void DebugVars::matlabDebugActivate()
{
  MatlabDebug::activate(matDebBaseFName);
}

void DebugVars::matlabDebugDeactivate()
{
  MatlabDebug::deactivate();
}
