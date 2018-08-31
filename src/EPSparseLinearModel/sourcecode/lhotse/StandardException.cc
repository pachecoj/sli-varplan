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
 * Desc.:  Definition class StandardException
 * ------------------------------------------------------------------- */

#include "lhotse/global.h"
#include "lhotse/StandardException.h"
#ifdef HAVE_DEBUG
#include "lhotse/DebugVars.h"
#endif

// Static members

string StandardException::name;

// Constructors

/*
 * Supporting the messing up LHOTSE exceptions debug facility (see
 * 'DebugVars'. If this is activated, we execute some code in the constructor
 * leading to a segmentation fault.
 */
StandardException::StandardException(const char* nam,const char* mess,
				     const char* file,int line) :
  message() {
  name=nam;
  if (mess==0 || strlen(mess)==0) {
    message=name+": unspecified";
  } else
    message=mess;
  if (file!=0) {
    message+="\nFile: "; message+=file;
    char buff[40];
    sprintf(buff," (line %d)",line);
    message+=buff;
  }
#ifdef HAVE_DEBUG
  if (DebugVars::doWePrintExcEarly()) {
    cout << "DEBUG: Exception created and thrown. Message:" << endl
	 << message << endl;
  }
  if (DebugVars::doWeMessUpExceptions(name)) {
    // Causes a segmentation fault
    int* breakIt=0;
    *breakIt=1234;
  }
#endif
}
