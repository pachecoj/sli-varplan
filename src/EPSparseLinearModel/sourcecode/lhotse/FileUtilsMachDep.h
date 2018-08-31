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
 * Desc.:  Header class FileUtilsMachDep
 * ------------------------------------------------------------------- */

#ifndef FILEUTILSMACHDEP_H
#define FILEUTILSMACHDEP_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/global.h" // global header

/**
 * Static methods for fiel handling services.
 * As opposed to 'FileUtils', methods here are almost certainly architecture
 * dependent. In the moment, they are implemented for the GNU C standard
 * library only.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
class FileUtilsMachDep
{
public:
  // Public static methods

  /**
   * Given a list of filenames in 'fnameLst', we check for each whether
   * the file exists, and if so, which is the last modification time.
   * We return the index of the last recently modified file, or -1 if none
   * of the files in 'fnameLst' exist.
   *
   * @param fnameLst List of filenames
   * @return         S.a.
   */
  static int lastModified(const ArrayHandle<string>& fnameLst);

  /**
   * Removes file with name 'fname'.
   *
   * @param fname Name of file to remove
   * @return      0: success, -1: failure
   */
  static int removeFile(const char* fname);
};

#endif
