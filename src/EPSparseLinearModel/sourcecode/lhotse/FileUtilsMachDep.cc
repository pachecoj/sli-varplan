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
 * Desc.:  Definition class FileUtilsMachDep
 * ------------------------------------------------------------------- */

#include "lhotse/FileUtilsMachDep.h"

#include <sys/stat.h>
#include <unistd.h>

// Public static methods

int FileUtilsMachDep::lastModified(const ArrayHandle<string>& fnameLst)
{
  int i,n=fnameLst.size(),ret=-1;
  time_t max_mtime;

  for (i=0; i<n; i++) {
    const char* fname=fnameLst[i].c_str();
    struct stat attr;
    if (stat(fname,&attr)==-1) continue;
    if (ret==-1 || attr.st_mtime>max_mtime) {
      max_mtime=attr.st_mtime; ret=i;
    }
  }

  return ret;
}

int FileUtilsMachDep::removeFile(const char* fname)
{
  return unlink(fname);
}
