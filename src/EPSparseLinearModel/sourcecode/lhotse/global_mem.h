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
 * Desc.:  Declarations for new/delete overload
 * ------------------------------------------------------------------- */

#ifndef GLOBALMEM_H
#define GLOBALMEM_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#if defined(USE_OWN_MEMMAN) || (defined(MATLAB_MEX) && defined(USE_MATLAB_MM))

extern void* operator new(size_t n);
extern void operator delete(void* ptr);
extern void* operator new[](size_t n);
extern void operator delete[](void* ptr);

#endif

#endif
