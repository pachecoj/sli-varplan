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
 * Module: matif
 * Desc.:  Standard exceptions
 * ------------------------------------------------------------------- */

#ifndef EX_MATIF_H
#define EX_MATIF_H

//BEGINNS(matif)
  /**
   * Field name does not exist (see 'MatIF_baseclass::getMember')
   */
  class FieldNotExistException : public StandardException {
  public:
    FieldNotExistException(const char* mess=0,const char* file=0,int line=0) :
      StandardException("FieldNotExistException",mess,file,line) {}
  };
//ENDNS

#endif
