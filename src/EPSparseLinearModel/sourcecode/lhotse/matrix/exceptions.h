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
 * Module: matrix
 * Desc.:  Standard exceptions
 * ------------------------------------------------------------------- */

#ifndef EX_MATRIX_H
#define EX_MATRIX_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

//BEGINNS(matrix)
  /**
   * Wrong dimensions passed as arguments
   */
  class WrongDimensionException : public StandardException {
  public:
    WrongDimensionException(const char* mess=0,const char* file=0,int line=0) :
      StandardException("WrongDimensionException",mess,file,line) {}
  };

  /**
   * Mismatch between the dimensions of argument matrices
   */
  class DimMismatchException : public StandardException {
  public:
    DimMismatchException(const char* mess=0,const char* file=0,int line=0) :
      StandardException("DimMismatchException",mess,file,line) {}
  };

  /**
   * Two or more identical arguments (or arguments ident. to the target
   * object) have been passed to an operation who
   * permits that (ex.: multiplication)
   */
  class EqualArgsException : public StandardException {
  public:
    EqualArgsException(const char* mess=0,const char* file=0,int line=0) :
      StandardException("EqualArgsException",mess,file,line) {}
  };

  /**
   * Error while creating matrix object
   */
  class MatrixCreationException : public StandardException {
  public:
    MatrixCreationException(const char* mess=0,const char* file=0,int line=0) :
      StandardException("MatrixCreationException",mess,file,line) {}
  };

  /**
   * Invalid use of a mask vector or a mask matrix
   */
  class MaskObjectException : public StandardException {
  public:
    MaskObjectException(const char* mess=0,const char* file=0,int line=0) :
      StandardException("MaskObjectException",mess,file,line) {}
  };

  /**
   * Method not implemented for indexed mask objects
   */
  class MaskNotImplException : public StandardException {
  public:
    MaskNotImplException(const char* mess=0,const char* file=0,int line=0) :
      StandardException("MaskNotImplException",mess,file,line) {}
  };

  /**
   * Elements have wrong type, violating object type constraints
   */
  class WrongTypeException : public StandardException {
  public:
    WrongTypeException(const char* mess=0,const char* file=0,int line=0) :
      StandardException("WrongTypeException",mess,file,line) {}
  };
//ENDNS

#endif
