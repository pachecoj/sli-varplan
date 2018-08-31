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
 * Desc.:  Header abstract class MatDefMembers
 * ------------------------------------------------------------------- */

#ifndef MATDEFMEMBERS_H
#define MATDEFMEMBERS_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/matrix/default.h"
#include "lhotse/NumberFormats.h"
#include "lhotse/TypeCode.h"

//BEGINNS(matrix)
  /*
   * Template function used to determine default fill value for type T.
   * This requires template specialization for the elementary types.
   */
  template<class T> void MatDefMembers_setDefFill(T& a)
  {
    // Default implementation, may not work for concrete T!
    T def;
    a=def;
  }

  template<> void MatDefMembers_setDefFill(char& a);
  template<> void MatDefMembers_setDefFill(uchar& a);
  template<> void MatDefMembers_setDefFill(int& a);
  template<> void MatDefMembers_setDefFill(uint& a);
  template<> void MatDefMembers_setDefFill(float& a);
  template<> void MatDefMembers_setDefFill(double& a);
  template<> void MatDefMembers_setDefFill(bool& a);
  template<> void MatDefMembers_setDefFill(long& a);
  template<> void MatDefMembers_setDefFill(ulong& a);

  /**
   * Defines members/methods for default values in matrix/vector classes.
   * - def. fill value 'defFill'
   * - variables for buffer increase strategy: 'incrBuffOffset',
   *   'incrBuffDouble'
   * <p>
   * Default constructors must call 'setDefValues' to init. the members
   * to def. init. values. Copy constructors should call 'copyDefValues'
   * to copy the members.
   * <p>
   * Default fill value:
   * The member 'defFill' is used as fill value in methods like
   * 'fill','expand', etc. It is initialised using
   * 'MatDefMembers_setDefFill', which is implemented using template
   * specialization.
   * This seems to work properly now, but we had some problems before. To
   * be sure, we overwrite 'setDefValues' in subclasses where T has a
   * definite value, calling 'setDefFillValue'. This ensures 'defFill' is
   * OK for sure in these cases.
   * <p>
   * Virtual methods:
   * These have to be implemented by all subclasses.
   * - setDefValues:
   *   Set def. values here to appr. values. Most important is the def.
   *   fill value which has to be set by all nonvirtual subclasses
   * - isValidElement:
   *   Returns true iff argument is valid as element of this matrix/vector.
   *   This method is called by 'setDefFillValue' as well.
   * <p>
   * Buffer increase strategy:
   * NOTE: These have to be implemented by subclasses, we only manage the
   * members here!
   * If 'size_n' is the old buffer size, 'num' the number of req. entries, the
   * strategies are:
   * - if 'incrBuffDouble'==true, the new buffer size is max(2*size_n,num),
   *   the size is at least doubled
   * - if 'incrBuffDouble'==false, the new buffer size is num+incrBuffOffset,
   *   allocating incrBuffOffset elements more than required. Otherwise,
   *   'incrBuffOffset' is ignored
   * <p>
   * NOTE: Neither 'defFill' nor the 'incrBuffXXX' members are static
   * variables, each object has its own copy (if they were static, all
   * subclasses would share the values).
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  template<class T> class MatDefMembers
  {
  protected:

    T defFill;                // Default fill value
    int incrBuffOffset;       // See header comment
    bool incrBuffDouble;

    /**
     * To be called by default constructors. The def. implementation
     * init. the 'incrBuffXXX' members to a conservative setting, and
     * uses 'MatDefMembers_setDefFill' to set 'defFill'.
     */
    virtual void setDefValues() {
      // Def. buffer increase strategy: double the size
      incrBuffOffset=10; // not used
      incrBuffDouble=true;
      // Def. fill value
      MatDefMembers_setDefFill(defFill);
    }

    /**
     * To be called by copy constructors.
     * Copies default value members from 'arg' to this object.
     *
     * @param arg Source object
     */
    void copyDefValues(const MatDefMembers<T>& arg) {
      defFill=arg.defFill; incrBuffOffset=arg.incrBuffOffset;
      incrBuffDouble=arg.incrBuffDouble;
    }

  public:
    // Public methods

    /**
     * @param a New def. fill value
     */
    void setDefFillValue(const T& a) {
      if (!isValidElement(a))
	throw InvalidParameterException(EXCEPT_MSG("Invalid element"));
      defFill=a;
    }

    /**
     * Is 'elem' valid as element of this matrix.
     * NOTE: Has to be implemented by subclasses!
     *
     * @param elem Element for matrix
     * @return     Is 'elem' valid as an element in this matrix?
     */
    virtual bool isValidElement(const T& elem) const = 0;

    /**
     * @return Def. fill value
     */
    const T& getDefFillValue() const {
      return defFill;
    }

    /**
     * @param New value for 'incrBuffOffset'
     */
    void setIncrBuffOffset(int off) {
      if (off<0) throw InvalidParameterException("'off' must be non-neg.");
      incrBuffOffset=off;
    }

    int getIncrBuffOffset() const {
      return incrBuffOffset;
    }

    void setIncrBuffDouble(bool dbl) {
      incrBuffDouble=dbl;
    }

    bool getIncrBuffDouble() const {
      return incrBuffDouble;
    }

    /**
     * Returns new buffer size, given that 'currSz' is the current size and
     * 'num' is the number of req. entries. Dep. on 'incrBuffXXX' members.
     *
     * @param num    Number of req. entries
     * @param currSz Current buffer size
     * @return       New buffer size (>= 'currSz')
     */
    int getNewBuffSize(int num,int currSz) const {
      if (num<=currSz) return currSz; // buffer large enough
      else if (incrBuffDouble) return std::max(2*currSz,num);
      else return num+incrBuffOffset;
    }

    /**
     * Loads int value from file stream 'is' and compares against type code
     * for T. If they are different, 'FileFormatException' is thrown.
     * <p>
     * DOWNW. COMPAT.: Type codes in files written by earlier LHOTSE versions
     * were wrong, in that they were 'typeOther' for all elem. types. If
     * 'dwComp'==true and the stored type code is 'typeOther' for elementary
     * T, a warning message is printed instead of throwing an exception.
     *
     * @param is     Input file stream
     * @param dwComp S.a. Def.: false
     */
    void checkTypeCode(ifstream& is,bool dwComp=false) {
      int i;
      int tCode=getTypeCode(defFill); // Type code for T
      NumberFormats<int>::load(is,&i,1,1,0,4);
      if (i!=tCode) {
	if (dwComp && i==TypeCodeConsts::typeOther)
	  cout << "WARNING: File has wrong type code! Try to read it with type code of class..." << endl;
	else
	  throw FileFormatException("File has wrong type code");
      }
    }
  };
//ENDNS

#endif
