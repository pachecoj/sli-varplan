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
 * Desc.:  Header class ArgBlock
 * ------------------------------------------------------------------- */

#ifndef ARGBLOCK_H
#define ARGBLOCK_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/global.h" // global header
#include "lhotse/CommandParser.h"
#include "lhotse/matrix/predecl.h"

/**
 * Helper class for 'ArgBlock'. Description for one parameter of the
 * block (the descriptions are stored in a static list in 'ArgBlock').
 * A parameter is ident. by its position in the description list.
 * Its attributes are:
 * - name (string)
 * - value type: see 'CommandParser::typeXXX'
 * - def. value (maint. as void*)
 * - interval for range checking (optional, maint. as void*)
 * - presence check type (see 'ArgBlock::presXXX')
 * - dependent parameter (optional, given as its number)
 * See 'ArgBlock' header for details.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
class ArgBlockType
{
public:
  // Constants
  static const int presReqUncond  =0;
  static const int presOptUncond  =1;
  static const int presReqIfDep   =2;
  static const int presReqIfNotDep=3;

  // Members
  my_string name;
  int type;
  void* defVal;
  void* rngIv;
  int pcheck;
  int depPar;

  ArgBlockType(const my_string& nameP,int typeP,void* defValP,void* rngIvP=0,
	       int pcheckP=presOptUncond,int depParP=-1);

  ~ArgBlockType();

  /**
   * Checks whether value behind 'obj' falls into interval given by 'rngIv'.
   * Returns true if 'rngIv'==0 or if the type is bool or string.
   *
   * @param obj S.a.
   * @return     S.a.
   */
  bool checkIntVal(const void* obj) const;
};

/**
 * An argument blocks reads and maintains a set of parameters used to
 * control some entity like an algorithm. The parameter values are read
 * by processing a 'CommandParser' object, which in turn is constructed
 * from a control file.
 * The schema for an argument block is fixed. It consists of a list of
 * descriptions, one for each parameter. A parameter is ident. by its
 * name (string) or its number (pos. in the description list) and has
 * attributes 'ArgBlockType'.
 * <p>
 * Construction:
 * The descr. table has to be setup before an object can be created.
 * Once this is done, 'init' can be called to fill the value list from
 * a 'CommandParser' object. A subclass constructor has to:
 * - setup descr. table
 * - call 'init'
 * - optional: do add. consistency checks
 * NOTE: Declaring a virtual method for the descr. table setup leads to
 * trouble, because the virtual mech. does not work as long as the
 * object is not fully constructed!
 * NOTE: The descr. table is NOT maintained as static member, but rather
 * duplicated for each object. This is because it is defined here in
 * 'ArgBlock', but has to be different in different subclasses (a static
 * member of 'ArgBlock' would be shared by all subclasses!).
 *
 * In 'init, for each parameter there are two checks:
 * - presence check
 * - type and value check
 * The type check is a bit fragile(!) The type in the descr. is used to
 * read the value from 'pars'. If this causes a 'ParseException', there is
 * a type error. But type problems may remain undetected here!
 * If the descr. contains an interval in 'rngIv', we check whether the value
 * is in the interval. For vector types, each element is checked.
 * NOTE: The def. value for a param. need NOT lie in its interval!
 * The presence check depends on the 'pcheck' attribute:
 * - presReqUncond:   Is required, unconditional
 * - presOptUncond:   Need not be given, unconditional
 * - presReqIfDep:    Required if the dep. param. 'depPar' is given
 * - presReqIfNotDep: Required if the dep. param 'depPar' not given
 * NOTE: If 'depPar' is given, it has to be of "pure" 'pcheck':
 * 'presXXXUncond'.
 * If both checks are met, the param. value is read and stored in 'value'
 * if given, or the 'value' entry points to the 'defVal' field in the
 * description (def. value) if not given. The 'given' entry is true iff
 * the param. value is given.
 * <p>
 * Parameter names:
 * The basic name of each parameter is fixed in the description table. It
 * is used by 'getXXX', 'setXXX' to identify the param. In contrast, the
 * parameter may have a different full name in the 'CommandParser' object
 * (control file). 'getFullName' maps basic to full names. In the moment,
 * this is done by simply appending 'prefName' and 'suffName' at beginning
 * and end. Full names are used in the constructor AND in the messages of
 * exceptions.
 * <p>
 * Additional checks:
 * These should be implemented at the end of the subclass constructor,
 * after the 'ArgBlock' constructor has returned.
 * <p>
 * Adding values after the control file has been read:
 * - include parameter in the description, using 'presOptUncond' for
 *   'pcheck'
 * - use 'set' for set the value. This can also be used to modify values.
 *   The interval condition is checked in 'set'
 *   NOTE: It is not possible to un-set a value!
 * <p>
 * Polymorphism:
 * Most methods here are not virtual, because they may be required during
 * construction/destruction. Polymorphic use should be avoided!
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
class ArgBlock
{
  friend class ArgBlockType;

protected:
  // Types

  typedef MAP_CONSTITER(my_string,int) LstIter;

  // Members

  MAP_TYPE(my_string,int) nameToNum;   // Maps names to numbers
  ArrayPtrHandle<ArgBlockType> descr;  // List of descriptions
  ArrayHandle<void*> value;            // List of values
  ArrayHandle<bool> given;             // Has entry been given?
  my_string prefName,suffName;         // See 'getFullName'

public:
  // Methods

  /**
   * The subclass constructor has to call this method before doing
   * anything else!
   * <p>
   * 'init' calls 'setupDescription' before anything else.
   * 'pars' is the 'CommandParser' object from which the values are read.
   * 'prefNam' and 'suffNam' are optional. If given, 'prefNam' is appended
   * as prefix, 'suffNam' as suffix to each (basic) attribute name to generate
   * a full name (see 'getFullName'). The full name is used only when reading
   * the par. values here from 'pars', and in error messages. The access
   * methods 'getXXX' use the basic names!
   * If any of the checks fail, a descriptive error message is generated and
   * a 'ArgBlockException' is thrown.
   * <p>
   * If 'block' is given, it ref. to another object of the same type as
   * this (concrete class not checked here!). In this case, for every value
   * not present in 'pars', we copy the corr. value from 'block'.
   *
   * @param pars    S.a.
   * @param prefNam S.a. Def.: ""
   * @param suffNam S.a. Def.: ""
   * @param block   S.a. Def.: 0
   */
  void init(const CommandParser& pars,const my_string& prefNam="",
	    const my_string& suffNam="",const ArgBlock* block=0);

  virtual ~ArgBlock();

  // Public methods

  /*
   * The 'getXXX' methods return a ref. to the value for a parameter
   * of a given name. They throw an 'ArgBlockException' if the
   * parameter cannot be ident. or the return type is wrong.
   */

  int getInt(const my_string& name) const;
  double getDouble(const my_string& name) const;
  const my_string& getString(const my_string& name) const;
  bool getBool(const my_string& name) const;
  long getLong(const my_string& name) const;
  const BaseVector<int>& getListInt(const my_string& name) const;
  const StVector& getListDouble(const my_string& name) const;

  /**
   * Returns true if the corr. value has been given or was supplied later
   * by 'set'.
   *
   * @param name Parameter name
   * @return     Has corr. value been given in the param. file?
   */
  bool wasGiven(const my_string& name) const {
    int num=findPar(name);
    return given[num];
  }

  /**
   * Allows to set the value of a parameter. The new value is pointed to by
   * 'valP' (a copy of the object is drawn).
   * It is not possible to un-set the value of a parameter,
   * so the "presence" consistency cannot be affected. 'val' has to fulfil
   * interval consistency if an interval is defined. The param. is identified
   * by its name.
   * NOTE: After the generic checks, but before setting the value, the method
   * calls 'checkNewValue' passing the value (a void*) and the number of the
   * param.
   *
   * @param name Param. name
   * @param valP S.a.
   */
  void set(const my_string& name,const void* valP);

  /**
   * NOTE: Can be overwritten by subclasses if there is a different
   * convention.
   *
   * @param num Parameter number
   * @return    Full name of param. (see header)
   */
  my_string getFullName(int num) const {
    my_string nm(prefName); nm+=descr[num]->name; nm+=suffName;
    return nm;
  }

  /**
   * Version for users.
   *
   * @param name Parameter name (basic)
   * @return     Full name of param. (see header)
   */
  my_string getFullName(const my_string& name) const {
    return getFullName(findPar(name));
  }

protected:
  // Internal methods

  /**
   * Implement in subclasses!
   * Called by 'setXXX' method before a parameter is assigned a new value.
   * A void* to the new value is in 'valP', its number passed in 'num'.
   * The method should throw a descriptive 'ArgBlockException' if there's
   * something wrong.
   *
   * @param valP S.a.
   * @param num  S.a.
   */
  virtual void checkNewValue(const void* valP,int num) const {}

  /**
   * Searches 'nameToNum' for param. 'name' and returns its number.
   * If not found, a 'ArgBlockException' is thrown. Otherwise, its
   * type is compared against 'typ'. If not identical, an
   * 'ArgBlockException' is thrown.
   * NOTE: If 'typ'==-1 (def.), the type comparison is not done.
   *
   * @param name S.a.
   * @param typ  S.a. Def.: -1
   * @return     Number
   */
  int findPar(const my_string& name,int typ=-1) const;

  /**
   * Does constructor's job for a single parameter with number 'num',
   * given the 'CommandParser' object 'pars'. When checking presence
   * for pars. with 'pcheck'=='presReqIfXXX', we assume that all
   * pars. without dependent par. ("pure" ones) have been read before.
   *
   * @param pars  S.a.
   * @param num   S.a.
   * @param block See constr.
   */
  void readValue(const CommandParser& pars,int num,const ArgBlock* block=0);

  /**
   * Generates object of type 'typ' on the heap and returns pointer
   *
   * @param typ S.a.
   * @return    Pointer to object
   */
  static void* genObj(int typ);

  /**
   * Generates copy of 'obj' on the heap, returns pointer
   *
   * @param obj Obj. to copy
   * @param typ Type of object
   * @return    S.a.
   */
  static void* copyObj(const void* obj,int typ);

  /**
   * Destroys object 'obj' (given by void*) of type 'typ'.
   *
   * @param obj S.a
   * @param typ S.a.
   */
  static void freeObj(void* obj,int typ);

  /**
   * @param num Param. number
   * @return    Is corr. 'pcheck' "pure" type?
   */
  bool isPure(int num) const {
    int pc=descr[num]->pcheck;
    return (pc==ArgBlockType::presReqUncond ||
	    pc==ArgBlockType::presOptUncond);
  }

  /**
   * Helper for 'setupDescription'. Checks whether 'pcheck' and 'depPar'
   * values in descr. table are consistent.
   */
  bool presConsistent() const;

  /**
   * Helper for 'checkExtConsistency'. If param. with number 'num' is not
   * given, an exception is thrown.
   *
   * @param num S.a.
   */
  void checkGiven(int num) const;
};

#endif
