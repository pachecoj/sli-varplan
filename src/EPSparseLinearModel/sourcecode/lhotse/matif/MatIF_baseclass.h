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
 * Desc.:  Header class MatIF_baseclass
 * ------------------------------------------------------------------- */

#ifndef MATIF_BASECLASS_H
#define MATIF_BASECLASS_H

#include "lhotse/matif/default.h"
#include "lhotse/matif/MatlabMatrix.h"

//BEGINNS(matif)
  /**
   * Wrapper for class 'baseclass'
   * <p>
   * 'baseclass' is the root of the Matlab class hierarchy and defines all
   * virtual methods. The class is virtual itself. The following methods
   * have to be implemented by subclasses:
   * - create (static): construct object of this class and return pointer.
   *   NOTE: Virtual classes such as 'baseclass' do not implement 'create'.
   *   The MEX function returns with error if 'create' is called for them
   * - run: called with method name and args. Search for method name in local
   *   list. Call method if found, otherwise call superclass 'run'
   *   NOTE: Wanted to do this with pointer to methods, but the concept seems
   *   pretty broken (cannot store array of those as static member)
   * <p>
   * List of method names:
   * Stored in static 'numMeth', 'methNames'. These are setup iff
   * 'methInit'==true. The setup is done by the static method 'setup' which
   * also sets the 'methInit' flag. 'setup' does nothing if the flag is
   * already set.
   * NOTE: 'setup' must be implemented by every subclass which exports new
   * methods. The class constructor has to call 'setup' before anything
   * else. It also has to call the constructor of the superclass. It is
   * important that the 'MatIF_baseclass' constructor is called for every
   * construction of a subclass object.
   * <p>
   * Subobjects and ref. counting:
   * Every object has a ref. counter 'refCount' which is init. with 0 at
   * construction. 'refCount' must be increased whenever a ref. to this object
   * is created in another object (i.e. this object becomes a subobject). It
   * must be decreased whenever such a ref. is destroyed.
   * ==> Reponsibility of classes which use subobjects!
   * ==> ATTENTION: No cyclic structures!
   * Each destructor must check whether 'refCount'==0 and throw exception
   * if not. Amounts to calling 'isRefCountZero' before anything else.
   * NOTE: Has to be done in every subclass, because destructors are called
   * starting from the child class upwards.
   * To be done in classes with subobjects:
   * - call 'incrRefCount' of subobj. when ref. is created
   * - call 'decrRefCount' of subobj. when ref. is destroyed
   * <p>
   * Wrapper classes:
   * Some 'MatIF_XXX' classes simply wrap STATSIM class 'XXX' for use within
   * Matlab. In this case, the wrapped object should be kept as
   * 'Handle<XXX>'. If the wrapped object is used as subobject in other
   * objects, pass the handle. Then, even if the wrapper object is
   * destroyed, the 'XXX' object is not.
   * <p>
   * Static methods:
   * If the class has static methods, 'runstatic' has to be implemented and
   * listed in the class list in 'statsimif_call'. The implem. is up to
   * subclasses, but typically a map from method name to ptr to functions
   * is kept as static member, init. upon first call of 'runstatic'.
   * ==> See 'MatIF_ivmmulticlass' for example.
   * <p>
   * Unique identifier:
   * Each existing object has a handle different from all other existing
   * objects (stored in 'hand' here), but handles of deleted objects can be
   * re-used. In general, new handle values are given by the Matlab side. For
   * purposes such as cleanup, each object which exists at some time during a
   * session has a unique identifier 'ident' which is the larger, the more
   * recently the object has been created. Values are not re-used. This is
   * managed by the static 'identCount' which is increased in the constructor
   * here.
   * <p>
   * Matlab debugging ('MatlabDebug'):
   * If this is done in MEX based mode (see 'MatlabDebug'), it can be
   * controlled from Matlab. Methods:
   * - MATLABDEBUG_ACTIVATE:     Activates/deactives facility
   * - MATLABDEBUG_SETDIFFTHRES: Sets 'diffThres' threshold
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class MatIF_baseclass
  {
  public:
    /* Types
     *
     * - PtrToStaticMeth: Pointer to static exported method
     * - PtrToCreateMeth: Pointer to 'create' method
     * - StatMethList:    Type for static method list, for 'runstatic'
     * - SMCIter:         Iterator for static method list
     */

    typedef void (*PtrToStaticMeth)(int,mxArray**,int,const mxArray**);
    typedef MatIF_baseclass* (*PtrToCreateMeth)(int,const mxArray**);
    typedef MAP_TYPE(my_string,PtrToStaticMeth) StatMethList;
    typedef MAP_CONSTITER(my_string,PtrToStaticMeth) SMCIter;

  protected:
    // Static members

    static long identCount;     // s.a.
    static bool methInit;       // method info init. already?
    static int numMeth;         // number of local methods
    static char** methName;     // method names

    // Members

    int hand;     // object handle
    int refCount; // ref. counter (s.a.)
    long ident;   // unique identifier (s.a.)

  public:
    // Constructors

    /**
     * Default constructor
     *
     * @param h Object handle
     */
    MatIF_baseclass(int h) : hand(h),refCount(0),ident(identCount) {
      setup();
      identCount++;
    }

    virtual ~MatIF_baseclass() {
      isRefCountZero();
    }

    // Public static methods

    /**
     * Must be called by constructor BEFORE anything else. If 'methInit' is
     * false, the static structures are init. and 'methInit' is set.
     * Otherwise, the method does nothing.
     * by 'run' to find methods.
     * NOTE: Names in 'methName' must be sorted in asc. order!
     */
    static void setup();

    // Public methods

    /**
     * @return Object handle
     */
    int getHand() const {
      return hand;
    }

    /**
     * @return Unique identifier
     */
    long getIdent() const {
      return ident;
    }

    /**
     * Executes method. The method name is given in the third input argument.
     * Input arguments 4,... (if any) are passed to the method.
     * If method is not found locally, superclass 'run' is called, except for
     * 'baseclass' which throws an exception.
     *
     * @param nlhs Matlab return args
     * @param plhs "
     * @param nrhs Matlab input args
     * @param prhs "
     */
    virtual void run(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]);

    /**
     * See header comment. Locally, 'name' is compared with Matlab class name
     * XXX (for C++ class name MatIF_XXX), return true if identical,
     * otherwise call superclass method.
     *
     * @param name S.a.
     * @return     S.a.
     */
    virtual bool isInstanceOf(const char* name) const {
      // This is the base class:
      return (strcmp(name,"baseclass")==0);
    }

    /** NOT TESTED!!
     * Provides access to members of this object maintained on the Matlab
     * side. These members are simply fields of
     *   OBJECTS{GETPOS4HAND('hand')},
     * where OBJECTS is a global cell array and GETPOS4HAND is a function
     * mapping handles to positions in OBJECTS. The entries of OBJECTS
     * are structures which can have arbitrary, uncontrolled fields.
     *
     * @param name Member name
     * @return     Handle to MatlabMatrix object
     */
    virtual Handle<MatlabMatrix> getMember(const char* name) const;

    /** NOT TESTED!!
     * Provides write access to members of this object maintained on the
     * Matlab side. See 'getMember'.
     *
     * @param name Member name
     * @param val  New value
     */
    virtual void setMember(const char* name,const MatlabMatrix& val) const;

    /**
     * Increase 'refCount' by 1
     */
    void incrRefCount() {
      refCount++;
    }

    /**
     * Decrease 'refCount' by 1
     */
    void decrRefCount() {
      if (--refCount<0) {
	char msg[100];
	sprintf(msg,"INTERNAL ERROR: Object with handle %d has negative ref. counter",hand);
	throw MatIFException(msg);
      }
    }

    // Public exported methods

    /**
     * Activates/deactivates Matlab debugging (see 'MatlabDebug').
     * Works only in MEX based mode.
     *
     * Input:
     * - ACT: Activate feature? [scalar bool]
     */
    void matlabdebug_activate(int nlhs,mxArray *plhs[],int nrhs,
			      const mxArray *prhs[]);

    /**
     * Sets 'diffThres' threshold for Matlab debugging (see 'MatlabDebug').
     * Works only in MEX based mode.
     *
     * Input:
     * - THRES: Threshold value [scalar posit.]
     */
    void matlabdebug_setdiffthres(int nlhs,mxArray *plhs[],int nrhs,
				  const mxArray *prhs[]);

  protected:
    // Internal static methods

    /**
     * Helper method for 'setup', can be used by subclasses.
     * Allocates array of strings 'tarr' with 'num' entries and copies strings
     * from 'sarr'. All allocations via Matlab and persistent.
     * NOTE: If 'num'==0, do nothing.
     *
     * @param sarr Source array
     * @param num  Number of 'sarr' entries
     * @param tarr Target array ret. here
     */
    static void setupNames(const char** sarr,int num,char**& tarr);

    /**
     * Helper method for 'run', can be used in subclasses.
     * Search for entry 'name' in string array 'arr', size 'num'. If found,
     * the position is returned, otherwise -1.
     * NOTE: 'arr' must be sorted in increasing order. Binary search is used.
     *
     * @param name S.a.
     * @param arr  String array
     * @param num  Number of 'arr' entries
     * @return     S.a.
     */
    static int findName(const char* name,char* arr[],int num);

    // Internal methods

    /**
     * Tests whether 'refCount'==0, throws 'MatIFException' if not.
     * NOTE: Must be called by every destructor!
     */
    void isRefCountZero() const {
      if (refCount>0) {
	char msg[100];
	sprintf(msg,"Cannot destroy object with handle %d, ref. counter is %d (references exist)",hand,refCount);
	throw MatIFException(msg);
      }
    }
  };
//ENDNS

#endif
