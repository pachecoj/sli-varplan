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
 * Desc.:  Header class TempBaseVector
 * ------------------------------------------------------------------- */

#ifndef TEMPBASEVECTOR_H
#define TEMPBASEVECTOR_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/matrix/default.h"
#include "lhotse/matrix/TempMatMethods.h"
#include "lhotse/matrix/TempMatMethods_AssignElement.h"

//BEGINNS(matrix)
  /**
   * Helper class for 'BaseVector'.
   * Used to maintain temp. objects returned by 'operator()' and 'mask'.
   * In contrast to a normal return-by-value, this class encaps. a 'BaseVector'
   * object which can also be used as l-value in expressions.
   * Internally, 'TempBaseVector' works like 'Handle' but:
   * - cannot be copied or def. constructed
   * - no ref. counting
   * - conversion operator to 'BaseVector&'
   * - mirrors some methods of 'BaseVector' for convenience
   * <p>
   * The problem is that we cannot know 'BaseVector' so cannot even call the
   * destructor to clean up 'rep'. The solution is to keep two pointers:
   * 'rep' for virtual calls, and 'repCorr' for explicit 'BaseVector'
   * calls.
   *
   * Methods of 'BaseVector' are called through the '->' operator, or
   * through explicit cast to 'BaseVector&'. Examples:
   *   vec(Range(1,3))->fill(...);
   *   ((BaseVector<T>&) vec(Range(1,3))).fill(...);
   *   a=((BaseVector<T>&) vec(Range(1,3)))[1];
   *   if (((BaseVector<T>&) a(Range(...)))==b(Range(...))) ...
   * Since the explicit cast is awkward, some elementary 'BaseVector'
   * operators are mirrored here, using 'TempMatMethods' support:
   * - operator= (vector and scalar variant)
   * Instead of
   *   ((BaseVector<T>&) vec(Range(1,3)))=vec2;
   * you can use
   *   vec(Range(1,3))=vec2;
   * NOTE: The variant using an explicit cast is still valid, and is doing
   * exactly the same thing.
   * ==> TODO: Mirror more operators here
   * See doc/system/masking.txt for details.
   * <p>
   * NOTE: This class has to be derived in parallel for each subclass of
   * 'BaseVector'. F.ex., if B is a subclass of 'BaseVector', then 'TempB' has
   * to be defined. The 'TempXXX' classes do not form a hierarchy, they do
   * not know of each other.
   * <p>
   * NOTE: This is an internal class with some hacks:
   * - To allow for copying through temp. versions, the copy constructor is
   *   implemented, but it resets the repr. pointer of its argument. Otherwise,
   *   the arg. object would destroy the repres. when going out of scope!
   * - Need two repres. pointers 'repCorr', 'rep' to the same object. 'rep'
   *   is not enough, because we cannot cast it to 'BaseVector'
   * - 'operator=' with 'const TempBaseVector&' argument has to implemented
   *   here, otherwise the compiler invents a def. variant which just copies
   *   all members across!
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  template<class T> class TempBaseVector
  {
    friend class BaseVector<T>;
    friend class BaseMatrix<T>;

  protected:
    // Members

    BaseVector<T>* repCorr; // pointer to representation
    mutable TempMatMethods* rep;

    // Methods
  protected:
    /**
     * Must only be invoked by methods of 'BaseVector', 'BaseMatrix'!
     *
     * @param pp Pointer to newly dyn. created 'BaseVector' object
     */
    TempBaseVector(BaseVector<T>* ppCorr,TempMatMethods* pp) : repCorr(ppCorr),
    rep(pp) {}

  public:
    /**
     * Copy constructor.
     * NOTE: Not a CC in the usual sense! The repr. pointer of the 'src'
     * object is set to 0, so after the copying 'src' can only be destroyed.
     * ==> just to allow for creation and copying of temp. versions
     */
    TempBaseVector(const TempBaseVector<T>& src) : rep(src.rep),
    repCorr(src.repCorr) {
      src.rep=0; // 'src' is not usable from now on!
    }

    ~TempBaseVector() {
      if (rep!=0) delete rep;
    }

    /**
     * Allows to call 'BaseVector' methods via dereferencing.
     *
     * @return Pointer to underlying reference
     */
    BaseVector<T>* operator->() const {
      return repCorr;
    }

    /**
     * Drives automatic conversion to 'BaseVector&' (non-const), so objects
     * can be used as l- or r-value arguments.
     * NOTE: Compiler does not invoke automatic conversion in many cases,
     * f.ex. for l-values, or when target type is not clear. Use explicit
     * cast then.
     *
     * @return Ref. to repres.
     */
    operator BaseVector<T>&() const {
      return *repCorr;
    }

    /**
     * Implements '=' for user convenience. It allows you to write
     *   vec(Range(...))=vec2;
     * instead of
     *   ((BaseVector<T>&) vec(Range(...)))=vec2;
     * NOTE: The latter is still valid, and does the same thing.
     *
     * @param arg Source argument
     */
    BaseVector<T>& operator=(const BaseVector<T>& arg) {
      //cout << "TempBaseVector::operator=" << endl;
      rep->assignVirtual(&arg);

      return *repCorr;
    }

    /**
     * We need this one as well, otherwise
     *   vec1(rng1)=vec2(rng2);
     * breaks! Why? In this case, the automatic conversion to
     * BaseVector<T>& is somehow not used, instead the compiler invents a
     * '=' operator which just copies members (which are pointers) across!
     *
     * @param arg Source argument
     */
    BaseVector<T>& operator=(const TempBaseVector<T>& arg) {
      //cout << "TempBaseVector-=TV: rep=" << rep << ",repc=" << repCorr << ",arg.repc=" << arg.repCorr << endl;
      rep->assignVirtual(arg.repCorr);

      return *repCorr;
    }

    /**
     * Implements '=' with scalar r.h.s. for user convenience. It allows
     * you to write
     *   vec(Range(...))=scalar;
     * instead of
     *   ((BaseVector<T>&) vec(Range(...)))=scalar;
     * NOTE: The latter is still valid, and does the same thing.
     *
     * @param arg Source argument
     */
    BaseVector<T>& operator=(T arg) {
      //cout << "TempBaseVector::operator=" << endl;
      TempMatMethods_AssignElement<T> a(arg);
      rep->assignVirtual(&a);

      return *repCorr;
    }
  };
//ENDNS

#endif
