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
 * Desc.:  Header class BaseVector
 * ------------------------------------------------------------------- */

#ifndef BASEVECTOR_H
#define BASEVECTOR_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

/*
 * TODO:
 * - methods based on matching (sorting): 'isDisjointFrom', 'isSubset',
 *   'containsDupl', ... These should first assume that vectors are
 *   sorted. If it is found they are not, go back and sort them. Draw
 *   copies of non-flat vectors. Also: merging of sorted vecs.
 *   NOTE: STL has much of this!
 * - STL bidirec. iterator! How difficult is that?? Quite tough!
 */

/*
 * Internal storage format:
 * A vector can be normal or a mask vector.
 * The repres. of a normal n vector is given by a flat buffer 'buff' of size
 * 'size_n', where 'size_n'>=n. The last 'size_n'-n elements are not used. The
 * buffer exists iff 'size_n'>0. An empty vector may have an existing buffer.
 *
 * A mask vector does not own its buffer. It can be hooked to a matrix/vector
 * base object via 'baseObj'. In this case, we check the validity of the buffer
 * of the base object using the timestamp mechanism (see doc/masking.txt and
 * doc/memManage.txt). The member 'buff_n' is not used for a mask vector.
 * A mask vector can be flat, linear or indexed. If indexed, it uses a position
 * index 'mindex'.
 * NOTE: In this case, the size of 'mindex' is >= n. If it is > n, then only
 * the initial positions are used.
 * A linear mask vector can have an assoc. step size 'step'. A mask is flat iff
 * 'step'==1. Note that 'step' < 0 is allowed.
 */

#include "lhotse/FileUtils.h"
#include "lhotse/NumberFormats.h"
#include <functional>
#include "lhotse/matrix/default.h"
#include "lhotse/matrix/Vector.h"
#include "lhotse/matrix/ArrayUtils.h"
#include "lhotse/matrix/MatTimeStamp.h"
#include "lhotse/matrix/TempMatMethods.h"
#include "lhotse/matrix/TempBaseVector.h"
#include "lhotse/matrix/MatDefMembers.h"
#include "lhotse/matrix/WriteBackVec.h"
#include "lhotse/matrix/BaseLinVec.h"
#ifdef MATLAB_MEX
#include "lhotse/matif/predecl.h"
#include "lhotse/matif/MatlabMatrix.h"
#endif

// Local macro for args to 'checkTS'
#define TSARG __FILE__,__LINE__

//BEGINNS(matrix)
  /**
   * Helper class to allow 'BaseVector' methods (e.g., 'find') to return
   * temp. created vector objects by value without copying. See 'BaseVector'
   * header comment for details.
   * NOTE: A return-by-value is automatically restricted to be an r-value
   * in the enclosing expression. If this is not appropriate, use the
   * 'TempBaseVector' mechanism (e.g., 'operator()').
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  template<class T> class BaseVecWrapper
  {
  protected:
    // Members

    BaseVector<T>* repr; // representation

  public:
    // Methods

    /**
     * Default constructor
     *
     * @param rep Repres.
     */
    BaseVecWrapper(BaseVector<T>* rep) : repr(rep) {}

    /**
     * @return Representation
     */
    BaseVector<T>* getRep() const {
      return repr;
    }

    /**
     * Reset repres. to 0.
     */
    void reset() {
      repr=0;
    }
  };

  /**
   * Base (non-abstract) class for vectors which use a flat buffer, whose
   * elements can be copied and assigned.
   * <p>
   * Mask Vectors:
   * See comment above (storage format) and doc/masking.txt. In general, a
   * flat or linear mask vector is as efficient as a normal one, while an
   * indexed mask vector is much less efficient to use. It is not allowed to
   * change the size of a mask vector ('MaskObjectException').
   * If a mask vector is hooked to an underlying matrix/vector object, we
   * use a timestamp mech. ('MatTimeStamp') to find out about size changes of
   * this object. Such invalidate the mask vector (a 'MaskObjectException' is
   * thrown if it is used).
   * <p>
   * Indexed mask vectors:
   * Most methods, esp. the arithmetic ones in subclasses, do not deal with
   * indexed masks directly, but draw flat copies. On the other hand,
   * assignment or initialisation methods typically work with indexed masks,
   * and some O(n) arithm. methods may have special implementations which
   * do not copy (see comments).
   * The index vector 'mindex' can be shared between different mask vectors
   * (see 'reassign'), but it is NEVER shared between an indexed range passed
   * to a reassign method and a mask vector.
   * ==> OK, because the content of 'mindex' is never changed explicitly.
   * <p>
   * Dealing with non-flat vectors:
   * For many arithmetic methods, it is inefficient or cumbersome to support
   * indexed or non-flat vectors directy. In these cases, flat copies are
   * drawn and written back at the end.
   * - 'getLinVec' creates 'BaseLinVec' object which represents this vector.
   *   Features:
   *   - if this vector is flat or linear, the object repr. it as such without
   *     copying. Write-back is not required, the 'BaseLinVec' object acts
   *     like a mask
   *   - if this vector is indexed, a flat copy is drawn. If this vector is
   *     used as l-value, an automatic write-back upon destruction of the
   *     flat copy is supported.
   *   - pointer conversion for difference between our and Fortran's way of
   *     repres. linear vectors with neg. step size
   *   ==> Preferred way to deal with non-flat or non-linear vectors, supports
   *       Fortran interface directly
   * - 'getFlatBuff' returns flat repres. of this vector (as 'ArrayHandle').
   *   If this vector is flat, the repres. is like a mask and uses the same
   *   mem. watcher. Otherwise, a flat copy is drawn.
   *   ==> Simpler than 'getLinVec', but no write-back support. Also, linear
   *       non-flat mask vectors are copied
   * <p>
   * Buffer control:
   * Every vector is defined relative to a buffer 'buff'. Several vectors or
   * even other objects can share access to a memory region, esp. once mask
   * vectors are used. Some control is imposed by trying to ensure that for
   * every dyn. alloc. memory region, a 'MemWatchBase' object is maintained.
   * This objects counts the number of "users" of this region. It allows
   * deallocation only if this ref. counter drops to 0. Every user has to
   * keep a ref. to the watcher. The 'buffWatch' member is defined in the
   * 'MatTimeStamp' superclass (see comments there).
   * NOTE: Buffer control works only if fully supported by all methods here
   * and in subclasses!
   * - if the object de-assoc. from its buffer, 'deassocBuff' has to
   *   be called. This decr. the ref. counter and destroys the mem. watcher
   *   if the counter drops to 0
   * - if the object assoc. with a new buffer, 'assocBuff' has to be called
   *   passing the corr. mem. watcher. 'assocBuff' increases the ref.
   *   counter iff the object is a mask vector
   * - if the vector size changes, 'incrTS' has to be called (see
   *   'MatTimeStamp')
   * - call 'checkTS' for this object and every vector argument (see
   *   'MatTimeStamp') in each method. 'checkTS' does nothing if this object
   *   is not a mask, for a mask it checks whether the timestamp of the
   *   underlying normal object is the same as ours. If not, a
   *   'MaskObjectException' is thrown (the mask is invalid)
   * NOTE: For compat. with other code, mask vectors may have no mem.
   * watcher at all. Such masks are NOT protected and unsafe!
   * <p>
   * Size of memory buffer:
   * Typically, the size of the memory buffer never shrinks. An exception are
   * 'resize' and 'resizeSave' which adjust the buffer size explicitely.
   * Buffer size is increased on demand only (again, 'resizeXXX' are
   * exceptions). There are two different situations we distinguish:
   * - buffer increases to a fixed target size to be able to hold a result,
   *   previous content not copied (will be overwritten). This is typ. done
   *   by 'ensureCapacity'. The new buffer size is exactly what is required.
   * - buffer increases because elements are inserted or are appended. The
   *   prev. content is copied. See 'expand', 'insert'. Here, the new buffer
   *   size is det. by a strategy (see 'MatDefMembers') and it may well be
   *   larger than abs. required. The idea is to avoid extensive copying if
   *   f. ex. a vector is build up by inserting single elements or small
   *   patches.
   * <p>
   * Methods returning vector objects by value:
   * Two different situations where return by value is required:
   * (1) Object created locally, to be assigned to a non-temp. vector. The
   *     returned object can be directly used as r-value only. It can be
   *     assigned to a new or existing vector object which can then be
   *     modified. The feature is that this does NOT involve copying and
   *     re-alloc. buffers. It also allows to return mask vectors by value.
   * (2) Mask to a non-temp. vector to be returned as inherently temp.
   *     object, on which l-value operations can be performed. Method (1)
   *     cannot be used because the compiler allows r-values only to be
   *     returned as in (1).
   *
   * (1): To avoid copying the locally created object, we use the wrapper
   * class 'BaseVecWrapper'. This works as follows:
   * - method in question (say: f) has ret. type 'BaseVector<T>'
   * - the return statement returns a 'BaseVecWrapper<T>' object (locally
   *   created, wraps dyn. alloc. 'BaseVector<T>' (or subclass) object)
   * - 'BaseVector' and all subclasses implement "copy" constructor with
   *   'BaseVecWrapper' argument. This uses 'convertWrapped' to "weed out" the
   *   argument, use its fields for the new vector and dispose of the "empty
   *   hull". This avoids copying
   * NOTE: Avoiding copying is not only imp. for efficiency, also to retain
   * mask vector status, because the normal CC creates a normal vector from
   * a mask vector!
   *
   * (2): Requires more work because of the r-value restriction of the
   * compiler (meaning that temp. created objects cannot be used as
   * l-values). It is important for methods like 'operator()' and 'mask' which
   * allow to access part of a given vector/matrix/mem. region as if it would
   * be a genuine 'BaseVector' object itself, also as l-value. This is done
   * WITHOUT creating a copy. Instead, the methods create temp. mask objects.
   * The methods creating temp. vectors return 'TempBaseVector' objects which
   * wrap the dyn. alloc. vector. 'TempBaseVector' has a conversion operator
   * to 'BaseVector&' (non-const) and the deref. operator '->' to call
   * 'BaseVector' methods. It also implements 'operator=' which calls the
   * variant here through 'assignVirtual'. See 'TempBaseVector' and
   * doc/system/masking.txt for details.
   *
   * NOTE: (2) is messy and should only be used if (1) does not apply bec.
   * of the r-value restriction! If possible, use of (2) should be restricted
   * to 'operator()' and 'mask'.
   * On the other hand, our mechanism for (2) works just as well for (1),
   * so we could save the additional mechanism for (1). We keep it, because
   * it is simpler and more stable (in some situations, the (2) mechanism
   * requires explicit casts to 'BaseVector&', where (1) does not).
   *
   * For every subclass A of this one, a corr. subclass TempA of
   * 'TempBaseVector' has to be defined. This is to avoid inconveniences
   * such as dynamic casts. Also, the methods ret. vector objects have to
   * be overwritten in A to return the correct type TempA.
   * <p>
   * Default fill value:
   * See superclass 'MatDefMembers'
   * <p>
   * Persistence:
   * 'save' and 'load' use type codes def. in class 'TypeCodeConsts'. Types
   * which are not recognized there, have type code 'typeOther'.
   * <p>
   * Restrictions on elements (in subclasses):
   * The elements of a 'BaseVector' are by def. unrestricted members of T.
   * In subclasses, there may be restrictions. Assume the restrictions are
   * placed on the elements: a vector is valid if all its elements are.
   * In this case, a lot of overwriting of methods can be avoided by just:
   * - overwrite 'doCheckValid' to return true
   * - overwrite 'isValidElement'
   * - overwrite 'isValidVecType'
   * The methods 'isValidVec' and 'isValidConvVec' have def. implem. based on
   * the others.
   *
   * If a subvector is of a type which enforces the element restriction, then
   * its elements do not have to be checked. This is realised via
   * 'isValidVecType': given a vector pointer, this method uses a dynamic
   * cast to find out whether the req. restriction is already enforced.
   * The def. implem. of 'isValidVec' calls 'isValidVecType' before it checks
   * the single elements.
   *
   * All methods which change elements will call 'doCheckValid' and check
   * validity if this ret. true. Some short ones will use 'isValidElement'
   * without calling 'doCheckValid'.
   * The def. implement. are:
   * - doCheckValid: returns false
   * - isValidElement: returns true
   * - isValidVecType: returns true
   * NOTE: If the restrictions are more complicated, or if the internal repres.
   * of a vector changes, then (at least) all methods changing elements of
   * the vector have to be overwritten.
   * See doc/basevector.txt for more details.
   *
   * Example: 'IndexVector' requires non-negativity. Its implement. has a
   * 'doCheckValid' which returns true, a 'isValidElement' which checks for
   * non-negat., and a 'isValidVecType' which does a dyn. cast to
   * 'IndexVector'.
   *
   * NOTE: Also the methods which reassign masks do validity checking. This
   * is because in subclasses, these methods can be used to do type
   * conversions. Internally, if 'doCheckValid' ret. true and a check is
   * required, the 'reassignInt' methods (which don't do validity checking)
   * are used to create a temp. mask which is assigned to this object if
   * the valid. checks are OK.
   *
   * Holes in the checking:
   * operator[] is defined even if there are checks in place, so an assignment
   *   vec[0]=a
   * is possible even if 'a' is not a valid element.
   * ==> In earlier version, the l-value operator[] threw exception if
   *     'doCheckValid' ret. true. Problem: l-value version is sometimes
   *     used if the r-value version would be appropriate! SOLUTION??
   * <p>
   * Implementation of subclasses:
   * As much of the code as possible is written generically. Specialised
   * issues:
   * - methods which change entries of the vector (also 'reassign') have to
   *   check whether the new entries confirm to invariants of the concrete
   *   type. This is done via the 'isValidXXX' methods
   * - methods which create temp. mask vectors ('operator()', 'mask') are
   *   specific in several respects:
   *   - their return type changes (to 'TempA' if the class is 'A')
   *   - they have to create temp. mask vectors of the correct type, and wrap
   *     them into 'TempA'
   * - if 'A' is the class, 'TempA' has to be def. in the 'TempBaseVector'
   *   hierarchy and several methods have to be overwritten (see
   *   'TempBaseVector')
   *
   * The following methods have to be overwritten in subclasses:
   * - if restrictions on elements: 'doCheckValid', etc. (see above)
   * - 'convertWrapped': for 'BaseVecWrapper' constructor
   *   (only if new members/restrictions are added)
   * - 'newEmpty': call default constructor (for empty vector)
   * - 'copy': call copy constructor
   * - "copy" constructor from 'BaseVecWrapper': just copy the implem. here,
   *   but provide 'convertWrapped'
   * - normal copy constructor: copy code here.
   *   NOTE: Do not simply call copy constructor of superclass. If the
   *   subclass adds element restrictions, these are checked in 'assignInt',
   *   but they are not checked if the superclass constructor is called
   * - constructor with expl. fill element: copy code here ('init' is virtual).
   *   NOTE: Do not simply call superclass constructor, this would not invoke
   *   element checks of the subclass.
   * - 'mask': create vector of correct type (static): just copy the implem.
   *   here, but change return type to 'TempA'
   * - 'operator()', r-value and -lvalue version: just copy implem. here,
   *   but dyn. cast pointer ret. from internal method to A* before creating
   *   'TempA' object. Return type: 'TempA'
   *   NOTE: Since 'operator()' is not virtual, overloading does not work, so
   *   also the generic versions (with int argument) have to be re-def. in
   *   every subclass (just copy implem. here)
   * - operator=: Not virtual, due to changing return type. Just copy code
   *   here
   * - 'assignVirtual': This must do exactly the same as 'operator=', except
   *   returning *this. If the subclass 'operator=' just calls the
   *   superclass variant, 'assignVirtual' does not have to be overwritten
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  template<class T> class BaseVector : public Vector,public MatTimeStamp,
				       public TempMatMethods,
				       public MatDefMembers<T>,
				       public WriteBackVec<T>
  {
    friend class BaseMatrix<T>;
    friend class BaseLinVec<T>;

  public:
    // Static members/constants
    // - mask vector iff != 'statNormal'
    // - uses 'mindex' iff == 'statMaskInd'
    static const unsigned char statNormal     =0;
    static const unsigned char statMaskFlat   =1;
    static const unsigned char statMaskLin    =2;
    static const unsigned char statMaskInd    =3;

  protected:
    // Additional members

    T* buff;                  // Buffer pointer
    unsigned char maskStat;   // Mask status (see 'statXXX' constants)
    union {
      int size_n;             // Number buffer entries (normal vector only)
      const MatTimeStamp* baseObj; // Underlying object (mask vector only)
    };
    int step;                 // Step size (not for indexed masks)
    ArrayHandle<int> mindex;  // Position index (indexed masks only)

  public:

    // Constructors

    /**
     * Constructs empty vector.
     */
    BaseVector() : Vector(),MatTimeStamp(),TempMatMethods(),size_n(0),buff(0),
		   maskStat(statNormal),step(1) {
      this->setDefValues();
    }

    /**
     * Constructor.
     * If 'a' is not given, the def. fill value is used. NOTE: If 'a' is given,
     * the def. fill value is NOT set to 'a', but init. by def. as usual.
     *
     * @param num Number of entries
     * @param a   Fill value. Optional
     */
    explicit BaseVector(int num) : Vector(),MatTimeStamp(),TempMatMethods(),
				   size_n(0),buff(0),maskStat(statNormal),
				   step(1) {
      this->setDefValues();
      init(num);
      if (num>0) incrTS();
    }

    BaseVector(int num,const T& a) : Vector(),MatTimeStamp(),TempMatMethods(),
    size_n(0),buff(0),maskStat(statNormal),step(1) {
      this->setDefValues();
      init(num,a);
      if (num>0) incrTS();
    }

    /**
     * Copy constructor.
     * Creates a normal vector, no matter what 'vec' is.
     *
     * @param vec Vector to copy
     */
    BaseVector(const BaseVector<T>& vec) : Vector(),MatTimeStamp(),
    TempMatMethods(),size_n(0),buff(0),maskStat(statNormal),step(1) {
      vec.checkTS(TSARG);
      this->copyDefValues(vec); // def. value members
      assignInt(vec); // copy
    }

    /**
     * Special "copy" constructor from 'BaseVecWrapper'. See header comment
     * for details. This constructor has to be implemented in all subclasses!
     * It creates a new vector by casting the repres. of the wrapper to the
     * class type, using all its fields (it is a temp. mask object) for the
     * new object (method 'convertWrapped') ==> no copying involved.
     *
     * @param arg Wrapped argument
     */
    BaseVector(const BaseVecWrapper<T>& arg) : Vector(),MatTimeStamp(),
    TempMatMethods() {
      // use 'arg' fields to init. ours, then dispose of 'arg'
      convertWrapped((BaseVecWrapper<T>&) arg);
    }

    /**
     * Copy constructor from STL vector
     *
     * @param arg STL vector
     */
    BaseVector(const VEC_TYPE(T)& arg) :
      Vector(),MatTimeStamp(),TempMatMethods(),size_n(0),buff(0),
      maskStat(statNormal),step(1) {
      this->setDefValues();
      int an=arg.size();
      if (an>0) {
	init(an);
	incrTS();
	typename VEC_CONSTITER(T) it=arg.begin();
	for (int i=0; i<an; i++,++it) buff[i]=*it;
      }
    }

#ifdef MATLAB_MEX
    /**
     * Special copy constructor from 'MatlabMatrix' object. This uses
     * static cast to convert from double to T. For T==double, it is
     * more economic to mask the object using 'StVector::maskMatlab'.
     * NOTE: 'mat' is read column-order.
     * NOTE: Cast double -> T must work as expected!
     *
     * @param mat Matlab matrix
     */
    BaseVector(const MatlabMatrix& mat);
#endif

    /**
     * Destructor
     */
    virtual ~BaseVector() {
      //cout << "this=" << this << endl;
      //cout << "~BaseVector" << endl;
      dealloc();
    }

    // Information methods

    /**
     * @return Is this vector a mask?
     */
    bool isMask() const {
      return (maskStat!=statNormal);
    }

    /**
     * A normal empty vector can be transformed into a mask vector if it has
     * been empty since creation. This method returns true iff this vector
     * is a mask or is empty and has been since creation. The latter is
     * checked by comp. its timestamp value with 0.
     * A 'reassign' method works for this object only if this method returns
     * true.
     *
     * @return Is this vector a mask or the empty vector (without a buffer)?
     */
    bool isMaskOrEmpty() const {
      return (maskStat!=statNormal || (size_n==0 && getTS()==0));
    }

    /**
     * @return Is this a mask vector with a position index?
     */
    bool isMaskInd() const {
      return (maskStat==statMaskInd);
    }

    /**
     * @return Is this a flat mask vector (step size 1)?
     */
    bool isMaskFlat() const {
      return (maskStat==statMaskFlat);
    }

    /**
     * If the vector is "flat", then 'getFlatBuff' will return a handle which
     * coincides with the vector buffer, i.e. modifications will be done
     * directly on the buffer. Otherwise, 'getFlatBuff' returns a flat copy.
     *
     * @return Is this vector flat?
     */
    virtual bool isFlat() const {
      return (maskStat==statNormal || maskStat==statMaskFlat);
    }

    /**
     * @return Mask status (see constants 'statXXX')
     */
    unsigned char getMaskStatus() const {
      return maskStat;
    }

    /**
     * @return Size of vector buffer, or 0 for a mask object
     */
    virtual int getBufferSize() const {
      return (!isMask())?size_n:0;
    }

    /**
     * Returns pos. index for indexed mask vector, 0 handle otherwise.
     *
     * @return S.a.
     */
    const ArrayHandle<int>& getMaskIndex() const {
      return isMaskInd()?mindex:ArrayHandleZero<int>::get();
    }

    /*
     * Mask vector methods
     * - reassign:   reassigns mask vector to memory region
     * - mask:       static version of 'reassign', returns temp. mask vector
     * - operator(): subselection, returns temp. mask vector
     */

    /**
     * Reassigns a mask vector (or the empty vector) to point to some buffer
     * region 'pbuff'. 'pbuff' must span at least 0,...,'pn'-1. If 'rng' is
     * not given, the mask vector will be flat with range 0,...,'pn'-1.
     * Otherwise, it will use the range 'rng'. In this case, the mask vector
     * is indexed iff 'rng' is an index. 'rng' is checked against 0,...,'pn'-1
     * for buffer violations.
     * NOTE: If 'rng' is indexed, a copy of its index vector is drawn here.
     * <p>
     * A mem. watcher for the buffer region can be passed via 'watch'. If so,
     * it is used as 'buffWatch' here (ref. counter is increased, call to
     * 'assocBuff'). If not, there is no buffer control (dangerous!).
     * ==> Not recommended!
     * If 'noCp'==true and 'rng' is indexed, its index vector is re-used here
     * if possible ==> only for internal methods!
     * <p>
     * NOTE: This is just to interface foreign code! Use 'ArrayHandle' version
     * below whenever possible!
     * <p>
     * NOTE: An empty vector can be made a mask vector here only if it has
     * been empty since creation!
     *
     * @param pbuff See above
     * @param pn    "
     * @param rng   ". Optional.
     * @param watch "
     * @param noCp  ". Def.: false
     */
    virtual void reassign(const T* pbuff,int pn,
			  const Range& rng=RangeFull::get(),
			  MemWatchBase* watch=0,bool noCp=false) {
      if (!isMaskOrEmpty())
	throw MaskObjectException("Needs to be mask or empty!");
      if (rng.checkRange(pn)) throw OutOfRangeException(EXCEPT_MSG(""));
      if (doCheckValid() && !isValidVec(pbuff,pn,rng))
	throw InvalidParameterException("Invalid vector elements");
      reassignInt(pbuff,pn,rng,watch,noCp);
    }

    /**
     * Same as above, but 'pbuff' is an ArrayHandle. We use the mem. watcher
     * of this handle.
     *
     * @param pbuff See above
     * @param rng   "
     */
    virtual void reassign(const ArrayHandle<T>& pbuff,
			  const Range& rng=RangeFull::get()) {
      reassign(pbuff.p(),pbuff.size(),rng,pbuff.getMemWatch());
    }

    /**
     * Reassigns a mask vector (or the empty vector) to point to (part of)
     * another vector.
     * The semantics is the same as for 'reassign' with a flat buffer, i.e.
     * 'vec' is regarded as flat vector. The difference is that here, 'vec'
     * (or its parent) is registered as parent (in 'baseObj').
     * If 'vec' is itself an indexed mask vector, then the position index
     * is mapped through the range given by 'rng'. The def. for 'rng' is the
     * full range.
     * NOTE: If 'vec' is an indexed mask and 'rng' the full range, we do not
     * draw a copy of the pos. index, but just copy the handle 'mindex'.
     * <p>
     * NOTE: An empty vector can be made a mask vector here only if it has
     * been empty since creation!
     * <p>
     * NOTE: 'reassign' can be used to draw a copy of a mask vector. F.ex. if
     * a is a mask vector, b is empty or a mask, then
     *   b.reassign(a);
     * draws a copy. If c is a normal vector and ind a position index, then
     *   b.reassign(c(Range(ind)));
     * creates a temp. mask which is then reassigned to b. Note that in this
     * example, the temp. mask c(Range(ind)) and b share the same pos. index.
     *
     * @param vec   Source vector (to refer to)
     * @param rng   Range. Def.: full range
     */
    virtual void reassign(const BaseVector<T>& vec,
			  const Range& rng=RangeFull::get());

    /**
     * Creates a mask vector picking out range 'rng' from the memory region
     * starting at 'pbuff', length 'pn'. 'rng' must not leave this region.
     * The def. range is the empty one (whole region).
     * NOTE: As opposed to 'operator()', there is no notion of "expanding
     * the buffer" here: 'pbuff' must be large enough.
     * If 'watch' is given, it is used as mem. watcher 'buffWatch' for the
     * mem. region, otherwise no watcher is used (dangerous!).
     * <p>
     * NOTE: This is just to interface foreign code! Use 'ArrayHandle' version
     * below whenever possible!
     * <p>
     * Has to be implemented in every subclass with the corr. ret. type.
     *
     * @param pbuff See above
     * @param pn    "
     * @param rng   ". Optional.
     * @param watch "
     */
    static TempBaseVector<T> mask(const T* pbuff,int pn,
				  const Range& rng=RangeFull::get(),
				  MemWatchBase* watch=0) {
      // In subclasses: copy code, but replace following line by one where the
      // empty vector is created using the target type.
      BaseVector<T>* vec=new BaseVector<T>();
      vec->reassign(pbuff,pn,rng,watch);
      return TempBaseVector<T>(vec,vec);
    }

    /**
     * Same as above, but with 'ArrayHandle' argument 'parr'. We use the mem.
     * watcher of this handle.
     *
     * @param parr Array to be masked as vector
     * @param rng  Range to be picked out. Optional. Def.: full
     */
    static TempBaseVector<T> mask(const ArrayHandle<T>& parr,
				  const Range& rng=RangeFull::get()) {
      // In subclasses: copy code, but replace following line by one where the
      // empty vector is created using the target type.
      BaseVector<T>* vec=new BaseVector<T>();
      vec->reassign(parr.p(),parr.size(),rng,parr.getMemWatch());
      return TempBaseVector<T>(vec,vec);
    }

    /**
     * The purpose for this static method is type conversion without copying.
     * It checks whether 'vec' can be seen as vector of this class. In general
     * this is true if:
     * - 'vec' passes the element checks for this class
     * - if this class has additional members as comp. to 'BaseVector', 'vec'
     *   has to be at least of the most general type which has all the same
     *   members as this class
     * NOTE: 'vec' may be of a type more specific than this class, in which
     * case information can get lost in the conversion.
     *
     * @param vec S.a.
     * @param rng Range to pick out. Def.: full
     * @return    Temp. mask vector
     */
    static TempBaseVector<T> mask(const BaseVector<T>& vec,
				  const Range& rng=RangeFull::get()) {
      // In subclasses: copy code, but replace following line by one where the
      // empty vector is created using the target type.
      BaseVector<T>* mvec=new BaseVector<T>();
      mvec->reassign(vec,rng);
      return TempBaseVector<T>(mvec,mvec);
    }

    /**
     * RValue version:
     * Creates a mask vector picking out range 'rng' from this vector. The
     * resulting mask vector can be used as rvalue, e.g. with the assignment
     * operator. 'rng' must fit within this vector.
     * <p>
     * Has to be implemented in every subclass with the corr. ret. type.
     * The implementation here just has to be copied.
     *
     * @param rng Range to pick out (can be a single position)
     */
    const TempBaseVector<T> operator()(const Range& rng=RangeFull::get())
    const {
      // In subclasses: dyn. cast pointer to target type!
      BaseVector<T>* vec=subrvalInt(rng);
      return TempBaseVector<T>(vec,vec);
    }

    const T& operator()(int rng) const {
      return operator[](rng);
    }

    /**
     * LValue version:
     * Creates a mask vector picking out range 'rng' from this vector. The
     * resulting mask vector can be used as lvalue, e.g. as left-hand side
     * in an assignment.
     * If 'rng' does not fit within this vector, then this vector is expanded
     * to the smallest size which just fits 'rng'. The new entries are init.
     * with the def. fill value.
     * <p>
     * ATTENTION: Strictly avoid usings masks to the same underlying object
     * as rvalue and lvalue in an expression. This typically leads to a
     * 'MaskObjectException'.
     *   a(rng1)=a(rng2); // ERROR!
     * <p>
     * NOTE: Careful when using ranges without right end. In this case, the
     * right end is set to the maximum of 'n' and the left end of 'rng', so if
     * the vector has to be enlarged, the res. mask vector will have size 1
     * only! Example: if 'a' has length 5, then
     *   a(Range(5))=b;
     * works only if 'b' has length 1, because the l.h.s. will be a temp.
     * mask of length 1 ref. to element 5 of 'a' (which is 0 and has been
     * appended to 'a' by this method).
     * NOTE: In the example above, the size of 'a' cannot be adjusted based
     * on the size of 'b'. The mask vector for 'a(...)' is created before it
     * enters the assignment. The correct way would be:
     *   a(Range(5,4+b.size()))=b;
     * <p>
     * Has to be implemented in every subclass with the corr. ret. type.
     * The implementation here just has to be copied.
     *
     * @param rng Range to pick out. Vector expanded if necessary. Can be
     *            single pos.
     */
    TempBaseVector<T> operator()(const Range& rng=RangeFull::get()) {
      // In subclasses: dyn. cast to target type!
      BaseVector<T>* vec=sublvalInt(rng);
      return TempBaseVector<T>(vec,vec);
    }

    /*
     * NOTE: Access NOT protected, even if 'doCheckValid' is true!
     */
    T& operator()(int rng) {
      if (rng>=n)
	expand(rng+1);
      return operator[](rng);
    }

    // Access to flat/linear vector (copy if necessary)

    /**
     * Returns 'ArrayHandle' containing a flat representation of this vector.
     * If this vector is flat (normal or flat mask with step size 1), we return
     * a handle which uses the same mem. watcher 'buffWatch'. Otherwise, a flat
     * copy of this vector is drawn and wrapped into a new handle.
     * NOTE: If this vector is a flat mask, but does not have a mem. watcher,
     * the handle returned does not have a mem. watcher either.
     * ==> Not recommended, there is no protection against premature
     *     destruction!
     * NOTE: If this vector is empty, the zero handle is returned.
     * NOTE: Use to interface code which requires explicit flat arrays only!
     * ==> most arithmetic methods in subclasses will use this method if the
     *     mask is indexed, while they usually can deal with linear (non-flat)
     *     masks (see 'BaseLinVec'). Use 'getLinVec' instead!
     *
     * @return See above
     */
    virtual ArrayHandle<T> getFlatBuff() const;

    /**
     * The 'BaseLinVec' object 'blVec' (whose entries must be empty) is init.
     * with a flat or linear representation of this vector. This is used to
     * interface with code which does not allow for indexing.
     * A flat copy of this vector is drawn only if it is an indexed mask.
     * However, if 'doFlat'==true (def.: false) and this is a linear non-flat
     * mask, a flat copy is drawn as well.
     * NOTE: Method does not work for empty vector.
     * <p>
     * Write-back:
     * If 'writeBack'==true, the 'BaseLinVec' object is conf. s.t. upon
     * its destruction the content is written back into this vector. See
     * 'BaseLinVec' for details.
     * <p>
     * Negative step size:
     * See 'BaseLinVec' and member 'negFort' for details. 'negFort' here
     * supplies the value. Def.: true (Fortran convention).
     *
     * @param blVec     BaseLinVec object ret. here
     * @param writeBack S.a.
     * @param negFort   S.a. Def.: true
     * @param doFlat    S.a. Def.: false
     */
    virtual void getLinVec(BaseLinVec<T>& blVec,bool writeBack,
			   bool negFort=true,bool doFlat=false) const;

    // Subscripting

    /**
     * Subscripting operator
     * NOTE: Even though this is virtual, the compiler should inline it if used
     * in a non-polymorphic context(?!).
     * NOTE: Access NOT protected, even if 'doCheckValid' is true.
     *
     * @param pos Position
     */
    virtual const T& operator[](int pos) const;
    virtual T& operator[](int pos);

    /**
     * Safe version of l-value operator[]. If 'doCheckValid' returns true, i.e.
     * if there are restrictions on the elements, only this method works.
     *
     * @param pos  Position
     * @param elem New element value
     */
    virtual void set(int pos,T elem);

    // Changing vector and buffer size

    /**
     * Resizes vector to given size and sets vector entries to def. fill
     * value.
     * <p>
     * Note: To simply ensure that the buffer can hold an n vector, it's
     * better to use 'fill', the latter will NOT reallocate the buffer
     * if it's too large
     *
     * @param num Desired size
     */
    virtual void resize(int num);

    /**
     * Resizes vector to given size, but keeps the old entries. A size smaller
     * than the current is refused. The new entries are set to the def. fill
     * value. The buffer size is equal to 'num', i.e. the vector size, after
     * a call of this method. See also 'expand'.
     *
     * @param num Desired size
     */
    virtual void resizeSave(int num);

    /**
     * Expands vector to given size, but keeps the old entries. A size that
     * would lead to loss of entries is refused. The new entries are set to
     * the def. fill value.
     * This method behaves different from 'resizeSave'. While 'resizeSave'
     * reallocates the buffer if the buffer size is != 'num', 'expand'
     * reallocates only if 'num' is larger than the buffer size. In case of
     * realloc., the new buffer size is det. by the buffer increase strat.
     * (see 'MatDefMembers'), but it is >= 'num'.
     * <p>
     * ATTENTION: Old code uses 'expand(num,offset)'. Replace!
     *
     * @param num    Desired vector size
     */
    virtual void expand(int num);

    /**
     * "Unsafe" version of 'expand'. If 'num' > 'n', 'expand' is called. If
     * 'num' < 'n', the vector is shrunk to size 'num', retaining the first
     * 'num' entries, the remaining ones are lost. In the latter case, the
     * buffer size is not changed.
     *
     * @param num New vector size. Must be positive
     */
    virtual void shrink(int num);

    // Assignment, copy, convert

    /**
     * Assignment operator
     * <p>
     * Physically copies all elements of 'vec'.
     * NOTE: To obtain a mask, use 'reassign'. 'operator()' gives a temporary
     * mask.
     * <p>
     * NOTE: This vector may be a mask vector as well. In this case, if
     * 'vec' is shorter than this vector, instead of throwing a
     * 'MaskObjectException', we copy 'vec' to the first positions of this
     * vector. If this vector is too short and a mask, a 'MaskObjectException'
     * is thrown, though. If this vector is not a mask, its size will be the
     * same as of 'vec' afterwards.
     * ==> Allows more intuitive use together with operator().
     *
     * @param vec Vector to copy
     */
    BaseVector<T>& operator=(const BaseVector<T>& vec) {
      //cout << "BaseVector-=: this=" << this << ",arg=" << &vec << endl;
      assignInt(vec);
      return *this;
    }

    /**
     * Assignment operator, r.h.s. a scalar
     * <p>
     * Sets all elements of this vector to 'elem'. Same as 'fill'. The length
     * of this vector is not changed, except if it is empty (grows to size 1
     * in this case).
     *
     * @param elem Fill element
     */
    BaseVector<T>& operator=(T elem) {
      assignInt(elem);
      return *this;
    }

    /**
     * Exchanges content of this vector and 'a' (must have same size).
     * <p>
     * NOTE: This is done from left to right, using the element swap
     * (this->temp,a->this,temp->a) (important if they overlap).
     *
     * @param a S.a.
     */
    virtual void exchange(BaseVector<T>& a);

    /**
     * Conversion from vector with different entry type. This is a version
     * of the assignment operator, thus elements are copied physically.
     * Casts from T2 to T must be possible and result in valid conversions.
     * <p>
     * If this vector is a mask and longer than 'src', we copy 'src' into
     * the initial positions of this vector.
     *
     * @param src Source object, element type T2
     */
    template<class T2> void convert(const BaseVector<T2>& src) {
      if (!isValidConvVec(src))
	throw InvalidParameterException(EXCEPT_MSG("Invalid vector elements after type casts"));
      checkTS(TSARG); src.checkTS(TSARG);
      int limit=src.size();
      if (n<limit || !isMask()) ensureCapacity(limit);
      
      ArrayUtils<T>::convert(buff,src.getBuffPtr_INT(),limit,getStep_INT(),
			     src.getStep_INT(),getMindex_INT(),
			     src.getMindex_INT());
    }

    // Insert, remove, swap

    /**
     * Insert element 'a' 'num' times into vector at position 'pos'. The
     * elements at and after 'pos' are moved to the right. If 'pos'=='n' or
     * ==-1, 'a' is appended to the end. The buffer size is increased if
     * required.
     *
     * @param a   Element to insert
     * @param num Number of times. Def.: 1
     * @param pos Insert position. Def.: Append to end
     */
    virtual void insert(T a,int num=1,int pos=-1);

    /**
     * Insert vector 'vec' into vector at position 'pos'. The elements at and
     * after 'pos' are moved to the right. If 'pos'=='n' or ==-1, 'vec' is
     * appended to the end. The buffer size is increased if required.
     *
     * @param vec Vector to insert
     * @param pos Insert position. Def.: Append to end
     */
    virtual void insert(const BaseVector<T>& vec,int pos=-1);

    /**
     * Removes element at position 'pos'. The elements right of 'pos' are
     * moved to the left. The buffer size remains unchanged.
     * If 'num' is given, 'num' elements from 'pos' to the right are removed.
     * If 'pos'+'num' is beyond the vector end, all elements from 'pos' to the
     * end are removed. If 'num'==-1, the same happens.
     *
     * @param pos Removal position
     * @param num See above. Def.: 1
     */
    virtual void remove(int pos,int num=1);

    /**
     * Moves block of size 'len' starting at position 'oldpos' to new position
     * 'newpos'. The elements in between are moved towards 'oldpos'. Specif.
     * the outcome is the same as for
     *   // copy part from oldpos,...,oldpos+len-1 -> vec
     *   remove(oldpos,len);
     *   insert(vec,newpos);
     * but this method is usually more efficient than this combination. The
     * vector must not change size.
     * <p>
     * NOTE: Use this method instead of s.th. like
     *   a(rng1)=a(rng2),
     * which is not allowed!
     *
     * @param oldpos Old position
     * @param newpos New position
     * @param len    Optional. Def.: 1
     */
    virtual void move(int oldpos,int newpos,int len=1);

    /**
     * Swaps elements at positions 'apos','bpos'
     *
     * @param apos S.a.
     * @param bpos S.a.
     */
    virtual void swap(int apos,int bpos);

    /**
     * Permutes elements of this vector (in place)
     * <p>
     * Element at pos. i is written to pos. perm[i]. If 'inv'==true, the
     * inverse permutation is used, i.e. elem. at pos. perm[i] is written to
     * pos. i.
     * ATTENTION: 'perm' must be a permutation of the range of the vector,
     * otherwise the method fails. Failure may or may not be detected.
     * Uses temp. bool array of size of this vector.
     *
     * @param perm   S.a.
     * @param inv    S.a. Def.: false
     */
    virtual void permute(const BaseVector<int>& perm,bool inv=false);

    /*
     * Vectorisation methods (apply func. object, accumulators)
     * NOTE: If 'doCheckValid' returns true, a result vector is first
     * written into a temp. 'BaseVector<T>' before being written back to
     * this vector using 'assignInt' (which forces element checks).
     */

    /**
     * Computes nullary function once for each position in this vector,
     * outputs overwrite vector elements (start to end).
     * NOTE: Makes sense only if some internal state exists which makes the
     * function non-constant (example: PRN generators).
     *
     * @param func Nullary function object
     * @param sz   New vector size. Optional
     */
    template<class T2> void apply0(const NullaryFunc<T2>& func,int sz=-1) {
      checkTS(TSARG);
      if (sz==-1) sz=n;
      else if (sz<=0)
	throw InvalidParameterException(EXCEPT_MSG("sz"));
      if (sz>0) {
	if (!doCheckValid()) {
	  ensureCapacity(sz);
	  ArrayUtils<T>::applyNulFunc(buff,n,func,getStep_INT(),
				      getMindex_INT());
	} else {
	  // To force an element check during the final assignment, we have to
	  // build the vector in a 'BaseVector' object (which does not come
	  // with any element guarantees!)
	  BaseVector<T> tempVec(sz);
	  ArrayUtils<T>::applyNulFunc(tempVec.buff,sz,func);
	  assignInt(tempVec); // checks elements
	}
      }
    }

    /**
     * Applies a unary function T2 -> T (given by 'func') to each
     * element of 'a'. The result is returned in this vector.
     * NOTE: 'a' and this vector may be the same object.
     *
     * @param a    Source vector (elements in T2). Def.: this vector
     * @param func Unary T2 -> T function
     */
    template<class UnOp,class T2> void
    apply1(const BaseVector<T2>& a,const UnOp& func) {
      checkTS(TSARG); a.checkTS(TSARG);
      if (!doCheckValid()) {
	ensureCapacity(a.size());
	ArrayUtils<T>::applyFunc(buff,a.getBuffPtr_INT(),n,func,
				 getStep_INT(),a.getStep_INT(),getMindex_INT(),
				 a.getMindex_INT());
      } else {
	// To force an element check during the final assignment, we have to
	// build the vector in a 'BaseVector' object (which does not come with
	// any element guarantees!)
	BaseVector<T> tempVec(a.size());
	ArrayUtils<T>::applyFunc(tempVec.buff,a.getBuffPtr_INT(),
				 a.size(),func,1,a.getStep_INT(),0,
				 a.getMindex_INT());
	assignInt(tempVec); // checks elements
      }
    }

    template<class UnOp> void apply1(const UnOp& func) {
      checkTS(TSARG);
      ArrayUtils<T>::applyFunc(buff,buff,n,func,getStep_INT(),getStep_INT(),
			       getMindex_INT(),getMindex_INT());
      if (doCheckValid() && !isValidVec(*this))
	throw InvalidParameterException(EXCEPT_MSG("Invalid elements in vector after applic. of 'func'!"));
    }

    /**
     * Applies a binary function (T2,T3) -> T (given by 'func') to each
     * pair of elements from 'a','b' (must have same length). The result is
     * returned in this vector.
     * NOTE: 'a','b' and this vector may be the same object.
     *
     * @param a    Source vector. Def.: this vector
     * @param b    Source vector
     * @param func Function (T2,T3) -> T
     */
    template<class BinOp,class T2,class T3> void
    apply2(const BaseVector<T2>& a,const BaseVector<T3>& b,const BinOp& func) {
      checkTS(TSARG); a.checkTS(TSARG); b.checkTS(TSARG);
      int newn=a.size();
      if (newn!=b.size())
	throw DimMismatchException(EXCEPT_MSG("a,b must have same length"));
      if (!doCheckValid()) {
	ensureCapacity(newn);
	ArrayUtils<T>::applyBinFunc(buff,a.getBuffPtr_INT(),
				    b.getBuffPtr_INT(),n,func,
				    getStep_INT(),a.getStep_INT(),
				    b.getStep_INT(),getMindex_INT(),
				    a.getMindex_INT(),b.getMindex_INT());
      } else {
	// To force an element check during the final assignment, we have to
	// build the vector in a 'BaseVector' object (which does not come with
	// any element guarantees!)
	BaseVector<T> tempVec(newn);
	ArrayUtils<T>::applyBinFunc(tempVec.buff,a.getBuffPtr_INT(),
				    b.getBuffPtr_INT(),newn,func,1,
				    a.getStep_INT(),b.getStep_INT(),0,
				    a.getMindex_INT(),b.getMindex_INT());
	assignInt(tempVec); // checks elements
      }
    }

    template<class BinOp,class T2> void
    apply2(const BaseVector<T2>& b,
	   const BinOp& func) {
      checkTS(TSARG); b.checkTS(TSARG);
      if (n!=b.size())
	throw DimMismatchException(EXCEPT_MSG("a,b must have same length"));
      ArrayUtils<T>::applyBinFunc(buff,buff,b.getBuffPtr_INT(),n,func,
				  getStep_INT(),getStep_INT(),b.getStep_INT(),
				  getMindex_INT(),getMindex_INT(),
				  b.getMindex_INT());
      if (doCheckValid() && !isValidVec(*this))
	throw InvalidParameterException(EXCEPT_MSG("Invalid elements in vector after applic. of 'func'!"));
    }

    /**
     * Applies accumulator object 'acc' to each element (start to end) of
     * this vector, then returns result. If 'func' is given, 'acc' is applied
     * to each element of 'func' applied to this vector.
     * NOTE: 'acc' is NOT reset.
     * NOTE: 'func' is somewhat redundant, in that we could configure 'acc'
     * with 'compose22(g,func)' instead of g.
     *
     * @param acc  Accumulator object
     * @param func Optional. Def.: identity
     */
    template<class T1,class T2> T2 accumulate(const AccumulFunc<T1,T2>& acc)
      const {
      checkTS(TSARG);
      ArrayUtils<T>::applyAcc(buff,n,acc,getStep_INT(),getMindex_INT());

      return acc.get();
    }

    template<class UnOp,class T1,class T2> T2
    accumulate(const AccumulFunc<T1,T2>& acc,const UnOp& func) const {
      checkTS(TSARG);
      ArrayUtils<T>::applyAccFunc(buff,n,func,acc,getStep_INT(),
				  getMindex_INT());

      return acc.get();
    }

    /**
     * Same as 'accumulate' above, but 'func' is binary and fed by elements
     * of this vector and 'vec' (must have same size).
     * NOTE: 'acc' is NOT reset.
     *
     * @param acc  Accumulator object
     * @param func Binary function
     * @param vec  Argument vector
     */
    template<class BinOp,class T1,class T2,class T3> T2
    accumulate(const AccumulFunc<T1,T2>& acc,const BinOp& func,
	       const BaseVector<T3>& vec)
      const {
      checkTS(TSARG);
      if (vec.size()!=n)
	throw WrongDimensionException(EXCEPT_MSG(""));
      ArrayUtils<T>::applyAccBinFunc(buff,vec.getBuffPtr_INT(),n,func,acc,
				     getStep_INT(),vec.getStep_INT(),
				     getMindex_INT(),vec.getMindex_INT());
      return acc.get();
    }

    // Methods requiring order predicates (find, sort, min, ...)

    /**
     * If 'pos' is given, the position of the maximizer is returned there (if
     * there are more than one maximizer, the position of the leftmost is
     * returned).
     * NOTE: T must have the '>' operator implemented.
     *
     * @param pos Optional. See above. Def.: 0
     * @return    Maximum element
     */
    virtual T max(int* pos=0) const {
      if (pos==0) {
	AccumMaxVal<T> accu;
	return accumulate(accu);
      } else {
	AccumMax<T> accu;
	pair<T,int> res(accumulate(accu));
	*pos=res.second;
	return res.first;
      }
    }

    /**
     * If 'pos' is given, the position of the minimizer is returned there (if
     * there are more than one minimizer, the position of the leftmost is
     * returned).
     * NOTE: T must implement the '<' operator.
     *
     * @param pos Optional. See above. Def.: 0
     * @return    Minimum element
     */
    virtual T min(int* pos=0) const {
      if (pos==0) {
	AccumMinVal<T> accu;
	return accumulate(accu);
      } else {
	AccumMin<T> accu;
	pair<T,int> res(accumulate(accu));
	*pos=res.second;
	return res.first;
      }
    }

    /**
     * Sort elements in vector. If 'asc'==true, we sort in ascending, otherwise
     * in descending order. If 'index'!=0, 'index' returns the original
     * positions of the elements occur. in the sorted vector, i.e.
     *   sorted[i]=orig[index[i]].
     * <p>
     * NOTE: If this is an indexed or linear (non-flat) mask vector, we first
     * draw a flat copy.
     *
     * @param asc   Sort in ascending order (default)? False: desc. order
     * @param index Optional. See above
     */
    virtual void sort(bool asc=true,BaseVector<int>* index=0);

    /**
     * Same as 'sort', but here, 'index' is mandatory and has to have the
     * same length as this vector (or be longer). Each time elements in this
     * vector are exchanged, the same operation is done on 'index'.
     * NOTE: 'sort' (with 'index') is a special case of this method with
     * T2==int, namely 'index' is first init. with the identity.
     * <p>
     * NOTE: Uses 'ArrayUtils<T>::sortInd', which draws a flat copy to pair
     * elements of this vector with elements of 'index'.
     *
     * @param asc   Sort in ascending order?
     * @param index See above
     */
    template<class T2> void sortAlt(bool asc,BaseVector<T2>& index) {
      if (index.size()<n)
	throw InvalidParameterException("'index' has wrong size");
      if (n>0) {
	checkTS(TSARG); index.checkTS(TSARG);
	// Do not have to check 'index' elems: changed just by a permutation
	ArrayUtils<T>::sortInd(buff,(T2*) index.getBuffPtr_INT(),n,
			       getStep_INT(),index.getStep_INT(),
			       getMindex_INT(),index.getMindex_INT(),asc);
      }
    }

    /**
     * Finds pos. of first elem. for which 'pred' is true. If such an elem. is
     * not found, we return -1.
     *
     * @param pred Unary predicate
     * @return     S.a.
     */
    template<class UnOp> int findfst(const UnOp& pred) const {
      typedef typename UnOp::argument_type Arg;
      int i;
      checkTS(TSARG);
      if (!isMaskInd()) {
	const T* tempr=buff;
	for (i=0; i<n; i++,tempr+=step)
	  if (pred((Arg) *tempr)) return i;
      } else {
	const int* iB=mindex.p();
	for (i=0; i<n; i++)
	  if (pred((Arg) buff[*(iB++)])) return i;
      }
      return -1;
    }

    /**
     * Same as 'findfst' above, but 'pred' is binary p(x,y), with x from this,
     * y from 'vec' (same size).
     *
     * @param pred Binary predicate
     * @param vec  Argument vector
     * @return     S.a.
     */
    template<class BinOp,class T2> int
    findfst(const BinOp& pred,const BaseVector<T2>& vec) const {
      typedef typename BinOp::first_argument_type A1;
      typedef typename BinOp::second_argument_type A2;
      int i;
      checkTS(TSARG); vec.checkTS(TSARG);
      if (!isMaskInd()) {
	const T* tP=buff;
	const T2* vP=vec.getBuffPtr_INT();
	for (i=0; i<n; i++,tP+=step,vP+=getStep_INT())
	  if (pred((A1) *tP,(A2) *vP)) return i;
      } else {
	for (i=0; i<n; i++)
	  if (pred((A1) operator[](i),(A2) vec[i])) return i;
      }
      return -1;
    }

    /**
     * 'findfst' for the indicator function [X == 'elem'].
     *
     * @param elem See above
     * @return     "
     */
    int findfst(const T& elem) const {
      return findfst(bind2nd(equal_to<T>(),elem));
    }

    /**
     * Special case of 'findfst' in which this vector is assumed to be
     * ordered ('asc'==true: ascending order; desc. otherwise). The pos.
     * returned is one pos. for which the elem. is =='elem', or -1 if all
     * elem. !='elem'. If there are several, one of them is picked.
     * NOTE: If the vector is not sorted, the outcome is undefined.
     * <p>
     * 'high' is optional and ret. only if 'elem' is not found. It is the
     * smallest pos. for which elem. > 'elem' ('asc'==true), or the smallest
     * for which elem. < 'elem' ('asc'==false); -1 if this ineq. does not
     * hold for any elem.
     * NOTE: Passing '*high' to 'insert' keeps the vector sorted.
     *
     * @param elem Element value to find
     * @param asc  S.a. Def.: true
     * @param high S.a. Optional. Only if ret. is -1
     * @return     S.a.
     */
    virtual int findbin(const T& elem,bool asc=true,int* high=0) const;

    /**
     * Finds positions of all elements for which 'pred' is true, returns
     * the corr. position index (newly created here). The main usage is to
     * access these elements via masking:
     *   vec(Range(vec.find(...)))
     * The returned vector is sorted in asc. order. If 'off' is given,
     * positions are processed starting from 'off'.
     *
     * @param pred Unary pred.
     * @param off  Optional. Nonneg. Def.: 0
     * @return     S.a.
     */
    template<class UnOp> BaseVector<int> find(const UnOp& pred,int off=0)
      const {
      typedef typename UnOp::argument_type Arg;
      int i;
      checkTS(TSARG);
      if (off<0) throw InvalidParameterException(EXCEPT_MSG("off"));
      // Use 'BaseVecWrapper' to avoid copying upon return
      BaseVecWrapper<int> resWrap(new BaseVector<int>());
      BaseVector<int>* res=(BaseVector<int>*) resWrap.getRep();
      if (!isMaskInd()) {
	const T* tempr=buff;
	for (i=off; i<n; i++,tempr+=step)
	  if (pred((Arg) *tempr)) res->insert(i);
      } else {
	for (i=off; i<n; i++)
	  if (pred((Arg) operator[](i))) res->insert(i);
      }

      return resWrap;
    }

    /**
     * Same as 'find' above, but 'pred' is binary p(x,y), with x from this,
     * y from 'vec' (same size).
     *
     * @param pred Binary pred.
     * @param vec  Arg. vector
     * @param off  Optional. Nonneg. Def.: 0
     * @return     S.a.
     */
    template<class BinOp,class T2> BaseVector<int>
    find(const BinOp& pred,const BaseVector<T2>& vec,int off=0) const {
      typedef typename BinOp::first_argument_type A1;
      typedef typename BinOp::second_argument_type A2;
      int i;
      checkTS(TSARG); vec.checkTS(TSARG);
      if (n!=vec.size())
	throw DimMismatchException(EXCEPT_MSG(""));
      if (off<0) throw InvalidParameterException(EXCEPT_MSG("off"));
      // Use 'BaseVecWrapper' to avoid copying upon return
      BaseVecWrapper<int> resWrap(new BaseVector<int>());
      BaseVector<int>* res=(BaseVector<int>*) resWrap.getRep();
      if (!isMaskInd() && !vec.isMaskInd()) {
	const T* tP=buff;
	const T2* vP=vec.getBuffPtr_INT();
	for (i=off; i<n; i++,tP+=step,vP+=vec.getStep_INT())
	  if (pred((A1) *tP,(A2) *vP)) res->insert(i);
      } else {
	for (i=off; i<n; i++)
	  if (pred((A1) operator[](i),(A2) vec[i])) res->insert(i);
      }

      return resWrap;
    }

    /**
     * Counts how many elem. of this vector make unary predicate 'pred' true.
     *
     * @param pred Unary pred.
     * @return     S.a.
     */
    template<class UnOp> int count(const UnOp& pred) const {
      // NOTE: Cast bool -> int done automatically
      return accumulate(accum_fun(std::plus<int>(),0),pred);
    }

    /**
     * Version of 'count' with binary predicate. Entries of this vector and
     * of 'vec' become 1st and 2nd args.
     *
     * @param pred Binary pred.
     * @param vec  Argument vector
     * @return     S.a.
     */
    template<class BinOp,class T2> int
    count(const BinOp& pred,const BaseVector<T2>& vec) const {
      // NOTE: Cast bool -> int done automatically
      return accumulate(accum_fun(std::plus<int>(),0),pred,vec);
    }

    /**
     * Special case of 'findfst' for intervals. Checks whether all elements
     * fall within the interval given by 'ival' (see class 'Interval'). If
     * not, then the pos. of the first element outside of 'ival' can be
     * returned via 'pos', and a status can be ret. via 'stat': 1 -> elem.
     * too small, 2 -> elem. too large.
     * <p>
     * NOTE: Req. '<' and '==' operators.
     * NOTE: If 'stat' is not required, 'findfst(ival,...)' can be used
     * instead.
     *
     * @param ival Interval range
     * @param pos  Optional. See above
     * @param stat "
     * @return     All elements within interval?
     */
    virtual bool checkBounds(const Interval<T>& ival,int* pos=0,int* stat=0)
      const;

    /**
     * Compares this vector with 'a'. Position of first discrep. returned.
     * If the vectors are identical and have the same length, -1 is returned.
     * If one is the prefix of the other, the length of the shorter is
     * returned.
     *
     * @param a Other vector
     * @return  S.a.
     */
    virtual int compare(const BaseVector<T>& a) const;

    /**
     * Comparison operator
     *
     * @param a Other vector
     * @return  Are this vector and 'a' identical?
     */
    virtual bool operator==(const BaseVector<T>& a) const;

    // Methods operating on sorted vectors

    /**
     * Merges sorted vectors (asc. order if 'asc'==true, desc. otherwise)
     * 'a', 'b' into sorted vector written into this.
     * ATTENTION: 'a', 'b' must be different object from this vector.
     * NOTE: Does not check whether args are sorted!
     *
     * @param a   Source vector (sorted)
     * @param b   Source vector (sorted)
     * @param asc S.a. Def.: true
     */
    virtual void merge(const BaseVector<T>& a,const BaseVector<T>&b,
		       bool asc=true);

    /**
     * Writes elements occuring in 'a', but not in 'b', into this vector.
     * 'a', 'b' must be sorted (asc. order if 'asc'==true, desc. otherwise).
     * This vector is sorted in the same way.
     * ATTENTION: 'a', 'b', this vector must be different objects!
     * If 'saveMem'==false, the buffer of this vector is enlarged to the
     * size of 'a' at the beginning, otherwise single elements are inserted
     * using 'insert'.
     * NOTE: Does not check whether args are sorted!
     *
     * @param a       Source vector (sorted)
     * @param b       Source vector (sorted)
     * @param asc     S.a. Def.: true
     * @param saveMem S.a. Def.: false
     */
    virtual void difference(const BaseVector<T>& a,const BaseVector<T>&b,
			    bool asc=true,bool saveMem=false);

    /**
     * Checks whether all elements of 'a' are contained in this vector. Both
     * 'a' and this vector must be sorted (asc. order if 'asc'==true, desc.
     * otherwise).
     * NOTE: Does not check whether args are sorted!
     *
     * @param a   Source vector
     * @param asc S.a. Def.: true
     * @return    Is 'a' a subset?
     */
    virtual bool isSubset(const BaseVector<T>& a,bool asc=true) const;

    /**
     * Checks whether this vector contains any duplicate entries.
     * Requires comparison operators for T.
     * NOTE: We first assume that the entries are already sorted in ascending
     * order. If they are, the check is O(n). If they are not, a copy is drawn
     * and sorted, then the check repeated on there: O(n log(n)).
     * NOTE: Given that the elements are unique and the method returns true,
     * additional information can be returned: min. value via 'minVal', max.
     * value via 'maxVal'. These are not written if the method returns false.
     *
     * @param minVal See above. Optional
     * @param maxVal "
     * @return       Is every element unique within the vector?
     */
    virtual bool areElemUnique(T* minVal=0,T* maxVal=0) const;

    // Initialization methods

    /**
     * Initializes vector with element 'a'. With no further argument, the
     * original dimensions are kept
     *
     * @param num Optional. Number of entries
     * @param a   Optional. Fill element. Def.: def. fill value
     */
    virtual void fill() {
      checkTS(TSARG);
      ArrayUtils<T>::fill(buff,this->defFill,n,getStep_INT(),getMindex_INT());
    }

    virtual void fill(int num,T a);

    virtual void fill(int num) {
      fill(num,this->defFill);
    }

    // Persistance, IO

    /**
     * Writes repres. of object into (binary) file stream 'os'.
     * Generic format (FF version 1):
     * - tag @BaseVector (not written if 'noTag'==true)
     * - FF version [int(4)]
     * - type code [int(4)]
     * - vector [using 'saveInt']
     * Change this by specialisation if not appropriate!
     * Files in FFV 0 do not have the byte size of T stored (see
     * 'saveInt'). For these files, the system byte size is used,
     * which can lead to mistakes when loading.
     *
     * @param os    Output file stream
     * @param noTag If true, then no tag is written
     */
    virtual ofstream& save(ofstream& os,bool noTag=false) const;

    /**
     * Read repres. of object from (binary) file stream 'is'.
     * NOTE: Can also load files from old 'OpenVector'
     * format, same as this one, but has no type field (type is
     * assumed to be correct here). The tag is @OpenVector.
     * If 'noTag'==true, no tag is read (format must be
     * 'BaseVector' then).
     * <p>
     * DOWNW. COMPAT.: Files written by earlier LHOTSE versions had
     * wrong type code fields ('typeOther' for elementary T). If
     * this is encountered, a warning message is printed, rather
     * than throwing an exception (see 'checkTypeCode').
     *
     * @param is    Input file stream
     * @param noTag If true, no tag is read
     */
    virtual ifstream& load(ifstream& is,bool noTag=false);

    /**
     * If for a special T, the vector cannot be printed, overwrite this
     * method and also 'print'.
     *
     * @return Does 'print' work for this class?
     */
    virtual bool canPrint() const {
      return true;
    }

    /**
     * Prints vector (as row vector) into output stream 'os'. There is no
     * newline after the row.
     * NOTE: To print as column vector, mask it as matrix:
     *   BaseMatrix<T>::mask(*this)->print(os);
     * <p>
     * NOTE: Requires the '<<' operator to be implemented for T.
     * Override if the def. implementation does not make sense!
     */
    virtual void print(ostream& os) const;

#ifdef MATLAB_MEX
    /**
     * Print to stdout as row vector (MEX file).
     * Requires that T can be cast to double.
     */
    virtual void matlabPrint() const {
      for (int i=0; i<n-1; i++)
	mexPrintf("%f ",(double) operator[](i));
      mexPrintf("%f",(double) operator[](n-1));
    }
#endif

    // Other (specialized) methods

    /**
     * Specialized method for selection and accumulation.
     * 'imap' has 2*k-1 int entries, ie. k consec. tupels i(s),j(s), with
     * j(k-1) not given. Works by running a counter pos(s), with pos(0)=0.
     * For s=0,...,k-1:
     *   this[pos(s)] += a[i(s)],
     *   pos(s+1) = pos(s) + j(s)
     * This vector must have suff. size, and be init. properly.
     * 'a' and this vector can be the same object.
     * NOTE: Cannot use 'BaseVector<int>' for 'imap'. Use 'getFlatBuff'.
     *
     * @param a    Source vector
     * @param imap S.a.
     */
    virtual void accuMap(const BaseVector<T>& a,const ArrayHandle<int>& imap);

    // Methods for 'TempMatMethods' superclass

    /*
     * 'assignVirtual' must do the same as 'operator=' !
     */

    void assignVirtual(const TempMatMethods* arg) {
      const BaseVector<T>* aptr=DYNCAST(const BaseVector<T>,arg);
      if (aptr==0) throw InternalException(EXCEPT_MSG(""));
      //cout << "BaseVector-AV: this=" << this << ",arg=" << arg << endl;
      assignInt(*aptr);
    }

    void assignVirtual(const TempMatMethods_AssignElement_Generic* arg) {
      //cout << "BaseVector::assignVirtual" << endl;
      const TempMatMethods_AssignElement<T>* aptr=
	DYNCAST(const TempMatMethods_AssignElement<T>,arg);
      if (aptr==0) throw InternalException(EXCEPT_MSG(""));
      //cout << "BaseVector-AV: this=" << this << ",arg=" << aptr->elem
      //   << endl;
      assignInt(aptr->elem);
    }

    // Methods for 'WriteBackVec' superclass

    void writeBackBuff(const ArrayHandle<T>& buff) {
      Handle<BaseVector<T> > temp(newEmpty());
      temp->reassign(buff);
      assignInt(*temp);
    }

    /*
     * "Internal" methods which are declared public:
     * ATTENTION: These methods must only be called by methods of this class
     * or subclasses!
     * Why not protected? Template methods have to call them for
     * 'BaseVector<T2>' with T2!=T, and it seems we cannot be friend of
     * 'BaseVector<T2>' for arb. T2 (SOLUTION???).
     */

    /**
     * INTERNAL METHOD! DO NOT USE!
     *
     * @return Buffer pointer
     */
    T* getBuffPtr_INT() const {
      return buff;
    }

    /**
     * INTERNAL METHOD! DO NOT USE!
     *
     * @return Value of 'mindex' if this is an indexed mask vector; 0 otherwise
     */
    const ArrayHandle<int>& getMindex_INT() const {
      return isMaskInd()?mindex:ArrayHandleZero<int>::get();
    }

    /**
     * INTERNAL METHOD! DO NOT USE!
     *
     * @return Value of step size 'step' if this is not an indexed mask; 0
     *         otherwise
     */
    int getStep_INT() const {
      return isMaskInd()?0:step;
    }

  protected:
    /*
     * Internal methods
     * NOTE: Some internal methods do NOT call 'checkTS' and do not call
     *       'incrTS' when changing the vector size. Check comments
     *       ("TS not served").
     */

    /**
     * Called by "copy" constructor from 'BaseVecWrapper'. If things change
     * in subclasses, this method has to be overridden!
     * 'arg' has wrapped a temp. vector. We initialise all fields of this
     * object with the fields of this vector. Then, the fields of the temp.
     * vector are reset, so that 'arg' is disposed of.
     * <p>
     * NOTE: This method must be called only for an "empty" vector (either
     * from a constructor or after 'dealloc' has been used)!
     * <p>
     * NOTE: The difference to a normal CC is not just that a physical copy
     * is not done, but also that the mask status of the object wrapped
     * by 'arg' is transferred to this vector, while a normal CC always
     * creates a normal vector.
     *
     * @param arg Wrapped "source" object
     */
    virtual void convertWrapped(BaseVecWrapper<T>& arg);

    /**
     * Draws a copy of this object. The return object is always normal.
     * Same as the copy constructor, but polymorphic.
     * NOTE: Dyn. created on local heap, wrap into handle!
     * ==> Subclasses have to overwrite this method!
     *
     * @return Copy
     */
    virtual BaseVector<T>* copy() const {
      return new BaseVector<T>(*this);
    }

    /**
     * Creates an empty vector of the type of this class. Used in virtual
     * methods. NOTE: Dyn. created on local heap, wrap into handle!
     * ==> Subclasses have to overwrite this method!
     *
     * @return See above
     */
    virtual BaseVector<T>* newEmpty() const {
      return new BaseVector<T>();
    }

    /**
     * Does job of assignment operator =. In contrast to 'operator=',
     * this is a virtual method.
     * <p>
     * NOTE: Overlapping buffers. If both this vector and 'vec' are flat,
     * this is guaranteed to work, otherwise there can be mistakes!
     *
     * @param vec Source vector
     */
    virtual void assignInt(const BaseVector<T>& vec);

    /**
     * Does job of assignment operator =. In contrast to 'operator=',
     * this is a virtual method.
     *
     * @param elem Element to fill vector with
     */
    virtual void assignInt(T elem) {
      if (!isValidElement(elem))
	throw InvalidParameterException(EXCEPT_MSG("'elem' invalid element"));
      fill(std::max(n,1),elem);
    }

    /**
     * See header comment. If this method returns true, then all methods which
     * modify the vector, will use the 'isValidXXX' methods to check
     * validity.
     * <p>
     * ATTENTION: This def. implementation for unconstrained element types
     * has to be overwritten in subclasses with constraints!
     *
     * @return See above
     */
    virtual bool doCheckValid() const {
      return false;
    }

    /**
     * See header comment and 'MatDefMembers'.
     * <p>
     * ATTENTION: This def. implementation for unconstrained element types
     * has to be overwritten in subclasses with constraints!
     */
    bool isValidElement(const T& elem) const {
      return true;
    }

    /**
     * See header comment. Returns true if 'vecP' points to an object whose
     * type enforces the required validity constraints. In this case, checking
     * each element is not necessary for 'vecP'.
     * NOTE: Shortcut for elementwise checking, which is done if this method
     * returns false.
     *
     * @param vecP See above
     * @return     "
     */
    virtual bool isValidVecType(const BaseVector<T>* vecP) const {
      return true;
    }

    /**
     * See header comment. Returns true iff the given vector is s.t. it
     * could be inserted into this vector. The vector to be tested is given
     * via buffer point 'vbuff', vector size 'vn' and range 'vrng'.
     * <p>
     * NOTE: The def. impl. just calls 'isValidElement' for each element and
     * is probably sufficient in all cases.
     * NOTE: Should be called AFTER the validity of the range 'vrng' has been
     * checked!
     *
     * @param vbuff See above
     * @param vn    ". Used only if 'vrng' is open
     * @param vrng  ". Def.: full
     * @return      "
     */
    virtual bool isValidVec(const T* vbuff,int vn,
			    const Range& vrng=RangeFull::get()) const;

    /**
     * See other 'isValidVec' version.
     * NOTE: If 'doCheck' is false (def.), we call 'isValidVecType' first.
     * If this returns true, the vector is autom. valid.: no element checks
     * are done. If 'doCheck'==true, element checks are always done.
     * <p>
     * TS not served.
     *
     * @param vec     Vector to check
     * @param doCheck S.a. Def.: false
     * @return        See above
     */
    virtual bool isValidVec(const BaseVector<T>& vec,bool doCheck=false) const;

    /**
     * Same as 'isValidVec' above, but here every element is first cast from
     * T2 to T before it is checked.
     *
     * @param vbuff See above
     * @param vn    ". Used only if 'vrng' is open
     * @param vrng  ". Def.: full
     * @return      "
     */
    template<class T2> bool
    isValidConvVec(const T2* vbuff,int vn,const Range& vrng=RangeFull::get())
    const {
      int i,sz=vrng.size(vn);
      for (i=0; i<sz; i++)
	if (!isValidElement((T) vbuff[vrng[i]])) return false;

      return true;
    }

    /**
     * See other 'isValidConvVec' version.
     * TS not served.
     *
     * @param vec Vector to check
     * @return    See above
     */
    template<class T2> bool isValidConvVec(const BaseVector<T2>& vec) const {
      int pn=vec.size();
      if (pn==0)
	return true;
      else if (vec.isMaskInd())
	return isValidConvVec(vec.getBuffPtr_INT(),pn,
			      Range(vec.getMaskIndex()));
      else {
	int vstep=vec.getStep_INT();
	if (vstep>0)
	  return isValidConvVec(vec.getBuffPtr_INT(),pn,
				Range(0,(pn-1)*vstep,vstep));
	else
	  return isValidConvVec(vec.getBuffPtr_INT()+(pn-1)*vstep,pn,
				Range((1-pn)*vstep,0,vstep));
      }
    }

    /*
     * Required by 'MatTimeStamp' superclass.
     */
    const MatTimeStamp* getBaseObj() const {
      return isMask()?baseObj:0;
    }

    /**
     * Generic code for 'operator()', r-value version.
     *
     * @param rng See operator()
     */
    virtual BaseVector<T>* subrvalInt(const Range& rng) const;

    /**
     * Generic code for 'operator()', l-value version. This vector is
     * expanded if necessary.
     *
     * @param rng See operator()
     */
    virtual BaseVector<T>* sublvalInt(const Range& rng);

    /**
     * Writes repres. of object into (binary) file stream 'os'.
     * Format of default implementation:
     * - length n [int(4)]
     * - byte size of T [int(4)]: Only if 'storeBSize'==true!
     * - content [using 'NumberFormats<T>::save']
     * Override this by specialisation if not appropriate!
     * <p>
     * NOTE: No tag is written here, and no type info.
     *
     * @param os         Output file stream
     * @param storeBSize S.a. Def.: true
     */
    virtual ofstream& saveInt(ofstream& os,bool storeBSize=true) const;

    /**
     * Read repres. of object from (binary) file stream 'is'.
     * NOTE: No tag is read here, and no type info.
     * NOTE: If 'doCheckValid' is true, the elements are checked after
     * loading. This is because many subclasses use the 'BaseVector'
     * format even though they implement element checks.
     * If 'loadBSize'==true, the byte size of T is loaded from file (see
     * 'saveInt'), otherwise the system byte size is used, which can lead
     * to mistakes if the file was produced on a diff. system.
     *
     * @param is        Input file stream
     * @param loadBSize S.a. Def.: true
     */
    virtual ifstream& loadInt(ifstream& is,bool loadBSize=true);

    /**
     * Allocates memory for vector and sets all entries to 'a'. If buffer
     * is used, it is deallocated first. The buffer size is 'num' exactly.
     * <p>
     * NOTE: TS not served. 'a' is checked for validity.
     *
     * @param num Number of elements
     * @param a   Fill element. Def.: def. fill value
     */
    virtual void init(int num) {
      init(num,this->defFill);
    }

    virtual void init(int num,T a);

    /**
     * Allocates memory for vector. If buffer is used, it is deallocated
     * first. The entries are not initialized.
     * If 'num'==0, we only dealloc. the buffer
     * <p>
     * NOTE: TS not served.
     *
     * @param num Number of elements. Must be non-negative
     */
    virtual void alloc(int num,uchar debStat);

    /**
     * Ensures that the buffer has a given size (or is larger). If the buffer
     * is less than the given size, it's reallocated with the given size. In
     * this case, the old elements are lost. The elements are not initialized.
     * The vector size is set to 'num' exactly.
     * <p>
     * NOTE: The timestamp is increased iff 'num' is different from the current
     * length.
     *
     * @param num Minimum number of entries. Must be non-negative
     */
    virtual void ensureCapacity(int num);

    /**
     * Deallocates resources of a vector.
     * If this vector is assoc. with a buffer, the watcher's ref. counter
     * is decreased. This leads to de-alloc. of the vector buffer iff the
     * counter drops to 0. This is done calling 'MatTimeStamp::deassocBuff'.
     * <p>
     * NOTE: TS not served.
     *
     * @param noReset If true, we just call 'deassocBuff' for buffer de-assoc.,
     *                but do not change any other members. Def.: false
     */
    virtual void dealloc(bool noReset=false);

    /**
     * Does job of 'reassign', but without checking validity of the range,
     * of the arguments, or of the elements in 'pbuff'.
     *
     * @param pbuff See above
     * @param pn    "
     * @param rng   ". Optional.
     * @param watch "
     * @param noCp  "
     */
    virtual void reassignInt(const T* pbuff,int pn,
			     const Range& rng=RangeFull::get(),
			     MemWatchBase* watch=0,bool noCp=false);

    /**
     * Does job of 'reassign', but without checking validity of the range,
     * of the arguments, or of the elements in 'vec'.
     * The parent object ('baseObj') will be 'vec', or (if 'vec' is a mask)
     * the parent object of 'vec'.
     *
     * @param vec   Source vector (to refer to)
     * @param rng   Range. Def.: full range
     */
    virtual void reassignInt(const BaseVector<T>& vec,
			     const Range& rng=RangeFull::get());

    /**
     * Helper for resizing methods. Reallocates buffer to size 'num' while
     * keeping the vector content.
     *
     * @param num New buffer size
     */
    void reallocBuff(int num,uchar debStat);
  };

  // Public methods

#ifdef MATLAB_MEX
  template<class T> BaseVector<T>::BaseVector(const MatlabMatrix& mat) :
    Vector(),MatTimeStamp(),TempMatMethods(),size_n(0),buff(0),
    maskStat(statNormal),step(1)
  {
    this->setDefValues();
    if (mat.m()>0 && mat.n()>0) {
      ensureCapacity(mat.m()*mat.n());
      ArrayUtils<T>::convert(buff,mat.buff(),n);
    }
  }
#endif

  /*
   * Either the validity of each element of vec(rng) has to be checked or
   * not. The latter case is true iff
   * - doCheckValid returns false, or
   * - isValidVecType(&vec) returns true (restriction is already enforced
   *   by 'vec's type)
   * In the former case, we first create the desired mask temporarily using
   * 'reassignInt', then check it using 'isValidVec'. Only if this is OK
   * will this vector be changed.
   */
  template<class T> void BaseVector<T>::reassign(const BaseVector<T>& vec,
						 const Range& rng)
  {
    if (!isMaskOrEmpty())
      throw MaskObjectException("Has to be mask or empty");
    if (rng.checkRange(vec.n)) throw OutOfRangeException(EXCEPT_MSG(""));
    if (vec.n==0)
      dealloc(); // empty vector
    else {
      if (!doCheckValid() || isValidVecType(&vec)) {
	// No elem. checks necessary
	reassignInt(vec,rng);
      } else {
	// Have to check elems:
	// We use 'BaseVecWrapper' here to first create the mask vector to be
	// assigned to this vector. If the element check fails, this vector is
	// not changed.
	// 'tempVec' has type 'BaseVector'. 'isValidVec' calls
	// 'isValidVecType' internally, but this will result in false, since
	// 'BaseVector' does not enforce any restrictions.
	BaseVecWrapper<T> tvWrap(new BaseVector<T>());
	BaseVector<T>* tempVec=tvWrap.getRep();
	tempVec->reassignInt(vec,rng); // no elem. check here
	// Check elements:
	if (!isValidVec(*tempVec))
	  throw InvalidParameterException("'vec' contains invalid elements");
	// Use fields of 'tempVec' for this vector
	dealloc(); // remove old stuff
	convertWrapped(tvWrap); // take over fields, clean-up
      }
    }
  }

  template<class T> inline ArrayHandle<T> BaseVector<T>::getFlatBuff() const
  {
    checkTS(TSARG);
    if (n==0)
      return ArrayHandleZero<T>::get(); // zero handle
    else if (isFlat()) {
      // Handle with same mem. watcher (internal constructor)
      // If 'buffWatch'==0, the handle does not own its buffer and will
      // not touch it
      return ArrayHandle<T>(buff,n,buffWatch);
    } else {
      // Draw flat copy
      ArrayHandle<T> flatCp(n);
      ArrayUtils<T>::copy(flatCp,buff,n,1,step,0,getMindex_INT());
      return flatCp;
    }
  }

  template<class T> inline void BaseVector<T>::getLinVec(BaseLinVec<T>& blVec,
							 bool writeBack,
							 bool negFort,
							 bool doFlat) const
  {
    if (n==0) throw WrongDimensionException("Vector must not be empty");
    if (!isMaskInd() && (!doFlat || step==1)) {
      // Flat or linear vector
      if (step>0) {
	ArrayHandle<T> hand(buff,(n-1)*step+1,buffWatch); // int. constr.
	blVec.init(n,hand,step,negFort); // write-back automatic
      } else {
	// Buffer pointer in BLAS (Fortran) style
	ArrayHandle<T> hand(buff+((n-1)*step),(1-n)*step+1,buffWatch);
	blVec.init(n,hand,step,negFort); // write-back automatic
      }
    } else {
      // Draw flat copy, step size = 1. Write-back dep. on 'writeBack'
      blVec.init(n,getFlatBuff(),1,negFort,
		 writeBack?((BaseVector<T>*) this):0);
    }
  }

  template<class T> inline const T& BaseVector<T>::operator[](int pos) const
  {
    if (pos<0 || pos>=n) throw OutOfRangeException(EXCEPT_MSG(""));
    checkTS(TSARG);
    if (maskStat==statNormal || maskStat==statMaskFlat) return buff[pos];
    else if (maskStat==statMaskLin) return *(buff+(pos*step));
    else return buff[mindex[pos]];
  }

  template<class T> inline T& BaseVector<T>::operator[](int pos)
  {
    if (pos<0 || pos>=n) throw OutOfRangeException(EXCEPT_MSG(""));
    checkTS(TSARG);
    if (maskStat==statNormal || maskStat==statMaskFlat) return buff[pos];
    else if (maskStat==statMaskLin) return *(buff+(pos*step));
    else return buff[mindex[pos]];
  }

  template<class T> inline void BaseVector<T>::set(int pos,T elem)
  {
    if (pos<0 || pos>=n) throw WrongDimensionException(EXCEPT_MSG(""));
    if (doCheckValid() && !isValidElement(elem))
      throw InvalidParameterException("Invalid vector element");
    checkTS(TSARG);
    if (maskStat==statNormal || maskStat==statMaskFlat) buff[pos]=elem;
    else if (maskStat==statMaskLin) *(buff+(pos*step))=elem;
    else buff[mindex[pos]]=elem;
  }

  template<class T> inline void BaseVector<T>::resize(int num)
  {
    if (!isMask() && size_n!=num) {
      init(num);
      incrTS(); // buffer changed
    } else
      fill(num,this->defFill);
  }

  template<class T> inline void BaseVector<T>::resizeSave(int num)
  {
    if (num!=n) {
      if (isMask()) throw MaskObjectException(EXCEPT_MSG(""));
      if (n>num) throw WrongDimensionException(EXCEPT_MSG(""));
      if (size_n!=num)
	reallocBuff(num,1); // re-allocate buffer to size 'num'
      // Set "margin entries" to 'this->defFill'
      ArrayUtils<T>::fill(buff+n,this->defFill,num-n);
      n=num;
      incrTS(); // size has changed
    }
  }

  template<class T> inline void BaseVector<T>::expand(int num)
  {
    if (num<n) throw WrongDimensionException(EXCEPT_MSG("num"));
    if (num!=n) {
      if (isMask()) throw MaskObjectException(EXCEPT_MSG(""));
      if (size_n<num)
	// Re-alloc. buffer
	reallocBuff(this->getNewBuffSize(num,size_n),2);
      // Set "margin entries" to def. fill value
      ArrayUtils<T>::fill(buff+n,this->defFill,num-n);
      n=num;
      incrTS(); // size has changed
    }
  }

  template<class T> inline void BaseVector<T>::shrink(int num)
  {
    if (num<=0) throw InvalidParameterException(EXCEPT_MSG("num"));
    if (num>n) expand(num); // expand instead of shrink!
    else if (num<n) {
      if (isMask()) throw MaskObjectException(EXCEPT_MSG(""));
      n=num; // shrink size
      incrTS(); // size change
    }
  }

  template<class T> inline void BaseVector<T>::exchange(BaseVector<T>& a)
  {
    checkTS(TSARG); a.checkTS(TSARG);
    if (n!=a.n) throw WrongDimensionException(EXCEPT_MSG(""));
    if (this!=&a)
      ArrayUtils<T>::swap(buff,(T*) a.getBuffPtr_INT(),n,getStep_INT(),
			  a.getStep_INT(),getMindex_INT(),a.getMindex_INT());
  }

  // If the buffer has to be enlarged, we do not use 'expand' in order to
  // save copying twice
  template<class T> inline void BaseVector<T>::insert(T a,int num,
						      int pos)
  {
    if (num<0) throw InvalidParameterException("num");
    else if (num>0) {
      if (!isValidElement(a))
	throw InvalidParameterException(EXCEPT_MSG("Invalid vector element"));
      if (isMask()) throw MaskObjectException(EXCEPT_MSG(""));
      if (pos==-1) pos=n; // append at end
      else if (pos<0 || pos>n)
	throw OutOfRangeException(EXCEPT_MSG("pos"));
      // Make space
      if (n+num<=size_n) {
	// Just copy
	ArrayUtils<T>::copy(buff+(pos+num),buff+pos,n-pos);
      } else {
	// Expand and copy (see 'reallocBuff')
	int newSize=this->getNewBuffSize(n+num,size_n);
	T* oldBuff=buff; buff=0;
	int oldN=n; n=0;
	MemWatchBase* oldWatch=buffWatch; buffWatch=0;
	size_n=0; // looks like empty to 'alloc'
	alloc(newSize,3); // expand buffer size
	n=oldN;
	if (oldBuff!=0) {
	  // Copy
	  ArrayUtils<T>::copy(buff,oldBuff,pos);
	  ArrayUtils<T>::copy(buff+(pos+num),oldBuff+pos,n-pos);
	  // De-assoc. with old buffer
	  MemWatchBase* newWatch=buffWatch;
	  buffWatch=oldWatch; // looks like old to 'dealloc'
	  dealloc(true);
	  buffWatch=newWatch; // restore watcher
	}
      }
      n+=num; // new vector size
      ArrayUtils<T>::fill(buff+pos,a,num);
      incrTS(); // size changed
    }
  }

  // Do not use 'expand' to avoid copying twice
  template<class T> void BaseVector<T>::insert(const BaseVector<T>& vec,
					       int pos)
  {
    int num=vec.size();
    if (num>0) {
      vec.checkTS(TSARG);
      if (doCheckValid() && !isValidVec(vec))
	throw InvalidParameterException(EXCEPT_MSG("'vec' has invalid elements"));
      if (isMask()) throw MaskObjectException(EXCEPT_MSG(""));
      if (pos==-1) pos=n; // append at end
      if (pos<0 || pos>n) throw OutOfRangeException(EXCEPT_MSG("pos"));
      // Make space
      if (n+num<=size_n) {
	// Just copy
	ArrayUtils<T>::copy(buff+(pos+num),buff+pos,n-pos);
      } else {
	// Expand and copy
	int newSize=this->getNewBuffSize(n+num,size_n);
	T* oldBuff=buff; buff=0;
	int oldN=n; n=0;
	MemWatchBase* oldWatch=buffWatch; buffWatch=0;
	size_n=0; // looks like empty to 'alloc'
	alloc(newSize,4); // expand buffer
	n=oldN;
	if (oldBuff!=0) {
	  // Copy
	  ArrayUtils<T>::copy(buff,oldBuff,pos);
	  ArrayUtils<T>::copy(buff+(pos+num),oldBuff+pos,(n-pos));
	  // De-assoc. with old buffer
	  MemWatchBase* newWatch=buffWatch;
	  buffWatch=oldWatch; // looks like old to 'dealloc'
	  dealloc(true);
	  buffWatch=newWatch; // restore watcher
	}
      }
      n+=num; // new vector size
      Handle<BaseVector<T> > maskVec(newEmpty());
      maskVec->reassignInt(buff+pos,num,RangeFull::get(),buffWatch);
      maskVec->assignInt(vec); // copy 'vec'
      incrTS(); // size changed
    }
  }

  template<class T> inline void BaseVector<T>::remove(int pos,int num)
  {
    if (isMask()) throw MaskObjectException(EXCEPT_MSG(""));
    if (pos<0 || pos>=n) throw OutOfRangeException(EXCEPT_MSG("pos"));
    if (num==-1 || num+pos>=n)
      n=pos;
    else if (num>0) {
      ArrayUtils<T>::copy(buff+pos,buff+(pos+num),n-pos-num);
      n-=num;
    }
    incrTS(); // size has changed
  }

  template<class T> inline void BaseVector<T>::move(int oldpos,int newpos,
						    int len)
  {
    if (oldpos==newpos) return;
    if (len<0) throw InvalidParameterException("len");
    else if (len>0) {
      if (oldpos<0 || oldpos+len>n || newpos<0 || newpos+len>n)
	throw OutOfRangeException(EXCEPT_MSG("oldpos,newpos"));
      checkTS(TSARG);
      // Copy of part to move
      Handle<BaseVector<T> > maskVec(newEmpty());
      maskVec->reassign(*this,Range(oldpos,oldpos+len-1));
      Handle<BaseVector<T> > tempVec(maskVec->copy());
      int bLen=newpos-oldpos;
      if (!isMaskInd()) {
	if (bLen>0)
	  // Start to end
	  ArrayUtils<T>::copy(buff+oldpos*step,buff+(oldpos+len)*step,bLen,
			      step,step);
	else
	  // End to start
	  ArrayUtils<T>::copy(buff+(oldpos+len-1)*step,buff+(oldpos-1)*step,
			      -bLen,-step,-step);
      } else {
	int i;
	if (bLen>0) {
	  for (i=oldpos; i<newpos; i++)
	    operator[](i)=operator[](i+len);
	} else {
	  for (i=oldpos-1; i>=newpos; i--)
	    operator[](i+len)=operator[](i);
	}
      }
      maskVec->reassign(*this,Range(newpos,newpos+len-1));
      maskVec->assignInt(*tempVec); // copy back
    }
  }

  template<class T> inline void BaseVector<T>::swap(int apos,int bpos)
  {
    checkTS(TSARG);
    if (apos<0 || apos>=n || bpos<0 || bpos>=n)
      throw OutOfRangeException(EXCEPT_MSG(""));
    T* ap,*bp;
    if (!isMaskInd()) {
      ap=buff+(apos*step); bp=buff+(bpos*step);
    } else {
      ap=buff+mindex[apos]; bp=buff+mindex[bpos];
    }
    T temp=*ap;
    *ap=*bp; *bp=temp;
  }

  template<class T> inline void
  BaseVector<T>::permute(const BaseVector<int>& perm,bool inv)
  {
    int numDone,targ,start,prev;
    T val,swp;
    T* ptr;

    if (n==0) return; // empty vector
    checkTS(TSARG); perm.checkTS(TSARG);
    if (perm.size()!=n)
      throw InvalidParameterException(EXCEPT_MSG("perm"));
    if (!perm.checkBounds(Interval<int>(0,n-1,IntVal::ivClosed,
					IntVal::ivClosed)))
      throw InvalidParameterException(EXCEPT_MSG("perm"));
    BaseVector<bool> tickoff(n,false);
    start=0; numDone=0;
    for (;;) {
      // Start new cycle with 'start'
      val=operator[](start); numDone++;
      targ=perm[start];
      if (!inv) {
	while (targ!=start) {
	  swp=operator[](targ);
	  set(targ,val); val=swp;
	  tickoff[targ]=true; numDone++;
	  targ=perm[targ];
	  // Sanity check (avoids infinite cycles if 'perm' is not a perm.)
	  if (numDone>n)
	    throw InvalidParameterException(EXCEPT_MSG("'perm' not a permutation"));
	}
	set(start,val); // close cycle
      } else {
	prev=start;
	while (targ!=start) {
	  set(prev,operator[](targ));
	  tickoff[targ]=true; numDone++;
	  prev=targ; targ=perm[targ];
	  // Sanity check (avoids infinite cycles if 'perm' is not a perm.)
	  if (numDone>n)
	    throw InvalidParameterException(EXCEPT_MSG("'perm' not a permutation"));
	}
	set(prev,val); // close cycle
      }
      if (numDone==n) break; // leave loop
      // Find next starting point
      for (start++; start<n && tickoff[start]; start++);
      // Sanity check
      if (start==n)
	throw InvalidParameterException(EXCEPT_MSG("'perm' not a permutation"));
    }
  }

  template<class T> inline void
  BaseVector<T>::sort(bool asc,BaseVector<int>* index)
  {
    if (n>0) {
      checkTS(TSARG);
      if (index!=0) {
	index->checkTS(TSARG);
	index->fill(n);
	for (int i=0; i<n; i++) (*index)[i]=i;
	sortAlt(asc,*index);
      } else {
	BaseLinVec<T> flatb;
	getLinVec(flatb,true,true,true); // need flat vector
	ArrayUtils<T>::sort(flatb.getBuff(),n,asc); // sort
      }
    }
  }

  template<class T> inline int BaseVector<T>::findbin(const T& elem,bool asc,
						      int* high) const
  {
    int up,lo,i;
    T act;

    checkTS(TSARG);
    if (n==0) {
      if (high!=0) *high=-1;
      return -1;
    }
    act=operator[](0);
    if (elem==act) return 0;
    else if (asc^(elem>act)) {
      if (high!=0) *high=0;
      return -1;
    }
    act=operator[](n-1);
    if (elem==act) return n-1;
    else if (asc^(elem<act)) {
      if (high!=0) *high=-1;
      return -1;
    }
    up=n-1; lo=0;
    while (lo<up-1) {
      // Loop invariant: lo<up-1, x[lo] < elem < x[up] (for 'asc'==true).
      // For the chosen i: lo<i<up, so there is progress in each iter.
      i=(up+lo)/2; act=operator[](i);
      if (act==elem) return i;
      else if (asc^(act>elem)) lo=i;
      else up=i;
    }

    if (high!=0) *high=up;
    return -1;
  }

  template<class T> inline bool
  BaseVector<T>::checkBounds(const Interval<T>& ival,int* pos,int* stat) const
  {
    int i,ret;

    checkTS(TSARG);
    if (!isMaskInd()) {
      const T* tempr=buff;
      for (i=0; i<n; i++,tempr+=step)
	if ((ret=ival.check(*tempr))!=0) {
	  if (pos!=0) *pos=i;
	  if (stat!=0) *stat=ret;
	  return false;
	}
    } else {
      const int* iB=mindex;
      for (i=0; i<n; i++)
	if ((ret=ival.check(buff[*(iB++)]))!=0) {
	  if (pos!=0) *pos=i;
	  if (stat!=0) *stat=ret;
	  return false;
	}
    }
    return true;
  }

  template<class T> inline
  int BaseVector<T>::compare(const BaseVector<T>& a) const
  {
    checkTS(TSARG); a.checkTS(TSARG);
    int i,bnd=std::min(n,a.size());
    if (!isMaskInd() && !a.isMaskInd()) {
      const T* tempr=buff,*tempa=a.getBuffPtr_INT();
      for (i=0; i<bnd; i++,tempr+=step,tempa+=a.getStep_INT())
	if (*tempr!=*tempa) return i;
    } else {
      for (i=0; i<bnd; i++)
	if (operator[](i)!=a[i]) return i;
    }

    return (n==a.size())?(-1):bnd;
  }

  template<class T> inline bool
  BaseVector<T>::operator==(const BaseVector<T>& a) const
  {
    checkTS(TSARG); a.checkTS(TSARG);
    return (n!=a.size())?false:(compare(a)==-1);
  }

  template<class T> inline void
  BaseVector<T>::accuMap(const BaseVector<T>& a,const ArrayHandle<int>& imap)
  {
    if (imap.size()>0) {
      checkTS(TSARG); a.checkTS(TSARG);
      if (imap.size()%2!=1)
	throw InvalidParameterException(EXCEPT_MSG(""));
      if (!isValidVec(a))
	throw InvalidParameterException(EXCEPT_MSG(""));
      ArrayUtils<T>::accuMap(buff,a.getBuffPtr_INT(),imap.p(),
			     (imap.size()+1)/2,step,a.getStep_INT(),
			     getMindex_INT(),a.getMindex_INT());
    }
  }

  template<class T> inline void
  BaseVector<T>::merge(const BaseVector<T>& a,const BaseVector<T>&b,bool asc)
  {
    int apos=0,bpos=0,an=a.size(),bn=b.size(),pos=0;
    T aelem=a[0],belem=b[0];

    if (isMask()) throw MaskObjectException(EXCEPT_MSG(""));
    a.checkTS(TSARG); b.checkTS(TSARG);
    if (doCheckValid() && (!isValidVec(a) || !isValidVec(b)))
      throw InvalidParameterException(EXCEPT_MSG(""));
    ensureCapacity(an+bn);
    for (;;) {
      if (asc^(aelem>belem)) {
	set(pos++,aelem);
	if (++apos<an) aelem=a[apos];
	else break;
      } else {
	set(pos++,belem);
	if (++bpos<bn) belem=b[bpos];
	else break;
      }
    }
    for (; apos<an; apos++) set(pos++,a[apos]);
    for (; bpos<bn; bpos++) set(pos++,b[bpos]);
    n=pos;
  }

  template<class T> inline void
  BaseVector<T>::difference(const BaseVector<T>& a,const BaseVector<T>&b,
			    bool asc,bool saveMem)
  {
    if (isMask())
      throw MaskObjectException(EXCEPT_MSG(""));
    a.checkTS(TSARG); b.checkTS(TSARG);
    if (a.size()==0)
      ensureCapacity(0);
    else if (b.size()==0)
      operator=(a);
    else {
      int apos=0,bpos=0,an=a.size(),bn=b.size(),pos=0;
      T aelem=a[0],belem=b[0];
      if (!saveMem)
	ensureCapacity(an);
      else
	ensureCapacity(0);
      for (;;) {
	if (aelem==belem) {
	  if (++apos==an) break;
	  else aelem=a[apos];
	} else if (asc^(aelem>belem)) {
	  if (!saveMem) set(pos,aelem);
	  else insert(aelem);
	  pos++;
	  if (++apos==an) break;
	  else aelem=a[apos];
	} else {
	  if (++bpos==bn) break;
	  else belem=b[bpos];
	}
      }
      if (apos<an) {
	Handle<BaseVector<T> > smsk(a.subrvalInt(Range(apos))); // rest of 'a'
	if (!saveMem) {
	  Handle<BaseVector<T> > tmsk(sublvalInt(Range(pos,pos+an-apos-1)));
	  *tmsk=*smsk;
	} else
	  insert(*smsk);
	pos+=(an-apos);
      }
      n=pos;
    }
  }

  template<class T> inline bool
  BaseVector<T>::isSubset(const BaseVector<T>& a,bool asc) const
  {
    checkTS(TSARG); a.checkTS(TSARG);
    if (n==0) return false;
    else if (a.size()==0) return true;
    else {
      int apos=0,pos=0,an=a.size();
      T aelem=a[0],elem=operator[](0);
      // Invariant: 'elem' <= 'aelem' (for ascending)
      for (;;) {
	if (elem==aelem) {
	  if (++apos==an) return true; // is subset
	  else aelem=a[apos];
	} else if (asc^(elem>aelem)) {
	  if (++pos==n) return false; // not a subset
	  else elem=operator[](pos);
	} else
	  return false; // not a subset
      }
    }
  }

  // In the first run, we assume the vector is sorted in asc. order. If not,
  // we bail out, draw a copy, sort it and run again.
  template<class T> inline bool BaseVector<T>::areElemUnique(T* minVal,
							     T* maxVal) const
  {
    int i;
    bool runFlag=false;
    Handle<BaseVector<T> > maskV(newEmpty());
    Handle<BaseVector<T> > flatVec(newEmpty());
    T elem,prev;

    checkTS(TSARG);
    if (n==0) return true; // empty vector
    maskV->reassign(*this);
    for (;;) {
      prev=(*maskV)[0];
      for (i=1; i<n; i++) {
	elem=(*maskV)[i];
	if (elem==prev) return false; // duplicate
	else if (elem<prev) break; // not sorted
      }
      if (i==n) {
	// OK, no duplicates
	if (minVal!=0) *minVal=(*maskV)[0];
	if (maxVal!=0) *maxVal=(*maskV)[n-1];
	return true;
      } else {
	if (runFlag)
	  throw InternalException("BaseVector::areElemUnique: BUG!");
	// Draw copy, sort in asc. order
	*flatVec=*this;
	flatVec->sort();
	maskV->reassign(*flatVec);
	runFlag=true;
      }
    }
  }

  template<class T> inline void BaseVector<T>::fill(int num,T a)
  {
    if (num<0) throw WrongDimensionException(EXCEPT_MSG(""));
    if (num>0 && !isValidElement(a))
      throw InvalidParameterException("Invalid vector element");
    if (num!=n) ensureCapacity(num); // deals with masks, timestamps
    else checkTS(TSARG);
    ArrayUtils<T>::fill(buff,a,num,getStep_INT(),getMindex_INT());
  }

  template<class T> inline ofstream& BaseVector<T>::save(ofstream& os,
							 bool noTag) const
  {
    int i;
    checkTS(TSARG);
    if (!noTag)
      FileUtils::saveHeader(os,"@BaseVector",1);
    else {
      i=1; NumberFormats<int>::save(os,&i,1,1,0,4);
    }
    i=getTypeCode(this->defFill);
    NumberFormats<int>::save(os,&i,1,1,0,4);
    return saveInt(os);
  }

  template<class T> inline ifstream& BaseVector<T>::load(ifstream& is,
							 bool noTag)
  {
    int i,ffv;
    bool isBV=true,loadBSize=false;
    checkTS(TSARG);
    if (!noTag) {
      ArrayHandle<string> tags(2);
      tags[0]="@BaseVector";
      tags[1]="@OpenVector";
      ffv=FileUtils::loadHeaderMulti(is,tags,i);
      isBV=(i==0);
    } else
      NumberFormats<int>::load(is,&ffv,1,1,0,4);
    if (ffv!=1 && (!isBV || ffv!=0))
      throw FileFormatException(EXCEPT_MSG("Unknown FF version number"));
    if (isBV) {
      // Read type code. If it is 'typeOther', a warning is displayed inst.
      // of throwing an exception (so that old files can be read)
      NumberFormats<int>::load(is,&i,1,1,0,4);
      if (i!=getTypeCode(this->defFill)) {
	ArrayHandle<char> msg(101);
	if (i==TypeCodeConsts::typeOther) {
	  sprintf(msg.p(),"BaseVector::load: File has wrong type code %d (required: %d). Continuing...",i,getTypeCode(this->defFill));
	  printMsgStdout(msg.p());
	} else {
	  sprintf(msg.p(),"BaseVector::load: File has wrong type code %d (required: %d)",i,getTypeCode(this->defFill));
	  throw WrongTypeException(msg.p());
	}
      }
      loadBSize=(ffv==1);
    }
    return loadInt(is,loadBSize);
  }

  template<class T> inline void BaseVector<T>::print(ostream& os) const
  {
    if (n>0) {
      checkTS(TSARG);
      for (int i=0; i<n-1; i++)
	os << operator[](i) << " ";
      os << operator[](n-1);
    }
  }

  // Internal methods

  template<class T> inline void
  BaseVector<T>::convertWrapped(BaseVecWrapper<T>& arg) {
    // Take over fields from 'arg'
    BaseVector<T>* vec=arg.getRep(); // dyn. cast in subclasses!
    n=vec->n;
    copyMembers(vec); // copy 'MatTimeStamp' stuff, deassoc. 'vec'
    buff=vec->buff;
    maskStat=vec->maskStat;
    if (!vec->isMask()) size_n=vec->size_n;
    else baseObj=vec->baseObj;
    if (vec->isMaskInd()) mindex=vec->mindex;
    else step=vec->step;
    this->copyDefValues(*vec); // def. values
    // Dispose of 'arg' via 'dealloc'. This is safe because 'copyMembers'
    // above set 'vec->buffWatch' to 0, so de-assoc. in 'dealloc' does not
    // have effect
    vec->dealloc();
    delete vec; // safe now: vector is empty
    arg.reset(); // remove repres. in 'arg'
  }

  template<class T> inline void
  BaseVector<T>::assignInt(const BaseVector<T>& vec)
  {
    if (this!=&vec) {
      checkTS(TSARG); vec.checkTS(TSARG);
      if (doCheckValid() && !isValidVec(vec))
	throw InvalidParameterException(EXCEPT_MSG("'vec' contains invalid elements"));
      if (n<vec.n || !isMask()) ensureCapacity(vec.n);
      ArrayUtils<T>::copy(buff,vec.getBuffPtr_INT(),vec.n,step,
			  vec.getStep_INT(),getMindex_INT(),
			  vec.getMindex_INT());
    }
  }

  template<class T> inline bool
  BaseVector<T>::isValidVec(const T* vbuff,int vn,const Range& vrng) const {
    int i,sz=vrng.size(vn);
    for (i=0; i<sz; i++)
      if (!isValidElement(vbuff[vrng[i]])) return false;

    return true;
  }

  template<class T> inline bool
  BaseVector<T>::isValidVec(const BaseVector<T>& vec,bool doCheck) const
  {
    if ((!doCheck && isValidVecType(&vec)) || vec.n==0)
      return true;
    else if (vec.isMaskInd())
      return isValidVec(vec.getBuffPtr_INT(),vec.size(),
			Range(vec.getMindex_INT()));
    else {
      int vstep=vec.getStep_INT();
      if (vstep>0)
	return isValidVec(vec.getBuffPtr_INT(),vec.size(),
			  Range(0,(vec.size()-1)*vstep,vstep));
      else
	return isValidVec(vec.getBuffPtr_INT()+(vec.size()-1)*vstep,vec.size(),
			  Range((1-vec.size())*vstep,0,vstep));
    }
  }

  template<class T> inline BaseVector<T>*
  BaseVector<T>::subrvalInt(const Range& rng) const {
    checkTS(TSARG);
    BaseVector<T>* vec=newEmpty();
    vec->reassign(*this,rng);

    return vec;
  }

  template<class T> inline BaseVector<T>*
  BaseVector<T>::sublvalInt(const Range& rng) {
    int newl;

    checkTS(TSARG);
    if (rng.isOpen())
      // open range:
      newl=std::max(rng.getStart()+1,n);
    else
      newl=rng.getMaxPos(0)+1;
    if (newl>n) expand(newl); // expand to new size
    BaseVector<T>* vec=newEmpty();
    vec->reassignInt(*this,rng);

    return vec;
  }

  template<class T> inline ofstream&
  BaseVector<T>::saveInt(ofstream& os,bool storeBSize) const {
    checkTS(TSARG);
    NumberFormats<int>::save(os,&n,1,1,0,4);
    if (storeBSize) {
      int i=sizeof(T);
      NumberFormats<int>::save(os,&i,1,1,0,4);
    }
    NumberFormats<T>::save(os,buff,n,getStep_INT(),getMindex_INT());
    return os;
  }

  template<class T> inline ifstream&
  BaseVector<T>::loadInt(ifstream& is,bool loadBSize) {
    int newn,fsize;

    checkTS(TSARG);
    NumberFormats<int>::load(is,&newn,1,1,0,4);
    if (newn<0)
      throw FileFormatException(EXCEPT_MSG("Negative vector length"));
    ensureCapacity(newn);
    if (!loadBSize) fsize=sizeof(T);
    else {
      NumberFormats<int>::load(is,&fsize,1,1,0,4);
      if (fsize<1)
	throw FileFormatException(EXCEPT_MSG("Invalid byte size"));
    }
    NumberFormats<T>::load(is,buff,n,getStep_INT(),getMindex_INT(),fsize);
    // Force element checks
    if (doCheckValid() && !isValidVec(*this,true))
      throw FileFormatException(EXCEPT_MSG("Entries in file violate vector element constraints!"));

    return is;
  }

  template<class T> inline void BaseVector<T>::init(int num,T a)
  {
    if (num>0 && !isValidElement(a))
      throw InvalidParameterException("'a' invalid element");
    alloc(num,5);
    ArrayUtils<T>::fill(buff,a,num);
  }

  /*
   * The memory alloc. is done in the 'MemWatcher' constructor. Makes
   * sure that the alloc. mech. is compatible with the de-alloc. mech.
   * used in the destructor there!
   */
  template<class T> inline void BaseVector<T>::alloc(int num,uchar debStat)
  {
    if (num<0) throw WrongDimensionException(EXCEPT_MSG(""));
    if (isMask()) throw MaskObjectException(EXCEPT_MSG(""));
    dealloc();
    if (num>0) {
      n=size_n=num;
      MemWatcher<T>* watch=new MemWatcher<T>(n,0,debStat); // allocation
      assocBuff(watch); // association
      buff=watch->getBuff(); // buffer pointer
    }
  }

  template<class T> inline void BaseVector<T>::ensureCapacity(int num)
  {
    if (num!=n) {
      if (isMask()) throw MaskObjectException(EXCEPT_MSG(""));
      if (size_n<num)
	alloc(num,6); // re-alloc
      else {
	n=num;
	incrTS(); // size change
      }
    }
  }

  template<class T> inline void BaseVector<T>::dealloc(bool noReset)
  {
    deassocBuff(); // remove assoc. with buffer
    if (!noReset) {
      mindex.changeRep(0);
      n=size_n=0;
      buff=0; maskStat=statNormal; step=1;
    }
  }

  template<class T> inline void
  BaseVector<T>::reassignInt(const T* pbuff,int pn,const Range& rng,
			     MemWatchBase* watch,bool noCp)
  {
    dealloc(); // transform to empty vector
    if (pn>0) {
      if (rng.isLiteralRange()) {
	// Literal range
	maskStat=rng.isFlatRange()?statMaskFlat:statMaskLin;
	n=rng.size(pn);
	buff=((T*) pbuff)+rng.getStart();
	step=rng.getStep();
      } else {
	// Indexed range
	maskStat=statMaskInd;
	n=rng.size(pn);
	buff=(T*) pbuff;
	if (!noCp)
	  mindex.copy(rng.getIndex()); // draw copy
	else
	  mindex=rng.getIndex(); // do not copy, share
      }
      assocBuff(watch); // assoc. with buffer
      baseObj=0; // refers to a memory region, not a matrix/vector object
    }
  }

  template<class T> inline void
  BaseVector<T>::reassignInt(const BaseVector<T>& vec,const Range& rng)
  {
    int i,vn=vec.size();

    dealloc(); // transform to empty vector
    if (vn>0) {
#ifdef DEBUG_TRACKHANDLES
      debugCause=6; // DEBUG
#endif
      if (vec.isFlat())
	// 'vec' is flat
	// 'baseObj', TS is dealt with below
	reassignInt(vec.getBuffPtr_INT(),vn,rng,vec.buffWatch);
      else {
	if (vec.maskStat==statMaskLin) {
	  // 'vec' is linear, but not flat
	  int vstep=vec.getStep_INT();
	  if (rng.isLiteralRange()) {
	    // Mask vector will be linear
	    maskStat=statMaskLin;
	    buff=vec.getBuffPtr_INT()+rng.getStart()*vstep;
	    step=vstep*rng.getStep();
	    if (step==1) maskStat=statMaskFlat; // can happen, e.g. both -1
	    n=rng.size(vn);
	  } else {
	    // Range is indexed: will be indexed mask vector
	    const ArrayHandle<int>& rind=rng.getIndex();
	    maskStat=statMaskInd;
	    n=rng.size(vn);
	    mindex.changeRep(n); // alloc. new index
	    if (vstep>0) {
	      // 'vec' has pos. step size, but >1
	      buff=vec.getBuffPtr_INT();
	      for (i=0; i<n; i++) mindex[i]=rind[i]*vstep;
	    } else {
	      // 'vec' has neg. step size. Have to set 'buff' s.t. the new
	      // 'mindex' does not have negative entries!
	      buff=vec.getBuffPtr_INT()+vstep*(vn-1);
	      for (i=0; i<n; i++)
		mindex[i]=(rind[i]-vn+1)*vec.step;
	    }
	  }
	} else {
	  // 'vec' is an indexed mask vector. This means that this vector will
	  // be one as well. We try to use 'vec.mindex', which we can do iff
	  // 'rng' is flat and starts with 0 (if 'rng' is the full range, this
	  // is true).
	  maskStat=statMaskInd;
	  n=rng.size(vn);
	  if (rng.isFlatRange() && rng.getStart()==0) {
	    // Just copy handle for 'mindex'
	    mindex=vec.mindex;
	  } else {
	    // Construct new 'mindex' by mapping
	    mindex.changeRep(n);
	    rng.mapIndex(vec.getMindex_INT(),vn,mindex,false);
	  }
	  buff=vec.getBuffPtr_INT();
	}
	// Assoc. with 'vec' buffer
	assocBuff(vec.buffWatch); 
      }
      // Base object pointer and timestamp: copy from 'vec'
      if (vec.isMask())
	// NOTE: It is important that we copy 'vec's timestamp value, not the
	// one of the object underlying 'vec', because if 'vec' is an invalid
	// mask, so should be this mask vector.
	baseObj=vec.baseObj; // same base object as 'vec'
      else
	baseObj=&vec; // hook onto 'vec' (which is normal)
      setTS(vec.getTS()); // 'vec's current timestamp value
    }
  }

  template<class T> inline void
  BaseVector<T>::reallocBuff(int num,uchar debStat)
  {
    if (num<n) throw WrongDimensionException(EXCEPT_MSG(""));
    // Store old members, make vector look empty
    if (buff==0 && buffWatch!=0) // sanity check
      throw InternalException(EXCEPT_MSG(""));
    T* oldBuff=buff; buff=0;
    int oldN=n; n=0;
    MemWatchBase* oldWatch=buffWatch; buffWatch=0;
    int oldSizeN=size_n; size_n=0; // looks empty to 'alloc' now
    alloc(num,debStat);
    n=oldN;
    if (oldBuff!=0) {
      ArrayUtils<T>::copy(buff,oldBuff,n); // copy
      // De-assoc. with old buffer
      T* newBuff=buff; buff=oldBuff;
      MemWatchBase* newWatch=buffWatch;
      buffWatch=oldWatch; // looks like old to 'dealloc'
      size_n=oldSizeN;
      dealloc(true);
      buffWatch=newWatch; // restore watcher
      buff=newBuff;
      n=oldN; size_n=num;
    }
  }

  // Definition of constants

  template<class T> const unsigned char BaseVector<T>::statNormal     ;
  template<class T> const unsigned char BaseVector<T>::statMaskFlat   ;
  template<class T> const unsigned char BaseVector<T>::statMaskLin    ;
  template<class T> const unsigned char BaseVector<T>::statMaskInd    ;

  // Global functions

  template<class T> inline ostream&
  operator<<(ostream& os,const BaseVector<T>& vec)
  {
    if (!vec.canPrint())
      throw NotImplemException(EXCEPT_MSG("'print' not implemented"));
    vec.print(os);

    return os;
  }
//ENDNS

#undef TSARG

#endif
