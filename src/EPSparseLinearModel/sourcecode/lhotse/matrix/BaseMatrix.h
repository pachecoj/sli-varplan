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
 * Desc.:  Header class BaseMatrix
 * ------------------------------------------------------------------- */

#ifndef BASEMATRIX_H
#define BASEMATRIX_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

/*
 * Storage format. Mask matrices
 *
 * In order to be able to use Fortran numerical libraries, the storage
 * format has been changed drastically from what has been used in earlier
 * LHOTSE versions. Earlier versions used row-major ordering, a C-style
 * array-of-arrays structure and some hacks to support the Numerical
 * Recipes distribution. The new format is column-major, uses a flat buffer
 * and "striding" (i.e. addition/subtraction of a fixed offset to move to the
 * next/prev. column).
 *
 * Internal storage format:
 * Distinguish between normal and mask matrices.
 * A normal m*n matrix consists of a flat data buffer of buffSz T entries.
 * Here, m<=stride, n<=buffSz/stride. If one of n,m are zero, both must be
 * zero. If one of stride, buffSz are zero, both must be zero.
 * We use column-major ordering: element (i,j) is stored in
 *   buff[j*stride + i].
 * 'stride' is also referred to as "striding constant" or LDA (leading
 * dimension) in the context of the Fortran libraries (BLAS,LAPACK).
 *
 * Empty matrix:
 * Must have m==0,n==0, but may have a buffer. A matrix without buffer must
 * have buffSz==0,stride==0.
 * If a matrix has a buffer, then buffSz>0, and also stride>0, even if the
 * matrix is empty.
 *
 * Buffer size:
 * Under normal usage of a matrix, its buffer never shrinks and is grown on
 * demand. The striding value 'stride' is kept constant if possible, even if
 * the matrix size changes.
 * The only public methods which can shrink the buffer are:
 * - 'resize': Destroys content and shrinks buffer and matrix size to given
 *   values. Elements set to given fill value.
 * - 'resizeSave': Retains old entries and brings matrix and buffer to given
 *   size, filling remaining entries with given fill value. Shrinkage below
 *   the matrix size is not allowed
 * - 'shrinkBuffer': If the buffer size is smaller than the given size, this
 *   method does nothing. Otherwise, the buffer is reallocated to the given
 *   size. The content and size of the matrix is not changed. Shrinkage below
 *   the matrix size is not permitted.
 * - 'reshape', 'reshapeSave' do not change the buffer size, but modify the
 *   striding constant 'stride'
 * NOTE: In contrast to 'resizeSave', the method 'expand' will NOT shrink the
 * buffer, but only expand it if necessary.
 * NOTE: 'resize' can be used to completely deallocate the buffer of a matrix.
 * NOTE: None of these methods can be used for mask matrices.
 *
 * Mask matrices:
 * A matrix can be normal, a flat mask or an indexed mask. A mask matrix does
 * not own its flat buffer, but uses the buffer of another normal matrix or
 * simply some flat memory region. In general, it puts a matrix structure on
 * a flat memory region relative to 'buff'. For mask matrices, 'buffSz' is
 * not used. If the mask uses the mem. region of some matrix or vector, this
 * is called the parent object. 'baseObj' stores a ptr. to the parent object.
 * Since 'buffSz', 'baseObj' are never used both, they form a union.
 * NOTE: A mask need not have a parent object. For example, it can refer to
 * some memory region.
 *
 * Flat mask matrices:
 * In a flat mask matrix, columns are still contiguous, and we can move to
 * the next/prev. column by adding/subtracting 'stride' to the pointer.
 * 'buff' points to the first element of the first column. Since all methods
 * use 'stride' to move from column to column, they do not need to distinguish
 * between a normal and a flat mask matrix.
 * NOTE: In acc. with the Fortran libraries, 'stride' cannot be negative,
 * in fact 'stride'>='m'.
 *
 * Indexed mask matrices:
 * Additionally, an indexed mask matrix has a row and column position index.
 * Element (i,j) is stored at buff[colInd[j]*stride+rowInd[i]]. 'colInd' and
 * 'rowInd' can be 0 (repres. identity), but not both at the same time.
 * Both must be non-neg., 1-1. 'rowInd' must map into 0,...,stride-1, while
 * 'colInd' can be arbitrary. If M is the max. of all 'colInd' entries, the
 * contig. region up to M*stride+stride-1 must be accessible.
 */

#include "lhotse/TypeCode.h"
#include "lhotse/FileUtils.h"
#include "lhotse/NumberFormats.h"
#include <algorithm>
#include "lhotse/matrix/default.h"
#include "lhotse/matrix/Matrix.h"
#include "lhotse/matrix/ArrayUtils.h"
#include "lhotse/matrix/BaseVector.h"
#include "lhotse/matrix/MatTimeStamp.h"
#include "lhotse/matrix/TempMatMethods.h"
#include "lhotse/matrix/WriteBackMat.h"
#include "lhotse/matrix/BaseLinMat.h"
#include "lhotse/matrix/TempBaseMatrix.h"
#include "lhotse/matrix/TransInPlace.h"
#ifdef MATLAB_MEX
#include "lhotse/matif/mex_for_cpp.h"
#endif

// Local macro for args to 'checkTS'
#define TSARG __FILE__,__LINE__

// Macros for indexed masks
#define COLPOS(i) ((colInd==0)?(i):colInd[i])
#define ROWPOS(i) ((rowInd==0)?(i):rowInd[i])

//BEGINNS(matrix)
  /**
   * Helper class to allow 'BaseMatrix' methods to return temp. created
   * vector objects by value without copying. Same as 'BaseVecWrapper' for
   * 'BaseVector', see comments there.
   * NOTE: A return-by-value is automatically restricted to be an r-value
   * in the enclosing expression. If this is not appropriate, use the
   * 'TempBaseMatrix' mechanism (e.g., 'operator()').
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  template<class T> class BaseMatWrapper
  {
  protected:
    // Members

    BaseMatrix<T>* repr; // representation

  public:
    // Methods

    /**
     * Default constructor
     *
     * @param rep Repres.
     */
    BaseMatWrapper(BaseMatrix<T>* rep) : repr(rep) {}

    /**
     * @return Representation
     */
    BaseMatrix<T>* getRep() const {
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
   * Provides standard rectangular matrix for elementary types T with support
   * for masking and default methods which can be implemented for general T
   * (numerical methods are not implemented here).
   * <p>
   * A 'BaseMatrix' can be normal or a mask. A mask matrix does not own the
   * data buffer it uses. See doc/masking.txt and the header comment above for
   * details.
   * <p>
   * General:
   * Many mechanisms are the same or similar to the ones for 'BaseVector',
   * they are not mentioned here again, see 'BaseVector' header comments.
   * <p>
   * Implementation of subclasses:
   * The following methods have to be overwritten:
   * - constructors, copy constr.
   * - copy constr. from 'BaseMatWrapper'
   * - 'mask': create temp. with target type, return corr. TempXXX
   * - 'operator()': return corr. TempXXX, dyn. cast
   * - 'diag': return corr. TempXXX, dyn. cast
   * - 'operator[]': like 'operator()'
   * - 'operator=': return type
   * - 'max','min' ('BaseVector' return): return type
   * - 'save','load'
   * - 'convertWrapped': if new attribs.
   * - 'copy', 'newEmpty', 'newEmptyVector': use corr. type
   * - 'doCheckValid', 'isValidElement', 'isValidMatType': if elem. constraints
   * - 'saveInt', 'loadInt': if new attribs.
   * - 'assignVirtual': This has to do exactly the same thing as the corr.
   *   'operator=' variant. If 'operator=' does the same thing as the
   *   superclass (except returning a diff. ref. type), 'assignVirtual' does
   *   not have to be overwritten.
   * Furthermore:
   * - if this is 'XXX', 'TempXXX' has to be implem. (see 'TempBaseMatrix')
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  template<class T> class BaseMatrix : public Matrix,public MatTimeStamp,
				       public TempMatMethods,
				       public MatDefMembers<T>,
				       public WriteBackMat<T>
  {
  public:
    /*
     * Constants
     * statNormal:   Normal matrix
     * statMaskFlat: Flat mask matrix
     * statMaskInd:  Indexed mask matrix
     */
    static const unsigned char statNormal  =0;
    static const unsigned char statMaskFlat=1;
    static const unsigned char statMaskInd =2;

  protected:
    // Additional members

    T* buff;                 // Buffer pointer
    unsigned char maskStat;  // Mask status (see 'statXXX')
    int stride;              // Striding constant (LDA). Number of buffer
                             // rows for normal matrix
    union {
      int buffSz;            // Buffer size (normal matrix)
      const MatTimeStamp* baseObj; // Parent object (mask matrix)
    };
    ArrayHandle<int> rowInd; // Row position index (indexed mask, can be 0)
    ArrayHandle<int> colInd; // Column position index (indexed mask, can be 0)

  public:
    // Constructors

    /**
     * Constructs empty matrix, no buffers.
     */
    BaseMatrix() : Matrix(),MatTimeStamp(),TempMatMethods(),buff(0),
    maskStat(statNormal),stride(0),buffSz(0) {
      this->setDefValues();
    }

    /**
     * Constructor. If 'a' is not given, the matrix is init. with the
     * def. fill value. If 'cols' is not given, 'cols'=='rows'.
     *
     * @param rows Number of rows
     * @param cols Optional. Number of columns
     * @param a    Optional. Fill value
     */
    explicit BaseMatrix(int rows,int cols=-1) :
      Matrix(),MatTimeStamp(),TempMatMethods(),buff(0),stride(0),
      maskStat(statNormal),buffSz(0) {
      this->setDefValues();
      if (cols==-1) cols=rows;
      init(rows,cols,7);
      if (buff!=0) incrTS();
    }

    BaseMatrix(int rows,int cols,const T& a) : Matrix(),MatTimeStamp(),
    TempMatMethods(),buff(0),stride(0),maskStat(statNormal),buffSz(0) {
      this->setDefValues();
      init(rows,cols,12,a);
      if (buff!=0) incrTS();
    }

    /**
     * Copy constructor.
     * Creates a normal matrix, no matter what status 'mat' has.
     * To create a mask or to copy a mask, use 'reassign'.
     *
     * @param mat Source matrix
     */
    BaseMatrix(const BaseMatrix<T>& mat) : Matrix(),MatTimeStamp(),
    TempMatMethods(),buff(0),stride(0),maskStat(statNormal),buffSz(0) {
      mat.checkTS(TSARG);
      this->copyDefValues(mat); // def. value members
      assignInt(mat); // copy
    }

    /**
     * Special "copy" constructor from 'BaseMatWrapper'. See header comment
     * for details. This constructor has to be implemented in all subclasses!
     * It creates a new matrix by casting the repres. of the wrapper to the
     * class type, using all its fields (it is a temp. mask object) for the
     * new object (method 'convertWrapped') ==> no copying involved.
     *
     * @param arg Wrapped argument
     */
    BaseMatrix(const BaseMatWrapper<T>& arg) : Matrix(),MatTimeStamp(),
    TempMatMethods() {
      // use 'arg' fields to init. ours, then dispose of 'arg'
      convertWrapped((BaseMatWrapper<T>&) arg);
    }

    /**
     * Destructor
     */
    virtual ~BaseMatrix() {
      dealloc();
    }

    // Information methods

    /**
     * @return Is this a mask matrix?
     */
    bool isMask() const {
      return (maskStat!=statNormal);
    }

    /**
     * A normal empty matrix can be transformed into a mask iff it has
     * been empty since creation. This method returns true iff this object
     * is a mask or is empty and has been since creation. The latter is
     * checked by comp. its timestamp value with 0.
     * A 'reassign' method works for this object only if this method returns
     * true.
     *
     * @return Is this object a mask or the empty matrix (without a buffer)?
     */
    bool isMaskOrEmpty() const {
      return (maskStat!=statNormal || (buff==0 && getTS()==0));
    }

    /**
     * @return Is this an indexed mask matrix?
     */
    bool isMaskInd() const {
      return (maskStat==statMaskInd);
    }

    /**
     * @return Mask status (see constants 'statXXX')
     */
    unsigned char getMaskStatus() const {
      return maskStat;
    }

    /**
     * Returns column pos. index for indexed mask matrix, 0 handle otherwise.
     *
     * @return S.a.
     */
    const ArrayHandle<int>& getMaskColIndex() const {
      return isMaskInd()?colInd:ArrayHandleZero<int>::get();
    }

    /**
     * Returns row pos. index for indexed mask matrix, 0 handle otherwise.
     *
     * @return S.a.
     */
    const ArrayHandle<int>& getMaskRowIndex() const {
      return isMaskInd()?rowInd:ArrayHandleZero<int>::get();
    }

    /**
     * This matrix is flat iff not an indexed mask and the striding value is
     * equal to m. 'getFlatBuff' draws a buffer copy iff the matrix is not
     * flat.
     *
     * @return Is this matrix flat?
     */
    virtual bool isFlat() const {
      return (maskStat!=statMaskInd && m==stride);
    }

    /**
     * Returns buffer size 'buffSz' for a normal matrix. For a mask, an exc.
     * is thrown.
     *
     * @return Buffer size
     */
    int getBufferSize() const {
      if (isMask()) throw MaskObjectException(EXCEPT_MSG(""));
      return buffSz;
    }

    // Mask matrix methods

    /**
     * Reassigns a mask matrix (or the empty matrix) to point to some buffer
     * region 'pbuff'. The size of the matrix is passed via 'pm' (rows),
     * 'pn' (cols), the striding constant as 'strid'. We must have
     * strid >= pm. The mask is a flat one.
     * <p>
     * A mem. watcher for the buffer region can be passed via 'watch'. If so,
     * it is used as 'buffWatch' here (ref. counter is increased, call to
     * 'assocBuff'). If not, there is no buffer control (dangerous!).
     * ==> Not recommended!
     * <p>
     * NOTE: This is just to interface foreign code! Use 'ArrayHandle' version
     * below whenever possible!
     * NOTE: An empty matrix can be made a mask here only if it has been empty
     * since creation!
     *
     * @param pbuff Buffer
     * @param pm    Number of rows
     * @param pn    Number of cols
     * @param strid S.a.
     * @param watch Optional. S.a.
     */
    virtual void reassign(const T* pbuff,int pm,int pn,int strid,
			  MemWatchBase* watch=0) {
      if (!isMaskOrEmpty())
	throw MaskObjectException(EXCEPT_MSG("Needs to be mask or empty!"));
      if (pm<=0 || pn<=0 || strid<pm)
	throw InvalidParameterException(EXCEPT_MSG(""));
      if (doCheckValid() && !isValidMat(pbuff,pm,pn,strid))
	throw WrongTypeException(EXCEPT_MSG("Invalid matrix elements"));
      reassignInt(pbuff,pm,pn,strid,watch);
    }

    /**
     * Same as above, but 'pbuff' is an ArrayHandle. We use the mem. watcher
     * of this handle.
     *
     * @param pbuff See above
     * @param pm    Number of rows
     * @param pn    Number of cols
     * @param strid S.a.
     */
    virtual void reassign(const ArrayHandle<T>& pbuff,int pm,int pn,
			  int strid) {
      reassign(pbuff.p(),pm,pn,strid,pbuff.getMemWatch());
    }

    /**
     * Reassigns a mask matrix (or the empty matrix) to point to some buffer
     * region 'pbuff'. The striding constant is passed as 'strid', ranges as
     * 'rngR' (row positions) and 'rngC' (col. positions). Both must be 1-1,
     * 'rngR' values in 0,...,strid-1. 'rngR' can be open (applied to
     * 0,...,strid-1), 'rngC' must be closed.
     * The ranges are applied to the hypothetical matrix with 'strid' rows
     * whose buffer starts at 'pbuff'.
     * <p>
     * The method tries to create a flat mask. This is possible if 'rngR' is
     * a flat range (step size +1) and 'rngC' is a flat or linear range with
     * pos. step size (the striding factor for the flat mask is 'strid' times
     * the step size).
     * If a flat mask cannot be created, but 'rngR' is flat, the row index
     * 'rowInd' is set to 0.
     * NOTE: If the mask is not flat, linear ranges with step size !=1 are
     * converted to indexes.
     * <p>
     * NOTE: This is just to interface foreign code! Use 'ArrayHandle' version
     * below whenever possible!
     * NOTE: An empty matrix can be made a mask here only if it has been empty
     * since creation!
     *
     * @param pbuff Buffer
     * @param strid See above
     * @param rngR  "
     * @param rngC  "
     * @param watch Optional. S.a.
     */
    virtual void reassign(const T* pbuff,int strid,const Range& rngR,
			  const Range& rngC,MemWatchBase* watch=0) {
      if (!isMaskOrEmpty())
	throw MaskObjectException(EXCEPT_MSG("Needs to be mask or empty!"));
      if (rngC.isOpen())
	throw InvalidParameterException(EXCEPT_MSG("'rngC' must not be open"));
      if (strid<=0 || rngR.checkRange(strid))
	throw InvalidParameterException(EXCEPT_MSG("'strid','rngR' not compatible"));
      if (!rngC.isUniqueMap() || !rngR.isUniqueMap())
	throw InvalidParameterException(EXCEPT_MSG("Ranges must be 1-1"));
      const Range* rngrP=&rngR;
      Handle<Range> rngHand;
      if (rngR.isOpen()) {
	// Close open range
	rngHand.changeRep(new Range(rngR.getStart(),strid-1));
	rngrP=rngHand;
      }
      if (doCheckValid() && !isValidMat(pbuff,strid,*rngrP,rngC))
	throw WrongTypeException(EXCEPT_MSG("Invalid matrix elements"));
      reassignInt(pbuff,strid,*rngrP,rngC,watch);
    }

    /**
     * Same as above, but 'pbuff' is 'ArrayHandle' (use mem. watcher from
     * there).
     * NOTE: If 'rngC' is open, then the size of 'pbuff' must be a multiple
     * of 'strid'.
     *
     * @param pbuff Buffer array
     * @param strid See above
     * @param rngR  ". Def.: full
     * @param rngC  ". Def: full (see comm. for open range)
     */
    virtual void reassign(const ArrayHandle<T>& pbuff,int strid,
			  const Range& rngR=RangeFull::get(),
			  const Range& rngC=RangeFull::get()) {
      const Range* rngcP;
      Handle<Range> hrng;
      if (rngC.isOpen()) {
	int sz=pbuff.size();
	if (sz%strid!=0)
	  throw InvalidParameterException(EXCEPT_MSG("'rngC' is open"));
	if (rngC.getStart()>=(sz/=strid))
	  throw InvalidParameterException(EXCEPT_MSG("rngC"));
	hrng.changeRep(new Range(rngC.getStart(),sz-1));
	rngcP=hrng.p();
      } else
	rngcP=&rngC;
      reassign(pbuff.p(),strid,rngR,*rngcP,pbuff.getMemWatch());
    }

    /**
     * Reassigns a mask matrix (or the empty matrix) to point to a part of
     * another matrix 'mat'. The semantics is similar to 'reassign' with a
     * flat buffer, where the content of the flat buffer is given by 'mat',
     * read in column-major ordering. Here, the ranges 'rngR' (row pos.) and
     * 'rngC' (col. pos.) can be open. The def. for either is the full range.
     * <p>
     * If 'mat' is itself an indexed mask, the pos. indexes are mapped through
     * the indexes given by the ranges.
     * The method tries to avoid creating pos. indexes. This is possible if
     * 'mat' is normal or a flat mask.
     * If 'mat' is an indexed mask, this object will be an indexed mask as
     * well. In this case, we try to use row/col. index of 'mat' without
     * copying. However, index vectors within the ranges are always copied.
     * <p>
     * NOTE: An empty matrix can be made a mask here only if it has been empty
     * since creation!
     *
     * @param mat  Source matrix
     * @param rngR See above
     * @param rngC "
     */
    virtual void reassign(const BaseMatrix<T>& mat,
			  const Range& rngR=RangeFull::get(),
			  const Range& rngC=RangeFull::get());

    /**
     * Reassigns a mask matrix (or the empty matrix) to point to vector
     * 'vec'. If 'col'==true (def.), this is a column vector, oth. a row.
     * <p>
     * NOTE: An empty matrix can be made a mask here only if it has been empty
     * since creation!
     *
     * @param vec  Source vector
     * @param col  Column vector? Def.: true
     */
    virtual void reassign(const BaseVector<T>& vec,bool col=true);

    /**
     * Static wrapper for 'reassign' above, returns temp. mask matrix.
     * NOTE: As opposed to 'operator()', there is no notion of "expanding
     * the buffer" here: 'pbuff' must be large enough.
     * If 'watch' is given, it is used as mem. watcher 'buffWatch' for the
     * mem. region, otherwise no watcher is used (dangerous!).
     * <p>
     * NOTE: This is just to interface foreign code! Use 'ArrayHandle' version
     * below whenever possible!
     * Has to be implemented in every subclass with the corr. ret. type.
     *
     * @param pbuff See 'reassign'
     * @param pm    Number of rows
     * @param pn    Number of cols
     * @param strid S.a.
     * @param watch Optional. S.a.
     * @return      S.a.
     */
    static TempBaseMatrix<T> mask(const T* pbuff,int pm,int pn,int strid,
				  MemWatchBase* watch=0) {
      // In subclasses: copy code, but replace following line by one where the
      // empty matrix is created using the target type.
      BaseMatrix<T>* mat=new BaseMatrix<T>();
      mat->reassign(pbuff,pm,pn,strid,watch);
      return TempBaseMatrix<T>(mat,mat);
    }

    /**
     * Same as above, but with 'ArrayHandle' argument 'parr'. We use the mem.
     * watcher of this handle.
     * Has to be implemented in every subclass with the corr. ret. type.
     *
     * @param parr Array to be masked
     * @param pm    Number of rows
     * @param pn    Number of cols
     * @param strid S.a.
     */
    static TempBaseMatrix<T> mask(const ArrayHandle<T>& parr,int pm,int pn,
				  int strid) {
      // In subclasses: copy code, but replace following line by one where the
      // empty matrix is created using the target type.
      BaseMatrix<T>* mat=new BaseMatrix<T>();
      mat->reassign(parr.p(),pm,pn,strid,parr.getMemWatch());
      return TempBaseMatrix<T>(mat,mat);
    }

    /**
     * Static wrapper for 'reassign' above. Returns temp. mask matrix.
     * <p>
     * NOTE: This is just to interface foreign code! Use 'ArrayHandle' version
     * below whenever possible!
     * Has to be implemented in every subclass with the corr. ret. type.
     *
     * @param pbuff Buffer
     * @param strid See above
     * @param rngR  "
     * @param rngC  "
     * @param watch Optional. S.a.
     * @return      S.a.
     */
    static TempBaseMatrix<T> mask(const T* pbuff,int strid,const Range& rngR,
				  const Range& rngC,MemWatchBase* watch=0) {
      // In subclasses: copy code, but replace following line by one where the
      // empty matrix is created using the target type.
      BaseMatrix<T>* mat=new BaseMatrix<T>();
      mat->reassign(pbuff,strid,rngR,rngC,watch);
      return TempBaseMatrix<T>(mat,mat);
    }

    /**
     * Same as above, but uses 'ArrayHandle'.
     *
     * @param pbuff Buffer
     * @param strid See above
     * @param rngR  ". Def.: full
     * @param rngC  ". Def.: full
     * @return      S.a.
     */
    static TempBaseMatrix<T> mask(const ArrayHandle<T>& pbuff,int strid,
				  const Range& rngR=RangeFull::get(),
				  const Range& rngC=RangeFull::get()) {
      // In subclasses: copy code, but replace following line by one where the
      // empty matrix is created using the target type.
      BaseMatrix<T>* mat=new BaseMatrix<T>();
      mat->reassign(pbuff,strid,rngR,rngC);
      return TempBaseMatrix<T>(mat,mat);
    }

    /**
     * Static wrapper for 'reassign' above.
     * The purpose is type conversion without copying ('mat' is masked iff
     * it can be seen as matrix of this class).
     * <p>
     * Has to be implemented in every subclass with the corr. ret. type.
     *
     * @param mat  S.a.
     * @param rngR "
     * @param rngC "
     * @return     Temp. mask
     */
    static TempBaseMatrix<T> mask(const BaseMatrix<T>& mat,
				  const Range& rngR=RangeFull::get(),
				  const Range& rngC=RangeFull::get()) {
      // In subclasses: copy code, but replace following line by one where the
      // empty matrix is created using the target type.
      BaseMatrix<T>* mmat=new BaseMatrix<T>();
      mmat->reassign(mat,rngR,rngC);
      return TempBaseMatrix<T>(mmat,mmat);
    }

    /**
     * Allows to mask vector as matrix (col. vector by default).
     * <p>
     * Has to be implemented in every subclass with the corr. ret. type.
     *
     * @param vec Vector to be masked
     * @param col Column vector? Otherwise: row vector
     * @return    Temp. mask
     */
    static TempBaseMatrix<T> mask(const BaseVector<T>& vec,bool col=true) {
      // In subclasses: copy code, but replace following line by one where the
      // empty matrix is created using the target type.
      BaseMatrix<T>* mmat=new BaseMatrix<T>();
      mmat->reassign(vec,col);
      return TempBaseMatrix<T>(mmat,mmat);
    }

    /**
     * RValue version:
     * Creates a mask matrix picking out the part given by the row/col. pos.
     * ranges 'rngR'/'rngC' from this matrix. The res. mask can be used as
     * rvalue, e.g. with the assignment operator. The ranges must fit within
     * this matrix.
     * <p>
     * Either of the ranges can be a single position (int value). In this case,
     * we create a 'BaseVector<T>' or a 'T' object instead.
     * NOTE: A 'Range' object repr. a range of size 1 is NOT recognized as
     * single pos. in this sense.
     * <p>
     * Has to be implemented in every subclass with the corr. ret. type.
     * The implementation here just has to be copied.
     *
     * @param rngR See above. Def.: full range. Can be a single position
     * @param rngC "
     */
    const TempBaseMatrix<T> operator()(const Range& rngR=RangeFull::get(),
				       const Range& rngC=RangeFull::get())
      const {
      // In subclasses: dyn. cast pointer to target type!
      BaseMatrix<T>* mat=subrvalInt(rngR,rngC);
      return TempBaseMatrix<T>(mat,mat);
    }

    const TempBaseVector<T> operator()(const Range& rngR,int rngC) const {
      // In subclasses: dyn. cast pointer to target type!
      BaseVector<T>* vec=subrvalInt(rngR,rngC);
      return TempBaseVector<T>(vec,vec);
    }

    const TempBaseVector<T> operator()(int rngR,const Range& rngC=
				       RangeFull::get()) const {
      // In subclasses: dyn. cast pointer to target type!
      BaseVector<T>* vec=subrvalInt(rngR,rngC);
      return TempBaseVector<T>(vec,vec);
    }

    T operator()(int rngR,int rngC) const {
      return get(rngR,rngC);
    }

    /**
     * ATTENTION: L-value versions often used instead of r-value ones!
     * ==> SOLUTION??
     *
     * LValue version:
     * The resulting mask matrix/vector can be used as l-value.
     * If the ranges do not fit within this matrix, it is expanded to the
     * smallest size for which the ranges fit. The new entries are init. with
     * the def. fill value.
     * <p>
     * Comments of 'BaseVector<T>::operator()' l-value apply here as well!
     * <p>
     * Has to be implemented in every subclass with the corr. ret. type.
     * The implementation here just has to be copied.
     *
     * @param rngR See above. Def.: full range. Can be a single position
     * @param rngC "
     */
    TempBaseMatrix<T> operator()(const Range& rngR=RangeFull::get(),
				 const Range& rngC=RangeFull::get()) {
      // In subclasses: dyn. cast pointer to target type!
      BaseMatrix<T>* mat=sublvalInt(rngR,rngC);
      return TempBaseMatrix<T>(mat,mat);
    }

    TempBaseVector<T> operator()(const Range& rngR,int rngC) {
      // In subclasses: dyn. cast pointer to target type!
      BaseVector<T>* vec=sublvalInt(rngR,rngC);
      return TempBaseVector<T>(vec,vec);
    }

    TempBaseVector<T> operator()(int rngR,const Range& rngC=RangeFull::get()) {
      // In subclasses: dyn. cast pointer to target type!
      BaseVector<T>* vec=sublvalInt(rngR,rngC);
      return TempBaseVector<T>(vec,vec);
    }

    T& operator()(int rngR,int rngC);

    /**
     * Returns mask vector picking out (part of) the diagonal or off-diagonal
     * of this matrix. For 'off'>=0, the area is (i,i+'off'), i>=0 (above
     * diagonal). For 'off'<0, the area is (i-'off',i), i>=0 (below diagonal).
     * The mask vector is obtained by applying 'rng' to this area (def.: full).
     * NOTE: As opposed to 'operator()' l-value, the matrix is not enlarged.
     * If 'rng' does not fit, an exc. is thrown.
     *
     * @param off S.a. Def.: 0
     * @param rng S.a. Def.: full
     * @return    S.a.
     */
    TempBaseVector<T> diag(int off=0,const Range& rng=RangeFull::get())
      const {
      // In subclasses: dyn. cast pointer to target type!
      BaseVector<T>* vec=diagInt(off,rng);
      return TempBaseVector<T>(vec,vec);
    }

    // Conversion to flat matrix

    /**
     * Returns flat contiguous representation of this matrix (i.e. the striding
     * value is the same as the number of rows). If this matrix is flat and
     * contiguous, its buffer is returned as 'ArrayHandle' with the same mem.
     * watcher (more specific, the initial part of size m*n). Otherwise, a flat
     * copy is drawn.
     * NOTE: Usage of 'getLinMat' is safer (write-back).
     * NOTE: A copy is drawn for flat, but non-contiguous matrices as well.
     *
     * @return S.a.
     */
    virtual ArrayHandle<T> getFlatBuff() const;

    /**
     * The 'BaseLinMat' object 'blMat' (whose entries must be empty) is init.
     * with a flat representation of this matrix. This is used to interface
     * with code which does not allow for indexing.
     * A flat copy of this matrix is drawn only if it is an indexed mask.
     * <p>
     * Write-back:
     * If 'writeBack'==true, the 'BaseLinMat' object is conf. s.t. upon
     * its destruction the content is written back into this matrix (if a
     * flat copy was drawn). In this case, a structure pattern can be def. via
     * 'strpatt' (def.: full matrix; see 'MatStrct::XXX'). Only
     * the part of the matrix def. by the pattern is written back then. See
     * 'BaseLinMat' for details.
     *
     * @param blMat     BaseLinMat object ret. here
     * @param writeBack S.a.
     * @param strpatt   S.a. Def.: full matrix
     */
    virtual void getLinMat(BaseLinMat<T>& blMat,bool writeBack,
			   unsigned char strpatt=MatStrct::normal)
      const;

    // Subscripting and element access

    /**
     * Subscripting operator
     * This is for convenience only, use 'get','set' or 'operator()' for a
     * more efficient access.
     * NOTE: As opposed to 'operator()', the matrix is not expanded. An
     * exception is thrown if 'row' is out of range.
     * NOTE: Unfortunately, s.th. like
     *   mat[i][j]
     * does not work, bec. the automatic conversion for 'mat[i]' cannot
     * det. the correct type. Could do
     *   ((BaseMatrix<T>) mat[i])[j]
     * <p>
     * Have to be overwritten in subclasses to return the correct type!
     *
     * @param row Row position
     * @return    S.a.
     */
    TempBaseVector<T> operator[](int row) const {
      checkTS(TSARG);
      // In subclasses: dyn. cast pointer to target type!
      BaseVector<T>* vec=subrvalInt(row,RangeFull::get());
      return TempBaseVector<T>(vec,vec);
    }

    /**
     * @param row Row number
     * @param col Column number
     * @return    Matrix element
     */
    T get(int row,int col) const {
      checkTS(TSARG);
      if (row<0 || row>=m || col<0 || col>=n)
	throw OutOfRangeException(EXCEPT_MSG(""));
      if (!isMaskInd())
	return *(buff+(col*stride+row));
      else
	return *(buff+(COLPOS(col)*stride+ROWPOS(row)));
    }

    /**
     * NOTE: As opposed to 'operator[]' l-value, we check the validity of
     * 'el' here.
     * NOTE: As opposed to 'operator()' l-value, the matrix is not expanded
     * if 'row','col' is out of range.
     *
     * @param row Row number
     * @param col Column number
     * @param el  New matrix element
     */
    void set(int row,int col,T el) {
      checkTS(TSARG);
      if (row<0 || row>=m || col<0 || col>=n)
	throw OutOfRangeException(EXCEPT_MSG(""));
      if (!isValidElement(el))
	throw WrongTypeException("Invalid matrix element");
      if (!isMaskInd())
	*(buff+(col*stride+row))=el;
      else
	*(buff+(COLPOS(col)*stride+ROWPOS(row)))=el;
    }

    // Methods to control matrix and buffer size

    /**
     * Resizes matrix buffer to given size and sets matrix entries to the def.
     * fill value.
     * The new matrix has size and buffer size 'rows'-by-'cols' and striding
     * factor =='rows'.
     * <p>
     * Note: To simply ensure that the buffer can hold an m*n matrix, it's
     * better to use zeros(m,n), the latter will NOT reallocate the buffer
     * if it's too large.
     * NOTE: For lossless shrinking, use 'resizeSave'.
     * NOTE: If either 'rows' or 'cols' is 0, the buffer will be deallocated
     * (empty matrix).
     *
     * @param rows Number of desired rows
     * @param cols Number of desired columns
     */
    virtual void resize(int rows,int cols);

    /**
     * Resizes matrix buffer to given size, but keeps the old entries.
     * The new matrix has size and buffer size 'rows'-by-'cols' and striding
     * constant 'rows'.
     * A size that would lead to loss of entries is refused. The new entries
     * are set to the def. fill value.
     * NOTE: If 'rows' and/or 'cols' is -1, it is substituted with the
     * actual numbers of rows/columns of the current matrix. F.ex., calling
     * the method without parameters will resize the buffer to the actual
     * matrix size.
     *
     * @param rows Number of desired rows. Def: -1
     * @param cols Number of desired columns. Def.: -1
     */
    virtual void resizeSave(int rows=-1,int cols=-1);

    /**
     * Resizes matrix to given size, but keeps the old entries. A size that
     * would lead to loss of entries is refused. The new entries are set
     * to the def. fill value.
     * <p>
     * New buffer size and dimensions:
     * - If buffer too small ('rows'*'cols' > 'buffSz'):
     *   Buffer reallocated to size
     *     getNewBuffSize(rows*cols,buffSz),
     *   then 'stride' set to 'rows'
     * - Otherwise: no realloc., 'buffSz' stays the same.
     *   - If 'stride' < 'rows':
     *     Buffer content copied internally, then 'stride' set to 'rows'
     *   - Oth.: 'stride' stays same as before if possible (i.e. if
     *     'stride'*'cols' <= 'buffSz'), oth. is reduced to 'rows' (this
     *     req. internal copying as well)
     *
     * @param rows   Number of desired rows
     * @param cols   Number of desired columns
     */
    virtual void expand(int rows,int cols);

    /**
     * Safe (public) version of 'ensureCapacity'. If the matrix has size
     * 'rows'-by-'cols', this method does nothing. Otherwise, the matrix
     * is brought to size 'rows'-by-'cols' and filled with 'a'.
     *
     * @param rows S.a.
     * @param cols S.a. Def.: same as 'rows'
     * @param a    Optional. Def.: Def. fill value
     */
    virtual void ensureSize(int rows,int cols,T a) {
      checkTS(TSARG);
      if (m!=rows || n!=cols) fill(rows,cols,a);
    }

    virtual void ensureSize(int rows,int cols=-1) {
      if (cols==-1) cols=rows;
      ensureSize(rows,cols,this->defFill);
    }

    /**
     * Here, 'rows'<='m', 'cols'<='n'. This matrix is replaced by its upper
     * left 'rows'-by-'cols' block. The buffer dimensions are not changed.
     * <p>
     * NOTE: Leads to loss of entries (as opposed to 'expand', 'resizeSave')!
     *
     * @param rows S.a.
     * @param cols S.a.
     */
    virtual void shrink(int rows,int cols) {
      if (rows>m || cols>n)
	throw InvalidParameterException(EXCEPT_MSG(""));
      if (rows==0 || cols==0) rows=cols=0;
      ensureCapacity(rows,cols);
    }

    /**
     * Does not change size or content of matrix, but makes sure that the
     * matrix buffer fits into the frame given by 'rows','cols'. We must
     * have m<='rows', n<='cols'. If the matrix buffer does not fit, it is
     * reallocated to have size 'rows'-by-'cols' exactly. Otherwise, it
     * keeps its old size.
     * <p>
     * Difference to 'resizeSave': buffer is not changed if it fits into
     * 'rows'-by-'cols'. Also, the matrix size does not change here.
     * NOTE: If 'rows' ('cols') is -1, it is substituted with m (n).
     * NOTE: May involve internal copying without re-alloc. if the buffer
     * fits, but 'stride'>'rows'. In this case, 'stride'=='rows' afterwards.
     *
     * @param rows See above. Def.: -1
     * @param cols "
     */
    void shrinkBuffer(int rows=-1,int cols=-1);

    /**
     * Specialised method.
     * Reshapes matrix to size WITHOUT changing the buffer size 'buffSz'. If
     * this is not possible, a 'WrongDimensionException' is thrown. Either of
     * 'rows','cols' can be -1 (but not both), in which case it is chosen
     * maximally s.t. the res. matrix fits into the given buffer.
     * By def., the new matrix is filled with the def. fill value, but if
     * 'noInit'==true the buffer entries are not changed.
     * NOTE: The new buffer must have at least one row and one column.
     * <p>
     * NOTE: If this method is called with both 'rows','cols'!=-1, then w.r.t.
     * matrix buffer dimensions we try to max. the number of rows, i.e. max.
     * 'stride' s.t. the buffer has 'cols' columns.
     * However, if 'strict'==true, this is NOT done, the new buffer will have
     * 'stride'=='rows' (max. number of cols), the matrix will be flat.
     * <p>
     * NOTE: Differs from 'ensureCapacity'! Here, the buffer dimensions are
     * changed also if the new size 'rows'-by-'cols' would fit into the
     * current dimensions.
     *
     * @param rows   S.a.
     * @param cols   S.a. Def.: -1
     * @param noInit S.a. Def.: false
     * @param strict S.a. Def.: false
     */
    void reshape(int rows,int cols=-1,bool noInit=false,bool strict=false);

    /**
     * Version of 'reshape' which retains elements of the matrix. Namely, the
     * first min(n,'cols') elements of the first min(m,'rows') rows of this
     * matrix are retained in the new one. The remaining entries are filled
     * with the def. fill value.
     * NOTE: As opposed to 'resizeSave', 'expand' we can have 'rows'<m and/or
     * 'cols'<n, res. in loss of elements.
     *
     * @param rows   S.a.
     * @param cols   S.a.
     * @param strict S.a. Def.: false
     */
    void reshapeSave(int rows,int cols=-1,bool strict=false);

    // Assignment, conversion, insert, remove, exchange

    /**
     * Assignment operator
     * <p>
     * Physically copies all elements of 'mat'. To obtain a mask, use
     * 'reassign'. 'operator()' gives a temporary mask.
     * <p>
     * NOTE: This matrix may be a mask. In this case, if 'mat' fits into this
     * matrix, instead of throwing a 'MaskObjectException', we copy 'mat' into
     * the upper left corner of this matrix, all other entries are not
     & changed. If 'mat' does not fit, an exc. is thrown though.
     * ==> Allows more intuitive use together with 'operator()'.
     * <p>
     * NOTE: Overwrite in subclasses to return correct type (chaining).
     *
     * @param mat Matrix to copy
     */
    BaseMatrix<T>& operator=(const BaseMatrix<T>& mat) {
      assignInt(mat);
      return *this;
    }

    /**
     * Assignment operator, r.h.s. a (column) vector
     * <p>
     * Physically copies all elements of 'vec', which is taken as column
     * vector. To obtain a mask, or to deal with row vectors, use
     * 'reassign'.
     * <p>
     * NOTE: This matrix may be a mask, if it is already a column vector
     * of size >= than 'vec'. In this case, 'vec' is copied to the upper
     * part of this column, the other entries are not changed.
     * ==> Allows more intuitive use together with 'operator()'.
     * <p>
     * NOTE: Overwrite in subclasses to return correct type (chaining).
     *
     * @param mat Matrix to copy
     */
    BaseMatrix<T>& operator=(const BaseVector<T>& vec) {
      assignInt(vec);
      return *this;
    }

    /**
     * Assignment operator, r.h.s. a scalar
     * <p>
     * Sets all elements of this matrix to 'elem' (equiv. to 'fill'). The size
     * of this matrix is not changed, except if it is empty (grows to size 1,1
     * in this case).
     * <p>
     * NOTE: Overwrite in subclasses to return correct type (chaining).
     *
     * @param elem Fill element
     */
    BaseMatrix<T>& operator=(T elem) {
      assignInt(elem);
      return *this;
    }

    /**
     * Conversion from matrix with different entry type. This is a version
     * of the assignment operator, thus elements are copied physically.
     * Casts from T2 to T must be possible and result in valid conversions.
     * <p>
     * If this matrix is a mask and 'src' fits into it, we copy 'src' into
     * the upper left corner.
     *
     * @param src Source object, element type T2
     */
    template<class T2> void convert(const BaseMatrix<T2>& src) {
      if (!isValidConvMat(src))
	throw WrongTypeException(EXCEPT_MSG("Invalid matrix elements after type casts"));
      checkTS(TSARG); src.checkTS(TSARG);
      int cpm=src.rows(),cpn=src.cols();
      if (m<cpm || n<cpn || !isMask()) ensureCapacity(cpm,cpn);
      Handle<BaseVector<T> > trgM;
      BaseVector<T2> srcM;
      Range rng(0,cpm-1);
      for (int i=0; i<cpn; i++) {
	trgM.changeRep(subrvalInt(rng,i));
	srcM.reassign(src(rng,i)); // cannot use internal meth.
	trgM->convert(srcM);
      }
    }

    /**
     * Matrix 'mat' is inserted into this matrix as follows: rows 'rpos' and
     * lower are moved down by 'mat.rows()'. If 'rngC' applied to rows of
     * 'mat' is larger than n, this matrix is expanded by adding a min.
     * number of cols. If 'rngC' is not open, its size must be 'mat.cols()'.
     * Then, row i of 'mat' is copied into row 'rpos'+i of
     * this matrix, applying 'rngC'. This means that if 'rngC' does not cover
     * the new rectangle, the remaining entries will have the def. fill value.
     * If 'rpos'==-1,the rows are appended at the bottom.
     * NOTE: 'rngC' can be open, but otherwise has to have the same size as
     * 'mat's number of cols.
     * NOTE: If 'rpos'>m, the new rows in between are filled with the def.
     * fill value.
     * <p>
     * NOTE: Use 'mask' to sell vector as row vector, in order to insert a
     * single row, f.ex.:
     *   a.insertRows(BaseMatrix<T>::mask(vec,false),rpos);
     *
     * @param mat  S.a.
     * @param rpos S.a. Def.: -1
     * @param rngC S.a. Def.: full
     */
    void insertRows(const BaseMatrix<T>& mat,int rpos=-1,
		    const Range& rngC=RangeFull::get());

    /**
     * Inserts 'num' rows into this matrix starting from row 'pos'. The
     * new rows are filled with 'val'.
     * NOTE: Here, 'pos'>'m' is not permitted.
     *
     * @param num Number of new rows
     * @param pos S.a. Def.: -1 (bottom)
     * @param val S.a. Def.: def. fill value
     */
    void insertRows(int num,int pos,T val);

    void insertRows(int num,int pos=-1) {
      insertRows(num,pos,this->defFill);
    }

    /**
     * Same as 'insertRows', but column- instead of row-oriented.
     *
     * @param mat  S.a.
     * @param cpos S.a.
     * @param rngR S.a.
     */
    void insertCols(const BaseMatrix<T>& mat,int cpos=-1,
		    const Range& rngR=RangeFull::get());

    /**
     * Same as 'insertRows', but column- instead of row-oriented.
     *
     * @param num S.a.
     * @param pos S.a. Def.: -1 (right)
     * @param val S.a. Def.: def. fill value
     */
    void insertCols(int num,int pos,T val);

    void insertCols(int num,int pos=-1) {
      insertCols(num,pos,this->defFill);
    }

    /**
     * Removes rows 'pos' to 'pos'+'num'-1 from this matrix. If 'num'==-1, it
     * is set to m-'pos' (all remaining rows).
     * NOTE: If 'pos'+'num'>'m', but 'pos' is valid (<'m'), all remaining
     * rows are removed.
     *
     * @param pos S.a.
     * @param num S.a. Def.: 1
     */
    void removeRows(int pos,int num=1);

    /**
     * Removes cols 'pos' to 'pos'+'num'-1 from this matrix. If 'num'==-1, it
     * is set to n-'pos'.
     *
     * @param pos S.a.
     * @param num S.a. Def.: 1
     */
    void removeCols(int pos,int num=1);

    /**
     * Exchanges rows/cols at pos. 'p1', 'p2' with each other.
     *
     * @param p1 S.a.
     * @param p2 s.a.
     */
    void exchangeCols(int p1,int p2) {
      if (p1<0 || p1>=n || p2<0 || p2>=n)
	throw OutOfRangeException(EXCEPT_MSG(""));
      if (p1!=p2) {
	Handle<BaseVector<T> > v1(subrvalInt(RangeFull::get(),p1));
	Handle<BaseVector<T> > v2(subrvalInt(RangeFull::get(),p2));
	v1->exchange(*v2);
      }
    }

    void exchangeRows(int p1,int p2) {
      if (p1<0 || p1>=m || p2<0 || p2>=m)
	throw OutOfRangeException(EXCEPT_MSG(""));
      if (p1!=p2) {
	Handle<BaseVector<T> > v1(subrvalInt(p1,RangeFull::get()));
	Handle<BaseVector<T> > v2(subrvalInt(p2,RangeFull::get()));
	v1->exchange(*v2);
      }
    }

    /**
     * Variant of assignment operator under structure constraints. A
     * structure pattern def. a part of a matrix (see 'MatStrct').
     * The source and target str. patt. 'srcstrct', 'trgstrct' must be
     * s.t. the parts have the same size, but they can be different
     * (upper/lower). The assignment is done for the parts only. This matrix
     * is enlarged if necessary and filled with def. value.
     * <p>
     * NOTE: Just a public wrapper for 'assignInt'.
     *
     * @param src Source matrix
     * @param srcstrct S.a. Def.: normal (full matrix)
     * @param trgstrct S.a. Def.: same as 'srcstrct'
     */
    virtual void assignStrct(const BaseMatrix<T>& src,
			     uchar srcstrct,uchar trgstrct) {
      assignInt(src,srcstrct,trgstrct);
    }

    virtual void assignStrct(const BaseMatrix<T>& src,
			     uchar srcstrct=MatStrct::normal) {
      assignInt(src,srcstrct,srcstrct);
    }

    /**
     * Builds column vector which contains the matrix element at
     * (rng[i],ind[rng[i]]) at pos. i. If 'rowMaj'==true, the same is done
     * on the transpose of this matrix.
     *
     * @param vec    Vector ret. here
     * @param ind    S.a.
     * @param rowMaj S.a. Def.: false
     * @param rng    S.a. Def.: full
     */
    virtual void select(BaseVector<T>& vec,const BaseVector<int>& ind,
			bool rowMaj=false,const Range& rng=RangeFull::get())
      const;

    virtual BaseVector<T> select(const BaseVector<int>& ind,
				 bool rowMaj=false,
				 const Range& rng=RangeFull::get()) const {
      BaseVecWrapper<T> resWrap(newEmptyVector());
      select(*(resWrap.getRep()),ind,rowMaj,rng);

      return resWrap;
    }

    /*
     * Vectorisation methods
     * All these methods process matrices in column-major ordering,
     * columns top-down, left to right.
     * <p>
     * NOTE: If 'doCheckValid' returns true, a result matrix is first
     * written into a temp. 'BaseMatrix<T>' before being written back to
     * this matrix using 'assignInt' (which forces element checks).
     */

    /**
     * Computes nullary function for each element. If 'rows', 'cols' are
     * given, the matrix size is adjusted.
     *
     * @param func Nullary T function
     * @param rows Optional. S.a.
     * @param cols Optional. Def.: same as 'rows'
     */
    template<class T1> void
    apply0(const NullaryFunc<T1>& func,int rows=-1,int cols=-1) {
      int i;

      checkTS(TSARG);
      Handle<BaseMatrix<T> > tempMat;
      BaseMatrix<T>* trgMat;
      if (!doCheckValid()) trgMat=this;
      else {
	// Build up in temp. matrix to allow for element check before
	// assignment
	tempMat.changeRep(new BaseMatrix<T>());
	trgMat=tempMat;
      }
      if (rows==-1) {
	rows=m; cols=n;
      } else if (rows<=0)
	throw InvalidParameterException(EXCEPT_MSG("rows"));
      else {
	if (cols==-1) cols=rows;
	else if (cols<=0)
	  throw InvalidParameterException(EXCEPT_MSG("cols"));
      }
      trgMat->ensureCapacity(rows,cols);
      if (!trgMat->isMaskInd()) {
	// Flat matrix
	if (trgMat->m==trgMat->stride) {
	  // Flat and contiguous
	  ArrayUtils<T>::applyNulFunc(trgMat->buff,trgMat->m*trgMat->n,func);
	} else {
	  T* trgB=trgMat->buff;
	  for (i=0; i<n; i++,trgB+=trgMat->stride)
	    ArrayUtils<T>::applyNulFunc(trgB,trgMat->m,func);
	}
      } else {
	// Use mask vectors for columns
	Handle<BaseVector<T> > trgM;
	Range rng(0,m-1);
	for (i=0; i<n; i++) {
	  trgM.changeRep(subrvalInt(rng,i));
	  trgM->apply0(func);
	}
      }
      if (doCheckValid())
	// copy back (forces element check)
	assignInt(*tempMat);
    }

    /**
     * Applies a unary function T2 -> T (given by 'func') to each
     * element of 'a'. The result is written into this matrix.
     * NOTE: 'a' and this matrix may be the same object.
     *
     * @param a    Source object. Def.: this matrix
     * @param func Unary T2 -> T function
     */
    template<class UnOp,class T2> void
    apply1(const BaseMatrix<T2>& a,const UnOp& func) {
      int i;

      checkTS(TSARG); a.checkTS(TSARG);
      Handle<BaseMatrix<T> > tempMat;
      BaseMatrix<T>* trgMat;
      if (!doCheckValid()) trgMat=this;
      else {
	// Build up in temp. matrix to allow for element check before
	// assignment
	tempMat.changeRep(new BaseMatrix<T>());
	trgMat=tempMat;
      }
      trgMat->ensureCapacity(a.rows(),a.cols());
      if (!trgMat->isMaskInd() && !a.isMaskInd()) {
	// Both matrices are flat
	if (trgMat->m==trgMat->stride && trgMat->m==a.getStride_INT()) {
	  // Both are flat and contiguous
	  ArrayUtils<T>::applyFunc(trgMat->buff,a.getBuffPtr_INT(),
				   trgMat->m*trgMat->n,func);
	} else {
	  T* trgB=trgMat->buff;
	  const T2* srcB=a.getBuffPtr_INT();
	  for (i=0; i<trgMat->n; i++,trgB+=trgMat->stride,
		 srcB+=a.getStride_INT())
	    ArrayUtils<T>::applyFunc(trgB,srcB,trgMat->m,func);
	}
      } else {
	// Use mask vectors for columns
	Handle<BaseVector<T> > trgM;
	BaseVector<T2> srcM;
	Range rng(0,trgMat->m-1);
	for (i=0; i<trgMat->n; i++) {
	  trgM.changeRep(trgMat->subrvalInt(rng,i));
	  srcM.reassign(a(rng,i));
	  trgM->apply1(srcM,func);
	}
      }
      if (doCheckValid())
	// copy back (forces element check)
	assignInt(*tempMat);
    }

    template<class UnOp> void apply1(const UnOp& func) {
      int i;

      if (m!=0) {
	checkTS(TSARG);
	Handle<BaseMatrix<T> > tempMat;
	BaseMatrix<T>* trgMat;
	if (!doCheckValid()) trgMat=this;
	else {
	  tempMat.changeRep(new BaseMatrix<T>(m,n));
	  trgMat=tempMat;
	}
	if (!isMaskInd()) {
	  // Flat
	  if (m==stride && m==trgMat->stride)
	    // Flat and contiguous
	    ArrayUtils<T>::applyFunc(trgMat->buff,buff,m*n,func);
	  else {
	    T* trgB=trgMat->buff;
	    const T* srcB=buff;
	    for (i=0; i<n; i++,srcB+=stride,trgB+=trgMat->stride)
	      ArrayUtils<T>::applyFunc(trgB,srcB,m,func);
	  }
	} else {
	  // Use mask vectors for columns
	  Handle<BaseVector<T> > trgM;
	  Handle<BaseVector<T> > srcM;
	  Range rng(0,m-1);
	  for (i=0; i<n; i++) {
	    trgM.changeRep(trgMat->subrvalInt(rng,i));
	    srcM.changeRep(subrvalInt(rng,i));
	    trgM->apply1(*srcM,func);
	  }
	}
	if (doCheckValid())
	  assignInt(*tempMat); // check elems., then assign
      }
    }

    /**
     * Applies a binary function (T2,T3) -> T (given by 'func') to each
     * pair of elements from 'a','b' (must have same size). The result is
     * returned in this matrix.
     * NOTE: 'a','b' and this matrix may be the same object.
     *
     * @param a    Source object. Def.: this matrix
     * @param b    Source object
     * @param func Binary function
     */
    template<class BinOp,class T2,class T3> void
    apply2(const BaseMatrix<T2>& a,const BaseMatrix<T3>& b,const BinOp& func) {
      int i;
 
      checkTS(TSARG); a.checkTS(TSARG); b.checkTS(TSARG);
      if (a.rows()!=b.rows() || a.cols()|=b.cols())
	throw DimMismatchException(EXCEPT_MSG("a,b must have same size"));
      Handle<BaseMatrix<T> > tempMat;
      BaseMatrix<T>* trgMat;
      if (!doCheckValid()) trgMat=this;
      else {
	// Build up in temp. matrix to allow for element check before
	// assignment
	tempMat.changeRep(new BaseMatrix<T>());
	trgMat=tempMat;
      }
      trgMat->ensureCapacity(a.rows(),a.cols());
      if (!trgMat->isMaskInd() && !a.isMaskInd() && !b.isMaskInd()) {
	// All matrices are flat
	if (trgMat->m==trgMat->stride && trgMat->m==a.getStride_INT() &&
	    trgMat->m==b.getStride_INT()) {
	  // All are flat and contiguous
	  ArrayUtils<T>::applyBinFunc(trgMat->buff,a.getBuffPtr_INT(),
				      b.getBuffPtr_INT(),trgMat->m*
				      trgMat->n,func);
	} else {
	  T* trgB=trgMat->buff;
	  const T2* aB=a.getBuffPtr_INT();
	  const T3* bB=b.getBuffPtr_INT();
	  for (i=0; i<trgMat->n; i++,trgB+=trgMat->stride,
		 aB+=a.getStride_INT(),bB+=b.getStride_INT())
	    ArrayUtils<T>::applyBinFunc(trgB,aB,bB,trgMat->m,func);
	}
      } else {
	// Use mask vectors for columns
	Handle<BaseVector<T> > trgM;
	BaseVector<T2> aM; // cannot use internal 'subrvalInt' of these
	BaseVector<T3> bM;
	Range rng(0,trgMat->m-1);
	for (i=0; i<trgMat->n; i++) {
	  trgM.changeRep(subrvalInt(rng,i));
	  aM.reassign(a(rng,i));
	  bM.reassign(b(rng,i));
	  trgM->apply2(aM,bM,func);
	}
      }
      if (doCheckValid())
	// copy back (forces element check)
	assignInt(*tempMat);
    }

    template<class BinOp,class T2> void
    apply2(const BaseMatrix<T2>& b,const BinOp& func) {
      int i;
 
      checkTS(TSARG); b.checkTS(TSARG);
      if (m!=b.rows() || n|=b.cols())
	throw DimMismatchException("a,b must have same size");
      Handle<BaseMatrix<T> > tempMat;
      BaseMatrix<T>* trgMat;
      if (!doCheckValid()) trgMat=this;
      else {
	// Build up in temp. matrix to allow for element check before
	// assignment
	tempMat.changeRep(new BaseMatrix<T>(m,n));
	trgMat=tempMat;
      }
      if (!trgMat->isMaskInd() && !b.isMaskInd()) {
	// All matrices are flat
	if (trgMat->m==trgMat->stride && trgMat->m==b.getStride_INT()) {
	  // All are flat and contiguous
	  ArrayUtils<T>::applyBinFunc(trgMat->buff,trgMat->buff,
				      b.getBuffPtr_INT(),trgMat->m*
				      trgMat->n,func);
	} else {
	  T* trgB=trgMat->buff;
	  const T2* bB=b.getBuffPtr_INT();
	  for (i=0; i<trgMat->n; i++,trgB+=trgMat->stride,
		 bB+=b.getStride_INT())
	    ArrayUtils<T>::applyBinFunc(trgB,trgB,bB,trgMat->m,func);
	}
      } else {
	// Use mask vectors for columns
	Handle<BaseVector<T> > trgM;
	BaseVector<T2> bM; // Cannot use internal 'subrvalInt' of this one
	Range rng(0,trgMat->m-1);
	for (i=0; i<trgMat->n; i++) {
	  trgM.changeRep(subrvalInt(rng,i));
	  bM.reassign(b(rng,i));
	  trgM->apply2(bM,func);
	}
      }
      if (doCheckValid())
	// copy back (forces element check)
	assignInt(*tempMat);
    }

    /**
     * For each column (left to right), 'acc' is applied to each element
     * (top-down). The results are ret. in vector 'res', which has one
     * entry for each column.
     * 'acc' is reset at the beginning of each column/row (also the first).
     * If 'func' is given, elements are mapped through 'func' before
     * passing them to 'acc'.
     * If 'rowMaj'==true, 'acc' is applied to rows (top-down, left-right).
     *
     * @param acc    Accumulator object
     * @param res    Result vector
     * @param rowMaj Optional. Def.: false
     * @param func   Optional. Def.: identity
     */
    template<class T2,class T3> void
    accumulate(const AccumulFunc<T2,T3>& acc,BaseVector<T3>& res,
	       bool rowMaj=false) const {
      int i;
      Handle<BaseVector<T> > srcM;

      checkTS(TSARG);
      if (m==0)
	throw WrongDimensionException(EXCEPT_MSG("Matrix must not be empty"));
      if (!rowMaj) {
	// Column-major
	res.fill(n);
	if (!isMaskInd()) {
	  // Flat matrix
	  T* srcB=buff;
	  for (i=0; i<n; i++,srcB+=stride) {
	    acc.reset();
	    ArrayUtils<T>::applyAcc(srcB,m,acc);
	    res.set(i,acc.get());
	  }
	} else {
	  // Use mask vectors for columns
	  Range rng(0,m-1);
	  for (i=0; i<n; i++) {
	    srcM.changeRep(subrvalInt(rng,i));
	    acc.reset();
	    srcM->accumulate(acc);
	    res.set(i,acc.get());
	  }
	}
      } else {
	// Row-major
	res.fill(m);
	if (!isMaskInd()) {
	  // Flat matrix
	  T* srcB=buff;
	  for (i=0; i<m; i++,srcB++) {
	    acc.reset();
	    ArrayUtils<T>::applyAcc(srcB,n,acc,stride);
	    res.set(i,acc.get());
	  }
	} else {
	  // Use mask vectors for rows
	  Range rng(0,n-1);
	  for (i=0; i<m; i++) {
	    srcM.changeRep(subrvalInt(i,rng));
	    acc.reset();
	    srcM->accumulate(acc);
	    res.set(i,acc.get());
	  }
	}
      }
    }

    template<class UnOp,class T2,class T3> void
    accumulate(const AccumulFunc<T2,T3>& acc,BaseVector<T3>& res,bool rowMaj,
	       const UnOp& func)
    const {
      int i;
      Handle<BaseVector<T> > srcM;

      checkTS(TSARG);
      if (m==0)
	throw WrongDimensionException(EXCEPT_MSG("Matrix must not be empty"));
      if (!rowMaj) {
	// Column-major
	res.fill(n);
	if (!isMaskInd()) {
	  // Flat matrix
	  T* srcB=buff;
	  for (i=0; i<n; i++,srcB+=stride) {
	    acc.reset();
	    ArrayUtils<T>::applyAccFunc(srcB,m,func,acc);
	    res.set(i,acc.get());
	  }
	} else {
	  // Use mask vectors for columns
	  Range rng(0,m-1);
	  for (i=0; i<n; i++) {
	    srcM.changeRep(subrvalInt(rng,i));
	    acc.reset();
	    srcM->accumulate(acc,func);
	    res.set(i,acc.get());
	  }
	}
      } else {
	// Row-major
	res.fill(m);
	if (!isMaskInd()) {
	  // Flat matrix
	  T* srcB=buff;
	  for (i=0; i<m; i++,srcB++) {
	    acc.reset();
	    ArrayUtils<T>::applyAccFunc(srcB,n,func,acc,stride);
	    res.set(i,acc.get());
	  }
	} else {
	  // Use mask vectors for rows
	  Range rng(0,m-1);
	  for (i=0; i<m; i++) {
	    srcM.changeRep(subrvalInt(i,rng));
	    acc.reset();
	    srcM->accumulate(acc,func);
	    res.set(i,acc.get());
	  }
	}
      }
    }

    /**
     * If a,b vectors in T2,T3, f : (T2,T3) -> T ('func'), this matrix is
     * set to (f(a_i,b_j))_{i,j}.
     *
     * @param a    Source vector
     * @param b    Source vector
     * @param func Binary (T2,T3) -> T function
     */
    template<class BinOp,class T2,class T3> void
    applyProd(const BaseVector<T2>& a,const BaseVector<T3>& b,
	      const BinOp& func) {
      typedef typename BinOp::second_argument_type Arg2;
      int i;

      checkTS(TSARG); a.checkTS(TSARG); b.checkTS(TSARG);
      if (a.size()==0 || b.size()==0)
	throw DimMismatchException(EXCEPT_MSG("a,b must not be empty"));
      Handle<BaseMatrix<T> > tempMat;
      BaseMatrix<T>* trgMat;
      if (!doCheckValid()) trgMat=this;
      else {
	// Build up in temp. matrix to allow for element check before
	// assignment
	tempMat.changeRep(new BaseMatrix<T>());
	trgMat=tempMat;
      }
      trgMat->ensureCapacity(a.size(),b.size());
      if (!trgMat->isMaskInd()) {
	// Flat matrix
	T* trgB=trgMat->buff;
	for (i=0; i<trgMat->n; i++,trgB+=trgMat->stride) {
	  ArrayUtils<T>::applyFunc(trgB,a.getBuffPtr_INT(),trgMat->m,
				   bind2nd(func,(Arg2) b[i]),1,a.getStep_INT(),
				   0,a.getMindex_INT());
	}
      } else {
	// Use mask vectors for columns
	Handle<BaseVector<T> > trgM;
	BaseVector<T2> aM; // cannot use int. 'subrvalInt' for this one
	Range rng(0,trgMat->m-1);
	for (i=0; i<trgMat->n; i++) {
	  trgM.changeRep(subrvalInt(rng,i));
	  aM.reassign(a(rng));
	  trgM->apply1(aM,bind2nd(func,(Arg2) b[i]));
	}
      }
      if (doCheckValid())
	// copy back (forces element check)
	assignInt(*tempMat);
    }

    /**
     * Specialized method for selection and accumulation.
     * Matrix version of 'BaseVector::accuMap'. Here, we have
     * two maps 'rmap' (rows), 'cmap' (cols). Either has an odd number of
     * entries describing a map of indexes (see 'BaseVector::accuMap').
     * We run over the columns as given by 'cmap' and apply
     * 'BaseVector::accuMap' with 'rmap' as map to the corr. columns.
     * NOTE: If either 'rmap', 'cmap' is empty, nothing is done here.
     * NOTE: Last entry of 'cmap' is not used, can be arbitrary.
     *
     * @param a    Source matrix
     * @param rmap S.a.
     * @param cmap S.a.
     */
    virtual void accuMap(const BaseMatrix<T>& a,const ArrayHandle<int>& rmap,
			 const ArrayHandle<int>& cmap);

    // Rearrangement methods

    /**
     * Transpose: this=A'
     * <p>
     * NOTE: If 'a' is the same as this matrix, 'trans()' is called.
     *
     * @param a A
     */
    virtual void trans(const BaseMatrix<T>& a);

    /**
     * Transpose (in-place): this=this^T
     * <p>
     * This is trivial for a square matrix. If this is not square, we use
     * the algorithm implemented in 'TransInPlace', which is ACM algorithm
     * 467.
     * NOTE: Since the algorithm requires a flat contiguous buffer without
     * gaps, we have to use 'reshapeSave' to first compress the buffer to just
     * contain this matrix, apply the algorithm, then "relax" the buffer
     * size again to its old values. In fact, the new 'stride' is the old
     * buffSz/stride.
     * <p>
     * NOTE: If 'notInPl'==true, the transposition is not done in-place for
     * non-square matrices, but a temp. copy is used.
     * ==> Use if somewhat unpred. running time for ACM algorithm is a
     *     problem (usually not)
     *
     * @param notInPl S.a. Def.: false
     */
    virtual void trans(bool notInPl=false);

    /**
     * Permutes columns/rows of this matrix (in place)
     * <p>
     * If 'rowMaj'==false (def.), columns are permuted, otherwise rows.
     * The operation moves col/row i to pos. 'perm[i]'. The method uses
     * an additional T and a bool vector of the col/row size.
     * <p>
     * If Pi is the permut. matrix for 'perm', for the row version the result
     * is Pi A (A this matrix), for the column version the result is
     * A Pi^T: row/col i in A moves to row/col pi(i) in the result.
     * <p>
     * ATTENTION: 'perm' must be a permutation of the range of cols/rows,
     * otherwise the method fails. Failure may or may not be detected.
     * <p>
     * NOTE: If p,ip are inverses of each other,
     *   b=a(RangeFull::get(),Range(p));
     * does the same as
     *   b=a; b.permute(ip);
     * The difference is that 'permute' can be used in place, while
     * something like
     *   a=a(RangeFull::get(),p);
     * fails with exception or leads to a wrong result! See 'operator='
     * comments (in lval=rval, lval and rval must NOT refer to the same
     * underlying object).
     *
     * @param perm   S.a.
     * @param rowMaj Permute rows? Def.: false
     */
    virtual void permute(const BaseVector<int>& perm,bool rowMaj=false);

    /**
     * Creates symmetric matrix by copying the lower triangular part of
     * this matrix onto the upper triangular part (which is overwritten).
     * <p>
     * This matrix must be square. Computations res. in symmetric matrices
     * can often be done in half the time by only computing the lower
     * triangular part and the diagonal, then using this method.
     * If 'low'==false, the upper triangle is copied onto the lower triangle.
     *
     * @param low S.a. Def.: true
     */
    virtual void makeSymm(bool low=true);

    // Methods requiring order predicates

    /**
     * Checks whether all elements fall within the interval given by 'ival'
     * (see class 'Interval'). If not, then the pos. of the first element
     * outside of 'ival' can be returned via 'pos', and a status can be ret.
     * via 'stat': 1 -> elem. too small, 2 -> elem. too large. Here, elements
     * are compared in column major ordering, left-right, top-down. Pos.
     * (i,j) is encoded as 'pos'==i+'m'*j.
     * <p>
     * NOTE: Req. '<' and '==' operators.
     *
     * @param ival Interval range
     * @param pos  Optional. See above
     * @param stat "
     * @return     All elements within interval?
     */
    virtual bool checkBounds(const Interval<T>& ival,int* pos=0,int* stat=0)
      const;

    /**
     * Comparison operator
     *
     * @param a Other matrix
     * @return  Are this matrix and 'a' identical?
     */
    virtual bool operator==(const BaseMatrix<T>& a) const;

    /**
     * Computes maximum element in each column and returns them as row
     * vector. If 'rowMaj'==true, the max. elem. along rows are comp.
     *
     * @param rowMaj S.a. Def.: false
     * @return       S.a.
     */
    BaseVector<T> max(bool rowMaj=false) const {
      BaseVecWrapper<T> resWrap(newEmptyVector());
      max(*(resWrap.getRep()),rowMaj);

      return resWrap;
    }

    /**
     * Computes maximum element in every col. and returns them in a row
     * vector. If 'rowMaj'==true, the max. elem. along rows are comp.
     *
     * @param vec    Vector to return max. elements in
     * @param rowMaj S.a. Def.: false
     */
    virtual void max(BaseVector<T>& vec,bool rowMaj=false) const {
      AccumMaxVal<T> accu;
      accumulate(accu,vec,rowMaj);
    }

    /**
     * Computes minimum element in each column and returns them as row
     * vector. If 'rowMaj'==true, the min. elem. along rows are comp.
     *
     * @param rowMaj S.a. Def.: false
     * @return       S.a.
     */
    BaseVector<T> min(bool rowMaj=false) const {
      BaseVecWrapper<T> resWrap(newEmptyVector());
      min(*(resWrap.getRep()),rowMaj);

      return resWrap;
    }

    /**
     * Computes minimum element in every col. and returns them in a row
     * vector. If 'rowMaj'==true, the min. elem. along rows are comp.
     *
     * @param vec    Vector to return min. elements in
     * @param rowMaj S.a. Def.: false
     */
    void min(BaseVector<T>& vec,bool rowMaj=false) const {
      AccumMinVal<T> accu;
      accumulate(accu,vec,rowMaj);
    }

    /**
     * Same as 'max', but the pos. of the max. elems. are returned instead
     * of the values.
     *
     * @param rowMaj S.a. Def.: false
     * @return       S.a.
     */
    virtual BaseVector<int> maxPos(bool rowMaj=false) const {
      BaseVecWrapper<int> resWrap(new BaseVector<int>());
      maxPos(*(resWrap.getRep()),rowMaj);

      return resWrap;
    }

    /**
     * Computes pos. of maximum element in every column and returns them in a
     * index vector. If 'rowMaj'==true, the max. pos. along rows are comp.
     *
     * @param vec    Vector to return pos. of max. elems. in
     * @param rowMaj S.a. Def.: false
     */
    virtual void maxPos(BaseVector<int>& vec,bool rowMaj=false) const {
      AccumMaxPos<T> accu;
      accumulate(accu,vec,rowMaj);
    }

    /**
     * Same as 'min', but the pos. of the min. elems. are returned instead
     * of the values.
     *
     * @param rowMaj S.a. Def.: false
     * @return       S.a.
     */
    virtual BaseVector<int> minPos(bool rowMaj=false) const {
      BaseVecWrapper<int> resWrap(new BaseVector<int>());
      minPos(*(resWrap.getRep()),rowMaj);

      return resWrap;
    }

    /**
     * Computes pos. of minimum element in every column and returns them in a
     * index vector. If 'rowMaj'==true, the min. pos. along rows are comp.
     *
     * @param vec    S.a.
     * @param rowMaj S.a. Def.: false
     */
    virtual void minPos(BaseVector<int>& vec,bool rowMaj=false) const {
      AccumMinPos<T> accu;
      accumulate(accu,vec,rowMaj);
    }

    /**
     * Finds pos. of first elem. for which 'pred' is true. The position
     * (row,col) is returned as pair. If such an elem. does not exist,
     * (-1,-1) is returned.
     * The matrix is traversed column-major if 'rowMaj'==false, row-major
     * otherwise.
     *
     * @param pred   Unary predicate
     * @param rowMaj S.a. Def.: false
     * @return       S.a.
     */
    template<class UnOp> pair<int,int> findfst(const UnOp& pred,
					       bool rowMaj=false) const {
      int i,j;
      Handle<BaseVector<T> > srcM;
      pair<int,int> ret(-1,-1);

      checkTS(TSARG);
      if (m==0)
	throw WrongDimensionException(EXCEPT_MSG("Matrix must not be empty"));
      if (!rowMaj) {
	// Column-major
	Range rng(0,m-1);
	for (i=0; i<n; i++) {
	  srcM.changeRep(subrvalInt(rng,i));
	  if ((j=srcM->findfst(pred))!=-1) {
	    ret.first=j; ret.second=i;
	    break;
	  }
	}
      } else {
	// Row-major
	Range rng(0,n-1);
	for (i=0; i<m; i++) {
	  srcM.changeRep(subrvalInt(i,rng));
	  if ((j=srcM->findfst(pred))!=-1) {
	    ret.first=i; ret.second=j;
	    break;
	  }
	}
      }

      return ret;
    }

    // Initialization methods

    /**
     * Initializes matrix with element 'a'. With no further argument, the
     * original dimensions are kept.
     *
     * @param rows Optional. Number of rows
     * @param cols Optional. Number of cols. Def.: same as rows
     * @param a    Optional. Fill element. Def.: def. fill value
     */
    virtual void fill() {
      fill(m,n,this->defFill);
    }

    virtual void fill(int rows,int cols,T a);

    virtual void fill(int rows,int cols=-1) {
      if (cols==-1)
	fill(rows,rows,this->defFill);
      else
	fill(rows,cols,this->defFill);
    }

    // IO methods

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
     * Prints matrix into output stream 'os', one row per line.
     * There is no newline after the last row.
     * <p>
     * NOTE: Requires the '<<' operator to be implemented for T.
     * Override if the def. implementation does not make sense!
     *
     * @param os Output stream
     */
    virtual void print(ostream& os) const;

#ifdef MATLAB_MEX
    /**
     * Print to stdout (MEX file). Requires that T can be cast to double.
     */
    virtual void matlabPrint() const {
      int i,j;

      for (i=0; i<m; i++) {
	for (j=0; j<n-1; j++)
	  mexPrintf("%f ",(double) get(i,j));
	if (i<m-1) mexPrintf("%f\n",(double) get(i,n-1));
	else mexPrintf("%f",(double) get(i,n-1));
      }
    }
#endif

    /**
     * Stores matrix into (binary) output file stream. Generic format
     * (FF version 1):
     * - tag @BaseMatrix (not written if 'noTag'==true)
     * - FF version [int(4)]
     * - type code [int(4)]
     * - matrix [using 'saveInt']
     * Files with FFV 0 did not store the byte size of T (see 'saveInt').
     *
     * @param os Binary output file stream
     * @param noTag Do not save tag?
     */
    virtual ofstream& save(ofstream& os,bool noTag=false) const;

    /**
     * Loads matrix from (binary) input file stream.
     * <p>
     * DOWNW. COMPAT.: Files written by earlier LHOTSE versions had
     * wrong type code fields ('typeOther' for elementary T). If
     * this is encountered, a warning message is printed, rather
     * than throwing an exception (see 'checkTypeCode').
     *
     * @param is Binary input file stream
     * @param noTag Do not read tag?
     */
    virtual ifstream& load(ifstream& is,bool noTag=false);

    // Methods for 'TempMatMethods' superclass

    /*
     * 'assignVirtual' must do the same as 'operator=' !
     */

    void assignVirtual(const TempMatMethods* arg) {
      //cout << "BaseMatrix::assignVirtual" << endl;
      // Is 'arg' a matrix?
      const BaseMatrix<T>* aptr=DYNCAST(const BaseMatrix<T>,arg);
      if (aptr==0) {
	// Is 'arg' a vector?
	const BaseVector<T>* avptr=DYNCAST(const BaseVector<T>,arg);
	if (avptr==0) throw InternalException(EXCEPT_MSG(""));
	assignInt(*avptr);
      } else
	assignInt(*aptr);
    }

    void assignVirtual(const TempMatMethods_AssignElement_Generic* arg) {
      //cout << "BaseMatrix::assignVirtual" << endl;
      const TempMatMethods_AssignElement<T>* aptr=
	DYNCAST(const TempMatMethods_AssignElement<T>,arg);
      if (aptr==0) throw InternalException(EXCEPT_MSG(""));
      assignInt(aptr->elem);
    }

    // Methods for superclass 'WriteBackMat'

    void writeBackBuff(const ArrayHandle<T>& buff,unsigned char strpatt) {
      BaseMatrix<T> temp;
      temp.reassign(buff,m,n,m); // contiguous matrix -> stride==m
      assignInt(temp,strpatt);
    }

    // Methods for structure pattern support (in subclasses)

    /**
     * This method must be overwritten by subclasses which support
     * structure patterns (such as 'StMatrix'). The default implementation
     * just returns 'normal'.
     *
     * @return Structure pattern (see 'MatStrct')
     */
    virtual uchar getStrctPatt() const {
      return MatStrct::normal;
    }

    /*
     * "Internal" methods which are declared public:
     * ATTENTION: These methods must only be called by methods of this class
     * or subclasses!
     * Why do we need them? Template methods have to call them for
     * 'BaseMatrix<T2>' with T2!=T, and it seems we cannot be friend of
     * 'BaseMatrix<T2>' for arb. T2 (SOLUTION???).
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
     * @return Value of striding constant 'stride'
     */
    int getStride_INT() const {
      return stride;
    }

    /*
     * Methods kept (in the moment) for downwards compatibility
     * ==> DO NOT USE THEM!
     */
    
    /** DOWNWARDS COMPAT. DO NOT USE
     * Reassigns mask vector 'mVec' (might also be an empty vector) to point
     * to row 'row' of this matrix (or a part of it).
     *
     * @param mVec  Mask vector
     * @param row   Row number
     * @param begin Optional. First column. Default is 0
     * @param end   Optional. Last column. Default is last col. of matrix
     */
    void maskRow(BaseVector<T>& mVec,int row,int begin=0,int end=-1) const {
      Handle<BaseVector<T> > tvec(subrvalInt(row,Range(begin,end)));
      mVec.reassign(*tvec);
    }

    /** DOWNWARDS COMPAT. DO NOT USE
     * Reassigns this mask matrix to point to matrix 'mat' (or a submatrix
     * thereof).
     *
     * @param mat    Matrix to point to
     * @param beginR First row
     * @param beginC Optional. First column. Default is 0
     * @param endR   Optional. Last row. Default is last row of 'mat'
     * @param endC   Optional. Last column. Default is last column of 'mat'
     */
    void reassign(const BaseMatrix<T>& mat,int beginR,int beginC=0,
		  int endR=-1,int endC=-1) {
      reassign(mat,Range(beginR,endR),Range(beginC,endC));
    }

    /** DOWNW. COMPAT. DO NOT USE
     * Stores copy of (part of a) row of this matrix into vector.
     *
     * @param vec   Target vector
     * @param row   Row number
     * @param begin Optional. First column. Default is 0
     * @param end   Optional. Last column. Default is last col. of matrix
     */
    void getRow(BaseVector<T>& vec,int row,int begin=0,int end=-1) const {
      Handle<BaseVector<T> > msk(subrvalInt(row,Range(begin,end)));
      vec=msk;
    }

    /** DOWNW. COMPAT. DO NOT USE
     * Stores copy of (part of a) column of this matrix into vector.
     *
     * @param vec   Target vector
     * @param col   Column number
     * @param begin Optional. First row. Default is 0
     * @param end   Optional. Last row. Default is last row of matrix
     */
    void getCol(BaseVector<T>& vec,int col,int begin=0,int end=-1) const {
      Handle<BaseVector<T> > msk(subrvalInt(Range(begin,end),col));
      vec=msk;
    }

    /** DOWNW. COMPAT. DO NOT USE
     * Assigns row of this matrix with elements of given vector. The matrix
     * is enlarged and filled with zeros if necessary. Old entries not covered
     * by the new row are not altered
     *
     * @param vec    Vector cont. elements for new row
     * @param row    Row number
     * @param beginC Column offset
     */
    void setRow(const BaseVector<T>& vec,int row,int beginC=0) {
      Handle<BaseVector<T> > msk(sublvalInt(row,Range(beginC,
						      beginC+vec.size()-1)));
      msk=vec;
    }

    /** DOWNW. COMPAT. DO NOT USE
     * Assigns col of this matrix with elements of given vector. The matrix
     * is enlarged and filled with zeros if necessary. Old entries not covered
     * by the new column are not altered
     *
     * @param vec    Vector cont. elements for new column
     * @param col    Column number
     * @param beginR Row offset
     */
    void setCol(const BaseVector<T>& vec,int col,int beginR=0) {
      Handle<BaseVector<T> > msk(sublvalInt(Range(beginR,beginR+vec.size()-1),
					    col));
      msk=vec;
    }

  protected:
    // Internal methods

    /*
     * NOTE: Some internal methods do NOT call 'checkTS', and they also do not
     *       in general call 'incrTS' when changing the vector size. Check
     *       comments ("TS not served").
     */

    /**
     * Called by "copy" constructor from 'BaseMatWrapper'. If things change
     * in subclasses, this method has to be overridden!
     * 'arg' has wrapped a temp. matrix. We initialise all fields of this
     * object with the fields of this vector. Then, the fields of the temp.
     * matrix are reset, so that 'arg' is disposed of.
     * 
     * NOTE: This method must be called only for an "empty" matrix (either
     * from a constructor or after 'dealloc' has been used).
     *
     * @param arg Wrapped "source" object
     */
    virtual void convertWrapped(BaseMatWrapper<T>& arg);

    /**
     * Draws a copy of this object. The return object is always normal.
     * Same as the copy constructor, but polymorphic.
     * NOTE: Dyn. created on local heap, wrap into handle!
     * ==> Subclasses have to overwrite this method!
     *
     * @return Copy
     */
    virtual BaseMatrix<T>* copy() const {
      return new BaseMatrix<T>(*this);
    }

    /**
     * Creates an empty matrix of the type of this class. Used in virtual
     * methods. NOTE: Dyn. created on local heap, wrap into handle!
     * ==> Subclasses have to overwrite this method!
     *
     * @return See above
     */
    virtual BaseMatrix<T>* newEmpty() const {
      return new BaseMatrix<T>();
    }

    /**
     * Creates empty vector of the type corr. to rows/cols of this matrix
     * class (operator() returns mask vectors of this type).
     * NOTE: If this class has element restrictions, the corr. vector class
     * must have the same or weaker ones!
     * ==> Subclasses have to overwrite this method!
     *
     * @return See above
     */
    virtual BaseVector<T>* newEmptyVector() const {
      return new BaseVector<T>();
    }

    /**
     * Does job of assignment operator =. In contrast to 'operator=',
     * this is a virtual method.
     * <p>
     * This matrix can be a mask. In this case, if 'mat' is smaller than
     * this one, it is copied into the upper left corner. If it is larger,
     * an exception is thrown. If this matrix is normal, its size will be
     * the same as 'mat' afterwards.
     * <p>
     * Structure patterns:
     * By def., the complete 'mat' is copied into this matrix. If 'srcstrct'
     * is used, a smaller pattern of 'mat' is copied only. If 'trgstrct' is
     * not given, it is the same as 'srcstrct'. Otherwise, 'trgstrct' and
     * 'srcstrct' must be compatible, i.e. describe structures of the same
     * size (but can be upper/lower). Ex.:
     *   a.assignInt(b,MatStrct::lower);
     * copies lower triangle, diagonal of 'b' into lower triangle, diagonal of
     * 'a'.
     *   a.assignInt(b,MatStrct::lower,
     *                 MatStrct::upper);
     * copies lower triangle, diagonal of 'b' into upper triangle, diagonal of
     * 'a'.
     * NOTE: If 'srcstrct' ref. to a pattern for square matrices, 'mat' must be
     * square. If this matrix is a mask and has place for 'mat', it can be
     * non-square.
     * NOTE: Overlapping buffers are NOT supported (not checked)!
     *
     * @param mat      Source matrix
     * @param srcstrct S.a. Def.: normal
     * @param trgstrct S.a. Def.: same as 'srcstrct'
     */
    virtual void assignInt(const BaseMatrix<T>& mat,
			   uchar srcstrct,uchar trgstrct);

    virtual void assignInt(const BaseMatrix<T>& mat,
			   uchar srcstrct=MatStrct::normal) {
      assignInt(mat,srcstrct,srcstrct);
    }

    /**
     * Does job of assignment operator =. In contrast to 'operator=',
     * this is a virtual method.
     * <p>
     * The argument 'vec' is a vector, which is taken to be a column vector.
     * If this matrix is a mask, it has to have a single column which is
     * >= the size of 'vec', otherwise an exception is thrown. 'vec' is copied
     * into the beginning of the column.
     * If this matrix is normal, it will be a column vector of the size of
     * 'vec' afterwards.
     * <p>
     * NOTE: Overlapping buffers are not supported!
     *
     * @param vec Source vector (taken as column)
     */
    virtual void assignInt(const BaseVector<T>& vec);

    /**
     * Does job of assignment operator =. In contrast to 'operator=',
     * this is a virtual method.
     * <p>
     * Structure pattern:
     * If 'trgstrct' is used, only the corr. structure pattern (see
     * 'strctXXX' constants) of this matrix is filled with 'elem'. If the
     * pattern is for a square matrix, this matrix must be square.
     *
     * @param elem     Element to fill matrix with
     * @param trgstrct S.a.
     */
    virtual void assignInt(T elem,
			   unsigned char trgstrct=
			   MatStrct::normal);

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
     * <p>
     * ATTENTION: This def. implementation for unconstrained element types
     * has to be overwritten in subclasses with constraints!
     *
     * @param vecP See above
     * @return     "
     */
    virtual bool isValidMatType(const BaseMatrix<T>* vecP) const {
      return true;
    }

    /**
     * See header comment. Returns true if 'vecP' points to an object whose
     * type enforces the required validity constraints. In this case, checking
     * each element is not necessary for 'vecP'.
     * <p>
     * ATTENTION: This def. implementation for unconstrained element types
     * has to be overwritten in subclasses with constraints!
     *
     * @param vecP See above
     * @return     "
     */
    virtual bool isValidVecType(const BaseVector<T>* vecP) const {
      return true;
    }

    /**
     * See header comment. Returns true iff the given matrix is s.t. it
     * could be inserted into this matrix. The matrix to be tested is given
     * via buffer point 'pbuff', size ('pm','pn'), striding const. 'pstride'.
     * <p>
     * NOTE: The def. impl. just calls 'isValidElement' for each element and
     * is probably sufficient in all cases.
     *
     * @param pbuff   S.a.
     * @param pm      S.a.
     * @param pn      S.a.
     * @param pstride S.a.
     * @return        Is it OK?
     */
    virtual bool isValidMat(const T* pbuff,int pm,int pn,int pstride) const;

    /**
     * Same as above, but the matrix is obtained by applying row/col range
     * 'rngR',rngC' to the frame with 'strid' rows starting at 'pbuff' (see
     * 'reassign').
     * NOTE: Ranges not checked for validity (must  not be open)!
     *
     * @param pbuff Buffer
     * @param strid See above
     * @param rngR  "
     * @param rngC  "
     * @return      Is it OK?
     */
    virtual bool isValidMat(const T* pbuff,int strid,const Range& rngR,
			    const Range& rngC) const;

    /**
     * See other 'isValidMat' version.
     * NOTE: If 'doCheck' is false (def.), we call 'isValidMatType' first.
     * If this returns true, the matrix is autom. valid.: no element checks
     * are done. If 'doCheck'==true, element checks are always done.
     * <p>
     * Structure pattern: If 'strct' is given, only the elements of the corr.
     & structure pattern (see 'MatStrct').
     * <p>
     * TS not served.
     *
     * @param mat     Matrix to check
     * @param doCheck S.a. Def.: false
     * @param strct   S.a. Def.: full matrix
     * @return        See above
     */
    virtual bool isValidMat(const BaseMatrix<T>& mat,bool doCheck=false,
			    unsigned char strct=MatStrct::normal)
      const;

    /**
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
    virtual bool isValidVec(const BaseVector<T>& vec,bool doCheck=false)
      const;

    /**
     * Same as 'isValidMat' above, but here every element is first cast from
     * T2 to T before it is checked.
     *
     * @param mat S.a.
     * @return    S.a.
     */
    template<class T2> bool
    isValidConvMat(const BaseMatrix<T2>& mat) const {
      int i,j;

      for (i=0; i<mat.cols(); i++)
	for (j=0; j<mat.rows(); j++)
	  if (!isValidElement((T) mat.get(j,i))) return false;
      return true;
    }

    /*
     * Required by 'MatTimeStamp' superclass.
     */
    const MatTimeStamp* getBaseObj() const {
      return isMask()?baseObj:0;
    }

    /**
     * Does job of 'reassign', but without checking validity of the elements
     * in 'pbuff'.
     *
     * @param pbuff See above
     * @param pm    "
     * @param pn    "
     * @param strid "
     * @param watch "
     */
    virtual void reassignInt(const T* pbuff,int pm,int pn,int strid,
			     MemWatchBase* watch=0);

    /**
     * Does job of 'reassign', but with argument, element or range checks.
     * Both ranges must be closed.
     *
     * @param pbuff S.a.
     * @param strid "
     * @param rngR  "
     * @param rngC  "
     * @param watch "
     */
    virtual void reassignInt(const T* pbuff,int strid,const Range& rngR,
			     const Range& rngC,MemWatchBase* watch=0);

    /**
     * Does the job of 'reassign', but without any checks.
     *
     * @param mat  S.a.
     * @param rngR "
     * @param rngC "
     */
    virtual void reassignInt(const BaseMatrix<T>& mat,
			     const Range& rngR=RangeFull::get(),
			     const Range& rngC=RangeFull::get());

    /**
     * Does the job of 'reassign', but without any checks.
     *
     * @param vec S.a.
     * @param col S.a.
     */
    virtual void reassignInt(const BaseVector<T>& vec,bool col);

    /**
     * Generic code for 'operator()', r-value version.
     *
     * @param rngR See operator()
     * @param rngC "
     */
    virtual BaseMatrix<T>* subrvalInt(const Range& rngR,const Range& rngC)
      const;
    virtual BaseVector<T>* subrvalInt(int rngR,const Range& rngC) const;
    virtual BaseVector<T>* subrvalInt(const Range& rngR,int rngC) const;

    /**
     * Generic code for 'operator()', l-value version.
     *
     * @param rngR See operator()
     * @param rngC "
     */
    virtual BaseMatrix<T>* sublvalInt(const Range& rngR,const Range& rngC);
    virtual BaseVector<T>* sublvalInt(int rngR,const Range& rngC);
    virtual BaseVector<T>* sublvalInt(const Range& rngR,int rngC);

    /**
     * Generic code for 'diag'.
     *
     * @param off S.a.
     * @param rng S.a.
     * @return    Dyn. alloc. mask vector
     */
    virtual BaseVector<T>* diagInt(int off,const Range& rng) const;

    /**
     * Writes repres. of object into (binary) file stream 'os'.
     * Format of default implementation:
     * - number of rows [int(4)]
     * - number of cols [int(4)]
     * - byte size of T [int(4)]: Only if 'storeBSize'==true!
     * - content (col-major) [using 'NumberFormats<T>::save']
     * Override this by specialisation if not appropriate!
     * <p>
     * NOTE: No tag is written here, and no type info.
     *
     * @param os         Output file stream
     * @param storeBSize Store byte size T in file? Def.: true
     */
    virtual ofstream& saveInt(ofstream& os,bool storeBSize=true) const;

    /**
     * Read repres. of object from (binary) file stream 'is'.
     * NOTE: No tag is read here, and no type info.
     * NOTE: If 'doCheckValid' is true, the elements are checked after
     * loading. This is because many subclasses use the 'BaseMatrix'
     * format even though they implement element checks.
     * ==> In this case, the content of this matrix is invalid!
     * For byte size of T, see 'saveInt'.
     *
     * @param is        Input file stream
     * @param loadBSize Load byte size of T? Def.: true
     */
    virtual ifstream& loadInt(ifstream& is,bool loadBSize=true);

    /**
     * Allocates memory for matrix and sets all entries to 'a'. An existing
     * buffer is dealloc. first.
     * NOTE: TS not served.
     *
     * @param rows    Number of rows. Must be >=0
     * @param cols    Number of columns. Must be >=0
     * @param debStat For debugging
     * @param a       Optional. Def.: def. fill value
     */
    virtual void init(int rows,int cols,uchar debStat,T a);

    virtual void init(int rows,int cols,uchar debStat) {
      init(rows,cols,debStat,this->defFill);
    }

    /**
     * Allocates memory for matrix. If a matrix buffer exists, it is
     * deallocated first. The matrix entries are not initialized. The matrix
     * size is the same as the buffer size, 'stride' is set to 'rows'.
     * NOTE: TS not served.
     *
     * @param rows Number of rows. Must be >= 0
     * @param cols Number of columns. Must be >= 0
     */
    virtual void alloc(int rows,int cols,uchar debStat);

    /**
     * Ensures that the buffer has a given size (or is larger). If the buffer
     * is less than the given size, it's reallocated with the given size.
     * The elements are not initialized. Old elements are not retained. The
     * matrix size is set to (rows,cols). If either 'rows'==0 or 'cols'==0, the
     * matrix becomes the empty one.
     * <p>
     * New buffer size and dimensions:
     * If 'rows'*'cols' > 'buffSz', the buffer is realloc. to have size
     * 'rows'-by-'cols'. Otherwise, the flat buffer is kept. In that case,
     * 'stride' is chosen as follows:
     * - If the new matrix is empty, 'stride' keeps its old value
     * - If 'stride' < 'rows': 'stride'='rows'
     * - Oth.: If 'stride'*'cols' <= 'buffSz': leave 'stride' as is; oth.
     *         choose 'stride' max. s.t. 'stride'*'cols' <= 'buffSz'
     * Thus, we try to keep the old 'stride' value if possible.
     * NOTE: If 'strict'==true, we always set 'stride'='rows'.
     * NOTE: If the flat buffer is not realloc., its entries are not changed
     * in any way. If it is realloc., the entries are not initialized.
     * <p>
     * NOTE: The timestamp is increased iff the matrix size changes.
     *
     * @param rows   Minimum number of rows
     * @param cols   Minimum number of columns
     * @param strict See above. Def.: false
     */
    virtual void ensureCapacity(int rows,int cols,bool strict=false);

    /**
     * Deallocates resources of matrix. The result is an empty matrix without
     * buffers.
     * NOTE: Works on mask matrices too.
     * NOTE: TS not served.
     *
     * @param noReset If true, we just call 'deassocBuff' for buffer de-assoc.,
     *                but do not change any other members. Def.: false
     */
    virtual void dealloc(bool noReset=false);

    /**
     * Helper for resizing methods. Copies matrix frame (srcBuff,srcM,srcN,
     * srcStride) to matrix frame (trgBuff,srcM,srcN,trgStride). Copying is
     * done column by column, top-down (starting from first column) if
     * 'top'==true (def.), bottom-up (starting from last column) oth.
     * The columns can overlap. For overlapping buffers, 'top' has to be
     * chosen correctly.
     *
     * @param srcBuff   S.a.
     * @param srcM      S.a.
     * @param srcN      S.a.
     * @param srcStride S.a.
     * @param trgBuff   S.a.
     * @param trgStride S.a.
     * @param top       S.a. Def.: true
     */
    void moveFrame(const T* srcBuff,int srcM,int srcN,int srcStride,
		   T* trgBuff,int trgStride,bool top=true) const;

    /**
     * Helper for resizing methods. Re-allocates the buffer of this matrix
     * to dimensions 'rows'-by-'cols' (size 'rows'*'cols') while keeping the
     * old buffer, then copies matrix content across and de-assoc. from old
     * buffer. Need 'rows'>=m, 'cols'>=n.
     * NOTE: TS not served.
     *
     * @param rows S.a.
     * @param cols S.a.
     */
    void reallocBuff(int rows,int cols,uchar debStat);
  };

  // Public methods

  template<class T>
  void BaseMatrix<T>::reassign(const BaseMatrix<T>& mat,const Range& rngR,
			       const Range& rngC)
  {
    if (!isMaskOrEmpty())
      throw MaskObjectException(EXCEPT_MSG("Has to be mask or empty"));
    if (rngR.checkRange(mat.m) || rngC.checkRange(mat.n))
      throw OutOfRangeException(EXCEPT_MSG(""));
    if (!rngR.isUniqueMap() || !rngC.isUniqueMap())
      throw InvalidParameterException(EXCEPT_MSG("Ranges must be 1-1"));
    mat.checkTS(TSARG);
    if (mat.n==0)
      dealloc(); // empty matrix
    else {
      if (!doCheckValid() || isValidMatType(&mat)) {
	// No elem. checks necessary
	reassignInt(mat,rngR,rngC);
      } else {
	// Have to check elems:
	// We use 'BaseMatWrapper' here to first create the mask to be
	// assigned to this matrix. If the element check fails, this matrix is
	// not changed.
	// 'tempMat' has type 'BaseMatrix'. 'isValidMat' calls
	// 'isValidMatType' internally, but this will result in false, since
	// 'BaseMatrix' does not enforce any restrictions.
	BaseMatWrapper<T> tvWrap(new BaseMatrix<T>());
	BaseMatrix<T>* tempMat=tvWrap.getRep();
	tempMat->reassignInt(mat,rngR,rngC);
	// Check elements:
	if (!isValidMat(*tempMat))
	  throw WrongTypeException(EXCEPT_MSG("'mat' contains invalid elements"));
	// Use fields of 'tempMat' for this vector
	dealloc(); // remove old stuff
	convertWrapped(tvWrap); // take over fields, clean-up
      }
    }
  }

  template<class T>
  void BaseMatrix<T>::reassign(const BaseVector<T>& vec,bool col)
  {
    if (!isMaskOrEmpty())
      throw MaskObjectException("Has to be mask or empty");
    vec.checkTS(TSARG);
    if (vec.size()==0)
      dealloc(); // empty matrix
    else {
      if (!doCheckValid()) {
	// No elem. checks necessary
	reassignInt(vec,col);
      } else {
	// Have to check elems:
	// We use 'BaseMatWrapper' here to first create the mask to be
	// assigned to this matrix. If the element check fails, this matrix is
	// not changed.
	// 'tempMat' has type 'BaseMatrix'. 'isValidMat' calls
	// 'isValidMatType' internally, but this will result in false, since
	// 'BaseMatrix' does not enforce any restrictions.
	BaseMatWrapper<T> tvWrap(new BaseMatrix<T>());
	BaseMatrix<T>* tempMat=tvWrap.getRep();
	tempMat->reassignInt(vec,col);
	// Check elements:
	if (!isValidMat(*tempMat))
	  throw WrongTypeException(EXCEPT_MSG("'vec' contains invalid elements"));
	// Use fields of 'tempMat' for this vector
	dealloc(); // remove old stuff
	convertWrapped(tvWrap); // take over fields, clean-up
      }
    }
  }

  template<class T> inline
  T& BaseMatrix<T>::operator()(int rngR,int rngC)
  {
    if (rngR<0 || rngR>=m || rngC<0 || rngC>=n)
      throw OutOfRangeException(EXCEPT_MSG(""));
    checkTS(TSARG);
    if (!isMaskInd())
      return *(buff+(rngC*stride+rngR));
    else
      return *(buff+(COLPOS(rngC)*stride+ROWPOS(rngR)));
  }

  template<class T> inline ArrayHandle<T> BaseMatrix<T>::getFlatBuff() const
  {
    checkTS(TSARG);
    if (m==0)
      return ArrayHandleZero<T>::get(); // zero
    else if (isFlat()) {
      // flat and contiguous
      return ArrayHandle<T>(buff,m*n,buffWatch); // int. constr.
    } else {
      // Draw flat copy
      ArrayHandle<T> flatCp(m*n);
      int i;
      T* trgB=flatCp;
      if (!isMaskInd())
	moveFrame(buff,m,n,stride,trgB,m);
      else {
	for (i=0; i<n; i++,trgB+=m)
	  ArrayUtils<T>::copy(trgB,buff+COLPOS(i)*stride,m,1,1,0,rowInd);
      }
      return flatCp;
    }
  }

  template<class T> inline void
  BaseMatrix<T>::getLinMat(BaseLinMat<T>& blMat,bool writeBack,
			   unsigned char strpatt) const
  {
    if (m==0)
      throw WrongDimensionException(EXCEPT_MSG("Matrix must not be empty"));
    checkTS(TSARG);
    if (isMaskInd() && writeBack) {
      if (strpatt>MatStrct::last)
	throw InvalidParameterException("'strpatt'");
      if (strpatt!=MatStrct::normal && m!=n)
	throw InvalidParameterException(EXCEPT_MSG("Matrix must be square for this 'strpatt'"));
    }
    if (!isMaskInd()) {
      // Flat or linear matrix
      ArrayHandle<T> hand(buff,(n-1)*stride+m,buffWatch); // int. constr.
      blMat.init(m,n,hand,stride,strpatt); // write-back automatic
    } else {
      // Indexed mask: draw flat copy
      blMat.init(m,n,getFlatBuff(),m,strpatt,
		 writeBack?((BaseMatrix<T>*) this):0);
    }
  }

  template<class T> inline void BaseMatrix<T>::resize(int rows,int cols)
  {
    if (isMask()) throw MaskObjectException(); // not allowed for mask matrices
    if (rows==0 || cols==0)
      // Deallocate only
      dealloc();
    else
      init(rows,cols,13);
    incrTS(); // size (or buffer) change
  }

  template<class T> inline void BaseMatrix<T>::resizeSave(int rows,int cols)
  {
    int i;

    if (isMask()) throw MaskObjectException();
    if (rows==-1) rows=m;
    if (cols==-1) cols=n;
    if (m>rows || n>cols) throw WrongDimensionException(EXCEPT_MSG(""));
    if (stride!=rows || buffSz!=rows*cols)
      reallocBuff(rows,cols,9); // re-alloc. buffer, retain matrix elems.
    // Set margin entries to def. fill value
    for (i=0; i<n; i++)
      ArrayUtils<T>::fill(buff+(i*stride+m),this->defFill,rows-m);
    ArrayUtils<T>::fill(buff+(n*stride),this->defFill,rows*(cols-n));
    m=rows; n=cols;
    incrTS();
  }

  template<class T> void BaseMatrix<T>::expand(int rows,int cols)
  {
    int i;

    if (rows<m || cols<n) throw WrongDimensionException(EXCEPT_MSG(""));
    if (rows!=m || cols!=n) {
      if (isMask()) throw MaskObjectException();
      if (buffSz<rows*cols) {
	// Have to re-alloc. buffer
	int newSz=this->getNewBuffSize(rows*cols,buffSz);
	reallocBuff(rows,newSz/rows,10); // realloc.
      } else if (stride<rows) {
	// Keep buffer, but change buffer dimensions.
	// Columns move down, have to copy bottom-up
	moveFrame(buff+stride,m,n-1,stride,buff+rows,rows,false);
	stride=rows;
      } else if (stride*cols>buffSz) {
	// Keep buffer, but change buffer dimensions.
	// Columns move up ('stride' >= 'rows'), have to copy top-down
	moveFrame(buff+stride,m,n-1,stride,buff+rows,rows,true);
	stride=rows;
      }
      // Set margin entries to def. fill value
      for (i=0; i<n; i++)
	ArrayUtils<T>::fill(buff+(i*stride+m),this->defFill,rows-m);
      for (i=n; i<cols; i++)
	ArrayUtils<T>::fill(buff+(i*stride),this->defFill,rows);
      m=rows; n=cols;
      incrTS();
    }
  }

  template<class T> inline void BaseMatrix<T>::shrinkBuffer(int rows,int cols)
  {
    int i;

    if (rows==-1) rows=m;
    if (cols==-1) cols=n;
    if (m>rows || n>cols) throw WrongDimensionException(EXCEPT_MSG(""));
    if (buffSz>rows*cols) {
      // Re-allocate buffer to size 'rows'*'cols'
      if (isMask()) throw MaskObjectException(); // not permitted for masks
      reallocBuff(rows,cols,11);
      incrTS();
    } else if (stride>rows) {
      // Keep buffer, but change dimensions.
      // Columns move up in buffer, so copy top-down
      if (isMask()) throw MaskObjectException(); // not permitted for masks
      moveFrame(buff+stride,m,n-1,stride,buff+rows,rows,true);
      stride=rows;
      incrTS();
    }
  }

  template<class T> inline void
  BaseMatrix<T>::reshape(int rows,int cols,bool noInit,bool strict)
  {
    if (isMask()) throw MaskObjectException(EXCEPT_MSG(""));
    if (rows==-1) {
      if (cols<1) throw InvalidParameterException(EXCEPT_MSG("rows,cols"));
      rows=buffSz/cols;
    } else if (cols==-1) {
      if (rows<1) throw InvalidParameterException(EXCEPT_MSG("rows,cols"));
      cols=buffSz/rows;
    } else if (buffSz<rows*cols)
      throw WrongDimensionException(EXCEPT_MSG("Buffer too small"));
    m=rows; n=cols;
    stride=strict?rows:(buffSz/cols);
    if (!noInit) fill(rows,cols);
    incrTS();
  }

  template<class T> void BaseMatrix<T>::reshapeSave(int rows,int cols,
						    bool strict)
  {
    int i;

    if (isMask()) throw MaskObjectException(EXCEPT_MSG(""));
    if (rows==-1) {
      if (cols<1) throw InvalidParameterException(EXCEPT_MSG("rows,cols"));
      rows=buffSz/cols;
    } else if (cols==-1) {
      if (rows<1) throw InvalidParameterException(EXCEPT_MSG("rows,cols"));
      cols=buffSz/rows;
    } else if (buffSz<rows*cols)
      throw InvalidParameterException(EXCEPT_MSG("Buffer too small"));
    // Keep buffer, but change dimensions
    int cpm=std::min(m,rows),cpn=std::min(n,cols);
    int oldStride=stride;
    stride=strict?rows:(buffSz/cols);
    // 'oldStride'>='stride' (<) ==> Columns move up (down) in buffer
    // ==> top-down (bottom-up)
    moveFrame(buff+oldStride,cpm,cpn-1,oldStride,buff+stride,stride,
	      oldStride>=stride);
    // Set margin entries to def. fill value
    for (i=0; i<cpn; i++)
      ArrayUtils<T>::fill(buff+(i*stride+cpm),this->defFill,rows-cpm);
    for (i=cpn; i<cols; i++)
      ArrayUtils<T>::fill(buff+(i*stride),this->defFill,rows);
    m=rows; n=cols;
    incrTS();
  }

  template<class T> void BaseMatrix<T>::insertRows(const BaseMatrix<T>& mat,
						   int rpos,const Range& rngC)
  {
    int i;

    if (mat.m==0) return; // empty matrix
    if (isMask()) throw MaskObjectException(EXCEPT_MSG(""));
    mat.checkTS(TSARG);
    if (doCheckValid() && !isValidMat(mat))
      throw WrongTypeException(EXCEPT_MSG("'mat' has invalid elements"));
    if (rpos==-1) rpos=m;
    else if (rpos<0) throw InvalidParameterException(EXCEPT_MSG("rpos"));
    // Need closed range 'rngP' from 'rngC'
    const Range* rngP;
    Handle<Range> rngHand;
    if (rngC.isOpen()) {
      rngHand.changeRep(new Range(rngC.getStart(),rngC.getStart()+mat.n-1));
      rngP=rngHand;
    } else {
      if (rngC.size(0)!=mat.n)
	throw InvalidParameterException(EXCEPT_MSG("'rngC' has wrong size"));
      rngP=&rngC;
    }
    int newm=std::max(m,rpos)+mat.m;
    int newn=std::max(n,rngP->getMaxPos(0)+1);
    int oldN=n,oldM=m;
    expand(newm,newn); // new entries filled
    Handle<BaseMatrix<T> > msk;
    if (rpos<oldM) {
      // Move elements down (bottom-up)
      moveFrame(buff+rpos,oldM-rpos,oldN,stride,buff+(rpos+mat.m),stride,
		false);
      // Fill new area
      msk.changeRep(subrvalInt(Range(rpos,rpos+mat.m-1),Range(0,oldN-1)));
      msk->fill();
    }
    // Copy 'mat' into new area
    msk.changeRep(subrvalInt(Range(rpos,rpos+mat.m-1),*rngP));
    *msk=mat;
    incrTS();
  }

  template<class T> inline void BaseMatrix<T>::insertRows(int num,int pos,
							  T val)
  {
    if (num==0) return;
    if (isMask()) throw MaskObjectException(EXCEPT_MSG(""));
    if (num<0) throw InvalidParameterException(EXCEPT_MSG(""));
    if (doCheckValid() && !isValidElement(val))
      throw WrongTypeException(EXCEPT_MSG(""));
    if (pos==-1) pos=m;
    else if (pos<0 || pos>m)
      throw OutOfRangeException(EXCEPT_MSG(""));
    int oldM=m;
    expand(m+num,n); // make space
    if (pos<m) {
      // Move elements down (bottom-up)
      moveFrame(buff+pos,oldM-pos,n,stride,buff+(pos+num),stride,false);
    }
    // Fill area
    if (pos<m || val!=this->defFill) {
      Handle<BaseMatrix<T> > msk(subrvalInt(Range(pos,pos+num-1),
					    Range(0,n-1)));
      msk->fill(num,n,val);
    }
    incrTS();
  }

  template<class T> void BaseMatrix<T>::insertCols(const BaseMatrix<T>& mat,
						   int cpos,const Range& rngR)
  {
    int i;

    if (mat.m==0) return; // empty matrix
    if (isMask()) throw MaskObjectException(EXCEPT_MSG(""));
    mat.checkTS(TSARG);
    if (doCheckValid() && !isValidMat(mat))
      throw WrongTypeException(EXCEPT_MSG("'mat' has invalid elements"));
    if (cpos==-1) cpos=n;
    else if (cpos<0) throw InvalidParameterException(EXCEPT_MSG("cpos"));
    // Need closed range 'rngP' from 'rngR'
    const Range* rngP;
    Handle<Range> rngHand;
    if (rngR.isOpen()) {
      rngHand.changeRep(new Range(rngR.getStart(),rngR.getStart()+mat.m-1));
      rngP=rngHand;
    } else {
      if (rngR.size(0)!=mat.m)
	throw InvalidParameterException(EXCEPT_MSG("'rngR' has wrong size"));
      rngP=&rngR;
    }
    int newn=std::max(n,cpos)+mat.n;
    int newm=std::max(m,rngP->getMaxPos(0)+1);
    int oldN=n,oldM=m;
    expand(newm,newn);
    Handle<BaseMatrix<T> > msk;
    if (cpos<oldN) {
      // Move elements to right (bottom-up)
      moveFrame(buff+(cpos*stride),oldM,oldN-cpos,stride,
		buff+((cpos+mat.n)*stride),stride,false);
      // Fill new area (needs to be done because 'rngR' might not cover all
      // of it!)
      msk.changeRep(subrvalInt(Range(0,oldM-1),Range(cpos,cpos+mat.n-1)));
      msk->fill();
    }
    // Copy 'mat' into new area
    msk.changeRep(subrvalInt(*rngP,Range(cpos,cpos+mat.n-1)));
    *msk=mat;
    incrTS();
  }

  template<class T> inline void BaseMatrix<T>::insertCols(int num,int pos,
							  T val)
  {
    if (num==0) return;
    if (isMask()) throw MaskObjectException(EXCEPT_MSG(""));
    if (num<0) throw InvalidParameterException(EXCEPT_MSG(""));
    if (doCheckValid() && !isValidElement(val))
      throw WrongTypeException(EXCEPT_MSG(""));
    if (pos==-1) pos=n;
    else if (pos<0 || pos>n) throw OutOfRangeException(EXCEPT_MSG(""));
    int oldN=n;
    expand(m,n+num); // make space
    if (pos<m) {
      // Move elements down (bottom-up)
      moveFrame(buff+pos*stride,m,oldN-pos,stride,buff+(pos+num)*stride,
		stride,false);
    }
    // Fill area
    if (pos<n || val!=this->defFill) {
      Handle<BaseMatrix<T> > msk(subrvalInt(Range(0,m-1),
					    Range(pos,pos+num-1)));
      msk->fill(m,num,val);
    }
    incrTS();
  }

  template<class T> inline void BaseMatrix<T>::removeRows(int pos,int num)
  {
    if (num==0) return;
    if (pos<0 || pos>=m) throw OutOfRangeException(EXCEPT_MSG("pos"));
    if (isMask()) throw MaskObjectException(EXCEPT_MSG(""));
    if (num==-1) num=m-pos;
    else if (num<0) throw InvalidParameterException(EXCEPT_MSG(""));
    if (pos+num<m) {
      // Copy (top-down)
      moveFrame(buff+(pos+num),m-pos-num,n,stride,buff+pos,stride);
      m-=num;
    } else
      m=pos;
    incrTS();
  }

  template<class T> inline void BaseMatrix<T>::removeCols(int pos,int num)
  {
    if (num==0) return;
    if (pos<0 || pos>=n) throw OutOfRangeException(EXCEPT_MSG("pos"));
    if (isMask()) throw MaskObjectException(EXCEPT_MSG(""));
    if (num==-1) num=n-pos;
    else if (num<0) throw InvalidParameterException(EXCEPT_MSG(""));
    if (pos+num<n) {
      // Copy (top-down)
      moveFrame(buff+(pos+num)*stride,m,n-pos-num,stride,buff+pos*stride,
		stride);
      n-=num;
    } else
      n=pos;
    incrTS();
  }

  template<class T> inline void
  BaseMatrix<T>::accuMap(const BaseMatrix<T>& a,const ArrayHandle<int>& rmap,
			 const ArrayHandle<int>& cmap)
  {
    int i,k,pos;
    const int* imP;
    Handle<BaseVector<T> > colA,colT;

    checkTS(TSARG); a.checkTS(TSARG);
    if (cmap.size()>0 && rmap.size()>0) {
      if (cmap.size()%2!=1)
	throw InvalidParameterException(EXCEPT_MSG(""));
      k=(cmap.size()+1)/2;
      for (i=pos=0,imP=cmap; i<k; i++) {
	colA.changeRep(a.subrvalInt(RangeFull::get(),*(imP++)));
	colT.changeRep(subrvalInt(RangeFull::get(),pos));
	colT->accuMap(*colA,rmap);
	if (i<k-1) pos+=(*(imP++));
      }
    }
  }

  template<class T> inline bool
  BaseMatrix<T>::checkBounds(const Interval<T>& ival,int* pos,int* stat) const
  {
    checkTS(TSARG);
    Handle<BaseVector<T> > colM;
    Range rng(0,m-1);
    for (int i=0; i<n; i++) {
      colM.changeRep(subrvalInt(rng,i));
      if (!colM->checkBounds(ival,pos,stat)) {
	if (pos!=0) (*pos)+=i*m;
	return false;
      }
    }
    return true;
  }

  template<class T> inline bool
  BaseMatrix<T>::operator==(const BaseMatrix<T>& a) const
  {
    checkTS(TSARG); a.checkTS(TSARG);
    if (m!=a.rows() || n!=a.cols()) return false;
    Handle<BaseVector<T> > colM,colA;
    Range rng(0,m-1);
    for (int i=0; i<n; i++) {
      colM.changeRep(subrvalInt(rng,i));
      colA.changeRep(a.subrvalInt(rng,i));
      if (colM->compare(*colA)!=-1) return false;
    }

    return true;
  }

  template<class T> inline void BaseMatrix<T>::trans(const BaseMatrix<T>& a)
  {
    checkTS(TSARG);
    if (this==&a) {
      trans(); // in-place variant
    } else {
      ensureCapacity(a.n,a.m);
      a.checkTS(TSARG);
      Handle<BaseVector<T> > srcM,trgM;
      Range rng(0,m-1);
      for (int i=0; i<n; i++) {
	trgM.changeRep(subrvalInt(rng,i));
	srcM.changeRep(a.subrvalInt(i,rng));
	trgM->assignInt(*srcM);
      }
    }
  }

  template<class T> void BaseMatrix<T>::trans(bool notInPl)
  {
    int i;

    if (m==0) return; // empty matrix
    checkTS(TSARG);
    if (m==n) {
      // Square matrix: simple special case
      if (!isMaskInd()) {
	T* uppP=buff+stride,*lowP=buff+1;
	for (i=0; i<n-1; i++) {
	  ArrayUtils<T>::swap(uppP,lowP,n-i-1,stride,1);
	  uppP+=(stride+1); lowP+=(stride+1);
	}
      } else {
	Handle<BaseVector<T> > uppM,lowM;
	for (i=0; i<n-1; i++) {
	  Range rng(i+1,n-1);
	  uppM.changeRep(subrvalInt(i,rng));
	  lowM.changeRep(subrvalInt(rng,i));
	  uppM->exchange(*lowM);
	}
      }
    } else {
      if (isMask()) throw MaskObjectException(EXCEPT_MSG("")); // not for masks
      if (notInPl) {
	// Use temp. copy
	BaseMatrix<T> tempMat(*this);
	trans(tempMat);
      } else {
	// "Compress" buffer to just contain this matrix.
	// This does not affect the true buffer size 'buffSz'
	int oldsize_n=buffSz/stride;
	reshapeSave(m,n,true);
	TransInPlace<T>::trans(buff,m,n);
	// Correct size, restore old buffer dim.
	// After 'trans', the content frame is flipped but m,n still as before.
	// First flip m,n, 'reshapeSave' restores buffer dim.
	// (involves moving the content frame), finally correct size.
	int newm=n; n=m; m=stride=newm; // flip (frame still packed)
	reshapeSave(oldsize_n,n,true); // sets stride to oldsize_n
	m=newm; // correct size
      }
    }
  }

  template<class T> void BaseMatrix<T>::permute(const BaseVector<int>& perm,
						bool rowMaj)
  {
    int numDone,targ,start;
    bool emptyCycle;

    if (m==0) return; // empty matrix
    checkTS(TSARG); perm.checkTS(TSARG);
    if (!rowMaj) {
      // Permute columns
      if (perm.size()!=n)
	throw InvalidParameterException(EXCEPT_MSG("perm"));
      if (!perm.checkBounds(Interval<int>(0,n-1,IntVal::ivClosed,
					  IntVal::ivClosed)))
	throw InvalidParameterException(EXCEPT_MSG("perm"));
      Handle<BaseVector<T> > tempVec(newEmptyVector());
      BaseVector<bool> tickoff(n,false);
      Handle<BaseVector<T> > msk;
      start=0; numDone=0;
      for (;;) {
	// Start new cycle with 'start'
	msk.changeRep(subrvalInt(RangeFull::get(),start));
	*tempVec=*msk; numDone++;
	targ=perm[start];
	emptyCycle=true;
	while (targ!=start) {
	  msk.changeRep(subrvalInt(RangeFull::get(),targ));
	  msk->exchange(*tempVec);
	  tickoff[targ]=true; numDone++; emptyCycle=false;
	  targ=perm[targ];
	  // Sanity check (avoids infinite cycles if 'perm' is not a perm.)
	  if (numDone>n)
	    throw InvalidParameterException(EXCEPT_MSG("'perm' not a permutation"));
	}
	if (!emptyCycle) {
	  // Close cycle
	  msk.changeRep(subrvalInt(RangeFull::get(),start));
	  *msk=*tempVec;
	}
	if (numDone==n) break; // leave loop
	// Find next starting point
	for (start++; start<n && tickoff[start]; start++);
	// Sanity check
	if (start==n)
	  throw InvalidParameterException(EXCEPT_MSG("'perm' not a permutation"));
      }
    } else {
      // Permute rows
      if (perm.size()!=m)
	throw InvalidParameterException(EXCEPT_MSG("perm"));
      if (!perm.checkBounds(Interval<int>(0,m-1,IntVal::ivClosed,
					  IntVal::ivClosed)))
	throw InvalidParameterException(EXCEPT_MSG("perm"));
      Handle<BaseVector<T> > tempVec(newEmptyVector());
      BaseVector<bool> tickoff(m,false);
      Handle<BaseVector<T> > msk;
      start=0; numDone=0;
      for (;;) {
	// Start new cycle with 'start'
	msk.changeRep(subrvalInt(start,RangeFull::get()));
	*tempVec=*msk; numDone++;
	targ=perm[start];
	emptyCycle=true;
	while (targ!=start) {
	  msk.changeRep(subrvalInt(targ,RangeFull::get()));
	  msk->exchange(*tempVec);
	  tickoff[targ]=true; numDone++; emptyCycle=false;
	  targ=perm[targ];
	  // Sanity check (avoids infinite cycles if 'perm' is not a perm.)
	  if (numDone>m)
	    throw InvalidParameterException(EXCEPT_MSG("'perm' not a permutation"));
	}
	if (!emptyCycle) {
	  // Close cycle
	  msk.changeRep(subrvalInt(start,RangeFull::get()));
	  *msk=*tempVec;
	}
	if (numDone==m) break; // leave loop
	// Find next starting point
	for (start++; start<m && tickoff[start]; start++);
	// Sanity check
	if (start==m)
	  throw InvalidParameterException(EXCEPT_MSG("'perm' not a permutation"));
      }
    }
  }

  template<class T> inline void BaseMatrix<T>::makeSymm(bool low)
  {
    int i;
    Handle<BaseVector<T> > uppM,lowM;

    if (m!=n)
      throw WrongDimensionException(EXCEPT_MSG("Matrix must be square"));
    checkTS(TSARG);
    if (!low) {
      if (!isMaskInd()) {
	T* uppP=buff+stride,*lowP=buff+1;
	for (i=0; i<n-1; i++) {
	  ArrayUtils<T>::copy(lowP,uppP,n-i-1,1,stride);
	  uppP+=(stride+1); lowP+=(stride+1);
	}
      } else {
	for (i=0; i<n-1; i++) {
	  Range rng(i+1,n-1);
	  uppM.changeRep(subrvalInt(i,rng));
	  lowM.changeRep(subrvalInt(rng,i));
	  lowM->assignInt(*uppM);
	}
      }
    } else {
      if (!isMaskInd()) {
	T* uppP=buff+stride,*lowP=buff+1;
	for (i=0; i<n-1; i++) {
	  ArrayUtils<T>::copy(uppP,lowP,n-i-1,stride,1);
	  uppP+=(stride+1); lowP+=(stride+1);
	}
      } else {
	for (i=0; i<n-1; i++) {
	  Range rng(i+1,n-1);
	  uppM.changeRep(subrvalInt(i,rng));
	  lowM.changeRep(subrvalInt(rng,i));
	  uppM->assignInt(*lowM);
	}
      }
    }
  }

  template<class T> inline void
  BaseMatrix<T>::select(BaseVector<T>& vec,const BaseVector<int>& ind,
			bool rowMaj,const Range& rng) const
  {
    int i,j,sz=rowMaj?n:m;

    checkTS(TSARG); vec.checkTS(TSARG); ind.checkTS(TSARG);
    i=rng.getMaxPos(sz);
    if (i>=ind.size() || i>=sz)
      throw OutOfRangeException(EXCEPT_MSG(""));
    sz=rng.size(sz);
    vec.fill(sz);
    if (!rowMaj) {
      for (i=0; i<sz; i++) {
	j=rng[i]; vec[i]=get(j,ind[j]);
      }
    } else {
      for (i=0; i<sz; i++) {
	j=rng[i]; vec[i]=get(ind[j],j);
      }
    }
  }

  template<class T> inline void
  BaseMatrix<T>::fill(int rows,int cols,T a)
  {
    int i;

    checkTS(TSARG);
    ensureCapacity(rows,cols);
    if (!isMaskInd()) {
      if (m==stride)
	ArrayUtils<T>::fill(buff,a,m*n);
      else {
	T* bP=buff;
	for (i=0; i<n; i++,bP+=stride)
	  ArrayUtils<T>::fill(bP,a,m);
      }
    } else {
      Handle<BaseVector<T> > msk;
      Range rng(0,m-1);
      for (i=0; i<n; i++) {
	msk.changeRep(subrvalInt(rng,i));
	msk->fill(m,a);
      }
    }
  }

  template<class T> inline void BaseMatrix<T>::print(ostream& os) const
  {
    int i,j;

    for (i=0; i<m; i++) {
      for (j=0; j<n-1; j++)
	os << get(i,j) << " ";
      os << get(i,n-1);
      if (i<m-1) os << endl;
    }
  }

  template<class T> inline ofstream&
  BaseMatrix<T>::save(ofstream& os,bool noTag) const
  {
    int i;
    checkTS(TSARG);
    if (!noTag)
      FileUtils::saveHeader(os,"@BaseMatrix",1);
    else {
      i=1; NumberFormats<int>::save(os,&i,1,1,0,4);
    }
    i=getTypeCode(this->defFill); // Type code for T
    NumberFormats<int>::save(os,&i,1,1,0,4);
    return saveInt(os);
  }

  template<class T> inline ifstream&
  BaseMatrix<T>::load(ifstream& is,bool noTag)
  {
    int i,ffv;
    checkTS(TSARG);
    if (!noTag)
      ffv=FileUtils::loadHeader(is,"@BaseMatrix");
    else
      NumberFormats<int>::load(is,&ffv,1,1,0,4);
    if (ffv!=0 && ffv!=1)
      throw FileFormatException("BaseMatrix::load: Unknown FF version number");
    // Check type code. If it is 'typeOther', a warning is displayed inst.
    // of throwing an exception (so that old files can be read)
    NumberFormats<int>::load(is,&i,1,1,0,4);
    if (i!=getTypeCode(this->defFill)) {
      ArrayHandle<char> msg(101);
      if (i==TypeCodeConsts::typeOther) {
	sprintf(msg.p(),"BaseMatrix::load: File has wrong type code %d (required: %d). Continuing...",i,getTypeCode(this->defFill));
	printMsgStdout(msg.p());
      } else {
	sprintf(msg.p(),"BaseMatrix::load: File has wrong type code %d (required: %d)",i,getTypeCode(this->defFill));
	throw WrongTypeException(msg.p());
      }
    }

    return loadInt(is,ffv==1);
  }

  // Internal methods

  template<class T> inline void
  BaseMatrix<T>::convertWrapped(BaseMatWrapper<T>& arg) {
    // Take over fields from 'arg'
    BaseMatrix<T>* mat=arg.getRep(); // dyn. cast in subclasses!
    m=mat->m; n=mat->n;
    copyMembers(mat); // copy 'MatTimeStamp' stuff, deassoc. 'mat' from buffer
    buff=mat->buff;
    maskStat=mat->maskStat;
    stride=mat->stride;
    if (!mat->isMask()) buffSz=mat->buffSz;
    else baseObj=mat->baseObj;
    if (mat->isMaskInd()) {
      colInd=mat->colInd;
      rowInd=mat->rowInd;
    }
    this->copyDefValues(*mat); // def. values
    // Dispose of object wrapped by 'arg' using 'dealloc'. Safe because
    // 'copyMembers' above set 'mat->buffWatch' to 0
    mat->dealloc();
    delete mat; // safe now: matrix is empty
    arg.reset(); // remove repres. in 'arg'
  }

  template<class T> void BaseMatrix<T>::assignInt(const BaseMatrix<T>& mat,
						  uchar srcstrct,
						  uchar trgstrct)
  {
    int i;
    Handle<BaseVector<T> > srcM,trgM;

    if (this==&mat && srcstrct==trgstrct) return;
    checkTS(TSARG); mat.checkTS(TSARG);
    if (srcstrct>MatStrct::last ||
	trgstrct>MatStrct::last)
      throw InvalidParameterException(EXCEPT_MSG(""));
    if (srcstrct!=MatStrct::normal && mat.m!=mat.n)
      throw WrongDimensionException(EXCEPT_MSG("'mat' must be square"));
    if (trgstrct!=srcstrct) {
      if (srcstrct==MatStrct::normal ||
	  (srcstrct==MatStrct::upper &&
	   trgstrct!=MatStrct::lower) ||
	  (srcstrct==MatStrct::lower &&
	   trgstrct!=MatStrct::upper) ||
	  (srcstrct==MatStrct::uppNDg &&
	   trgstrct!=MatStrct::lowNDg) ||
	  (srcstrct==MatStrct::lowNDg &&
	   trgstrct!=MatStrct::uppNDg))
	throw InvalidParameterException(EXCEPT_MSG("str. patt. incompat."));
    }
    if (doCheckValid() && !isValidMat(mat,false,srcstrct))
      throw WrongTypeException(EXCEPT_MSG(""));
    if (!isMask() || m<mat.m || n<mat.n) ensureSize(mat.m,mat.n);
    if (srcstrct==MatStrct::normal) {
      // Full rectangles
      if (!isMaskInd() && !mat.isMaskInd()) {
	// Both matrices are flat
	T* trgB=buff;
	const T* srcB=mat.buff;
	if (stride==mat.m && mat.stride==mat.m) {
	  // Both are contiguous (note that m==mat.m here, because
	  // stride >= m >= mat.m == stride
	  ArrayUtils<T>::copy(trgB,srcB,mat.m*mat.n);
	} else {
	  // Column by column
	  for (i=0; i<mat.n; i++,srcB+=mat.stride,trgB+=stride)
	    ArrayUtils<T>::copy(trgB,srcB,mat.m);
	}
      } else {
	// Use mask vectors for cols
	Range rng(0,mat.m-1);
	for (i=0; i<mat.n; i++) {
	  srcM.changeRep(mat.subrvalInt(rng,i));
	  trgM.changeRep(subrvalInt(rng,i));
	  (*trgM)=(*srcM);
	}
      }
    } else {
      // Structure pattern. Matrices both square
      bool srcUp=(srcstrct==MatStrct::upper ||
		  srcstrct==MatStrct::uppNDg),
	trgUp=(trgstrct==MatStrct::upper ||
	       trgstrct==MatStrct::uppNDg);
      int off=(srcstrct==MatStrct::upper ||
	       srcstrct==MatStrct::lower)?0:1;
      for (i=0; i<mat.n-off; i++) {
	Range rng(i+off,mat.n-1);
	// Row for upper, column for lower
	srcM.changeRep(srcUp?mat.subrvalInt(i,rng):mat.subrvalInt(rng,i));
	trgM.changeRep(trgUp?subrvalInt(i,rng):subrvalInt(rng,i));
	(*trgM)=(*srcM);
      }
    }
  }

  template<class T> void BaseMatrix<T>::assignInt(const BaseVector<T>& vec)
  {
    int i,sz=vec.size();
    Handle<BaseVector<T> > trgM;

    checkTS(TSARG); vec.checkTS(TSARG);
    if (doCheckValid() && !isValidVec(vec))
      throw WrongTypeException(EXCEPT_MSG(""));
    if (!isMask() || n>1 || m<sz) ensureSize(sz,1);
    trgM.changeRep(subrvalInt(Range(0,sz-1),0));
    (*trgM)=vec;
  }

  template<class T> inline void
  BaseMatrix<T>::assignInt(T elem,unsigned char trgstrct)
  {
    int i;
    Handle<BaseVector<T> > trgM;

    if (!isValidElement(elem))
      throw WrongTypeException(EXCEPT_MSG(""));
    checkTS(TSARG);
    if (trgstrct>MatStrct::last)
      throw InvalidParameterException(EXCEPT_MSG(""));
    if (trgstrct!=MatStrct::normal && m!=n)
      throw WrongDimensionException(EXCEPT_MSG("Matrix must be square"));
    if (m==0) ensureCapacity(1,1); // min. size
    if (trgstrct==MatStrct::normal) {
      fill(m,n,elem);
    } else if (trgstrct==MatStrct::upper ||
	       trgstrct==MatStrct::uppNDg) {
      int off=(trgstrct==MatStrct::upper)?0:1;
      for (i=off; i<n; i++) {
	trgM.changeRep(subrvalInt(Range(0,i-off),i));
	(*trgM)=elem;
      }
    } else if (trgstrct==MatStrct::lower ||
	       trgstrct==MatStrct::lowNDg) {
      int off=(trgstrct==MatStrct::lower)?0:1;
      for (i=0; i<n-off; i++) {
	trgM.changeRep(subrvalInt(Range(i+off,m-1),i));
	(*trgM)=elem;
      }
    }
  }

  template<class T> inline bool
  BaseMatrix<T>::isValidMat(const T* pbuff,int pm,int pn,int pstride) const {
    int i,j;

    const T* pB;
    for (i=0; i<pn; i++) {
      for (j=0,pB=pbuff+(i*pstride); j<pm; j++)
	if (!isValidElement(*(pB++))) return false;
    }

    return true;
  }

  template<class T> inline bool
  BaseMatrix<T>::isValidMat(const T* pbuff,int strid,const Range& rngR,
			    const Range& rngC) const {
    int i,j;

    if (rngR.isOpen() || rngC.isOpen())
      throw InvalidParameterException("rngR, rngC must not be open");
    const T* pB;
    for (i=0; i<rngC.size(0); i++) {
      for (j=0,pB=pbuff+(rngC[i]*strid); j<rngR.size(0); j++)
	if (!isValidElement(*(pB+rngR[j]))) return false;
    }

    return true;
  }

  template<class T> bool BaseMatrix<T>::isValidMat(const BaseMatrix<T>& mat,
						   bool doCheck,
						   unsigned char strct) const
  {
    int i,j;

    if (strct>MatStrct::last)
      throw InvalidParameterException(EXCEPT_MSG("strct"));
    if (strct!=MatStrct::normal && mat.m!=mat.n)
      throw InvalidParameterException(EXCEPT_MSG("Matrix must be square"));
    if ((!doCheck && isValidMatType(&mat)) || mat.m==0)
      return true;
    else if (strct==MatStrct::normal) {
      if (!mat.isMaskInd()) {
	const T* pB;
	for (i=0; i<mat.n; i++) {
	  for (j=0,pB=mat.buff+i*mat.stride; j<mat.m; j++)
	    if (!isValidElement(*(pB++))) return false;
	}
      } else {
	for (i=0; i<mat.n; i++)
	  for (j=0; j<mat.m; j++)
	    if (!isValidElement(mat.get(j,i))) return false;
      }
    } else if (strct==MatStrct::upper ||
	       strct==MatStrct::uppNDg) {
      // Upper triangle
      int off=(strct==MatStrct::upper)?0:1;
      if (!mat.isMaskInd()) {
	const T* pB;
	for (i=off; i<mat.n; i++) {
	  for (j=0,pB=mat.buff+i*mat.stride; j<=i-off; j++)
	    if (!isValidElement(*(pB++))) return false;
	}
      } else {
	for (i=off; i<mat.n; i++)
	  for (j=0; j<=i-off; j++)
	    if (!isValidElement(mat.get(j,i))) return false;
      }
    } else {
      // Lower triangle
      int off=(strct==MatStrct::lower)?0:1;
      if (!mat.isMaskInd()) {
	const T* pB;
	for (i=0; i<mat.n-off; i++) {
	  for (j=i+off,pB=mat.buff+i*(mat.stride+1)+off; j<mat.n; j++)
	    if (!isValidElement(*(pB++))) return false;
	}
      } else {
	for (i=0; i<mat.n-off; i++)
	  for (j=i+off; j<mat.n; j++)
	    if (!isValidElement(mat.get(j,i))) return false;
      }
    }

    return true;
  }

  template<class T> bool BaseMatrix<T>::isValidVec(const BaseVector<T>& vec,
						   bool doCheck) const
  {
    int i,sz=vec.size();

    if (!doCheck && isValidVecType(&vec))
      return true;
    for (i=0; i<sz; i++)
      if (!isValidElement(vec[i])) return false;

    return true;
  }

  template<class T> inline void
  BaseMatrix<T>::reassignInt(const T* pbuff,int pm,int pn,int strid,
			     MemWatchBase* watch)
  {
    dealloc(); // transform to empty matrix
    maskStat=statMaskFlat;
    buff=(T*) pbuff;
    m=pm; n=pn; stride=strid;
    baseObj=0; // ref. to buffer only
    rowInd.changeRep(0); colInd.changeRep(0);
    assocBuff(watch); // assoc. with buffer
  }

  template<class T> void
  BaseMatrix<T>::reassignInt(const T* pbuff,int strid,const Range& rngR,
			     const Range& rngC,MemWatchBase* watch)
  {
    int i,j;

    // Ranges are not open, not empty, 1-1 and checked
    dealloc(); // empty matrix
    buff=(T*) pbuff; // may be corrected below
    m=rngR.size(0); n=rngC.size(0);
    baseObj=0; // no assoc. object
    maskStat=statMaskInd; // may be changed below to flat mask
    // Row index
    if (rngR.isFlatRange()) {
      rowInd.changeRep(0); // row index is identity
      // Correct 'buff'
      buff+=rngR.getStart();
    } else if (rngR.isLiteralRange()) {
      // Have to create row index (step size != 1)
      rowInd.changeRep(m);
      for (i=0,j=rngR.getStart(); i<m; i++,j+=rngR.getStep()) rowInd[i]=j;
    } else {
      // Indexed row pos. range. Draw copy of index of 'rngR'
      rowInd.copy(rngR.getIndex());
    }
    // Column index
    if (rngC.isLiteralRange()) {
      if (rngC.getStep()>0) {
	// 'rngC' literal with pos. step size -> flat mask matrix
	if (rowInd==0) maskStat=statMaskFlat;
	stride=strid*rngC.getStep();
	colInd.changeRep(0); // no column pos. index
	// Correct 'buff'
	buff+=(rngC.getStart()*strid);
      } else {
	// Have to create column index (indexed mask)
	colInd.changeRep(n);
	for (i=0,j=rngC.getStart(); i<n; i++,j+=rngC.getStep())
	  colInd[i]=j;
	stride=strid;
      }
    } else {
      // Indexed col. pos. range
      stride=strid;
      colInd.copy(rngC.getIndex());
    }
    assocBuff(watch); // assoc. with buffer
  }

  template<class T> void
  BaseMatrix<T>::reassignInt(const BaseMatrix<T>& mat,const Range& rngR,
			     const Range& rngC)
  {
    dealloc(); // empty matrix
    if (mat.m>0) {
#ifdef DEBUG_TRACKHANDLES
      debugCause=2;
#endif
      if (!mat.isMaskInd()) {
	// 'mat' normal or flat mask: Can use flat buffer version of
	// 'reassignInt'. Need to close ranges if they are open. 'baseObj' and
	// timestamp are dealt with below.
	const Range* rpR=&rngR,*rpC=&rngC;
	Handle<Range> tempR,tempC;
	if (rngR.isOpen()) {
	  tempR.changeRep(new Range(rngR.getStart(),mat.m-1));
	  rpR=tempR;
	}
	if (rngC.isOpen()) {
	  tempC.changeRep(new Range(rngC.getStart(),mat.n-1));
	  rpC=tempC;
	}
	reassignInt(mat.buff,mat.stride,*rpR,*rpC,mat.buffWatch);
      } else {
	// 'mat' indexed mask -> so will be this object
	maskStat=statMaskInd;
	buff=mat.buff; // may be changed below
	m=rngR.size(mat.m); n=rngC.size(mat.n);
	stride=mat.stride; // may be changed below
	// Row pos. index
	if (mat.rowInd==0) {
	  // 'mat' has no row pos. index
	  if (rngR.isFlatRange()) {
	    rowInd.changeRep(0); // this row index is identity
	    buff+=rngR.getStart(); // modify 'buff'
	  } else if (rngR.isLiteralRange()) {
	    // Have to create row index (step size != 1)
	    rowInd.changeRep(m);
	    for (int i=0,j=rngR.getStart(); i<m; i++,j+=rngR.getStep())
	      rowInd[i]=j;
	  } else {
	    // Indexed row pos. range: draw copy
	    rowInd.copy(rngR.getIndex());
	  }
	} else {
	  // 'mat' has row pos. index. Can re-use this one iff 'rngR' is the
	  // full range
	  if (rngR.isFullRange(mat.m))
	    rowInd=mat.rowInd;
	  else {
	    rowInd.changeRep(m);
	    rngR.mapIndex(mat.rowInd,mat.m,rowInd,false); // map index
	  }
	}
	// Column pos. index. Can re-use 'mat' index iff 'rngC' is the full
	// range
	if (rngC.isFullRange(mat.n))
	  colInd=mat.colInd;
	else if (mat.colInd==0) {
	  // 'mat' has no column index
	  if (rngC.isLiteralRange() && rngC.getStep()>0) {
	    // Do not need column index
	    stride=mat.stride*rngC.getStep();
	    buff+=(mat.stride*rngC.getStart());
	  } else {
	    colInd.changeRep(n);
	    for (int i=0; i<n; i++) colInd[i]=rngC[i];
	  }
	} else {
	  // 'mat' has column index
	  colInd.changeRep(n);
	  rngC.mapIndex(mat.colInd,mat.n,colInd,false); // map index
	}
	assocBuff(mat.buffWatch);
      }
      // Base object pointer and timestamp: copy from 'mat'
      if (mat.isMask())
	// NOTE: It is important that we copy 'mat's timestamp value, not the
	// one of the object underlying 'mat', because if 'mat' is invalid, so
	// should be this mask matrix.
	baseObj=mat.baseObj; // same base object
      else
	baseObj=&mat; // hook onto 'mat' (which is normal)
      setTS(mat.getTS()); // 'mat's current timestamp value
    }
  }

  template<class T> void
  BaseMatrix<T>::reassignInt(const BaseVector<T>& vec,bool col)
  {
    int i;

    // 'baseObj' and timestamp dealt with below
    dealloc(); // empty matrix
    if (vec.n>0) {
      if (col) {
	// Column vector
	if (!vec.isMask() || vec.isMaskFlat()) {
	  // Flat vector -> flat mask matrix
	  reassignInt(vec.buff,vec.n,1,vec.n,vec.buffWatch);
	} else {
	  // Indexed mask matrix
	  // Although 'stride' is not really used, it must be set s.t. it
	  // is larger than any of 'rowInd' elem.
	  buff=vec.buff; // may be modified below
	  m=vec.n; n=1; maskStat=statMaskInd; // row index non-trivial
	  if (vec.isMaskInd()) {
	    rowInd=vec.mindex; // just use 'mindex'
	    stride=(*(std::max_element(rowInd.p(),rowInd.p()+rowInd.size())))+1;
	  } else {
	    // Linear (non-flat) mask vector. 'rowInd' has to be created.
	    // ATTENTION: 'rowInd' must have nonneg. entries.
	    rowInd.changeRep(m);
	    if (vec.step>0) {
	      for (i=0; i<m; i++) rowInd[i]=i*vec.step;
	      stride=(m-1)*vec.step+1;
	    } else {
	      buff+=((m-1)*vec.step); // correct 'buff'
	      for (i=0; i<m; i++) rowInd[i]=(m-i-1)*(-vec.step);
	      stride=(m-1)*(-vec.step)+1;
	    }
	  }
	  assocBuff(vec.buffWatch);
	}
      } else {
	// Row vector ('rowInd' not needed)
	buff=vec.buff; // may be modif. below
	m=1; n=vec.n;
	stride=1; maskStat=statMaskInd; // may be modif. below
	if (!vec.isMaskInd()) {
	  if (vec.step>0) {
	    // Not an indexed mask vector -> flat mask matrix
	    stride=vec.step; maskStat=statMaskFlat;
	  } else {
	    // Step size negative: have to create col. ind.
	    colInd.changeRep(n);
	    buff+=((n-1)*vec.step); // correct 'buff'
	    for (i=0; i<n; i++) colInd[i]=(n-i-1)*(-vec.step);
	  }
	} else {
	  // Indexed mask vector
	  colInd=vec.mindex; // re-use 'mindex'
	}
	assocBuff(vec.buffWatch);
      }
      // Base object pointer and timestamp: copy from 'vec'
      if (vec.isMask())
	// NOTE: It is important that we copy 'vec's timestamp value, not the
	// one of the object underlying 'vec', because if 'vec' is invalid, so
	// should be this mask matrix.
	baseObj=vec.baseObj; // same base object
      else
	baseObj=&vec; // hook onto 'vec' (which is normal)
      setTS(vec.getTS()); // 'vec's current timestamp value
    }
  }

  template<class T> inline BaseMatrix<T>*
  BaseMatrix<T>::subrvalInt(const Range& rngR,const Range& rngC) const {
    BaseMatrix<T>* mat=newEmpty();
    mat->reassign(*this,rngR,rngC);

    return mat;
  }

  template<class T> inline BaseVector<T>*
  BaseMatrix<T>::subrvalInt(const Range& rngR,int rngC) const {
    if (rngC<0 || rngC>=n) throw OutOfRangeException(EXCEPT_MSG("rngC"));
    if (rngR.checkRange(m)) throw OutOfRangeException(EXCEPT_MSG("rngR"));
    if (!rngR.isUniqueMap())
      throw InvalidParameterException(EXCEPT_MSG("Ranges must be 1-1"));
    BaseVector<T>* vec=newEmptyVector();
#ifdef DEBUG_TRACKHANDLES
    vec->debugCause=3; // DEBUG
#endif
    T* buffp=buff+stride*COLPOS(rngC);
    if (rowInd==0) {
      vec->reassignInt(buffp,m,rngR,buffWatch);
    } else {
      ArrayHandle<int> tempind(rngR.size(m));
      rngR.mapIndex(rowInd,m,tempind,false);
      vec->reassignInt(buffp,stride,Range(tempind),buffWatch,true);
    }

    return vec;
  }

  template<class T> inline BaseVector<T>*
  BaseMatrix<T>::subrvalInt(int rngR,const Range& rngC) const {
    if (rngR<0 || rngR>=m) throw OutOfRangeException(EXCEPT_MSG("rngR"));
    if (rngC.checkRange(n)) throw OutOfRangeException(EXCEPT_MSG("rngC"));
    if (!rngC.isUniqueMap())
      throw InvalidParameterException(EXCEPT_MSG("Ranges must be 1-1"));
    BaseVector<T>* vec=newEmptyVector();
#ifdef DEBUG_TRACKHANDLES
    vec->debugCause=4; // DEBUG
#endif
    T* buffp=buff+ROWPOS(rngR);
    if (colInd==0) {
      if (rngC.isLiteralRange()) {
	// Create linear mask vector
	vec->reassignInt(buffp,(n-1)*stride+1,
			 Range(rngC.getStart()*stride,rngC.getEnd(n)*stride,
			       rngC.getStep()*stride),buffWatch);
      } else {
	// Convert index
	ArrayHandle<int> tempind(rngC.size(n));
	for (int i=0; i<tempind.size(); i++) tempind[i]=rngC[i]*stride;
	vec->reassignInt(buffp,(n-1)*stride+1,Range(tempind),buffWatch,true);
      }
    } else {
      ArrayHandle<int> tempind(rngC.size(n));
      int mpos=colInd[rngC[0]]*stride,temp;
      tempind[0]=mpos;
      for (int i=1; i<tempind.size(); i++) {
	if ((temp=tempind[i]=colInd[rngC[i]]*stride)>mpos) mpos=temp;
      }
      vec->reassignInt(buffp,mpos+1,Range(tempind),buffWatch,true);
    }

    return vec;
  }

  template<class T> inline BaseMatrix<T>*
  BaseMatrix<T>::sublvalInt(const Range& rngR,const Range& rngC) {
    int newm,newn;

    if (rngR.isOpen())
      newm=std::max(rngR.getStart()+1,m); // open range
    else
      newm=rngR.getMaxPos(0)+1;
    if (rngC.isOpen())
      newn=std::max(rngC.getStart()+1,n); // open range
    else
      newn=rngC.getMaxPos(0)+1;
    if (newm>m || newn>n)
      expand(newm,newn); // expand to new size
    BaseMatrix<T>* mat=newEmpty();
    mat->reassignInt(*this,rngR,rngC);

    return mat;
  }

  template<class T> inline BaseVector<T>*
  BaseMatrix<T>::sublvalInt(const Range& rngR,int rngC) {
    int newm,newn;

    if (rngC<0) throw OutOfRangeException(EXCEPT_MSG("rngC"));
    if (!rngR.isUniqueMap())
      throw InvalidParameterException(EXCEPT_MSG("Ranges must be 1-1"));
    if (rngR.isOpen())
      newm=std::max(rngR.getStart()+1,m); // open range
    else
      newm=rngR.getMaxPos(0)+1;
    newn=std::max(rngC+1,n);
    if (newm>m || newn>n)
      expand(newm,newn); // expand to new size
    BaseVector<T>* vec=newEmptyVector();
#ifdef DEBUG_TRACKHANDLES
    vec->debugCause=3; // DEBUG
#endif
    T* buffp=buff+stride*COLPOS(rngC);
    if (rowInd==0) {
      vec->reassignInt(buffp,m,rngR,buffWatch);
    } else {
      ArrayHandle<int> tempind(rngR.size(m));
      rngR.mapIndex(rowInd,m,tempind,false);
      vec->reassignInt(buffp,stride,Range(tempind),buffWatch,true);
    }

    return vec;
  }

  template<class T> inline BaseVector<T>*
  BaseMatrix<T>::sublvalInt(int rngR,const Range& rngC) {
    int newm,newn;

    if (rngR<0) throw OutOfRangeException(EXCEPT_MSG("rngR"));
    if (!rngC.isUniqueMap())
      throw InvalidParameterException(EXCEPT_MSG("Ranges must be 1-1"));
    if (rngC.isOpen())
      newn=std::max(rngC.getStart()+1,n); // open range
    else
      newn=rngC.getMaxPos(0)+1;
    newm=std::max(rngR+1,m);
    if (newm>m || newn>n)
      expand(newm,newn); // expand to new size
    BaseVector<T>* vec=newEmptyVector();
#ifdef DEBUG_TRACKHANDLES
    vec->debugCause=4; // DEBUG
#endif
    T* buffp=buff+ROWPOS(rngR);
    if (colInd==0) {
      if (rngC.isLiteralRange()) {
	// Create linear mask vector
	vec->reassignInt(buffp,(n-1)*stride+1,
			 Range(rngC.getStart()*stride,rngC.getEnd(n)*stride,
			       rngC.getStep()*stride),buffWatch);
      } else {
	// Convert index
	ArrayHandle<int> tempind(rngC.size(n));
	for (int i=0; i<tempind.size(); i++) tempind[i]=rngC[i]*stride;
	vec->reassignInt(buffp,(n-1)*stride+1,Range(tempind),buffWatch,true);
      }
    } else {
      ArrayHandle<int> tempind(rngC.size(n));
      int mpos=colInd[rngC[0]]*stride,temp;
      tempind[0]=mpos;
      for (int i=1; i<tempind.size(); i++) {
	if ((temp=tempind[i]=colInd[rngC[i]]*stride)>mpos) mpos=temp;
      }
      vec->reassignInt(buffp,mpos+1,Range(tempind),buffWatch,true);
    }

    return vec;
  }

  /*
   * For a flat matrix, the elements are
   *   rng[i]*(stride+1)+start,
   * where start=off*stride for off>=0, start=-off for off<0.
   * This is linear if 'rng' is literal.
   * For an indexed mask, the elements are
   *   colInd[i+off]*stride+rowInd[i], off>=0
   *   colInd[i]*stride+rowInd[i-off], off<0
   */
  template<class T> inline BaseVector<T>*
  BaseMatrix<T>::diagInt(int off,const Range& rng) const
  {
    int len;

    if (off>=0)
      len=std::min(m,n-off);
    else
      len=std::min(m+off,n);
    if (len<=0)
      throw InvalidParameterException(EXCEPT_MSG("Diagonal empty"));
    if (rng.checkRange(len))
      throw OutOfRangeException(EXCEPT_MSG("rng"));
    BaseVector<T>* vec=newEmptyVector();
#ifdef DEBUG_TRACKHANDLES
    vec->debugCause=5; // DEBUG
#endif
    if (!isMaskInd()) {
      int start=(off>=0)?(off*stride):(-off);
      if (rng.isLiteralRange()) {
	// Diag. as linear mask vector
	int end=start+rng.getEnd(len)*(stride+1);
	start+=rng.getStart()*(stride+1);
	vec->reassignInt(buff,(n-1)*stride+m,
			 Range(start,end,rng.getStep()*(stride+1)),buffWatch);
      } else {
	ArrayHandle<int> tempind(rng.size(0));
	for (int i=0; i<tempind.size(); i++)
	  tempind[i]=rng[i]*(stride+1)+start;
	vec->reassignInt(buff,(n-1)*stride+m,Range(tempind),buffWatch,true);
      }
    } else {
      ArrayHandle<int> tempind(rng.size(len));
      int mpos,temp;
      if (off>=0) {
	tempind[0]=mpos=COLPOS(off)*stride+ROWPOS(0);
	for (int i=1; i<tempind.size(); i++)
	  if ((temp=tempind[i]=COLPOS(i+off)*stride+ROWPOS(i))>mpos)
	    mpos=temp;
      } else {
	tempind[0]=mpos=COLPOS(0)*stride+ROWPOS(-off);
	for (int i=1; i<tempind.size(); i++)
	  if ((temp=tempind[i]=COLPOS(i)*stride+ROWPOS(i-off))>mpos) mpos=temp;
      }
      vec->reassignInt(buff,mpos+1,Range(tempind),buffWatch,true);
    }

    return vec;
  }

  template<class T> inline ofstream&
  BaseMatrix<T>::saveInt(ofstream& os,bool storeBSize) const
  {
    int i;

    checkTS(TSARG);
    NumberFormats<int>::save(os,&m,1,1,0,4);
    NumberFormats<int>::save(os,&n,1,1,0,4);
    if (storeBSize) {
      i=sizeof(T);
      NumberFormats<int>::save(os,&i,1,1,0,4);
    }
    if (!isMaskInd()) {
      if (m==stride) {
	// Flat and contiguous
	NumberFormats<T>::save(os,buff,n*m);
      } else {
	for (i=0; i<n; i++)
	  NumberFormats<T>::save(os,buff+(i*stride),m);
      }
    } else {
      for (i=0; i<n; i++)
	NumberFormats<T>::save(os,buff+(COLPOS(i)*stride),m,1,rowInd);
    }

    return os;
  }

  template<class T> inline ifstream&
  BaseMatrix<T>::loadInt(ifstream& is,bool loadBSize)
  {
    int newm,newn,i,fsize;

    checkTS(TSARG);
    NumberFormats<int>::load(is,&newm,1,1,0,4);
    NumberFormats<int>::load(is,&newn,1,1,0,4);
    if (!loadBSize) fsize=sizeof(T);
    else {
      NumberFormats<int>::load(is,&fsize,1,1,0,4);
      if (fsize<1) throw FileFormatException(EXCEPT_MSG("File format error"));
    }
    if (newm<0 || newn <0)
      throw FileFormatException(EXCEPT_MSG("File format error"));
    if (newm==0 || newn==0)
      ensureCapacity(0,0);
    else {
      ensureCapacity(newm,newn);
      if (!isMaskInd()) {
	if (m==stride) {
	  // Flat and contiguous
	  NumberFormats<T>::load(is,buff,n*m,1,0,fsize);
	} else {
	  for (i=0; i<n; i++)
	    NumberFormats<T>::load(is,buff+(i*stride),m,1,0,fsize);
	}
      } else {
	for (i=0; i<n; i++)
	  NumberFormats<T>::load(is,buff+(COLPOS(i)*stride),m,1,rowInd,fsize);
      }
    }
    // Force element checks
    if (doCheckValid() && !isValidMat(*this,true))
      throw WrongTypeException(EXCEPT_MSG("Entries in file violate matrix element constraints!"));

    return is;
  }

  template<class T> inline void
  BaseMatrix<T>::init(int rows,int cols,uchar debStat,T a)
  {
    if (!isValidElement(a))
      throw WrongTypeException("'a' invalid element");
    alloc(rows,cols,debStat);
    ArrayUtils<T>::fill(buff,a,rows*cols);
  }

  template<class T> inline void
  BaseMatrix<T>::alloc(int rows,int cols,uchar debStat)
  {
    if (isMask()) throw MaskObjectException(EXCEPT_MSG(""));
    if (rows<0 || cols<0) throw WrongDimensionException(EXCEPT_MSG(""));
    dealloc(); // deallocate resources (if any)

    if (rows>0 && cols>0) {
      m=stride=rows; n=cols;
      buffSz=rows*cols;
      MemWatcher<T>* watch=new MemWatcher<T>(buffSz,0,debStat); // allocation
      assocBuff(watch); // association
      buff=watch->getBuff(); // buffer pointer
    }
  }

  template<class T> inline void
  BaseMatrix<T>::ensureCapacity(int rows,int cols,bool strict)
  {
    if (rows<0 || cols<0) throw WrongDimensionException(EXCEPT_MSG(""));
    if (rows==0 || cols==0)
      rows=cols=0;
    if (rows!=m || cols!=n) {
      if (isMask()) throw MaskObjectException(EXCEPT_MSG(""));
      if (rows!=0) {
	if (rows*cols>buffSz) {
	  alloc(rows,cols,8); // have to re-alloc. flat buffer
	} else if (stride<rows || strict) {
	  // Set 'stride'=='rows'
	  stride=rows;
	} else if (stride*cols>buffSz) {
	  // Need to reduce old 'stride' value
	  stride=buffSz/cols;
	}
      }
      m=rows; n=cols;
      incrTS(); // size changed
    }
  }

  template<class T> inline void BaseMatrix<T>::dealloc(bool noReset)
  {
    deassocBuff(); // remove assoc. with buffer
    if (!noReset) {
      m=n=0; buff=0;
      maskStat=statNormal;
      stride=buffSz=0;
      rowInd.changeRep(0); colInd.changeRep(0);
    }
  }

  template<class T> inline void
  BaseMatrix<T>::moveFrame(const T* srcBuff,int srcM,int srcN,int srcStride,
			   T* trgBuff,int trgStride,bool top) const
  {
    if (top) {
      for (int i=0; i<srcN; i++,srcBuff+=srcStride,trgBuff+=trgStride)
	ArrayUtils<T>::copy(trgBuff,srcBuff,srcM);
    } else {
      srcBuff+=(srcN-1)*srcStride; trgBuff+=(srcN-1)*trgStride;
      for (int i=0; i<srcN; i++,srcBuff-=srcStride,trgBuff-=trgStride)
	ArrayUtils<T>::copy(trgBuff,srcBuff,srcM);
    }
  }

  template<class T> inline void
  BaseMatrix<T>::reallocBuff(int rows,int cols,uchar debStat)
  {
    if (rows<m || cols<n) throw WrongDimensionException(EXCEPT_MSG(""));
    if (buff==0 && buffWatch!=0)
      throw InternalException(EXCEPT_MSG(""));
    T* oldBuff=buff; buff=0;
    int oldN=n; n=0;
    int oldM=m; m=0;
    int oldStride=stride; stride=0;
    int oldBuffSz=buffSz; buffSz=0;
    MemWatchBase* oldWatch=buffWatch; buffWatch=0; // looks empty now
    alloc(rows,cols,debStat); // re-alloc.
    m=oldM; n=oldN; // restore size
    if (oldBuff!=0) {
      // Copy
      moveFrame(oldBuff,m,n,oldStride,buff,rows);
      // De-assoc. with old buffer
      MemWatchBase* newWatch=buffWatch;
      buffWatch=oldWatch; // looks like old to 'dealloc'
      T* newBuff=buff; buff=oldBuff;
      stride=oldStride; buffSz=oldBuffSz;
      dealloc(true);
      buffWatch=newWatch; // restore watcher
      buff=newBuff;
      stride=rows; buffSz=rows*cols;
      m=oldM; n=oldN;
    }
  }

  // Definition of constants

  template<class T> const unsigned char BaseMatrix<T>::statNormal     ;
  template<class T> const unsigned char BaseMatrix<T>::statMaskFlat   ;
  template<class T> const unsigned char BaseMatrix<T>::statMaskInd    ;

  // Global functions

  template<class T> inline ostream&
  operator<<(ostream& os,const BaseMatrix<T>& mat)
  {
    if (!mat.canPrint())
      throw NotImplemException(EXCEPT_MSG("'print' not implemented"));
    mat.print(os);

    return os;
  }
//ENDNS

#undef TSARG
#undef COLPOS
#undef ROWPOS

#endif
