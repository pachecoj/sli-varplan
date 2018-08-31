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
 * Desc.:  Header class LogVector
 * ------------------------------------------------------------------- */

#ifndef LOGVECTOR_H
#define LOGVECTOR_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include <algorithm>

#include "lhotse/matrix/default.h"
#include "lhotse/matrix/Vector.h"
#include "lhotse/matrix/StVector.h"

//BEGINNS(matrix)
  /**
   * Manages a double vector of arbitrary size. Elements can only be appended
   * at the end. In normal mode, the vector keeps on growing, while in window
   * mode, once the vector has grown to a given (window) size, elements at the
   * beginning are removed so to keep the size constant.
   * <p>
   * The log vector can be deactivated at any time. For a deactivated log
   * vector, attempts to append elements are simply ignored, all other methods
   * work as usual. A deactivated log vector can be reactivated at any time.
   * By def., a log vector is activated.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class LogVector : public Vector
  {
  protected:
    // Constants

    static const int maxEmptySize=1000; // See 'reset'

    // Members

    ArrayHandle<double> buff; // Vector buffer
    int window;   // Window size (or -1 if in normal mode)
    bool active;  // Is log vector activated? See header comment

  public:
    // Constructors

    /**
     * Default constructor. Results in empty vector
     *
     * @param winSize Optional. Window size. Def.: Normal mode, no window
     * @param activeP Optional. Is log vector activated? Def.: true
     */
    LogVector(int winSize=-1,bool activeP=true) :
      Vector(),window(winSize),active(activeP) {
      if (winSize!=-1 && winSize<1)
	throw InvalidParameterException(EXCEPT_MSG(""));
      buff.changeRep(50); // init. size
    }

    LogVector(const LogVector& arg) : Vector(),window(-1),active(false) {
      operator=(arg);
    }

    // Public methods

    /**
     * Assignment operator.
     * NOTE: 'src' can also be a 'StVector' object. In this case,
     * this vector will retain its window size and activity status.
     *
     * @param src Source object
     */
    LogVector& operator=(const LogVector& src);

    LogVector& operator=(const StVector& src);

    /**
     * Activates log vector
     */
    void activate() {
      active=true;
    }

    /**
     * Deactivates log vector
     */
    void deactivate() {
      active=false;
    }

    /**
     * @return Is vector active?
     */
    bool isActive() const {
      return active;
    }

    /**
     * Adds element at end. Ignored if log vector deactivated
     *
     * @param elem Element to add
     */
    void add(double elem) {
      if (active) {
	if (window==-1 || n<window)
	  ensureCapacity(n+1);
	else
	  memmove(buff.p(),buff.p()+1,(n-1)*sizeof(double));
	buff[n-1]=elem;
      }
    }

    /**
     * Just for backwards compatibility. Use 'size'!
     *
     * @return Length of vector
     */
    int length() const {
      return n;
    }

    double operator[](int pos) const {
      if (pos<0 || pos>=n) throw OutOfRangeException("pos");
      return buff[pos];
    }

    /**
     * Resets log vector to be empty, keeping everything else (buffer size,
     * window size, activation status).
     * NOTE: If the buffer size is larger than 'maxEmptySize', the buffer is
     * deallocated. Otherwise, it is not changed.
     */
    void reset() {
      n=0;
      if (buff.size()>maxEmptySize)
	buff.changeRep(50);
    }

    /**
     * Copies (part of) the vector into a given double vector.
     *
     * @param vec Vector to return part in
     * @param rng Range for part. Def.: full
     */
    void getVec(StVector& vec,const Range& rng=RangeFull::get()) const {
      if (rng.checkRange(n)) throw OutOfRangeException(EXCEPT_MSG(""));
      if (!rng.isOpen())
	vec=(StVector&) StVector::mask(buff,rng);
      else
	vec=(StVector&) StVector::mask(buff,Range(rng.getStart(),n-1));
    }

    /**
     * Saves object into file stream 'os'. The object can be restored using
     * 'load'.
     * The file format is the current one used by 'StVector'.
     *
     * @param os Output file stream
     */
    ofstream& save(ofstream& os) const;

    /**
     * Loads object from file stream 'is'. The current object is overwritten.
     * The file format is the current one used by 'StVector'.
     * <p>
     * Loading does not affect the current setting of 'window' or 'active'.
     * If the stored vector is larger than 'window', only the right-end
     * elements are loaded.
     *
     * @param is Input file stream
     */
    ifstream& load(ifstream& is);

    // Internal methods

    /**
     * If the buffer is smaller than 'newLength', it is increased to size
     * max. of 'newLength' and twice its old size. If 'keepOld'
     * is true (def.), the old entries are retained, otherwise they are
     * lost. The new entries are not initialised. The vector size is set
     * to 'newLength'.
     *
     * @param newLength New vector size
     * @param keepOld   Optional. Def.: true
     */
    void ensureCapacity(int newLength,bool keepOld=true) {
      if (buff.size()<newLength) {
	int newSize=std::max(newLength,2*buff.size());
	if (keepOld) {
	  ArrayHandle<double> newBuff(newSize);
	  memmove(newBuff.p(),buff.p(),n*sizeof(double));
	  buff=newBuff;
	} else
	  buff.changeRep(newSize);
      }
      n=newLength;
    }
  };
//ENDNS

#endif
