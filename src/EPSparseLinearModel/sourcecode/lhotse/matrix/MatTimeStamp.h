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
 * Desc.:  Header abstract class MatTimeStamp
 * ------------------------------------------------------------------- */

#ifndef MATTIMESTAMP_H
#define MATTIMESTAMP_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/matrix/default.h"

//BEGINNS(matrix)
  /**
   * Both 'BaseMatrix' and 'BaseVector' inherit from this base class which
   * implements a timestamp member 'ts'. A normal object must increase the
   * timestamp each time its size changes. A mask object copies the
   * timestamp value of its base object upon association. A mask can then
   * check whether its base object has changed size by comparing its
   * timestamp value with the one of the base object. If they differ, the
   * mask association is corrupt, an exception is thrown.
   * <p>
   * The following methods have to implemented in subclasses:
   * - isMask: Is the object a mask? In this context, a object is a mask
   *   iff it is not the owner of the buffer that 'buffWatch' is watching.
   *   If 'buffWatch'==0, the object must be a mask (or be empty)
   * - getBaseObj: If the object is a mask which ref. to another
   *   'MatTimeStamp' object, ret. this one here. Otherwise, ret. 0.
   *   NOTE: A mask need not be assoc. with a 'MatTimeStamp' object. It can
   *   also just refer to a memory region (in which case the whole mechanism
   *   here does not apply)
   * <p>
   * How to use the services in subclass methods
   * NOTE: The mechanism here requires complete support from all subclass
   * methods!
   * - checkTS: 
   *   Must be called at the beginning of a method, for this object and all
   *   arguments. For a mask object, 'checkTS' checks for timestamp
   *   consistency, and also checks the 'isAssoc' flag (below).
   * - incrTS:
   *   Must be called within a method once this object changes size.
   *   NOTE: If the size change is done via another method, this one does
   *   call 'incrTS'.
   *   Failing to increase the timestamp for a size change means that
   *   masks to the object may be invalid, without the mechanisms here
   *   detecting this!
   *   NOTE: 'incrTS' is protected here, which means that it cannot always
   *   be called for an argument to a method (if the argument has type T2
   *   different from T, for example: below). In general, methods must be
   *   designed to never DIRECTLY change the size of an argument, but only
   *   through a method of the argument class.
   * <p>
   * Furthermore, the mem. watcher for the buffer of the object is managed
   * in 'buffWatch'.
   * NOTE: 'buffWatch' can be 0, namely if the object uses a buffer which does
   * not have a mem. watcher. This is used only in special situations for
   * downw. compatibility, or to wrap foreign code!
   * If the object is assoc. with a new buffer, 'assocBuff' must be called
   * passing the mem. watcher. If the object is de-assoc. from its buffer,
   * 'deassocBuff' must be called. Note that 'deassocBuff' will set the
   * 'isAssoc' flag in the mem. watcher to true if it is called for a non-mask
   * object. This signals mask objects later on in 'checkTS' to become invalid
   * if they have ref. to this buffer (see 'MemWatcher' and
   * doc/memManage.txt). 'deassocBuff' decr. the ref. counter of the mem.
   * watcher and destroys it if the counter drops to 0.
   * <p>
   * Special timestamp values:
   * - 0: Value of 'ts' after construction, before the first size change.
   *      This value is never used again (at a wrap-around, 'ts' is set to 1).
   *      ==> Since empty vectors/matrices are created with TS 0, this is a
   *          way to check whether an object has always been empty (see
   *          'reassign' methods)
   * <p>
   * NOTE: None of the methods here (with the exc. of 'isMask') are
   * be called by methods other than in subclasses, so why are not all
   * protected?
   * ==> Need to call them for 'BaseVector<T2>' (in template members of
   *     'BaseVector<T>') which cannot be our friend for all T2!
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class MatTimeStamp
  {
  protected:
    // Variables

    MemWatchBase* buffWatch; // See header comment
    uint ts;                 // timestamp value
#ifdef DEBUG_TRACKHANDLES
    int debugCause;
#endif

  public:
    /**
     * Default constructor. Timestamp init. to 0.
     *
     * @param buffW Initial 'buffWatch' value. Def.: 0
     */
    MatTimeStamp(MemWatchBase* buffW=0) : buffWatch(buffW),ts(0)
#ifdef DEBUG_TRACKHANDLES
      ,debugCause(0)
#endif
    {}

    /**
     * If this is a mask object with an assoc. base object (method 'getBaseObj'
     * is called to get the pointer) and our timestamp value is different from
     * the one of the base object, we return false. In any other case, we
     * return true. If 'noExc' is false (the default), then instead of
     * returning false, a 'MaskObjectException' is thrown.
     * If 'file' and 'line' are given, file name and line number can be passed
     * (from where 'checkTS' is called), which are then included in the
     * exception message (see 'BaseVector' for example).
     * <p>
     * If this is a normal object, or a mask without a mem. watcher
     * ('buffWatch'==0), we always return true.
     * <p>
     * NOTE: This method is public, for reasons noted in the header comments.
     *
     * @param file  S.a. Def.: 0
     * @param line  S.a. Must be given if 'file' is
     * @param noExc Optional. Def.: false
     * @return      See above
     */
    virtual bool checkTS(const char* file=0,int line=0,bool noExc=false)
      const {
      const MatTimeStamp* baseObj=getBaseObj(); // 0 if no mask
      if (baseObj!=0) {
	// NOTE: If 'baseObj'!=0:
	// - this object is a mask
	// - there is a mem. watcher ('buffWatch'), since the base object
	//   is normal and must have one
	if (buffWatch==0) // sanity check
	  throw InternalException(EXCEPT_MSG("MatTimeStamp::checkTS: buffWatch==0, baseObj!=0 (base object has no mem watcher)"));
	if (buffWatch->isAssoc) {
	  if (!noExc) {
	    char msg[100+strlen(file)];
	    if (file!=0)
	      sprintf(msg,"Parent object of mask has changed buffer association (file: %s, line: %d)",file,line);
	    else
	      strcpy(msg,
		     "Parent object of mask has changed buffer association");
	    throw MaskObjectException(msg);
	  } else
	    return false;
	}
	// Now we know that the base object still exists. Check timestamp
	// value
	if (ts!=baseObj->getTS()) {
	  if (!noExc) {
	    char msg[100+strlen(file)];
	    if (file!=0)
	      sprintf(msg,"Parent object of mask has changed size (file: %s, line: %d)",file,line);
	    else
	      strcpy(msg,"Parent object of mask has changed size");
	    throw MaskObjectException(msg);
	  } else
	    return false;
	}
      }

      return true;
    }

  protected:
    /**
     * Increases timestamp by 1.
     * Has to be called whenever object size changes!
     */
    virtual void incrTS() {
      ts++;
      if (ts==0) ts=1; // wrap around, but avoid 0
    }

    /**
     * @return Current timestamp
     */
    virtual uint getTS() const {
      return ts;
    }

    /**
     * Has to be called by superclass once it is assoc. with a new buffer.
     * The mem. watcher for the buffer has to be passed (can be 0 -> no
     * watcher).
     * NOTE: If there is a watcher and this object is a mask object, we
     * increase the watcher ref. counter by 1 here.
     * <p>
     * NOTE: The old value of 'buffWatch' is just overwritten here, there is
     * no de-assoc. Call 'deassocBuff' before to de-assoc. the old buffer!
     *
     * @param buffW New mem. watcher (or 0)
     */
    void assocBuff(MemWatchBase* buffW) {
      buffWatch=buffW;
      if (buffWatch!=0 && isMask()) {
#ifndef DEBUG_TRACKHANDLES
	buffWatch->incr(); // incr. ref. counter
#else
	buffWatch->incr(debugCause);
#endif
      }
    }

    /**
     * Has to be called by superclass once it deassoc. from its buffer. Two
     * things are done here:
     * - if this is not a mask, the 'isAssoc' flag in the mem. watcher is
     *   set to true (because this object is the buffer owner then)
     * - the ref. counter of the mem. watcher is decreased. If it drops to
     *   0, the watcher is destroyed (leads to buffer destruction).
     *   At end, 'buffWatch' is set to 0.
     */
    void deassocBuff() {
      if (buffWatch!=0) {
	if (!isMask()) buffWatch->isAssoc=true; // set flag
#ifdef DEBUG_TRACKHANDLES
	if (buffWatch->decr(debugCause)) {
#else
	if (buffWatch->decr()) {
#endif
	  delete buffWatch; // destroy watcher and buffer
	}
#ifdef DEBUG_TRACKHANDLES
	debugCause=0; // DEBUG: reset
#endif
	buffWatch=0;
      }
    }

  public:
    /**
     * Has to be implemented by subclasses!
     *
     * @return Is this a mask object?
     */
    virtual bool isMask() const = 0;

  protected:
    /**
     * Has to be implemented by subclasses!
     * NOTE: A non-mask object has to return 0! A mask may not have a
     * base object (it may just ref. to a memory region), in which case
     * 0 is returned.
     *
     * @return Pointer to base object (if we have one); or 0 otherwise
     */
    virtual const MatTimeStamp* getBaseObj() const = 0;

    /**
     * @param newts New timestamp value
     */
    void setTS(uint newts) {
      ts=newts;
    }

    /**
     * To be used by special subclass methods. Copies members managed here
     * from object 'src'. Then, the 'buffWatch' member of 'src' is set
     * to 0. The latter makes sure that subseq. actions on 'src' such as
     * its deletion do not affect the buffer watcher.
     *
     * @param src Source object
     */
    void copyMembers(MatTimeStamp* src) {
      ts=src->ts; buffWatch=src->buffWatch;
#ifdef DEBUG_TRACKHANDLES
      debugCause=src->debugCause;
#endif
      src->buffWatch=0;
    }
  };

//ENDNS

#endif
