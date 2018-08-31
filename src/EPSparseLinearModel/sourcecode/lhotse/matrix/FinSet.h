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
 * Desc.:  Header class FinSet
 * ------------------------------------------------------------------- */

#ifndef FINSET_H
#define FINSET_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

/* TODO!!
 * - other set services (union!)
 */

#include "lhotse/matrix/default.h"
#include "lhotse/matrix/BaseVector.h"

//BEGINNS(matrix)
  /**
   * Maintains a finite set of elements of T. T must implement comparison
   * operators. The set is kept as a sorted vector internally. New elements
   * are inserted one by one, det. their position by binary search. Larger
   * numbers of new elements are inserted by sorting the combined vectors,
   * then removing duplicates.
   * <p>
   * ATTENTION: The sorted set is kept as a handle and the copy constructor
   * just copies the handle, does NOT draw a copy of the set! If the set
   * is modified, all 'FinSet' copies are changed!
   * ==> Allows to use 'FinSet' as unary STL predicate (function objects
   *     are copied frequently)
   * <p>
   * NOTE: Cannot represent multi-sets, duplicates are not allowed.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  template<class T> class FinSet : public std::unary_function<T,bool>
  {
  protected:
    // Members

    Handle<BaseVector<T> > vec; // Vector repres. the set

  public:
    // Constructors

    /**
     * Default constructor. Creates empty set
     *
     * @param vecP Initial content. Def.: empty
     */
    FinSet() : vec(new BaseVector<T>()) {}

    /**
     * Copy constructor.
     * ATTENTION: Just the handle is copied. 'src' and new object share
     * the same set vector.
     */
    FinSet(const FinSet& src) : vec(src.vec) {}

    /**
     * @param vecP Init. set (need not be sorted)
     */
    FinSet(const BaseVector<T>& vecP) : vec(new BaseVector<T>()) {
      *vec=vecP;
      sortRemDuplicates();
    }

    // Public methods

    /**
     * Makes this class a unary predicate
     *
     * @return Is 'elem' contained?
     */
    bool operator()(const T& elem) const {
      int pos;
      return findPos(elem,pos);
    }

    /**
     * Inserts element into set if not already there. The insert position
     * is returned (-1 if the element is already contained).
     *
     * @param elem New element
     * @return     Insert position (-1 if 'elem' already contained)
     */
    virtual int insert(const T& elem) {
      int pos;
      if (!findPos(elem,pos)) {
	vec->insert(elem,1,pos);
	return pos;
      } else
	return -1;
    }

    /**
     * Resets to empty set
     */
    virtual void reset() {
      vec->fill(0);
    }

    /**
     * Inserts set of elements from 'src' into set. This is done either one
     * by one or by appending and calling 'sortRemDuplicates'.
     *
     * @param src Set to insert
     * @return    How many elements have been inserted?
     */
    virtual int insert(const BaseVector<T>& src);

    /**
     * @return Ref. to set as vector (sorted in asc. order)
     */
    virtual const BaseVector<T>& getVec() const {
      return *vec;
    }

    /**
     * @return Set size
     */
    virtual int size() const {
      return vec->size();
    }

  protected:
    // Internal methods

    /**
     * Finds position of 'elem' in the set vector 'vec'. We return true iff
     * the element is contained. The position is ret. via 'pos'. If 'elem'
     * is not contained, 'pos' returns the position where 'elem' should be
     * included to keep the new vector sorted.
     *
     * @param elem Element
     * @param pos  Position ret. here
     * @return     Is 'elem' contained?
     */
    bool findPos(const T& elem,int& pos) const;

    /**
     * Sorts set vector 'set' from scratch, then removes all duplicates.
     *
     * @return How many duplic. have been removed?
     */
    int sortRemDuplicates();
  };

  template<class T> inline int
  FinSet<T>::insert(const BaseVector<T>& src)
  {
    int sz=src.size(),n=vec->size(),num;
    if (sz>(int) log((double) (n+1))) {
      vec->insert(src);
      num=sz-sortRemDuplicates();
    } else {
      int num=0;
      for (int i=0; i<sz; i++)
	if (insert(src[i])!=-1) num++;
    }

    return num;
  }

  template<class T> inline bool
  FinSet<T>::findPos(const T& elem,int& pos) const
  {
    int ret;
    if (vec->size()==0) return -1;
    pos=vec->findbin(elem,true,&ret); // binary search
    if (pos!=-1) return true;
    else {
      pos=(ret!=-1)?ret:vec->size();
      return false;
    }
  }

  template<class T> inline int FinSet<T>::sortRemDuplicates()
  {
    vec->sort(); // sort in asc. order
    int i,n=vec->size(),pos;
    T lastEl=(*vec)[0],elem,numD=0;

    for (i=1; i<n; i++) {
      pos=i;
      while ((elem=(*vec)[pos])==lastEl) {
	if (++pos==n) break;
      }
      if ((pos-=i)>0) {
	vec->remove(i,pos); n-=pos; numD+=pos;
      }
      lastEl=elem;
    }

    return numD;
  }

//ENDNS

#endif
