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
 * Desc.:  Header class MatlabAlloc
 * ------------------------------------------------------------------- */

#ifndef MATLABALLOC_H
#define MATLABALLOC_H

#include "lhotse/matif/mex_for_cpp.h" // MEX interface

/**
 * Replaces the STL default allocator with one that uses Matlab's MM
 * for (de)allocation and makes all alloc. memory persistent.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class T> class MatlabAlloc
{
public:
  typedef size_t     size_type;
  typedef ptrdiff_t  difference_type;
  typedef T*         pointer;
  typedef const T*   const_pointer;
  typedef T&         reference;
  typedef const T&   const_reference;
  typedef T          value_type;

  template<typename T2> struct rebind {
    typedef MatlabAlloc<T2> other;
  };

  MatlabAlloc() throw() {}

  MatlabAlloc(const MatlabAlloc&) throw() {}

  template<typename T2> MatlabAlloc(const MatlabAlloc<T2>&) throw() {}

  ~MatlabAlloc() throw() {}

  pointer address(reference __x) const {
    return &__x;
  }

  const_pointer address(const_reference __x) const {
    return &__x;
  }

  T* allocate(size_type __n, const void* = 0) {
    void* ret=mxMalloc(__n*sizeof(T));
    mexMakeMemoryPersistent(ret);
    //mexPrintf("MatlabAlloc::allocate\n");
    return (T*) ret;
  }

  void deallocate(pointer __p, size_type __n) {
    //mexPrintf("MatlabAlloc::deallocate\n");
    mxFree((void*) __p);
  }

  size_type max_size() const throw() {
    return size_t(-1)/sizeof(T);
  }

  void construct(pointer __p, const T& __val) {
    new(__p) T(__val);
  }

  void destroy(pointer __p) {
    __p->~T();
  }
};

template<> class MatlabAlloc<void>
{
public:
  typedef size_t      size_type;
  typedef ptrdiff_t   difference_type;
  typedef void*       pointer;
  typedef const void* const_pointer;
  typedef void        value_type;

  template<typename T2> struct rebind {
    typedef MatlabAlloc<T2> other;
  };
};

template<typename _T1,typename _T2> inline bool
operator==(const MatlabAlloc<_T1>&,const MatlabAlloc<_T2>&)
{
  return true;
}

template<typename _T1,typename _T2> inline bool
operator!=(const MatlabAlloc<_T1>&,const MatlabAlloc<_T2>&)
{
  return false;
}

#endif
