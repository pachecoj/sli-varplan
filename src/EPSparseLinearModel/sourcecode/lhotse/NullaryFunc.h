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
 * Desc.:  Header class NullaryFunc, BinderOnly
 * ------------------------------------------------------------------- */

#ifndef NULLARYFUNC_H
#define NULLARYFUNC_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include <functional>

/**
 * Similar to STL 'unary_function', but for nullary function (no arguments).
 * Declares operator() (without argument).
 * NOTE: Makes sense if there is an internal state which makes the function
 * non-constant (example: PRN generator, see 'Generator').
 * <p>
 * 'mybindarg' extends STL 'bind1st' to map 'unary_function' ->
 * 'NullaryFunc'. 'ptr_0fun' extends STL 'ptr_fun' to 'NullaryFunc'.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class Res> class NullaryFunc
{
public:
  typedef Res result_type;

  /**
   * Function eval. operator
   */
  virtual result_type operator()() const = 0;
};

/**
 * Same as STL's 'bind1st', but maps 'unary_function' to 'NullaryFunc'.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class UnOp> class BinderOnly :
  public NullaryFunc<typename UnOp::result_type>
{
protected:
  UnOp op;
  typename UnOp::argument_type arg;

public:
  BinderOnly(const UnOp& x,const typename UnOp::argument_type& v) :
    op(x),arg(v) {}

  typename UnOp::result_type operator()() const {
    return op(arg);
  }
};

template<class UnOp,class T>
inline BinderOnly<UnOp> mybindarg(const UnOp& op,const T& v)
{
  typedef typename UnOp::argument_type ArgType;
  return BinderOnly<UnOp>(op,ArgType(v));
}

/**
 * Wraps pointer to nullary function
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class Res> class Ptr2NullFunc : public NullaryFunc<Res>
{
protected:
  typedef Res (*FuncType)();
  FuncType func;

public:
  Ptr2NullFunc(FuncType ptr) : func(ptr) {}

  Res operator()() const {
    return (*func)();
  }
};

// Was 'ptr2nullfunc' before
template<class Res>
inline Ptr2NullFunc<Res> ptr_0fun(Res (*func)())
{
  return Ptr2NullFunc<Res>(func);
}

#endif
