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
 * Module: optimize
 * Desc.:  Header class OneDimSolver
 * ------------------------------------------------------------------- */

#ifndef ONEDIMSOLVER_H
#define ONEDIMSOLVER_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

/*
 * TODO:
 * - Solver which does not need derivatives
 */

#include "lhotse/optimize/default.h"
#include "lhotse/optimize/FuncOneDim.h"

//BEGINNS(optimize)
  /**
   * Collects static methods which numerically solve f(x)=0 for x from a
   * subset of \R. f(x) is encoded as a instance of a subclass of 'FuncOneDim'.
   * The bracket-based methods will usually guarantee that a solution is found
   * within a given initial bracket, up to a given accuracy, given that f(x)
   * is continuously differentiable, and f(l)*f(r)<0 for the initial bracket
   * [l,r].
   * <p>
   * Typically, subclasses of 'FuncOneDim' require the derivative to be
   * computed alongside with the function value. Methods which do not need this
   * information, will simply ignore the derivative value.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class OneDimSolver
  {
  public:
    // Constants

    static const int brackRightRegular =0;
    static const int brackRightBound   =1;
    static const int brackRightInfinite=2;

    // Public static methods

    /**
     * One-dimensional Newton solver. The function object 'func' has to
     * return derivatives.
     * <p>
     * Requires initial bracket [l,r], l=='l',
     * r=='r', l<r, and f(l)*f(r)<0 (but see below). The function f, which must
     * be continuously differentiable, is given by the object 'func'. The
     * solution is accurate up to +/- 'acc', and is returned by the method.
     * If a solution is not found after a certain number of iterations (suff.
     * large for any reasonable function), we terminate and throw a
     * 'NumericalException'. This is only likely to happen if f violates the
     * assumptions.
     * <p>
     * Not passing the right bracket end r:
     * If 'brRight'!='brackRightRegular', we ignore the value passed via 'r'
     * and try to find the zero without an initial right bracket end. The
     * search is decomposed into a method to find an initial right bracket
     * end, and the actual search.
     * NOTE: While if an initial bracket is given and f is cont. there and
     * has a zero in the bracket, the algorithm will always find it, this is
     * NOT guaranteed if 'r' is not given. This mode should be used ONLY for
     * special f(x), e.g. if the sign of f'(x) is constant right of l and the
     * opposite of the sign of f(l). A convex f(x) (right of l) is even better.
     * There are two different poss. here. Either we know a R>l s.t. f(x)
     * becomes unbounded for x\to R from the left (with a sign opposite of
     * f(l)). R cannot be passed for 'r', since f(R) is undefined. In this
     * case, use 'brRight'=='brackRightBound' and pass R via 'boundR'. We
     * assume that f(x) can be eval. for x <= R-'acc'. If f(x) does not change
     * sign between l and R-'acc', the algorithm will fail. Or no such R is
     * known. In this case, pass 'brRight'=='brackRightInfinite'.
     * 'boundR' is ignored in this case.
     * We use the following method to find a right bracket end:
     * The first step is of size 'acc'. Say, f(l)<0. Up from then, we have two
     * positions l,x, with l<x and f(l)<0. We are done if f(x)>0. We now fit
     * the quadratic q with same deriv. as f at x, same value as f at x and
     * l. If this is convex, we set q = 0 and solve for the next x. l becomes
     * the old x. If q is concave, we use the line through x with same deriv.
     * as f at x, set it equal to 0 and solve for x (this is a Newton step). If
     * the eval. of f(x) fails, we instead do a step of 'acc' only. Note that
     * our method does smaller steps than pure Newton would do. This is to
     * avoid overshooting. The minimum step size is always 'acc'. We run this
     * method until a right bracket end is found.
     * NOTE: If 'brRight'=='brackRightBound', we never do steps with
     * x > 'boundR'-'acc'.
     * <p>
     * Note that (if 'brRight'=='brackRightRegular') the method is globally
     * convergent. Every time the current Newton step would move out of the
     * bracket or the derivative is very small, or the previous Newton step led
     * to a relative bracket shrinkage of <1/10, we do a bisection step
     * instead.
     * <p>
     * NOTE: Code similar to NR routine 'rtsafe', but corrected. Can be
     * improved upon!
     *
     * @param func      Function f
     * @param l         Initial left bracket end
     * @param r         Initial right bracket end (ignored if 'rightInf'==true)
     * @param acc       Accuracy
     * @param brRight   See above. Def: 'brackRightRegular'
     * @param boundR    Required iff 'brRight'=='brackRightBound'. See above
     * @return          Solution
     */
    static double newton(FuncOneDim* func,double l,double r,double acc,
			 int brRight=brackRightRegular,double boundR=0.0);
  };
//ENDNS

#endif
