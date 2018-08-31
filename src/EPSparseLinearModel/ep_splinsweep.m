%EP_SPLINSWEEP EP for linear model: sweep over (all) sites
%
%  [{DELTA,{NUMUPD}}] = EP_SPLINSWEEP(X,U,SIGSQ,TAU,SITEPI,SITEB,
%    UPDIND,L,GAMMA,PITHRES,REINIT,WKVEC1,WKVEC2,WKVEC3,{METH=0},
%    {REFRESH=1},{FRAC=1})
%
%  Does a single EP sweep over selected sites of a sparse linear model with
%  Laplace prior. The likelihood is N(U | X*A, SIGSQ*EYE), and each prior
%  factor has the form
%    (TT/2)*EXP( -TT*ABS(A(i)) ),
%  where TT = TAU/SQRT(SIGSQ). For all details, see 'ExpectPropLinear' in
%  LHOTSE project module 'eplin'.
%
%  Input:
%  - X:       Design matrix [m-by-n]
%  - U:       Vector of responses [m]
%  - SIGSQ:   Noise variance
%  - TAU:     Scale parameter of Laplace prior
%  - SITEPI:  Site parameters pi_i [n] (*)
%  - SITEB:   Site parameters b_i [n] (*)
%  - UPDIND:  Index of sites to be updated. Entries in 1:n. Can have any
%             size.
%             Can be empty, for example if only LOGZ, LOGTILZ are to
%             be computed.
%  - L:       Cholesky factor of representation [n-by-n / m-by-m] (*)
%  - GAMMA:   Vector of representation [n / m] (*)
%  - PITHRES: Threshold det. when updates are skipped because too small
%             (see 'ExpectPropLinear::updateSite')
%  - REINIT:  Boolean flag. If true, the representation is recomputed
%             from the site pars. before updates are done
%  - WKVEC1:  Scratch vector [n] (*)
%  - WKVEC2:  Scratch vector [n] (*)
%  - WKVEC3:  Scratch vector [n] (*)
%  - METH:    Approximate inference method to be used:
%             0: Standard EP [default]
%             2: Tipping SBL
%             3: Girolami variational, VMFB
%  - REFRESH: If true, the repres. is recomp. from scratch after
%             updates are done. Optional. Def.: true
%  - FRAC:    Value in (0,1]. If FRAC<1, a fractional EP update is done.
%             See below. Can only be used with METH==0. Def.: 1
%
%  Return:
%  - DELTA:   Optional. Maximum rel. marginal change over all sites
%             (see 'ExpectPropLinear::updateSite')
%  - NUMUPD:  Optional. Number of updates done (updates are not done if
%             change in pi_i is below PITHRES)
%
%  We use the undocumented feature that the content of input arguments
%  can be overwritten to get a "call by reference". This is done for the
%  (*) arguments, whose content is overwritten. We cannot change the
%  size that way, so these matrices must be of correct size, even if their
%  input content is not used here.
%
%  Different representations:
%  If m<n, you can use the non-degenerate representation (L n-by-n, GAMMA n)
%  or the degenerate one (L m-by-m, GAMMA m). If m>=n, only the non-degen.
%  can be used.
%  The kind of representation to be used is determined by the size of
%  L, GAMMA passed (even if REINIT==1).
%  The non-degen. is more numerically stable in general, and is more
%  efficient than the degen. one for about m >= n/2. If m << n, use the
%  degen. one to save time. Otherwise, or if max(n,m) is not large, we
%  recommend using the non-degen. one.
%
%  REFRESH and REINIT:
%  The low-rank updates of the representation after each site introduce
%  numerical errors. It is wise to "refresh" the representation by recomp.
%  it from scratch after each n updates.
%  If REFRESH==true, a refresh is done at the end. If REINIT==true, the
%  repres. is recomp. at the beginning. The default usage of this function
%  is to use a permutation of 1:n in UPDIND, together with REFRESH==true
%  and REINIT==false, but REINIT==true for the first sweep in order to
%  initialize the representation.
%  NOTE: A refresh after each sweep is fairly conservative, but
%  subdominant anyway. A complete EP sweep takes far longer than a
%  refresh, although they both have the same complexity.
%
%  Fractional EP updates:
%  If FRAC<1, a fractional EP update is done instead of a standard one.
%  This leads to different fixed points. It can improve robustness of the
%  method. In some cases, fractional EP leads to better solutions than
%  standard EP (if both converge properly).

%  Copyright (C) 2005 Matthias Seeger
%
%  This program is free software; you can redistribute it and/or
%  modify it under the terms of the GNU General Public License
%  as published by the Free Software Foundation; either version 2
%  of the License, or (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%  USA.
