function err = l1_diff_symmetric( pdf1, pdf2, x_vals )
% L1_DIFF_SYMMETRIC - Returns the symmetric L1 error between two
%   continuous PDFs computed on a discrete grid.  The error measured is the
%   symmetric "information theoretic" L1 error in the range [0,1].
%
% Jason L. Pacheco
% 12/15/11
%

  dx = x_vals(2) - x_vals(1);
  diff_L1 = abs( pdf1 - pdf2 );
  err = 0.5 * sum( diff_L1 ) * dx;

end

