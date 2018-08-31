function c=muldiag(a,b)
%MULDIAG Multiplies matrix with diagonal matrix
% C=MULDIAG(A,B) computes A*B where at least one of A,B must be a diagonal
%  matrix, repr. as a column vector.
%  NOTE: If one of A, B is a scalar, it has priority as diagonal matrix:
%  if A col. vec., B scalar, B is treated as diag. matrix.

[m,n]=size(a);
[p,q]=size(b);
if n==1
  if p==1 & q==1
    % Special case: B scalar
    c=b*a;
  elseif m~=p
    error('Wrong dimensions!');
  end
  c=a(:,ones(q,1)).*b;
elseif q==1
  if n~=p
    error('Wrong dimensions!');
  end
  vec=b';
  c=a.*vec(ones(m,1),:);
else
  error('One of A,B must be a diagonal matrix (column vector)!');
end
