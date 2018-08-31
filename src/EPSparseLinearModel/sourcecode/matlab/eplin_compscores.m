function scores = eplin_compscores(x,u,sigsq,tau,sitepi,siteb,l, ...
				   gamma,xcand,ucand)
%EPLIN_COMPSCORES Compute information gain scores for sparse linear model
%
%  SCORES = EPLIN_COMPSCORES(X,U,SIGSQ,TAU,SITEPI,SITEB,L,GAMMA,
%                            XCAND,UCAND)
%
%  The present linear model is given by the data X, U and the
%  parameters TAU, SIGSQ (see EP_SPLINSWEEP). The site parameters
%  are SITEPI, SITEB, and the representation is L, GAMMA. Depending
%  on the size of L, GAMMA, we use the non-degenerate or degenerate
%  representation. There are n variables and m datapoints, 
%  Here, we compute information gain scores, one for each candidate
%  (x_*,u_*). Here, x_* is a n-vector, and u_* a scalar. The x_*
%  are columns in XCAND.
%
%  NOTE: The gene network recovery application uses n independent
%  linear models with common X, but different response vectors
%  U. The candidates x_* are sampled based on candidate n-vectors
%  u_* and the current posterior, see EPLIN_SAMPXCAND. The score
%  for a pair (x_*,u_*) is the sum of scores computed here, passing
%  the different row models and entries of u_*.
%  The x_* for all u_* are sampled based on the same network matrix
%  drawn from the current posterior. Stable scores are obtained by
%  averaging over a few calls to this method (and EPLIN_SAMPXCAND).

[m,n]=size(x);
[d1,d2]=size(u);
if d1~=m || d2~=1
  error('U wrong size');
end
if sigsq<=0 || tau<=0
  error('SIGSQ, TAU wrong');
end
[d1,d2]=size(sitepi);
if d1~=n || d2~=1
  error('SITEPI wrong size');
end
[d1,d2]=size(siteb);
if d1~=n || d2~=1
  error('SITEB wrong size');
end
[d1,d2]=size(l);
if d1~=d2
  error('L wrong size');
end
if d1==n
  deg_repr=0;
elseif d1==m
  deg_repr=1;
else
  error('L wrong size');
end
[d2,d3]=size(gamma);
if d1~=d2 || d3~=1
  error ('GAMMA wrong size');
end
[d1,num]=size(xcand);
if d1~=n
  error('XCAND wrong size');
end
[d1,d2]=size(ucand);
if d1~=num || d2~=1
  error('UCAND wrong size');
end

% Compute ALPHA, BETA
if ~deg_repr
  % Non-degenerate representation
  vmat=l\xcand;
  alpha=sum(vmat.*vmat,1)'+1;
  beta=vmat'*gamma;
else
  % Degenerate representation
  vmat=muldiag(1./sitepi,xcand);
  alpha=sum(xcand.*vmat,1)'+1;
  beta=vmat'*(x'*u+siteb);
  vmat=l\(x*vmat);
  alpha=alpha-sum(vmat.*vmat,1)';
  beta=beta-vmat'*gamma;
end

% Compute SCORES
tvec=(ucand-beta)/sqrt(sigsq);
scores=0.5*(log(alpha)+(tvec.*tvec./alpha-1).*(alpha-1)./alpha);
