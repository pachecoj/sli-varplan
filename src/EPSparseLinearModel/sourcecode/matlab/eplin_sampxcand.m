function xcand = eplin_sampxcand(x,u,sigsq,tau,sitepi,siteb,l, ...
				 gamma.ucand)
%EPLIN_SAMPXCAND Samples x_* candidates, given u_* candidates
%
%  XCAND = EPLIN_SAMPXCAND(X,U,SIGSQ,TAU,SITEPI,SITEB,L,GAMMA,
%                          UCAND)
%
%  Gene network recovery application. n independent linear models,
%  one for each row of the network matrix A. Data is X, U (both
%  m-by-n). SITEPI, SITEB, L, GAMMA as in EP_SPLINSWEEP,
%  EPLIN_COMPSCORES, but they are cell arrays of size n, one entry
%  for each (row) model.
%  Here, for a set of candidates u_* (cols of UCAND), corresponding
%  candidates x_* are sampled and returned (cols of XCAND). This is
%  done by sampling a network matrix A from the current posterior.
%  Average over several calls (using different A's) to get stable
%  score values.

[m,n]=size(x);
[d1,d2]=size(u);
if d1~=m || d2~=n
  error('U wrong size');
end
if sigsq<=0 || tau<=0
  error('SIGSQ, TAU wrong');
end
[d1,num]=size(ucand);
if d1~=n
  error('UCAND wrong size');
end
d1=length(gamma{1});
if d1==n
  deg_repr=0;
elseif d1==m
  deg_repr=1;
else
  error('GAMMA wrong size');
end

% Sample A matrix
amat=zeros(n,n);
for i=1:n
  svec=sqrt(sigsq)*randn(n,1);
  if ~deg_repr
    % Non-degenerate representation
    amat(i,:)=(l{i}'\(svec+gamma{i}))';
  else
    % Degenerate representation
    [evecs,tmat]=eig(muldiag(x,1./sitepi{i})*x'+eye(m));
    tvec=sqrt(diag(tmat));
    svec=svec./sqrt(sitepi{i});
    rvec=(1./tvec)./(tvec+1);
    tvec=svec-(x'*muldiag(evecs,rvec)*evecs'*(x*svec))./sitepi{i};
    hvec=(x'*(u(:,i)-l{i}'\gamma{i})+siteb{i})./sitepi{i};
    amat(i,:)=(sqrt(sigsq)*tvec+hvec)';
  end
end

% Sample XCAND
xcand=amat\(ucand+sqrt(sigsq)*randn(n,num));
