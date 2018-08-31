function y = logdiffgausscdf(m,s,e)
%LOGDIFFGAUSSCDF Log of difference of Gaussian cdf's
%  Y=LOGDIFFGAUSSCDF(M,S,E) computes log of probability for A to lie
%  in [-E,E], where A is Gaussian with mean M, standard
%  deviation S. This is done in a stable way, using an asymptotic
%  expansion of the error function for large arguments.

n=size(m,1);
if size(m,2)~=1 || size(s,2)~=1 || size(s,1)~=n
  error('Wrong size');
end
y=zeros(n,1);
m=abs(m);
argp=(m+e)./s; argm=(m-e)./s;

ok=(argm<5);
okpos=find(ok);
% Compute y = F(argm), F(x) = log(1-Phi(x))
if ~isempty(okpos)
  % argm < 5
  y(okpos)=log(gausscdf(-argm(okpos)));
  temp=argp(okpos);
  okp=(temp>=5);
  ind=find(okp);
  prempos=okpos(ind); % argm < 5 <= argp
  ind=find(1-okp);
  pokpos=okpos(ind); % argp < 5
else
  prempos=[]; pokpos=[];
end
if length(okpos)<n
  % argm >= 5
  mrempos=find(1-ok); % 5 <= argm
  temp=argm(mrempos);
  temp2=(1./temp)./temp;
  gm=1-temp2.*(1-3*temp2.*(1-5*temp2.*(1-7*temp2)));
  y(mrempos)=-0.5*temp.*temp-log(temp)+log(gm)-0.5*log(2*pi);
else
  mrempos=[];
end
mrlen=length(mrempos);
if ~isempty(mrempos) || ~isempty(prempos)
  ind=[mrempos; prempos];
  temp=argp(ind);
  temp=(1./temp)./temp;
  gp=1-temp.*(1-3*temp.*(1-5*temp.*(1-7*temp)));
end
ratio=zeros(n,1);
if ~isempty(mrempos)
  % Tail approx. for both ratio parts
  mm=m(mrempos); ms=s(mrempos);
  ratio(mrempos)=exp(-(2*e*mm./ms)./ms).*(1-2*e./(mm+e)).* ...
      gp(1:mrlen)./gm;
end
if ~isempty(prempos)
  % Tail approx. for argp part of ratio
  temp=argp(prempos);
  ratio(prempos)=exp(-0.5*temp.*temp-log(temp)+ ...
		     log(gp((mrlen+1):end))-0.5*log(2*pi)- ...
		     y(prempos));
end
if ~isempty(pokpos)
  % No tail approx. for ratio
  ratio(pokpos)=exp(log(gausscdf(-argp(pokpos)))-y(pokpos));
end
y=y+log(1-ratio);
