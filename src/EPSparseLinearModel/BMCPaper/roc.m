function [fp,tp,auc,aucN,aucPR,score] = roc(truth,val)
% function [fp,tp,auc,aucN,aucPR,score] = roc(truth,val)

% some overhead
truth = truth(:) ~= 0;
Nyes = sum(truth);
Nno = sum(~truth);
val = val(:);
N = length(val)+1;

% compute the counts of tps and fps
[void,ind] = sort(val,'descend');
list = truth(ind)';
tp = [0 cumsum(list)];
fp = [0 cumsum(~list)];
% fn = Nyes - tp;
% tn = Nno - fp;

x = fp / Nno;
y = tp / Nyes;
auc = integrate(x,y);
aucN = integrate_upto(x,y,Nyes/Nno);


% precision recall 
ii = (tp + fp) ~= 0; y = zeros(size(tp));
% precision = prob. that classifieds are correct
x(ii) = tp(ii) ./ (tp(ii) + fp(ii)); 
% recall = percentage of correctly found trues
y = tp / Nyes; 
[x,ii] = sort(x);y = y(ii);
aucPR = integrate(x,y);


% score = # recovered edges / # true edges (TPR)  at a precision (see above) = #false
% predictions / # predictions = 10%
ii = (tp + fp) ~= 0; y = ones(size(tp));
x(ii) = fp(ii) ./ (tp(ii) + fp(ii)); 
y = tp / Nyes;
[x,ii] = sort(x);y = y(ii);
score = interpolate(x,y,0.1);


function ya = interpolate(x,y,a)
if a <= x(1)
    ya = y(1);
elseif a >= x(end)
    ya = y(end);
else
    dx = x(2:end)-x(1:end-1);
    ii = sum(x <= a);
    lambda = (a - x(ii)) /dx(ii);
    ya =  ((1-lambda/2)*y(ii)+lambda/2*y(ii+1));
end

function A = integrate(x,y)
dx = x(2:end)-x(1:end-1);
meanval = (y(2:end)+y(1:end-1)) / 2;
A = sum(dx .* meanval) + x(1)*y(1) + (1-x(end))*y(end);

function A = integrate_upto(x,y,a)
if a <= x(1)
    A = y(1)*a;
elseif a >= x(end)
    A = integrate(x,y) - y(end)*(1-a);
else
    dx = x(2:end)-x(1:end-1);
    meanval = (y(2:end)+y(1:end-1)) / 2;
    ii = sum(x <= a);
    lambda = (a - x(ii)) /dx(ii);
    mmeanval =  ((1-lambda/2)*y(ii)+lambda/2*y(ii+1));
    A = sum(dx(1:ii-1) .* meanval(1:ii-1)) + x(1)*y(1) + mmeanval * dx(ii)*lambda;
end
A = A / a;
