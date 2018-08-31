function run_example
% This is a demo example for sparse linear regression and experiment design
%
% Consider a function from R, so-called time, to R which is a linear combination of n sine waves of
% different freqeuncies. Assume further that only a few of these basis
% functions have non/zeros coefficients. We aim at finding a ranking of
% coefficients which puts the non-zero coeefficients at the top of the
% list.
%
% Furthermore, we can ask which at which time t should we measure the
% funciton value next in order to get the best ranking. Since we now the
% basis functions and we have a current belief of the coefficients
% including uncertainties, we could determine such t for which the
% uncertainty is reduced maximally after the next measurement.
%
% Notation: Coefficients alpha, basis function evaluations x(t) = [x_1(t) ... x_n(t)]', 
% collectively stored in the design matrix X = [x(t_1) ... x(t_m)],
% outcomes y.
% It then hold that y = X * alpha + eps, where the measurement noise eps can either be
% modelled as Gaussian or Laplace distributed.
%
% Note: In constrast to the paper on active design for network recovery, we
% here can choose X and measure y. It therefore holds that the active design
% criteria like information gain or expected entropy reduction are all
% equivalent to the maximisation of the predictive variance. This is NOT
% true in the case needed for network recovery, namely where y can be
% designed and X is measured.


% overhead 
addpath(fileparts(pwd));
rand('state',0);
randn('state',0);

% parameters defining the model
n = 40;                     % number of sine waves
n_nonzeros = round(n/7);    % number of used sine waves
spread = 2;                 % maximum for randomlz determined amplitude of sine waves
t = (0:0.094:2*pi)';        % possible time steps
m = length(t);
sigma_noise = 0.01;         % measurement noise
edge_threshold = 0.1;       % if alpha_i > edge_threshold, we consider it as non-zeros (needed for evaluation, see BMC paper)

rounds = 10;
% loop over independent experiements whose results are averaged
roc = zeros(m,4);
droc = zeros(m,4);
for r = 1:rounds
    fprintf('Round %i\n',r);
    [ROC,PROB,X,alpha,Q_GD,Q_LD] = one_round(n,m,t,edge_threshold,sigma_noise,n_nonzeros,spread);
    roc = roc + ROC;
    droc = droc + ROC.^2;
end
roc = roc / rounds;
droc = droc / rounds - roc.^2;

% plot the (averaged) results of different methods

% average ROC curves with error bars
subplot(221)
errorbar(repmat((1:m)',1,4),roc,droc,droc);
legend({'GR','GD','LR','LD'})
ylabel iAUC
xlabel('Experiment number');
title('ROC after each experiment');


% plot the probabilities of 
subplot(222)
imagesc([PROB(:,:,2);PROB(:,:,4)]);
caxis([0,1]);
colorbar
xlabel('Experiment number');
ylabel('Coefficient number');
title(sprintf('Probabilities for non-zero\nTop half: GD, bottom half: LD'));

% plot the basis functions that were used
subplot(223);
plot(1:m,X*diag(alpha));
xlabel time
title('Used Basis functions');

% mean and variance of the the different coefficients
subplot(224);
val = [alpha,Q_GD{m,1},Q_LD{m,1}];
dval = [zeros(n,1),sqrt(diag(Q_GD{m,2})),sqrt(diag(Q_LD{m,2}))];
cla
barerrorbar(val,dval);
legend({'true','GD','LD'})
xlabel('Coefficient number');
title(sprintf('Mean and Variance after %i Experiments',m));


% ******************************************************************
% This function performs one independent simulation including data
% generation, coefficient recovery using different models, and evaluation
% of these models
% ******************************************************************
function [ROC,Prob,X,alpha,Q_GD,Q_LD] = ...
    one_round(n,m,t,edge_threshold,sigma_noise,n_nonzeros,spread)

% create true parameters and linear model matrices
X = [];
for i = 1:n, X = [X, sin(0.45*i*t)]; end
alpha = zeros(n,1);
alpha(irand(n_nonzeros,1,n)) = (rand(n_nonzeros,1)-0.5)*2*spread;
n_nonzeros = sum(alpha ~= 0);
fprintf('%i nonzeros\n',n_nonzeros);
y = X*alpha;y = y + randn(size(y))*sigma_noise;

% do the reconstruction for each data point
dataorder = randperm(m)';
times = zeros(1,4);
tic;Q_GR = Gaussian_NoDesign(X,y,dataorder,sigma_noise,n_nonzeros,edge_threshold);times(1) = toc;
tic;Q_GD = Gaussian_Design(X,y,sigma_noise,n_nonzeros,edge_threshold);times(2) = toc;
tic;Q_LR = Laplace_NoDesign(X,y,dataorder,sigma_noise,n_nonzeros,edge_threshold);times(3) = toc;
tic;Q_LD = Laplace_Design(X,y,sigma_noise,n_nonzeros,edge_threshold);times(4) = toc;
fprintf('Times: GR = %g, GD = %g, LR = %g, LD = %g\n',times(1),times(2),times(3),times(4));

% evalutate the different methods
ROC = zeros(m,4);
Prob = zeros(n,m+1,4);
[ROC(:,1),Prob(:,:,1)] = evalMethod(Q_GR,alpha,edge_threshold);
[ROC(:,2),Prob(:,:,2)] = evalMethod(Q_GD,alpha,edge_threshold);
[ROC(:,3),Prob(:,:,3)] = evalMethod(Q_LR,alpha,edge_threshold);
[ROC(:,4),Prob(:,:,4)] = evalMethod(Q_LD,alpha,edge_threshold);


% ******************************************************************
% Helper functions
% ******************************************************************

function [ROC,Prob] = evalMethod(Q,alpha,edge_threshold)
ROC = zeros(length(Q),1);
true_coeffs = abs(alpha) > edge_threshold;

Prob = [];
for i = 1:length(Q)
    coeffs_score = -logdiffgausscdf(Q{i,1},sqrt(diag(Q{i,2})),edge_threshold);    
    Prob = [Prob, (1-exp(-coeffs_score))];
    [void,void,auc,iAUC] = roc(true_coeffs,coeffs_score);
    ROC(i) = iAUC;
end
Prob = [Prob, true_coeffs];



function Q = Gaussian_NoDesign(X,y,dataorder,sigma_noise,n_nonzeros,edge_threshold)
Q = Gaussian_inference(X,y,sigma_noise,n_nonzeros,edge_threshold,dataorder);

function [Q,order] = Gaussian_Design(X,y,sigma_noise,n_nonzeros,edge_threshold)
[Q,order] = Gaussian_inference(X,y,sigma_noise,n_nonzeros,edge_threshold,@select_next_exp);

function Q = Laplace_NoDesign(X,y,dataorder,sigma_noise,n_nonzeros,edge_threshold)
Q = Laplace_inference(X,y,sigma_noise,n_nonzeros,edge_threshold,dataorder);

function [Q,order] = Laplace_Design(X,y,sigma_noise,n_nonzeros,edge_threshold)
[Q,order] = Laplace_inference(X,y,sigma_noise,n_nonzeros,edge_threshold,@select_next_exp);


function ind = select_next_exp(i,m,cov,X,order)
m = size(X,1);n = size(X,2);

% find all that X that are not used yet
poss = true(m,1);
poss(order(1:i-1)) = false;
poss = find(poss);
poss_X = X(poss,:);

% choose the one where the predictive variance is the highest
pred_cov = poss_X*cov*poss_X';
[void,ind] = max(sqrt(diag(pred_cov)));
ind = poss(ind);

        
function [Q,order] = Gaussian_inference(X,y,sigma_noise,n_nonzeros,edge_threshold,orderfunction)
m = size(X,1);n = size(X,2);

opt = optimset('largescale','off','Display','off');
sigma_prior = fminunc(@(x)((2*mycdf(0,x,-edge_threshold)-n_nonzeros/n)^2),1,opt);
% fprintf('Gaussian sigma = %g, (P_edge = %g)\n',sigma_prior,2*mycdf(0,sigma_prior,-edge_threshold));
fprintf('Gaussian sigma = %g\n',sigma_prior);

order = zeros(m,1);
Q = cell(m,2); 
for i = 1:m
    % check for next candidate
    if isnumeric(orderfunction)
        order(i) = orderfunction(i);
    else
        if i == 1,
            order(i) = irand(1,1,m);
        else
            order(i) = orderfunction(i,Q{i-1,1},Q{i-1,2},X,order);
        end
    end
    
    % update posterior with true measurement
    XX = X(order(1:i),:);
    yy = y(order(1:i));
    Q{i,1} = (XX'*XX + eye(n).*(sigma_noise/sigma_prior)^2) \ (XX'*yy);
    Q{i,2} = inv(XX'*XX + eye(n).*(sigma_noise/sigma_prior)^2)*sigma_noise.^2;
%     Q{i,1} = (XX'*XX./sigma_noise.^2 + eye(n)./sigma_prior.^2) \ (XX'*yy./sigma_noise.^2);
%     Q{i,2} = inv(XX'*XX./sigma_noise.^2 + eye(n)./sigma_prior.^2);
end


function [Q,order] = Laplace_inference(X,y,sigma_noise,n_nonzeros,edge_threshold,orderfunction)
m = size(X,1);n = size(X,2);

tau = 1/edge_threshold * log(n/n_nonzeros);
fprintf('Laplace length = %g (tau = 1/length = %g)\n',1/tau,tau);

pithres = 1e-16;
sitepi = repmat(0.5*(tau*sigma_noise)^2,n,1);
siteb = zeros(n,1);
L = zeros(n);
gamma = zeros(n,1);
wkvec1 = zeros(n,1);wkvec2 = zeros(n,1);wkvec3 = zeros(n,1);

order = zeros(m,1);
Q = cell(m,2); 
for i = 1:m
    % check for next candidate
    if isnumeric(orderfunction)
        order(i) = orderfunction(i);
    else
        if i == 1,
            order(i) = irand(1,1,m);
        else
            order(i) = orderfunction(i,Q{i-1,1},Q{i-1,2},X,order);
        end
    end
    
    % update posterior with true measurement
    XX = X(order(1:i),:);
    yy = y(order(1:i));
    for j = 1:20
        sweep_ord = randperm(n);
        delta=ep_splinsweep(XX,yy,sigma_noise^2,tau*sigma_noise,sitepi,siteb,sweep_ord,L,gamma,pithres,double(j==1),wkvec1,wkvec2,wkvec3,0,1,0.5);
%         For checking whether the algotithm converges, delta should go to
%         zero
%         fprintf(1,'i = %i j = %i delta=%f\n',i,j,delta);
    end

    % compute mean and covariance
    Q{i,1} = L'\gamma;
    Q{i,2} = L'\(L\eye(n))*sigma_noise^2;
end


function res = irand(n,min,max)
% generate a random integer number between min and max
res = floor((rand(n,1)-1e-10)*(max-min)) + min;