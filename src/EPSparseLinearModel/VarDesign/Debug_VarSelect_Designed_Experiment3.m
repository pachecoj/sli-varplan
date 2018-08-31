function ind = Debug_VarSelect_Designed_Experiment3(i,X,U,sitepi,siteb,L,gamma,candU,sigma_noise,order,num_initial_rand)
  if i <= num_initial_rand
    ind = Select_Random_Experiment(i,X,U,sitepi,siteb,L,gamma,candU,sigma_noise,order);
  else
    
    n = size(X,1);    
    
    % find all that candU that are not used yet
    poss = true(size(candU,2),1);
%     poss(order(1:i-1)) = false;
    poss = find(poss);
    poss_candU = candU(:,poss);
    
    % precompute EP posterior natural params
    eta_a = cell(n,1);
    Lam_a = cell(n,1);
    for j = 1:n
      m_a = (L{j}'\gamma{j}); 
      Lam_a{j} = L{j} * L{j}' / sigma_noise^2;
      eta_a{j} = Lam_a{j} * m_a; % L{j} * gamma{j} / sigma_noise^2; %
    end
    
    % compute score for each U candidate    
    scores = zeros(size(poss_candU,2),1);
    scores_ADF = zeros(size(poss_candU,2),1);
    for ui = 1:size(poss_candU,2)
      fprintf('Scoring candidate %i of %i\n', ui, size(poss_candU,2));
      u = poss_candU(:,ui);
      
      % run EP/ADF
      max_iters = 1;
      conv_thresh = 1e-4;
      Lam_a_cond_x = EP_VarPlan( ...
        eta_a, Lam_a, u, sigma_noise, max_iters, conv_thresh );
      Lam_a_cond_x_ADF = EP_VarPlan( ...
        eta_a, Lam_a, u, sigma_noise, 1, conv_thresh );
        
      % compute mutual information
      for j=1:n
        scores(ui) = scores(ui) - 0.5 * sum( log( eig( Lam_a{j} ) ) ) ...
          + 0.5 * sum( log( eig( Lam_a_cond_x{j} ) ) );
        scores_ADF(ui) = scores_ADF(ui) - 0.5 * sum( log( eig( Lam_a{j} ) ) ) ...
          + 0.5 * sum( log( eig( Lam_a_cond_x_ADF{j} ) ) );
      end
    end
    
    [~,ind] = max( scores );
    ind = poss( ind );
    
    %% DEBUG: Plot EP / ADF estimates of MI
    figure('InvertHardcopy','off','Color',[1 1 1]);
    set(gca,'FontSize',14);
    set(gcf,'Position',[230   202   990   425]);
    hold on
    yval = NaN(size(candU,2),1);
    yval( poss ) = scores; % / n;
    yval_adf = NaN(size(candU,2),1);
    yval_adf( poss ) = scores_ADF;
    h_adf = plot( 1:numel(yval_adf), yval_adf, '-b', 'LineWidth', 2);
    h_ep = plot( 1:numel(yval), yval, '-k', 'LineWidth', 2);
    plot( 1:numel(yval), yval, '*b')
    plot( ind, scores(ind), '*r' );    
    xlabel('Intervention');
    ylabel('Mutual Information');    
    
    %% DEBUG: Plot sample-based MI estimates
    [scoresLV, scoresVI] = Debug_Select_Designed_Experiment(...
      10,500,i,X,U,sitepi,siteb,L,gamma,candU,sigma_noise,order,num_initial_rand);
    scoresLV = scoresLV .* n;    
    yval = NaN(size(candU,2),1);      
    yval( poss ) = mean( scoresLV, 2 );         
    ystd = NaN(size(candU,2),1);
    ystd( poss ) = std( scoresLV, 0, 2 );    
    h_mcmc = errorbar( 1:numel(yval), yval, ystd, '-s', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red');
    
    %% DEBUG: Plot variational MI bound with sample-based variance
    yval = NaN(size(candU,2),1);
    yval( poss ) = mean( scoresVI, 2 );
    ystd = NaN(size(candU,2),1);
    ystd( poss ) = std( scoresVI, 0, 2 );    
    xlim([0 (numel(yval)+1)])
    h_vi = errorbar( 1:numel(yval), yval, ystd, '-s', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red');
    xlabel('Intervention')
    ylabel('Mutual Information')
    title(sprintf('MI Estimates (Experiment %d)',size(X,2)));
    if max_iters>1
      legend([h_adf, h_ep, h_mcmc, h_vi], {'ADF','EP','MCMC','MCMCVI'});
    else
      legend([h_adf, h_ep, h_mcmc, h_vi], {'ADF','ADF','MCMC','MCMCVI'});
    end
    export_fig('~/Research/journal/121717/MI_est.pdf','-append');
    close;
  end
end


function [ Lam_a_cond_x, Lam_x, Lam_fac ] = EP_VarPlan( eta_a, Lam_a, u, sigma_noise, max_iters, conv_thresh )
  
  n = numel( Lam_a );
  Nruns = 20;
  MOMENTS_TYPE = 1;
  damp = 1.0;
  
  % optimization settings
  opt = optimoptions('fminunc','Display','off','GradObj','on','Hessian','off','DerivativeCheck','off',...
      'Algorithm','trust-region','TolFun',eps,'MaxIter',5000,'MaxFunEvals',5000);
    
  % declare stuff
  Lam_msg = cell(n,1);
  eta_msg = cell(n,1);
  Lam_x = zeros(n,n);  
  Lam_fac = cell(n,1);
  eta_fac = cell(n,1);   
  Lam_a_cond_x = cell(n,1);
  
  % initialize messages/marginal on X
  for i=1:n
    Lam_msg{i} = 1e-7/n * eye(n);
    eta_msg{i} = zeros(n,1);
    Lam_x = Lam_x + Lam_msg{i};    
  end
  eta_x = zeros(n,1);
    
  % do EP/ADF  
  for iter=1:max_iters
    if (max_iters>1), fprintf('EP Planning Iter %i of %i\n',iter,max_iters); end
    Lam_fac_old = Lam_fac;
    eta_fac_old = eta_fac;       
    
    % set update damping
    this_damp = 1.0;
    if iter>1, this_damp = damp; end

    % update factors
    for i=1:n
      
      % cavity "distribution" parameters
      Lam_cav = Lam_x - Lam_msg{i};
      eta_cav = eta_x - eta_msg{i};
%       Lam_cav = psdproj( Lam_cav, [] );
      
      % check conditioning of cavity distribution
      evals = eig( Lam_cav );
      max_e = max( evals );
      min_e = min( evals );
      if( (min_e < 0) || ( abs( min_e / max_e ) < eps ) )
        fprintf('Skipping node %d.\n', i);
        continue;
      end
                  
      %
      % JOINT OPTIMIZATION OVER (X,A)
      %      
      switch MOMENTS_TYPE
        
        % Laplace Approximation
        case 1
          [Lam_xa, eta_xa, V_xa, m_xa] = LaplaceMoments( ...
            Nruns, Lam_a{i}, eta_a{i}, u(i), Lam_cav, eta_cav, sigma_noise, opt );
          
        % Importance Sampling          
        case 2
          num_samples = 5000;
          [Lam_xa, eta_xa, V_xa, m_xa] = ISmoments( ...
            num_samples, Lam_a{i}, eta_a{i}, u(i), Lam_cav, eta_cav, sigma_noise );
          
        % Gibbs Sampling
        case 3
          num_samples = 1000;
          burnin = 1000;
          thinning = 10;
          [Lam_xa, eta_xa, V_xa, m_xa] = GibbsMoments( ...
            num_samples, burnin, thinning, Lam_a{i}, eta_a{i}, u(i), Lam_cav, eta_cav, sigma_noise );          
      end          
                  
      % update factor approx.
      Lam_fac{i} = Lam_xa - blkdiag( Lam_cav, Lam_a{i} );
      eta_fac{i} = eta_xa - [ eta_cav; eta_a{i} ];      
      Lam_a_cond_x{i} = Lam_xa((n+1):end,(n+1):end);
                        
      % marginal on X
      V_x = V_xa( 1:n, 1:n );
      m_x = m_xa(1:n);
      Lam_x = inv( V_x );      
      eta_x = V_x \ m_x;
      
      % update message from i^th factor      
      Lam_msg_new = Lam_x - Lam_cav;
      Lam_msg{i} = this_damp * Lam_msg_new + (1-this_damp) * Lam_msg{i};
      eta_msg_new = eta_x - eta_cav;
      eta_msg{i} = this_damp * eta_msg_new + (1-this_damp) * eta_msg{i};
      
      % reconstruct marginal from damped messages
      if this_damp<1.0
        thisLam_x = zeros( size( Lam_x ) );
        thisEta_x = zeros( size( eta_x ) );
        for j=1:n
          thisLam_x = thisLam_x + Lam_msg{j};
          thisEta_x = thisEta_x + eta_msg{j};
        end        
        Lam_x = thisLam_x;
        eta_x = thisEta_x;
      end
                            
    end        
  end
end


function [Lam_xa, eta_xa, V_xa, m_xa] = LaplaceMoments( ...
  Nruns, Lam_a, eta_a, u, Lam_cav, eta_cav, sigma_noise, opt )

  n = numel( eta_a );

  % maximize log-augmented distribution
  xaOpt = LaplaceOpt(Nruns, u, Lam_a, eta_a, Lam_cav, eta_cav, sigma_noise, opt);
  xOpt = xaOpt(1:n);
  aOpt = xaOpt((n+1):end);
  
  % compute Hessian and project to PSD / symmetric matrix
  Lam_xa = - LLHess(xOpt, aOpt, u, Lam_a, Lam_cav, sigma_noise);
  Lam_xa = psdproj( Lam_xa, [] );  
  eta_xa = Lam_xa*xaOpt;
  
  if nargout>2
    V_xa = inv( Lam_xa );
    m_xa = Lam_xa \ eta_xa;  %% TODO: Numerically unstable.  Should be equal to xaOpt
  end
  
end


function [Lam_xa, eta_xa, V_xa, m_xa] = ISmoments( ...
  num_samples, Lam_a, eta_a, u, Lam_cav, eta_cav, sigma_noise )

  n = numel( eta_a );

  % Compute moment parameters
  V_a = inv( Lam_a );
  m_a = Lam_a \ eta_a;
  m_x_cav = Lam_cav \ eta_cav;
  V_cav = inv(Lam_cav); 

  % draw importance samples     
  log_ws = zeros(num_samples,1);
  xasamp = zeros(2*n, num_samples);
  for s = 1:num_samples
    xsamp = mvnrnd(m_x_cav, V_cav)';
    asamp = mvnrnd(m_a, V_a)';
    xasamp(:,s) = [ xsamp; asamp ];
    log_ws(s) = - LLobj(xasamp(:,s), u, Lam_a, eta_a, Lam_cav, eta_cav, sigma_noise);
    log_ws(s) = log_ws(s) - 0.5 * (xsamp - m_x_cav)' * Lam_cav * (xsamp - m_x_cav);
    log_ws(s) = log_ws(s) - 0.5 * (asamp - m_a)' * Lam_a * (asamp - m_a);
  end
  
  % sample-based moment estimates
  ws = exp( log_ws - max( log_ws ) );
  ws = ws / sum( ws );
  m_xa = sum(bsxfun(@times, xasamp, ws'), 2);
  xx = bsxfun(@times, reshape(xasamp,2*n,1,[]), reshape(xasamp,1,2*n,[]));
  v_noncent_is = sum(bsxfun(@times, xx, reshape(ws,1,1,[])),3);
  V_xa = v_noncent_is - m_xa * m_xa';
  V_xa = psdproj( V_xa, [] );
  Lam_xa = V_xa \ eye(2*n);
  eta_xa = V_xa \ m_xa;
  
end


function [Lam_xa, eta_xa, V_xa, m_xa] = GibbsMoments( ...
  num_samples, burnin, thinning, Lam_a, eta_a, u, Lam_cav, eta_cav, sigma_noise )

  n = numel( eta_a );
  
  % Compute moment parameters
  m_x_cav = Lam_cav \ eta_cav;
  V_cav = inv(Lam_cav); 
  
  % run sampler
  xsamp = mvnrnd(m_x_cav, V_cav)';
  xasamp = zeros(2*n, num_samples);
  s = 1;
  for iter=1:(burnin+(thinning*num_samples))
        
    % sample a
    Lam_a_cond = Lam_a + xsamp * xsamp' / sigma_noise^2;    
    eta_a_cond = eta_a + xsamp * u / sigma_noise^2;    
    V_a_cond = Lam_a_cond \ eye(n);
    m_a_cond = Lam_a_cond \ eta_a_cond;    
    V_a_cond = psdproj( V_a_cond, [] );
    asamp = mvnrnd(m_a_cond', V_a_cond)';
    
    % sample x
    Lam_x_cond = Lam_cav + asamp * asamp' / sigma_noise^2;
    eta_x_cond = eta_cav + asamp * u / sigma_noise^2;
    V_x_cond = inv( Lam_x_cond );
    m_x_cond = Lam_x_cond \ eta_x_cond;
    V_x_cond = psdproj( V_x_cond, [] );
    xsamp = mvnrnd(m_x_cond', V_x_cond)';
    
    % keep sample?
    iter_adjusted = iter - burnin;
    if (iter>=burnin) && (mod(iter_adjusted, thinning)==1)
      xasamp(:,s) = [ xsamp; asamp ];
      s = s+1;
    end    
  end  
  
  % compute moments
  m_xa = mean( xasamp, 2 );
  xx = bsxfun(@times, reshape(xasamp,2*n,1,[]), reshape(xasamp,1,2*n,[]));
  Exx = mean( xx, 3 );  
  V_xa = Exx - m_xa * m_xa';
  V_xa = psdproj( V_xa, [] );
  Lam_xa = V_xa \ eye(2*n);
  eta_xa = V_xa \ m_xa;
  
end



function [xaOpt, Hopt] = LaplaceOpt(Nruns, u, Lam_a, eta_a, Lam_cav, eta_cav, sigma_noise, opt)

  fvalBest = Inf;
  xaOpt = [];
  Hopt = [];
  n = numel( eta_cav );
  
%   figure('InvertHardcopy','off','Color',[1 1 1]);
%   hold on
%   title('Argmax (x^{*},a^{*})')  
%   fval_all = zeros(Nruns,1);
%   eig_all = zeros(Nruns,n);
%   MI_all = zeros(Nruns,1);
  
  for i=1:Nruns
    xa_0 = -10 + 20 * rand( 2*n, 1 );
    [xaOptTmp, fval, ~, output, ~, H] = fminunc( ...
      @(X) LLobj(X, u, Lam_a, eta_a, Lam_cav, eta_cav, sigma_noise), xa_0, opt);
%     [xOptTmp, fval, ~, output] = fminunc( ...
%       @(X) LLobjX(X, u, Lam_a, eta_a, Lam_cav, eta_cav, sigma_noise), xa_0(1:n), opt);
%     aOptTmp = (Lam_a + xOptTmp*xOptTmp'./sigma_noise^2)\(eta_a + xOptTmp*u/sigma_noise^2);
%     xaOptTmp = [ xOptTmp; aOptTmp ];
    if ( fval < fvalBest )
      fvalBest = fval;
      xaOpt = xaOptTmp;
      Hopt = H;
    end
    
    %% DEBUG
%     plot(1:numel(xaOptTmp), xaOptTmp, '-b');
%     fval_all(i) = fval;
%     H = - LLHess(xaOptTmp(1:n), xaOptTmp((n+1):end), u, Lam_a, Lam_cav, sigma_noise);
%     Lam_a_cond_x = H((n+1):end,(n+1):end);
%     eig_all(i,:) = eig(Lam_a_cond_x);
%     MI_all(i) = sum( log( eig( Lam_a_cond_x ) ) );
  end 
  
%   m_a = Lam_a \ eta_a;
%   m_x_cav = Lam_cav \ eta_cav;
%   xa_0 = [ m_x_cav; m_a ];
%   [xaOptTmp, fval, ~, output] = fminunc( ...
%       @(X) LLobj(X, u, Lam_a, eta_a, Lam_cav, eta_cav, sigma_noise), xa_0, opt);
%   if ( fval < fvalBest )
%     fvalBest = fval;
%     xaOpt = xaOptTmp;
%   end
  
%   %% DEBUG
%   plot(1:numel(xaOptTmp), xaOptTmp, '-b');
%   fval_all(end) = fval;  
%   H = - LLHess(xaOptTmp(1:n), xaOptTmp((n+1):end), u, Lam_a, Lam_cav, sigma_noise);
%   Lam_a_cond_x = H((n+1):end,(n+1):end);
%   eig_all(end,:) = eig(Lam_a_cond_x);
%   MI_all(end) = sum( log( eig( Lam_a_cond_x ) ) );
%   figure('InvertHardcopy','off','Color',[1 1 1]);
%   plot(1:Nruns, fval_all, '-k');
%   title('Objective');
%   xlabel('Run');
%   figure('InvertHardcopy','off','Color',[1 1 1]);
%   hold on
%   for i=1:Nruns
%     plot(1:n, eig_all(i,:), '-b');    
%   end
%   title('Hessian Eigenvalues');
%   figure('InvertHardcopy','off','Color',[1 1 1]);
%   hold on
%   for i=1:Nruns
%     plot(1:Nruns, MI_all, '-b');    
%   end
%   title('Mutual Information (sort of)');
%   xlabel('Run');
  
end

function [negLogP, g, H] = LLobj(xa, u, Lam_a, eta_a, Lam_x, eta_x, sigma_noise)
  logP = 0;  
  n = length(Lam_x);
  x = xa(1:n);
  a = xa((n+1):end);

  % EP posterior
  logP = logP - 0.5 * a' * Lam_a * a + a'*eta_a;
  
  % cavity
  logP = logP - 0.5 * x' * Lam_x * x + x'*eta_x;
  
  % intervention likelihood
  logP = logP - 0.5 * (u-a'*x)^2/sigma_noise^2;

  % gradient and hessian
  negLogP = -logP;
  g = - LLgrad(x, a, u, Lam_a, eta_a, Lam_x, eta_x, sigma_noise);
  if nargout>2
    H = - LLHess(x, a, u, Lam_a, Lam_x, sigma_noise);
  end
  
end

function g = LLgrad(x, a, u, Lam_a, eta_a, Lam_x, eta_x, sigma_noise)
  gX = - Lam_x * x + eta_x - a*a'*x/sigma_noise^2 + a*u/sigma_noise^2;
  gA = - Lam_a * a + eta_a - x*x'*a/sigma_noise^2 + x*u/sigma_noise^2;
  g = [gX; gA];  
end

function H = LLHess(x, a, u, Lam_a, Lam_x, sigma_noise)
  n = numel(x);
  
  H_x = - Lam_x - a*a'/sigma_noise^2;
  H_a = - Lam_a - x*x'/sigma_noise^2;
  H_xa = - a*x'/sigma_noise^2;
  H_xa_diag = u/sigma_noise^2 - a'*x/sigma_noise^2;
  H_xa = H_xa + H_xa_diag * eye(n);
  
  % build hessian blockwise
  H = blkdiag(H_x, H_a);
  H(1:n,(n+1):end) = H_xa;
  H((n+1):end,1:n) = H_xa';
end


function [negLogP, g, H] = LLobjX(x, u, Lam_a, eta_a, Lam_x, eta_x, sigma_noise)
  logP = 0;  

  % EP posterior
  a = (Lam_a + x*x'./sigma_noise^2)\(eta_a + x*u/sigma_noise^2);
  logP = logP - 0.5 * a' * Lam_a * a + a'*eta_a;
  
  % cavity
  logP = logP - 0.5 * x' * Lam_x * x + x'*eta_x;
  
  % intervention likelihood
  logP = logP - 0.5 * (u-a'*x)^2/sigma_noise^2;
  
  % gradient and hessian
  negLogP = -logP;
  g = - LLgradX(x, a, u, Lam_x, eta_x, sigma_noise);
  H = - LLHessX(a, Lam_x, sigma_noise);
  
end


function gradX = LLgradX(x, a, u, Lam_x, eta_x, sigma_noise)
  gradX = - Lam_x * x + eta_x - a*a'*x/sigma_noise^2 + a*u/sigma_noise^2;
end

function H = LLHessX(a, Lam_x, sigma_noise)    
  H = - Lam_x - a*a'/sigma_noise^2;
end

