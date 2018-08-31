function ind = VarSelect_Designed_Experiment3(i,X,U,sitepi,siteb,L,gamma,candU,sigma_noise,order,num_initial_rand)
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
      eta_a{j} = L{j}*gamma{j};
      Lam_a{j} = L{j}*L{j}'/sigma_noise^2;
    end
    
    % compute score for each U candidate
    scores = zeros(size(poss_candU,2),1);
    for ui = 1:size(poss_candU,2)
      fprintf('Scoring candidate %i of %i\n', ui, size(poss_candU,2));
      u = poss_candU(:,ui);
      
      % run EP/ADF
      max_iters = 10;
      conv_thresh = 1e-4;
      [Lam_a_cond_x, Lam_x] = EP_VarPlan( ...
        eta_a, Lam_a, u, sigma_noise, max_iters, conv_thresh );
      
      % compute conditioinal entropy
%       H_x = sum( log( eig( Lam_x ) ) );
      for j=1:n
        scores(ui) = scores(ui) - 0.5 * sum( log( eig( Lam_a{j} ) ) ) + n/2 * log(2*pi*exp(1)) ...
          + 0.5 * sum( log( eig( Lam_a_cond_x{j} ) ) ) - n/2 * log(2*pi*exp(1));
      end
            
    end
    
    [~,ind] = max( scores );
    ind = poss( ind );
    
    %% DEBUG:
    figure('InvertHardcopy','off','Color',[1 1 1]);
    set(gca,'FontSize',14);
    set(gcf,'Position',[230   202   990   425]);
    yval = NaN(size(candU,2),1);
    yval( poss ) = scores; % / n;
    plot( 1:numel(yval), yval, 'LineWidth', 2)
    hold on
    plot( 1:numel(yval), yval, '*b')
    plot( ind, yval(ind), '*r' );
    xlabel('Intervention');
    ylabel('Mutual Information');    
    title(sprintf('Variational Estimate (Experiment %d)',size(X,2)));
%     export_fig('~/Research/journal/102917/var_MI.pdf','-append');
%     close;
    
    %% DEBUG
    scoresLV = Debug_Select_Designed_Experiment(...
      10,500,i,X,U,sitepi,siteb,L,gamma,candU,sigma_noise,order,num_initial_rand);
    scoresLV = scoresLV .* n;    
    yval = NaN(size(candU,2),1);      
    yval( poss ) = mean( scoresLV, 2 );         
    ystd = NaN(size(candU,2),1);
    ystd( poss ) = std( scoresLV, 0, 2 );    
    errorbar( 1:numel(yval), yval, ystd, '-s', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red')
    xlim([0 (numel(yval)+1)])
    xlabel('Intervention')
    ylabel('Mutual Information')
    title(sprintf('MC Estimate (Experiment %d)',size(X,2)));
    export_fig('~/Research/journal/102917/MI_est.pdf','-append');
    close;
  end
end


function [ Lam_a_cond_x, Lam_x, Lam_fac ] = EP_VarPlan( eta_a, Lam_a, u, sigma_noise, max_iters, conv_thresh )
  
  n = numel( Lam_a );
  min_eigval = 1e-2;
  
  % optimization settings
  opt = optimoptions('fminunc','Display','off','GradObj','on','Hessian','off','DerivativeCheck','off',...
      'FinDiffType','central','Algorithm','quasi-newton'); % ,'TolFun',eps); % ,'MaxIter',10,'MaxFunEvals',10);
    
  % declare stuff
  Lam_msg = cell(n,1);
  eta_msg = cell(n,1);
  Lam_x = zeros(n,n);
  eta_x = zeros(n,1);
  Lam_fac = cell(n,1);
  eta_fac = cell(n,1);   
  Lam_a_cond_x = cell(n,1);
  
  % initialize messages/marginal on X
  for i=1:n
    Lam_msg{i} = 1/n * eye(n);
    eta_msg{i} = zeros(n,1);
    Lam_x = Lam_x + Lam_msg{i};
    eta_x = eta_x + eta_msg{i};
  end
    
  % do EP/ADF  
  Aopt = zeros( n, n );
  Amean = zeros( n, n );
  Xopt = zeros( n, n );
  Xmean = zeros( n, n );
  for iter=1:max_iters
    if (max_iters>1) fprintf('EP Planning Iter %i of %i\n',iter,max_iters); end
    Lam_fac_old = Lam_fac;
    eta_fac_old = eta_fac;    
%     xa_0 = zeros(2*n,1);
    for i=1:n
      
      % cavity "distribution" parameters
      Lam_cav = Lam_x - Lam_msg{i};
      eta_cav = eta_x - eta_msg{i};
      Lam_cav = psdproj( Lam_cav, min_eigval );
                  
      %
      % JOINT OPTIMIZATION OVER (X,A)
      %      
      
      % maximize log-augmented distribution      
      m_a = Lam_a{i} \ eta_a{i};
      m_x_cav = Lam_cav \ eta_cav;
      if u(i)
        xa_0 = [ m_x_cav; m_a ];        
        %% DEBUG:
        xa_0 = rand( size( xa_0 ) );
        xaOpt = fminunc( ...
          @(X) LLobj(X, u(i), Lam_a{i}, eta_a{i}, Lam_cav, eta_cav, sigma_noise), xa_0, opt);        
      else
        xaOpt = [ m_x_cav; zeros(n,1) ];
      end
      xOpt = xaOpt(1:n);
      aOpt = xaOpt((n+1):end);
      
      %% DEBUG
      Xopt(i,:) = xOpt;
      Xmean(i,:) = m_x_cav;
      Aopt(i,:) = aOpt;      
      Amean(i,:) = m_a;      
            
      % compute Hessian and project to PSD / symmetric matrix
      H = - LLHess(xOpt, aOpt, u(i), Lam_a{i}, Lam_cav, sigma_noise);
      H = psdproj( H, min_eigval );      
      invH = inv( H );
      r = rcond(H);
      if ( r < eps )
        fprintf('WARNING: Hessian matrix poorly conditioned (rcond=%0.1e)\n', r);
      end
                  
      % update factor approx. using Laplace approx.
      Lam_xa = H; % L*L';
      eta_xa = Lam_xa*xaOpt;
      Lam_fac{i} = Lam_xa - blkdiag( Lam_cav, Lam_a{i} );
      eta_fac{i} = eta_xa - [ eta_cav; eta_a{i} ];      
      Lam_a_cond_x{i} = Lam_xa((n+1):end,(n+1):end);
                        
      % marginal on X
      V_xa = invH; % invL' * invL;
      V_x = V_xa( 1:n, 1:n );
      Lam_x = inv( V_x );
      m_xa = V_xa * eta_xa;
      eta_x = V_x \ m_xa(1:n);
      
      % update message from i^th factor
      Lam_msg{i} = Lam_x - Lam_cav;
      eta_msg{i} = eta_x - eta_cav;
                            
    end        
  end
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

