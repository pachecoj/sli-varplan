function ind = VarSelect_Designed_Experiment2(i,X,U,sitepi,siteb,L,gamma,candU,sigma_noise,order,num_initial_rand)
  if i <= num_initial_rand
    ind = Select_Random_Experiment(i,X,U,sitepi,siteb,L,gamma,candU,sigma_noise,order);
  else
    n = size(X,1);

    % find all that candU that are not used yet
    poss = true(size(candU,2),1);
    %     poss(order(1:i-1)) = false;
    poss = find(poss);
    poss_candU = candU(:,poss);

    % precompute EP posterior moments
    m = cell(n,1);
    Lam = cell(n,1);
    V = cell(n,1);
    for j = 1:n
      m{j} = (L{j}'\gamma{j});
      Lam{j} = L{j}*L{j}'/sigma_noise^2;  
      V{j} = inv( Lam{j} );
    end

    % compute score for each U candidate
    scores = zeros(size(poss_candU,2),1);
    [xSamp, aSamp] = eplin_sampxcand(X',U',sigma_noise^2,1,sitepi,siteb,L,gamma,poss_candU);
    aSampT = aSamp';
    for ui = 1:size(poss_candU,2)
      u = poss_candU(:,ui);
      t_start = tic();

      %
      % Optimize joint on (X,A)
      %
      
      % optimize local approximation
      xa_0 = [ xSamp(:,ui); aSampT(:) ];
      opt = optimoptions('fminunc','Display','off','GradObj','on','Hessian','on','DerivativeCheck','off',...
        'FinDiffType','central','Algorithm','trust-region'); % ,'TolFun',eps,'MaxIter',5000,'MaxFunEvals',5000);
      [xaOpt,fval,~,~,g,H] = fminunc( @(X) LLobj(X, u, m, Lam, sigma_noise), xa_0, opt);
      xOpt = xaOpt(1:n);
      AOpt = reshape( xaOpt((n+1):end), [n, n])';
      
      % compute Hessian and project to PSD / symmetric
      LamOpt = - LLhessXA(xOpt, AOpt, u, m, Lam, sigma_noise);
      LamOpt = psdproj(LamOpt, []);

      % compute Gaussian MI
      %     scores(ui) = 0.5 * log( sum( eig(VOpt) ) ) - 0.5 * log( sum( eig(Vx) ) ); %  - log( sum( eig(Va) ) );
      scores(ui) = 0;
      for j=1:n
        idx_start = 1 + n*j;
        idx_end = idx_start + n - 1;
        Lam_a_cond_x = LamOpt( idx_start:idx_end, idx_start:idx_end );
        scores(ui) = scores(ui) - 0.5 * sum( log( eig( Lam{j} ) ) ) ...
          + 0.5 * sum( log( eig( Lam_a_cond_x ) ) );
      end

      t_stop = toc(t_start);
      fprintf('%d (%0.1fs)  ', ui,  t_stop);
    end
    [~,ind] = max( scores );
    ind = poss( ind );


    %   %% DEBUG: Plot EP / ADF estimates of MI
    %   figure('InvertHardcopy','off','Color',[1 1 1]);
    %   set(gca,'FontSize',14);
    %   set(gcf,'Position',[230   202   990   425]);
    %   hold on
    %   yval = NaN(size(candU,2),1);
    %   yval( poss ) = scores; % / n;
    %   h_vp = plot( 1:numel(yval), yval, '-k', 'LineWidth', 2);
    %   plot( 1:numel(yval), yval, '*b')
    %   plot( ind, scores(ind), '*r' );
    %   xlabel('Intervention');
    %   ylabel('Mutual Information');
    %
    %   %% DEBUG: Plot sample-based MI estimates
    %   scoresLV = Debug_Select_Designed_Experiment(...
    %     10,500,i,X,U,sitepi,siteb,L,gamma,candU,sigma_noise,order,num_initial_rand);
    %   scoresLV = scoresLV .* n;
    %   yval = NaN(size(candU,2),1);
    %   yval( poss ) = mean( scoresLV, 2 );
    %   ystd = NaN(size(candU,2),1);
    %   ystd( poss ) = std( scoresLV, 0, 2 );
    %   h_mcmc = errorbar( 1:numel(yval), yval, ystd, '-s', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red');
    %
    %   %% DEBUG: Save plot
    %   title(sprintf('MI Estimates (Experiment %d)',size(X,2)));
    %   legend([h_vp, h_mcmc], {'Varplan','MCMC'});
    %   export_fig('~/Research/journal/010718/MI_est.pdf','-append');
    %   close;
  end
end

  

% function [negLogP, g, H] = LLobj2(x, u, m, Lam, sigma_noise)
%   n = length(m);
%   
%   % EP posterior  
%   logP = 0;
%   A = zeros(n,n);
%   for i = 1:n
%     h = Lam{i} * m{i};
%     a = ( Lam{i} + x*x'/sigma_noise^2 ) \ ( h + x*u(i)/sigma_noise^2);
%     logP = logP - 0.5 * (a - m{i})'*Lam{i}*(a - m{i}); % + -n/2 * log(2*pi) + 0.5 * log(det(Lam{i})) ...      
%     A(i,:) = a;
%   end
%   
%   % data term
%   logP = logP - 0.5/sigma_noise/sigma_noise * (u-A*x)' * (u-A*x); % -n/2 * log(2*pi*sigma_noise*sigma_noise) ...    
%   negLogP = - logP;
%   
%   % gradient
%   g = - LLgrad2(x, A, u, m, Lam, sigma_noise);
%   H = 0;% - LLhess2(x, A, u, m, Lam, sigma_noise);    
% end


  
function [negLogP, g, H] = LLobj(XA, u, m, Lam, sigma_noise)
  n = length(m);
  
  % unpack components
  x = XA(1:n);
  A = reshape( XA((n+1):end), [n,n] )';
    
  % EP posterior  
  logP = 0;
  for i = 1:n
    a = A(i,:)';    
    logP = logP + -n/2 * log(2*pi) + 0.5 * log(det(Lam{i})) ...
      - 0.5 * (a - m{i})'*Lam{i}*(a - m{i});
    logP = logP - 1/2 * log(2*pi*sigma_noise*sigma_noise) ...
      - 0.5 * (u(i)-a'*x)^2/sigma_noise/sigma_noise;
  end
  negLogP = - logP;
  
  % data term
%   logP = logP + -n/2 * log(2*pi*sigma_noise*sigma_noise) ...
%     - 0.5/sigma_noise/sigma_noise * (u-A*x)'*(u-A*x);  
  
  % gradient
  if nargout>1
    g = - LLgradXA(x, A, u, m, Lam, sigma_noise);
  end
  if nargout>2    
    H = - LLhessXA(x, A, u, m, Lam, sigma_noise);
  end
end

function H = LLhessSparsity(n)
  Hxx = ones(n,n);
  Haa = zeros(n^2,n^2);
  Hxa = ones(n,n^2);
    
  for i=1:n    
    idx_start = 1 + n*(i-1);
    idx_end = idx_start + n - 1;    
    Haa(idx_start:idx_end,idx_start:idx_end) = 1;    
  end
  
  % build hessian blockwise
  H = blkdiag(Hxx, Haa);
  H(1:n,(n+1):end) = Hxa;
  H((n+1):end,1:n) = Hxa';  
end