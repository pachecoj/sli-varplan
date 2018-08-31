function ind = VarSelect_Designed_Experiment4(i,X,U,sitepi,siteb,L,gamma,candU,sigma_noise,order,num_initial_rand)
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
    mCell = cell(n,1);
    LamCell = cell(n,1);
    LamMat = zeros(n,n,n);
    VMat = zeros(n,n,n);
    mMat = zeros(n,n);
    for j = 1:n
      mCell{j} = (L{j}'\gamma{j});
      mMat(:,j) = mCell{j};
      LamCell{j} = L{j}*L{j}'/sigma_noise^2;        
      LamMat(:,:,j) = LamCell{j};
      VMat(:,:,j) = inv( LamCell{j} );
    end    
    
    % compute score for each U candidate
    scores = zeros(size(poss_candU,2),1);
    [xSamp, aSamp] = eplin_sampxcand(X',U',sigma_noise^2,1,sitepi,siteb,L,gamma,poss_candU);
    for ui = 1:size(poss_candU,2)
      u = poss_candU(:,ui);
      t_start = tic();

      %
      % Optimize marginally on X
      %
      
%       x = adigatorCreateDerivInput( [ n, 1 ], 'x');
%       u = adigatorCreateAuxInput( size( u ) );
%       m = adigatorCreateAuxInput( size( mMat ) );
%       Lam = adigatorCreateAuxInput( [ numel( LamMat ), 1 ] );
%       V = adigatorCreateAuxInput( [ numel( VMat ), 1 ] );

      % optimize local approximation
      x_0 = xSamp(:,ui);
%       opt = optimoptions('fminunc','Display','off','GradObj','on','Hessian','off','DerivativeCheck','off',...
%         'FinDiffType','central','Algorithm','quasi-newton','MaxIter',Inf,'MaxFunEvals',Inf,'TolFun',1e-2);
%       [xOpt,fval,~,~,g,H] = fminunc( @(X) LLobjADwrapper(X, u, mMat, LamMat(:), VMat(:), sigma_noise), x_0, opt);
      xOpt = x_0;

      % compute variance
      AOpt = optA( xOpt, u, mCell, LamCell, sigma_noise );
      LamOpt = - LLhessXA(xOpt, AOpt, u, mCell, LamCell, sigma_noise);
      LamOpt = psdproj(LamOpt, []);
      g = LLgradXA(xOpt, AOpt, u, mCell, LamCell, sigma_noise);      

      % compute Gaussian MI
      %     scores(ui) = 0.5 * log( sum( eig(VOpt) ) ) - 0.5 * log( sum( eig(Vx) ) ); %  - log( sum( eig(Va) ) );
      scores(ui) = 0;
      for j=1:n
        idx_start = 1 + n*j;
        idx_end = idx_start + n - 1;
        Lam_a_cond_x = LamOpt( idx_start:idx_end, idx_start:idx_end );
        scores(ui) = scores(ui) - 0.5 * sum( log( eig( LamCell{j} ) ) ) ...
          + 0.5 * sum( log( eig( Lam_a_cond_x ) ) );
      end

      t_stop = toc(t_start);
      fprintf('%d (%0.1fs, g: %0.3f)  ', ui,  t_stop, norm(g));
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

function A = optA( x, u, m, Lam, sigma_noise )
  n = numel(x);
  A = zeros(n,n);
  for i = 1:n
    thisLam = Lam{i};
    this_m = m{i};    
    this_h = thisLam * this_m;    
%     invLamXX = V{i} - (1/sigma_noise^2) * V{i} * (x * x') * V{i} / (1 + x' * V{i} * x / sigma_noise^2 );
    a = ( thisLam + (x * x') / sigma_noise^2 ) \ ( this_h + x*u(i)/sigma_noise^2);
    A(i,:) = a;
  end
end