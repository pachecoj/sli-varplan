function [prior, trans, like] = eval_model(x_grid, loc, y_all, mu_0, sig_0, sig_u, alpha)
% EVAL_MODEL - Evaluate model probabilities.
%
% J. Pacheco, 2017
%
  T = size(y_all,1);
  Ngrid = numel(x_grid);
  N_sensors = numel( loc );

  % compute prior, likelihood and transition probabilities  
  trans = zeros(Ngrid,Ngrid,T);
  like = cell(1,T);
  for t=1:T
    for j=1:Ngrid
      mu_x = x_grid(j)/2 + 25*x_grid(j)/(1+x_grid(j)^2) + 8*cos(1.2*t);
      trans(j,:,t) = -0.5*log(2*pi*sig_u*sig_u) - 0.5*(x_grid-mu_x).^2./sig_u./sig_u;
    end
    
    % compute likelihood
    like{t} = NaN( Ngrid, N_sensors );
    for l=1:N_sensors
      sig_y = sqrt( alpha * ( loc(l) - x_grid + eps).^2 );
      mu_y = x_grid.^2./20;
      like{t}(:,l) = -0.5*log(2*pi.*sig_y.*sig_y) - 0.5*(y_all(t,l) - mu_y).^2./sig_y./sig_y;
%       like{t} = like{t}';
    end
  end
  prior = - 0.5 * log( 2*pi*sig_0*sig_0 ) - 0.5 * (x_grid' - mu_0).^2./sig_0./sig_0;
end