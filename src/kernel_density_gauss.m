function p = kernel_density_gauss( x_grid, x, w, sig )
% KERNEL_DENSITY_GAUSS - Gaussian (weighted) kernel density estimate.  
%
% INPUTS:
%
%   x_grid - Equally spaced linear array of x values to compute density on.
%
%   x - Unequally dispersed "particle" locations (e.g. kernel locations)
%
%   w - Weights
%
%   sig - Standard deviation (bandwidth parameter)
%
% OUTPUTS:
%
%   pdf - Numerical pdf
%
% J. Pacheco, 2013
%
  T = size(x,2);
  p = zeros([numel(x_grid),T]);
  dx = x_grid(2) - x_grid(1);
  for t=1:T
%     for i=1:size(x,1)
%       thisX = x(i,t);
%       thisPdf = w(i,t) * normpdf( x_grid, thisX, sig );
%       p(:,t) = p(:,t) + thisPdf';
%     end    
    p(:,t) = ksdensity( x(:,t), x_grid, 'weights', w(:,t), 'width', sig );
    p(:,t) = p(:,t) * dx;
  end  
end