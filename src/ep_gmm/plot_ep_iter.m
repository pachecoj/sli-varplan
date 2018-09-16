function plot_ep_iter( ...
  h_fig, ...      % Figure handle
  pdf_true, ...   % True posterior (discrete)
  pdf_est, ...    % Estimated posterior (discrete)
  pdf_gauss, ...  % Best single Gaussian posterior (discrete)
  msg_pdf, ...    % Message PDFs, (ith column is ith message)
  err, ...        % Vector of per-iteration errors
  y, ...          % Data vector
  x_vals ...      % Domain values of discrete PDFs
  )
% UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

  N = numel(y);

  %
  % PLOT ESTIMATES
  %
  figure( h_fig );
  subplot( 2, 1, 1 ), cla;
  hold on;    
    plot(x_vals, pdf_true, '-b', ...
      x_vals, pdf_est, '-r', ...
      x_vals, pdf_gauss, '-k', 'LineWidth', 2);
    plot( repmat(x_vals',1,N), real( msg_pdf ), '--k');
%     plot( x, 0.1, 'xk' );
    plot(y, 0.1, '.b');  
%     ylim([0 0.6]);
    xlim([-50, 50]);
    legend('True Posterior','EP','Best Gaussian','Messages');
  hold off;  

  %
  % PLOT ERROR
  %
  figure( h_fig );
  subplot( 2, 1, 2 ), cla;
  semilogy(err, '-k', 'LineWidth', 3);
  grid on;
  ylim([1.0e-03 2]);
%   xlim([0, numel(err)]);
  title('Iteration vs. L1 Error');
  drawnow;
    
end

