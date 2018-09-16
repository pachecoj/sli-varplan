function opt = init_ep_opt( ...
  sched, ...          % Message update schedule
  max_iters, ...      % Maximum # of iterations
  conv_thresh, ...    % Convergence threshold
  damp_type, ...      % Type of update damping (0=None, 1=Static, 2=Adaptive)
  min_damp_fact, ...  % Damp factor if damp_type=1, minimum damp if damp_type=2
  plot_flag ...       % Boolean, show plot
  )
% INIT_EP_OPT - Initializes options for EP inference.
%
% Jason L. Pacheco
% 12/14/11
%

  opt = struct('sched', sched, 'max_iters', max_iters, 'conv_thresh', conv_thresh, ...
    'damp_type', damp_type, 'min_damp_fact', min_damp_fact, 'plot_flag', plot_flag);

end

