function [x,y_all] = sample_gmm(gmm, N)
% SAMPLE_GMM - Samples N observations from a two-component Gaussian Mixture
%   Model parameterized by gmm as,
%
% p(x) = N(x; m0_x, V0_x)
% p(y|x) \propto w * N(y;Ax+a,V_0) + (1-w) N(y;Bx+b,V_1)
%
% Jason L. Pacheco
% 12/14/11
%

  x = gmm.m0_x + sqrt(gmm.V0_x) * randn();  
  pi_vec = rand(N,1);
  I_A = find(pi_vec <= gmm.w);
  I_B = find(pi_vec > gmm.w);
  y_all = zeros(numel(gmm.loc),N);  
  
  % generate observations
  for k=1:numel(gmm.loc)
    mu = gmm.loc(k);
    dist = abs( mu - x );
    var_tgt = gmm.c * dist + gmm.V_1;
    y_all(k,I_A) = (gmm.A*x + gmm.a) + sqrt(gmm.V_0) * randn(1,numel(I_A));
    y_all(k,I_B) = (gmm.B*x + gmm.b) + sqrt(var_tgt) * randn(1,numel(I_B));
  end
  
end