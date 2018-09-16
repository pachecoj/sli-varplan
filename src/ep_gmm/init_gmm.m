function model = init_gmm( A, B, a, b, w, c, V_0, V_1, m0_x, V0_x, Nsensors, sensorBounds )
% INIT_GMM - Constructs and returns a structure variable holding
%   parameters of a Gaussian Mixture Model.  The structure parameterizes
%   a two-component GMM of the form,
%
%   p(y|x) = w * N(y;Ax+a,V_0) + (1-w) N(y;Bx+b,V_1)
%
%   with a Gaussian prior on x as: 
%
%   p(x) = N(x; m0_x, V0_x)
%
% Jason L. Pacheco
% 12/14/11
%
  
  loc = linspace(sensorBounds(1), sensorBounds(end), Nsensors);
  model = struct( 'A', A, 'B', B, 'a', a, 'b', b, 'w', w, 'c', c, 'V_0', V_0, ...
    'V_1', V_1, 'm0_x', m0_x, 'V0_x', V0_x, 'loc', loc );

end

