function [ b_conv, delta ] = test_conv( lam, eta, lam_old, eta_old, thresh )
% TEST_CONV - Tests whether the message canonical parameters have 
%   converged.  The function returns true if the maximum absolute change
%   between the new parameters and old parameters is below thresh.
%
%   INPUTS:
%   lam - vector of message precisions (inverse variance)
%   
%   eta - vector of message scale parameters (inverse variance times mean)
%
%   lam_old, eta_old - vector of message parameters from previous iteration
%
%   thresh - convergence threshold
%
%   OUTPUTS:
%   b_conv - True if max( [ abs(lam - lam_old) , abs(eta - eta_old) ]) < thresh.
%            False otherwise.
%
%   delta - Max diff. in parameters
%
% Jason L. Pacheco
% 12/14/11
%

  delta_lam = abs( lam_old - lam );
  delta_eta = abs( eta_old - eta );
  delta = max( [ delta_lam, delta_eta ] );
  b_conv =  delta < thresh;

return
