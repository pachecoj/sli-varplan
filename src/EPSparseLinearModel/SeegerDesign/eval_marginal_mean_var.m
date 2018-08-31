function [Qm,Qv] = eval_marginal_mean_var(L,gamma,sigma_noise)
n = length(L);
Qm = zeros(n);
Qv = zeros(n);
for i = 1:n
  Qm(i,:) = (L{i}'\gamma{i})';
  Qv(i,:) = diag(L{i}'\(L{i}\eye(n))*sigma_noise^2)';
end