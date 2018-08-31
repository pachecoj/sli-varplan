function run_opt = get_run_opt(data_opt,num_inclus,method,initial_rand,tau,eps,sigma_method)
%  method means:
%  * - 1: Information gain, no EP update: Includes (x,u), but does not modify
%  *      EP site parameters, which remain at values before inclusion
%  * - 3: No score is used. u_* is randomly picked from the rem. candidates
%  *
%  * - 4: Gaussian prior with type 1 selection
%  * - 5: Gaussian prior with type 3 selection
    

if nargin ~=7,
    error('Wrong number of arguments to get_run_opt.');
end
run_opt.data_opt = data_opt;
run_opt.num_inclus = num_inclus;
run_opt.method = method;
run_opt.initial_rand = initial_rand;
run_opt.tau = tau;
run_opt.eps = eps;
run_opt.sigma_method = sigma_method;
