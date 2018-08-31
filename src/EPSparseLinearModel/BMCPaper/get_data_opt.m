function res = get_data_opt(name,N,type_input,num_input,num_nonzeros_input,num_initial,input_scale,linear,SDE,sigma_noise,net)
% prepare a structure containing all the necessarry information about one
% experiment

if nargin ~= 11,
    error('Wrong number of parameters');
end

res.name                = name;
res.N                   = N;
res.type_input          = type_input;
res.num_input           = num_input;
res.num_nonzeros_input  = num_nonzeros_input;
res.num_initial         = num_initial;
res.linear              = linear;
% res.SDE                 = SDE;
res.sigma_noise         = sigma_noise;
res.net                 = net;
res.input_scale         = input_scale;
