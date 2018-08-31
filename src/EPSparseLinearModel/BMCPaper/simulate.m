function [x,traj] = simulate(x0,G,linear,sigma_noise,u,time)
% function [x,traj] = simulate(x0,G,linear,sigma_noise,u,time)
% simulates system until it reaches steady state
%
% linear = 0 (non linear) / 1 (linearized)

if nargin < 6
    if isempty(G.lin)
        lambda = 1;
    else
        lambda = abs(real(eig(G.lin)));
    end
    T_long  = 1/min(lambda) * 10;
    T_short = 1/max(lambda) / 10;
    nstep = min(max(round(T_long/T_short),5e3),5e4);
    dt = T_long/nstep;
else
    nstep = time(1);
    dt = time(2);
end

N = size(G.C,1);
if ~linear
    noise = randn(N,nstep);
    p = nonlinear_prepare(G);
    [x,traj] = nonlinearsimulate(x0,p,G,nstep,noise,dt,u,sigma_noise);
else
    x = x0;
    traj = zeros(N,nstep);
    for j = 1:nstep
        f = G.lin * x;
        x = x + (f + u) * dt + sigma_noise*randn(N,1) * sqrt(dt);
        traj(:,j) = x;
    end
end


function p = nonlinear_prepare(G)
N = size(G.C,1);
p.Ai = cell(N,1);
p.Ii = cell(N,1);
p.AiIi = cell(N,1);
for i = 1:N
    p.Ai{i} = int32(find(G.C(i,:)>0) -1); %comvert to c standard
    p.Ii{i} = int32(find(G.C(i,:)<0) -1);
    p.AiIi{i} = [p.Ai{i} p.Ii{i}];
end

