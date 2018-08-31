function G = generate_network(N,model)
% return Graph object with approximately N nodes
% C                 = Conncetivity matrix ij

% kappa             = half-max value ij, kapa ii = d_i max decay rate
% n                 = Hill coefficient ij
% Vd                = degradation factor
% Vs                = Basal rate of expression
% AA                = Max Overexpression

% Glin              = linearized system matrix

% names, namesex    = (ext) names of the genes

G = generate_rand(N,model);
G = sample_nonlinear_dynamics(G);


function G = generate_rand(N,model);
% sample connections, undirected graph
% average degree 3
switch lower(model)
    case 'smallworld' % Watts-Strogatz
        % first connect the nearest neighors
        A = full(spdiags([rand(N,4) > 0.5],[-2 -1 1 2],N,N));
        for i = 1:round(0.3*N)
            for j = 1:10
                p = floor(rand(2,1)*N) + 1;
                if A(p(1),p(2)) ~= 0 %| p(1) >= N | G(p(1),p(1)+1) == 0
                    continue;
                end
                A(p(1),p(2)) = 1;
                %A(p(1),p(1)+1) = 0;
                break;
            end
        end
    otherwise
        error('Unknown model');
end

% sample wether edges are activating or repressing
ii = find(A(:)~=0);
A = double(A);
A(ii) = sign(rand(size(ii))-0.5);
A = A - diag(diag(A));

% Fill graph structure
G.C = A;
G.names = cell(N,1);
G.namesex = cell(N,1);
for i = 1:N
    G.names{i} = rand_str(4);
    G.namesex{i} = sprintf('%s-%i',G.names{i},i);
end


function G = sample_nonlinear_dynamics(G)
N = size(G.C,1);
r = 0;
while true
    r = r+1;
    if r > 100,
        error('Could not find random szstem with stable linearisation.');
    end

    G.kappa = zeros(N);
    G.n = zeros(N);
    ii = find(G.C(:)~=0);

    G.AA = zeros(N);
    G.AA(ii) = round(2 + rand(size(ii)) * 3);
    G.kappa(ii) = 0.2 + rand(size(ii)) * 1.5;
    G.kappa(sub2ind([N N],1:N,1:N)) = 20 + rand(N,1) * 50;
    G.n(ii) = round(1 + rand(size(ii)) * 1);

    G.Vs = 3 + rand(N,1)*2;
    G.Vd = round(150 + rand(N,1)*350);

    G.Vs = G.Vs / 10;
    G.Vd = G.Vd / 10;

    G.lin = [];
    [G.lin,G.x_steady] = linearize(G);
    if ~validate_deriv(G),
        fprintf('System derivative could not be validated (sample new).\n');
        continue;
    end;

    hh = max(real(eig(G.lin)));
    if hh > -0.03,
        fprintf('System not stable (sample new).\n');
        continue;
    end

    if max(G.x_steady) > 10,
        fprintf('System out of range (sample new).\n');
        continue;
    end
    break;
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

function fu = myfun(x)
global superG superp;
if ~all(x > 0),
    fu = -sum(x(x<0)) * 1e5;
else
    f = nonlinear_help(x,superp,superG);
    fu = f'*f;
end

function res = validate_deriv(G)
global superG superp;
superG = G;
superp = nonlinear_prepare(G);
N = size(G.C,1);

scale = 1e-5;
df = zeros(N);
for i = 1:N
    xh = zeros(N,1); xh(i) = scale;
    df(:,i) = nonlinear_help(G.x_steady+xh,superp,superG) - ...
        nonlinear_help(G.x_steady,superp,superG);
end
d = df - G.lin*scale;
res = sum(abs(d(:)) > scale*1e-2);
res = res == 0;

function [lin,xsteady] = linearize(G)
% compute linearisation around the steady state
N = size(G.C,1);

% compute steady state
global superG superp;
superG = G;
superp = nonlinear_prepare(G);
xsteady = simulate(ones(N,1)*0.2,G,0,0,zeros(N,1),[1e4 1e-2]);
% opt = optimset('Display','iter','TolX',1e-4,'Jacobian','off',...
%     'LargeScale','off','TolFun',1e-8,'TolCon',1e-5);
% xsteady = fmincon(@myfun,xsteady,[],[],[],[],2e-4*ones(N,1),[],[],opt);
lin = nonlinear_deriv(xsteady,superp,G);

function lin = nonlinear_deriv(x,p,G)
N = size(G.C,1);
lin = zeros(N);
f = nonlinear_help(x,p,G);
f = f + G.Vd .* x ./ (diag(G.kappa) + x);
for i = 1:N
    Ai = p.Ai{i}+1; Ii = p.Ii{i}+1; AiIi = [Ai Ii];
    lin(i,i) = -G.Vd(i) * (1. / (G.kappa(i,i) + x(i)) - x(i) / (G.kappa(i,i) + x(i))^2);
    for j = 1:length(Ai)
        lin(i,Ai(j)) = lin(i,Ai(j)) + f(i) / (1 + G.AA(i,Ai(j))*(x(Ai(j))/G.kappa(i,Ai(j)))^G.n(i,Ai(j))) * ...
            G.AA(i,Ai(j))*G.n(i,Ai(j))/G.kappa(i,Ai(j)) * (x(Ai(j))/G.kappa(i,Ai(j)))^(G.n(i,Ai(j))-1);
    end
    for j = 1:length(AiIi)
        lin(i,AiIi(j)) = lin(i,AiIi(j)) - f(i) / (1 + (x(AiIi(j))/G.kappa(i,AiIi(j)))^G.n(i,AiIi(j))) * ...
            G.n(i,AiIi(j))/G.kappa(i,AiIi(j)) * (x(AiIi(j))/G.kappa(i,AiIi(j)))^(G.n(i,AiIi(j))-1);
    end
end



