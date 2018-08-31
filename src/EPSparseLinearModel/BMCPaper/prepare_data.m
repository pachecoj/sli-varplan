function prepare_data(dir,data_opt,rounds)
if ~exist(dir,'dir'), mkdir(dir), end;

N = data_opt.N;
for r = 1:rounds
    fprintf('Sampling data for round %i. ... \n',r);
    % load graph *********************************
    load(sprintf('data/networks_%s/N=%i_r=%i.mat',data_opt.net,data_opt.N,r));
    A = G.lin;

    % create data  *********************************
    
    % create unit l2-norm data
    switch data_opt.type_input
        case 'sparse'
            Uh = unit_vecs(N,data_opt.num_nonzeros_input,data_opt.num_input);
        case 'non-sparse'
            Uh = randn(N,data_opt.num_input);
        otherwise
            error('Unkown input type');
    end
    % random up- or down-regulation
    ii = find(Uh(:) ~= 0);
    Uh(ii) = sign(rand(size(ii))-0.5);

    % normalize and random permutation
    Uh = Uh * diag(1./sqrt(sum(Uh.^2)));
    Uh = Uh(:,randperm(data_opt.num_input));
    
    % add a couple of zero inputs to determine the steady state of the system
    Nsteadystate = 10;  Uh = [zeros(N,Nsteadystate) Uh];

    % for nonlinear networks scale input (do not leave linearity region)
    if ~data_opt.linear,
        Uh      = data_opt.input_scale * Uh;
    end

    % generate data either for linear or nonlinear system
    if ~data_opt.linear
        % create x values by running the SDE, use -Uh because U appears on the lefthand side in U = AX
        [Xh,traj] = sample_data(G,-Uh,data_opt.sigma_noise,data_opt.linear);
        
        % plot
        figure('InvertHardcopy','off','Color',[1 1 1]);
        set(gca,'FontSize',14);
        plot(traj');
        xlabel('Time');
        ylabel('Expression');
        title('Simulated Gene Expressions');        
        if any(isnan(Xh(:))),
            error('Simulation lead to inconsistent values.');
        end
    else
        Xh = A \ (Uh + randn(size(Uh))*data_opt.sigma_noise);
    end


    % assign X,U, Xcand,Ucand
    XN = Xh(:,1:Nsteadystate);
    U = Uh(:,Nsteadystate+(1:data_opt.num_initial));
    X = Xh(:,Nsteadystate+(1:data_opt.num_initial));
    candU = Uh(:,Nsteadystate+data_opt.num_initial+1:end);
    candX = Xh(:,Nsteadystate+data_opt.num_initial+1:end);

    
    % save data  *********************************
    fname = sprintf('%s/%s-r=%i-data.mat',dir,data_opt.name,r);
    save(fname,'A','U','X','candU','candX','XN','G');
    fprintf('Saved: %s\n', fname);
end

function [X,example_trajectory] = sample_data(G,U,sigma_noise,linear)
% computes X for given U pairs for a network G by simulating stochastic model with input
N = size(U,1);
nu = size(U,2);
X = zeros(N,nu);
example_trajectory = [];
for iu = 1:nu
    x0 = G.x_steady;
    
    % get random state from equilibrium distribution
%     [x,t1] = simulate(x0,G,mode,sigma_noise,zeros(N,1));

    % get point from disturbed equilibrium
    [nx,t2] = simulate(x0,G,linear,sigma_noise,U(:,iu));

    X(:,iu) = nx(1:N)-x0(1:N);
    if isempty(example_trajectory) && max(abs(U(:,iu))) > 0,
%         example_trajectory = [t1;t2];
        example_trajectory = t2;
    end
end

