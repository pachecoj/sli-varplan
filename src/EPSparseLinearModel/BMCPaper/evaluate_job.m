function [aucN,aucPR,times] = evaluate_job(dir,run_opt,rounds,eps_multiplier)
% function [aucN,times] = evaluate_job(dir,run_opt,rounds,eps_multiplier)
%
% Returns the iAUC and running times
%
% Note: The quality for the non-linear case should judged relative to those
% edges which are present in the linearised version of the network. We have
% no chance to discover edges which are inactive around the null input steady state.
% in aucN we thus return these values.
%
% WARNING: Uses caching for the results: So if any change is done, delete the
% cache!

if ~exist(dir,'dir'),
    error('Could not open directory.');
end
file = sprintf('%s/%s-evaluation.mat',dir,get_run_opt_FN(run_opt));

% if exist(file,'file')
%     load(file,'aucN','times');
% else
    aucN = zeros(run_opt.num_inclus+1,rounds);
    aucPR = zeros(run_opt.num_inclus+1,rounds);
    ttimes = zeros(5,rounds);
    %     snapshot = struct;

    valid = true(rounds,1);
    for r = 1:rounds
        % load graph data *********************************
        load(sprintf('%s/%s-r=%i-data.mat',dir,run_opt.data_opt.name,r),'A','G','X','U','candX','candU');
        truth_eps = run_opt.eps*eps_multiplier;
        GG_act = A > truth_eps;
        GG_rep = A < -truth_eps;
        GG = GG_act | GG_rep;
        N = size(A,1);

        round_fn = sprintf('%s/%s-r=%i-',dir,get_run_opt_FN(run_opt),r);
        try
            load([round_fn 'results.mat'], 'Qm','Qv','times');
            % load easy stats *********************************
            %                 data.mm = reshape(loadbasematrix([round_fn 'mmeans.stm']),N,N,run_opt.num_inclus+1);
            %                 data.mv = reshape(loadbasematrix([round_fn 'mvars.stm']),N,N,run_opt.num_inclus+1);
            %                 data.mv = max(data.mv,1e-10*ones(size(data.mv)));
            %                 times(:,r) = loadbasevector([round_fn 'times.stv']);
        catch
            valid(r) = 0;
            continue;
        end

        for i = 0:run_opt.num_inclus
            Qv{i+1} = max(Qv{i+1},1e-10);
            rank = -logdiffgausscdf(reshape(Qm{i+1},N^2,1),reshape(sqrt(Qv{i+1}),N^2,1),run_opt.eps);
            rank = reshape(rank,N,N);
            [void,void,void,aucN(i+1,r),aucPR(i+1,r)] = roc(GG,rank);
        end
        ttimes(1:length(times),r) = times;

    end

    if ~any(valid),
        error(sprintf('Error evaluating: All runs crashed %s.',get_run_opt_FN(run_opt)));
    end

    % compute mean and variance of the stats
    aucN = prep(aucN,valid);
    aucPR = prep(aucPR,valid);
    times = prep(ttimes,valid);

    % save results in cache file
    save(file,'aucN','times');
% end

function a = prep(a,valid)
a = a(:,valid);
m = mean(a,2);
a = [m, m-min(a,[],2), max(a,[],2)-m]; % std(a,0,2)];
