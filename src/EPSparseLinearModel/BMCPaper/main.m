function main(action,exp_type,jobid)
% This is the main control script for all experiments.
% 
% The parameters are
%   action = 'perpare','run','eval'   Describes what should be done (see below)
%   exp_type                          Name of the experiments (see below)
%   jobid                             Job number within experiment.
% 
% 
%  For generating networks run main_generate_networks.m before calling this
%  function.

rng(0);

% some overhead
addpath(fileparts(pwd))
close all

if nargin < 2,    error('Not enough parameters.');end
if nargin < 3,    jobid = -1; end

labels = cell(100,1);
for i = 1:100,  labels{i} = sprintf('%i',i);end
props = {'',''};

% definition of experiements
switch exp_type       
    case 'design'
        rounds = 5; r = []; d = [];
        dir = 'data/dddesign';
        %         get_data_opt(name,N, type_input,num_input,nnz_input,num_initial,input_scale,linear,SDE,sigma_noise,  net)
        d = add(d,get_data_opt('N3n2',50,'sparse', 200,      3,         1,          1e-2,          0,     1,  1e-4,    'SW'));

        for i = 1:length(d)
            % get_run_opt(data_opt,num_inclus,method,initial_rand,tau,eps,sigma_method)
            r = add(r,get_run_opt(d(i), 50,      1,     0,         -1, 0.1,0.000184322));
            r = add(r,get_run_opt(d(i), 50,      1,    20,         -1, 0.1,0.000184322));
            r = add(r,get_run_opt(d(i), 50,      3,    0,         -1, 0.1,0.000184322));
            r = add(r,get_run_opt(d(i), 50,      4,     0,     3.1623, 0.1,0.000184322));
            r = add(r,get_run_opt(d(i), 50,      5,     0,     3.1623, 0.1,0.000184322));            
        end
        labels = {'LD','LM','LR','GD','GR'};
        
    case 'debug'
        rounds = 1; r = []; d = [];
        dir = 'data/debug';
        sigma_noise = 1e-4;
        
        % set data parameters
        %         get_data_opt(name,N, type_input,num_input,nnz_input,num_initial,input_scale,linear,SDE,sigma_noise,  net)        
        d = add(d,get_data_opt('N3n2',10,'sparse', 200,      3,         1,          1e-2,          0,     1,  sigma_noise,    'SW10'));

        for i = 1:length(d)
            % get_run_opt(data_opt,num_inclus,method,initial_rand,tau,eps,sigma_method)            
            r = add(r,get_run_opt(d(i), 20,      1,    1,         6.3096, 0.1, sigma_noise*1.2784));
            r = add(r,get_run_opt(d(i), 20,      3,    1,         6.3096, 0.1, sigma_noise*1.2784));
            r = add(r,get_run_opt(d(i), 20,      6,    1,         6.3096, 0.1, sigma_noise*1.2784));
        end
        labels = {'MCMC Design','Random','Variational Design'};
     
    case 'debug2'
        rounds = 5; r = []; d = [];
        dir = 'data/debug2';
        sigma_noise = 1e-4;
        %         get_data_opt(name,N, type_input,num_input,nnz_input,num_initial,input_scale,linear,SDE,sigma_noise,  net)
        d = add(d,get_data_opt('N3n2',50,'sparse', 200,      3,         1,          1e-2,          0,     1,  sigma_noise,    'SW'));

        for i = 1:length(d)
            % get_run_opt(data_opt,num_inclus,method,initial_rand,tau,eps,sigma_method)            
            r = add(r,get_run_opt(d(i), 50,      1,    20,         -1, 0.1, sigma_noise*1.84322));
            r = add(r,get_run_opt(d(i), 50,      3,    20,         -1, 0.1, sigma_noise*1.84322));
            r = add(r,get_run_opt(d(i), 50,      6,    20,         -1, 0.1, sigma_noise*1.84322));
        end
        labels = {'MCMC Design','Random','Variational Design'};

    case 'noise'
        rounds = 50; r = []; d = [];
        dir = 'data/ddnoise';

        noises = 10.^(-5:0.5:-2);
        for i = noises
            %         get_data_opt(name,           N, type_input,num_input,nnz_input,num_initial,input_scale,linear,SDE,sigma_noise,net)
            d = add(d,get_data_opt(sprintf('n%e',i),50,'sparse', 200,      3,         1,         1e-2,          0,     1,  i,    'SW'));
            if 1e-2/i < 1
                labels{length(d)} = sprintf('%.1g',1e-2/i); % SNR
            else
                labels{length(d)} = sprintf('%i',round(1e-2/i)); % SNR
            end
            if length(d) <= 5
                labels{length(d)} = sprintf('%s*',labels{length(d)});
            end
        end

        % determined a priori from networks of similar kind
        est_sigmas = [0.000123267, 0.000131226, 0.000184322, 0.000427546,...
            0.00127838, 0.00399288, 0.0123871, 0.0385018, 0.119991];
        taus = [-1, -1, -1, 12.5893, 6.3096, -1, -1, -1, -1];
        for i = 1:length(d)
            %   get_run_opt(data_opt,num_inclus,method,initial_rand,tau,eps,sigma_method)
            r = add(r,get_run_opt(d(i), 50,      1,     0,          taus(i), 0.1,est_sigmas(i)));
            r = add(r,get_run_opt(d(i), 50,      3,     0,          taus(i), 0.1,est_sigmas(i)));
        end
        props = {'Stochastic Noise','SNR'};

    case 'inputscale'
        rounds = 5; r = []; d = [];
        dir = 'data/ddinputscale';

        devs = [1e-3 5e-3 1e-2 5e-2 1e-1 2e-1 5e-1];
        for i = devs
            %         get_data_opt(name,           N, type_input,num_input,nnz_input,num_initial,input_scale,linear,SDE,sigma_noise,net)
            d = add(d,get_data_opt(sprintf('%e',i),50,'sparse', 200,      3,         1,          i,          0,     1,  1e-5,    'SW'));
            labels{length(d)} = sprintf('%g %%',i/1e-2);
            if length(d) <= 6
                labels{length(d)} = sprintf('%s*',labels{length(d)});
            end
        end

        est_sigmas = [1.29463e-05, 3.80607e-05, 0.000123267, 0.0021941,...
            0.00731039, 0.0224279, 0.0870597];
        taus = [-1, -1, -1, 12.5893, 6.3096, 6.3096, 6.3096];
        for i = 1:length(d)
            % get_run_opt(data_opt,num_inclus,method,initial_rand,tau,eps,sigma_method)
            r = add(r,get_run_opt(d(i), 50,      1,  0,           taus(i), 0.1,est_sigmas(i)));
            r = add(r,get_run_opt(d(i), 50,      3,  0,           taus(i), 0.1,est_sigmas(i)));
        end
        props = {'Pertubation Strength','change in steady state caused by perturbation'};


    case 'inputshape' % this should be the final version
        rounds = 50; r = []; d = [];
        dir = 'data/ddinputshape';

        shapes = [1 2 3 5 20];
        for i = shapes
            %         get_data_opt(name,           N, type_input,num_input,nnz_input,num_initial,input_scale,linear,SDE,sigma_noise,net)
            d = add(d,get_data_opt(sprintf('S%i',i),50,'sparse', 200,      i,         1,          1e-2,          0,     1, 1e-4,    'SW'));
            labels{length(d)} = sprintf('%i*',i);
        end
        %         get_data_opt(name,N, type_input,   num_input,nnz_input,num_initial,input_scale,linear,SDE,sigma_noise,net)
        d = add(d,get_data_opt('NS',50,'non-sparse', 200,      0,         1,          1e-2,          0,     1, 1e-4,    'SW'));
        labels{length(d)} = 'non-sparse*';

        for i = 1:length(d)
            %      get_run_opt(data_opt,num_inclus,method,initial_rand,tau,eps,sigma_method)
            r = add(r,get_run_opt(d(i),    50,        1,  0,           -1, 0.1,0.000184322));
            r = add(r,get_run_opt(d(i),    50,        3,  0,           -1, 0.1,0.000184322));
        end
        props = {'Type of Pertubations','Number of Pertubations per Experiment'};

    otherwise
        error('Unkown experiment.');
end

% definition of action that can be done with experiment
switch action
    case 'prepare'
        for i = 1:length(d)
            prepare_data(dir,d(i),rounds);
        end
    case 'preparedatasingle'
        prepare_data(dir,d(jobid),rounds);

    case 'run'
        for i = 1:length(r)
          run_job(dir,r(i),rounds);
        end
%      i = 3;
%      run_job(dir,r(i),rounds);

    case 'eval'
        if jobid == -1
            for i = 1:length(r)
                evaluate_job(dir,r(i),rounds,1);
            end
        else
            evaluate_job(dir,r(jobid),rounds,1);
        end
        
    case 'plot' 
        % produces the plot experiment number vs. iAUC score
        iAUC = zeros(r(1).num_inclus+1,length(r));
        diAUC = iAUC;
        lbAUC = iAUC;
        ubAUC = iAUC;
        iAUPRC = iAUC;
        diAUPRC = iAUC;
        lbAUPRC = iAUC;
        ubAUPRC = iAUC;
        for i = 1:length(r)
            [hh,hhPR] = evaluate_job(dir,r(i),rounds,1.0);
            iAUC(:,i) = hh(:,1);
            lbAUC(:,i) = hh(:,2);
            ubAUC(:,i) = hh(:,3);
            iAUPRC(:,i) = hhPR(:,1);
            lbAUPRC(:,i) = hhPR(:,2);
            ubAUPRC(:,i) = hhPR(:,3);
        end

        x = repmat(r(1).data_opt.num_initial+(0:r(1).num_inclus)',1,length(r));

        % ploto area under ROC
        figure('InvertHardcopy','off','Color',[1 1 1]);
%         subplot(1,1,1)
        cla
        set(gca,'Fontsize',14);
        set(gca,'Xscale','linear','Yscale','linear');
        hold on
        ww  = 1:size(iAUC,2);
        myplot(x(:,1),iAUC(:,ww),lbAUC(:,ww),ubAUC(:,ww));
%         errorbar(x,iAUC,diAUC,diAUC);

        set(gca,'XTickMode','auto');
        set(gca,'XTickLabelMode','auto');
        hold off
        axis tight
        axis([0 r(1).data_opt.N+1 0 1]);
        ylabel('iAUC');
        legend(labels(1:min(length(r),length(labels))),'Location','NorthWest');

        xlabel 'Experiment number'
        set(gca,'Box','on')
        xlim([1 r(1).num_inclus])
        
%         % ploto area under P-R
%         figure('InvertHardcopy','off','Color',[1 1 1]);
% %         subplot(1,1,1)
%         cla
%         set(gca,'Fontsize',14);
%         set(gca,'Xscale','linear','Yscale','linear');
%         hold on
%         ww  = 1:size(iAUPRC,2);
%         myplot(x(:,1),iAUPRC(:,ww),lbAUPRC(:,ww),ubAUPRC(:,ww));
% %         errorbar(x,iAUC,diAUC,diAUC);
% 
%         set(gca,'XTickMode','auto');
%         set(gca,'XTickLabelMode','auto');
%         hold off
%         axis tight
%         axis([0 r(1).data_opt.N+1 0 1]);
%         ylabel('iAUPRC');
%         legend(labels(1:min(length(r),length(labels))),'Location','NorthWest');
% 
%         xlabel 'Experiment number'
%         set(gca,'Box','on')

    case 'plotclass'
        % Produces the bar plots that compare different datasets
        % call this function with jobid equalling the number of groups in
        % the data (= 2 for the plots in the paper)
        
        % load all data
        iAUC = zeros(r(1).num_inclus+1,length(r));
        diAUC = zeros(r(1).num_inclus+1,length(r));
        for i = 1:length(r)
            hh = evaluate_job(dir,r(i),rounds,1);
            iAUC(:,i) = hh(:,1);
            diAUC(:,i) = hh(:,2);
        end
        % avarage over second half of experiements
        met = @(x)(sum(x(round(r(1).num_inclus/2):end,:))/(r(1).num_inclus+1));
        val = met(iAUC)';
        dval = met(diAUC)';

        % select different plots depending on which quantity is varyine
        x = (1:length(r))';
        dd = [r.data_opt];
        if std([dd.input_scale]) > 0
            x = [dd.input_scale]';
        elseif std([dd.sigma_noise]) > 0
            x = [dd.sigma_noise]';
        elseif std([dd.num_input]) > 0
            x = [dd.num_input]';
        end

        % How many experiments are there per x-group?
        num_compare = jobid;
        if num_compare ~= -1
            nv = []; dnv = [];
            for i = 1:num_compare
                nv = [nv,val(i:num_compare:end)];
                dnv = [dnv,dval(i:num_compare:end)];
            end
            val = nv; dval = dnv;
            x = x(1:num_compare:end);
        end

        cla
        subplot(1,1,1)
        set(gca,'Fontsize',25);
        barerrorbar(val,dval);
        set(gca,'XTickLabel',labels);
        axis tight
        ax = axis;
        axis([ax(1:2),0,1]);
        ylabel('iAUC averaged');
        xlabel(props{2});
        title(props{1});
        if strcmp(exp_type,'inputshape')
            legend({'LD','LR'},'Location','Northwest');
        else
            legend({'LD','LR'});
        end

    otherwise
        error('Unkown action.');
end




function arr = add(arr,dat)
if isempty(arr)
    arr = dat;
else
    arr(length(arr)+1) = dat;
end
