%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% Visualize dynamics in teclistamab model for individual virtual        %
% patients.                                                             %
% Updated: 3/21/2025                                                    %
%  - Code prompts user select dose of teclistamab                       %
%  - If select a dose of 0.72 or 1.5, it will include experimental      %
%    data in plots.                                                     %
%  - For any dose used in the protocol sweep code, and for a            %
%    frequency of either 1 or 2 (that the user must set), this code     %
%    will read in the parameters for the 1200 VPs generated in          %
%    protocol_sweep_VPs_subgroups.m                                     %
%  - Note: must download output_teclistamab_VCT_parallel.mat and save   %
%    in same directory as this code in order fo this to run for any     %
%    dose.                                                              %
%  - After plotting individual trajectories for plasma concentration    %
%    and helpful trimer (along with 95% confidence intervals), the      %
%    code will also separate out results by subgroups identified        %
%    in the protocol_sweep_VPs_subgroups.m code. If a dose is entered   %
%    that is not analyzed in that code, this code will throw an error   %
%    when it tries to display results by subgroups.                     %
%                                                                       %
% 80 VPs = 190 sec => 2850 sec for 1200 VPs                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; clc; close all; tic;
options=odeset('RelTol',1.0e-6); 
load output_teclistamab_VCT_parallel.mat
frequency = 1; % frequency to analyze for sensitivity
[p,ICs] = set_parameters_ICs();
baseline_p = p;
baseline_ICs = ICs; 
if (p.LD1 > 0) && (p.LD2 > 0) % loading dose
    time_add = p.intvlD*2;
else
    time_add = 0;
end
tmax = time_add +  frequency * 7 * (1+p.dosenumD); % changes based on frequency

dose_prompt = "Enter the dose (have experimental data for 0.72 and 1.5): ";
dose = input(dose_prompt);
[tDose, cDose, ti] = set_protocol(p, tmax, dose, frequency);
if dose == 0.72
    PK_SC_t    = [0  0.1  1  2	3  5  7]+ time_add;  % PK data
    PK_SC_drug = [1252.630175  1326.291716  3118.900845  4159.151398  ...
        4951.204559  3977.371187  4621.028168];									
elseif dose == 1.5
    PK_SC_t    = [0	 0.1  1	 2	3  7]+ time_add;	    % PK data
    PK_SC_drug = [1571.719325  1712.376561	5347.643161  7335.807348 ...
    	8488.477506	 7923.179002];
else % Generate VP trajectories for all VP parameterizations
    PK_SC_t    = [];
    PK_SC_drug = [];
end

%Nplot = 80;
Nplot = N_VPs;
VP_params_plot = VP_params(1:Nplot,:);
time_VPs = cell(Nplot,1);
Dp_VPs = cell(Nplot,1);
trimer_VPs = cell(Nplot,1);  
  
% npools = 4; % max 4 at home
% parpool('local', npools); % Open distributed processing 
% poolobj = gcp('nocreate'); % If no pool, don’t create 
% if isempty(poolobj)
%     poolsize = 0;
% else
%     poolsize = poolobj.NumWorkers;
% end

tic;
%% Equilibrate the system before dosing
p.equilibrate = 50;
equilibrated_ICs_all = zeros(Nplot, length(ICs)); % Store equilibrated ICs for all VPs
sBbm_values = zeros(Nplot, 1); % Store sBbm values for all VPs
mBbm_values = zeros(Nplot, 1); % Store mBbm values for all VPs
sBp_values = zeros(Nplot, 1); % Add this to store sBp values for all VPs
for vp = 1:Nplot
    p.Cl1 = VP_params_plot(vp, 1);
    p.kshed = VP_params_plot(vp, 2);
    ICs(10) = VP_params_plot(vp, 3); % mBbm0 
    p.k1bm = VP_params_plot(vp, 4); % Assign k1bm from the VP parameters
    ICs(12) = VP_params_plot(vp, 5); % Assign CD3bm0 from the VP parameters
    p.k10 = p.Cl1 / p.V1;
    p.ksyn_mB = p.kint_mB * ICs(10);
    p.ksyn_bm = ICs(12)*p.kint3;

    equilibrated_ICs = equilibrate_system(ICs, p, options);
    equilibrated_ICs_all(vp, :) = equilibrated_ICs;
    
    sBp_values(vp) = equilibrated_ICs(3);   % Store sBp after equilibration (new addition)
    mBbm_values(vp) = equilibrated_ICs(10); % mBbm after equilibration
    sBbm_values(vp) = equilibrated_ICs(11); % sBbm after equilibration

    fprintf('VP #%d equilibrated to sBP = %f, mBbm = %f, sBbM = %f\n',...
        vp,sBp_values(vp),mBbm_values(vp),sBbm_values(vp));
end


%par
for i = 1:Nplot
    ICs = zeros(size(equilibrated_ICs_all,2),1);
    for k = 1:size(equilibrated_ICs_all,2)
        ICs(k) = equilibrated_ICs_all(vp, k);
    end

    fprintf('#%d: Solving for this VP\n',i)
    [time, X] = solve_model(tDose, cDose, ICs, ti, p, options, VP_params_plot(i,:));
    time_VPs{i} = time;
    Dp_VPs{i} = X(:,6)* p.MW; % units scaled to ng/ml
    trimer_VPs{i} = replace_negatives(X(:, 18));
end
%delete(gcp('nocreate'));
toc

%% Model behavior at baseline/nominal parameterization
VP_params_baseline = [baseline_p.Cl1 baseline_p.kshed baseline_ICs(10) ...
    baseline_p.k1bm baseline_ICs(12)];
[time_baseline, X] = solve_model(tDose, cDose, ICs, ti, p, options, ...
    VP_params_baseline);
trimer_baseline = replace_negatives(X(:, 18)); % helpful trimer in the BM
Dp_baseline = X(:,6)* p.MW; % free drug in plasma, units scaled to ng/ml
AUC_deathterm_baseline = trapz(time_baseline, trimer_baseline);
fprintf('Baseline scenario has AUC(death_term) = %f\n',AUC_deathterm_baseline); 

figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05, 0.05, 0.45, 0.6]);
sgtitle(['Nominal Parameterization Dynamics: Dose of ' num2str(dose) ...
    'mg/kg (SC)'],'FontSize',18','FontWeight','bold')
subplot(1,2,1)
semilogy(time_baseline,Dp_baseline,'LineWidth',2,'Color','k');
xlabel('Time (days)','FontSize',16)
ylabel('Plasma concentration (ng/ml)','FontSize',16)
ylim([1e1,inf])

subplot(1,2,2)
semilogy(time_baseline,trimer_baseline,'LineWidth',2,'Color','k');
xlabel('Time (days)','FontSize',16)
ylabel('Concentration of active trimer in TME','FontSize',16)
ylim([1e-4,inf])

%% Model behavior at median parameterization
VP_params_median = median(VP_params)
[time_median, X_median] = solve_model(tDose, cDose, ICs, ti, p, options, ...
    VP_params_median);
trimer_median = replace_negatives(X_median(:, 18)); % helpful trimer in the BM
Dp_median = X_median(:,6)* p.MW; % free drug in plasma, units scaled to ng/ml
AUC_deathterm_median = trapz(time_median, trimer_median);
fprintf('Median scenario has AUC(death_term) = %f\n',AUC_deathterm_median); 

%% Now need 95% confidence interval
[find_min_Dp, find_max_Dp, min_Dp_all, max_Dp_all] = ...
    find_min_max(time_VPs,Dp_VPs,tDose);
[find_min_trimer, find_max_trimer, min_trimer_all, max_trimer_all] = ...
    find_min_max(time_VPs,trimer_VPs,tDose);

num_remove = floor(0.025*Nplot); 
[min_remove_Dp, max_remove_Dp] = find_95perc_CI(num_remove,min_Dp_all,...
    max_Dp_all);
[min_remove_trimer, max_remove_trimer] = find_95perc_CI(num_remove,...
    min_trimer_all,max_trimer_all);

% Plots each VP trajectory
f2 = figure('visible','off');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05, 0.05, 0.6, 0.65]);
sgtitle(['Virtual Population Dynamics at Dose of ' num2str(dose) ...
    ' mg/kg (SC)'], 'FontSize',18','FontWeight','bold')
subplot(1,2,1)
first_in_CI = 0; first_not_in_CI = 0; 
for i = 1:Nplot
    if isempty(find(min_remove_Dp == i,1)) && isempty(find(max_remove_Dp == i,1))
        if first_in_CI == 0
            h1 = plot(time_VPs{i},Dp_VPs{i},'LineWidth',2,'Color',...
                [0.7, 0.7, 0.7],'DisplayName','95% Confidence Interval');
            % h1 = semilogy(time_VPs{i},Dp_VPs{i},'LineWidth',2,'Color',...
            %    [0.7, 0.7, 0.7],'DisplayName','95% Confidence Interval');
        else 
            plot(time_VPs{i},Dp_VPs{i},'LineWidth',2,'Color',[0.7, 0.7, 0.7]);
            % semilogy(time_VPs{i},Dp_VPs{i},'LineWidth',2,'Color',[0.7, 0.7, 0.7]);
        end
        first_in_CI = first_in_CI+1; 
    else % Not in 95% confidence interval
        if first_not_in_CI == 0
            h2 = plot(time_VPs{i},Dp_VPs{i},'LineWidth',2,'Color',...
                [0.9, 0.9, 0.9],'DisplayName','VP Range');
            % h2 = semilogy(time_VPs{i},Dp_VPs{i},'LineWidth',2,'Color',...
            %     [0.9, 0.9, 0.9],'DisplayName','95% Confidence Interval');
        else 
            plot(time_VPs{i},Dp_VPs{i},'LineWidth',2,'Color',[0.9, 0.9, 0.9]);
            % semilogy(time_VPs{i},Dp_VPs{i},'LineWidth',2,'Color',[0.9, 0.9, 0.9]);
        end
        first_not_in_CI = first_not_in_CI+1; 
    end

    if i == 1
        hold on;
    end
end
%semilogy(time_baseline,Dp_baseline,'LineWidth',2,'Color','k');
h3 = plot(time_baseline,Dp_baseline,'--','LineWidth',2,'Color','k',...
    'DisplayName','Nominal Parameters');
h4 = plot(time_median,Dp_median,'LineWidth',2,'Color','k',...
    'DisplayName','Median Parameters');
if (dose == 1.5)||(dose==0.72)
    h5 = plot(PK_SC_t,PK_SC_drug,'xr','LineWidth',2,'DisplayName','Data');
end
hold off; 
xlim([0 inf])
xlabel('Time (days)','FontSize',16)
ylabel('Plasma concentration (ng/ml)','FontSize',16)
if (dose == 1.5)||(dose==0.72)
    legend([h1 h2 h3 h4 h5],'FontSize',16,'Location','NorthWest');
else
    legend([h1 h2 h3 h4],'FontSize',16,'Location','NorthWest');
end

subplot(1,2,2)
for i = 1:Nplot
    if isempty(find(min_remove_trimer == i,1)) && isempty(find(max_remove_trimer == i,1))
        plot(time_VPs{i},trimer_VPs{i},'LineWidth',2,'Color',...
            [0.7, 0.7, 0.7]);
        %semilogy(time_VPs{i},trimer_VPs{i},'LineWidth',2,'Color',[0, 0, 0.]);
    else% Not in 95% confidence interval
        %semilogy(time_VPs{i},trimer_VPs{i},'LineWidth',2,'Color',[0.2, 0.2,0.2]0);
        plot(time_VPs{i},trimer_VPs{i},'LineWidth',2,'Color',...
            [0.9, 0.9, 0.9]);
    end
    if i == 1
        hold on;
    end
end
%semilogy(time_baseline,trimer_baseline,'LineWidth',2,'Color','k');
plot(time_baseline,trimer_baseline,'--','LineWidth',2,'Color','k');
plot(time_median,trimer_median,'LineWidth',2,'Color','k');
hold off;
xlim([0 inf])
xlabel('Time (days)','FontSize',16)
ylabel('Concentration of active trimer in TME (nM)','FontSize',16)
saveas(f2, ['pop_dynamics_dose' num2str(dose) '.fig']);
saveas(f2, ['pop_dynamics_dose' num2str(dose) '.png']);

% Plots each VP trajectory by subpopulation
sBp_VPs = zeros(Nplot,1);
for i = 1:Nplot
    in_sBp_low = find(i ==subpopulations_sBp{1});
    in_sBp_high = find(i ==subpopulations_sBp{2});
    if isempty(in_sBp_low) && (in_sBp_high>=1) % in sBp high
        fprintf('Vpop #%d in sBp high group\n',i);
        sBp_VPs(i) = 2; % 2 for high
    elseif isempty(in_sBp_high) && (in_sBp_low>=1) % in sBp low
        fprintf('Vpop #%d in sBp low group\n',i);
        sBp_VPs(i) = 1; % 1 for low
    else
        fprintf('Problem: VPop %d not in either group\n',i);
        stop
    end
end

%% Same plot of individual trajectories, but instead of coloring by 
%% confidence interval, color by sBp subgroup
f3 = figure('visible','off');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05, 0.05, 0.55, 0.95]);
sgtitle(['Virtual Population Dynamics at Dose of ' num2str(dose) ...
    ' mg/kg (SC)'], 'FontSize',18','FontWeight','bold')
subplot(2,2,1)
first_in_sub1 = 0; first_in_sub2 = 0; 
for i = 1:Nplot
    if sBp_VPs(i) == 1 % in sBp low
        if first_in_sub1 == 0
            h1 = plot(time_VPs{i},Dp_VPs{i},'LineWidth',2,'Color',...
                [0., 0.45, 0.74],'DisplayName','sBp Low Subgroup');
        else 
            plot(time_VPs{i},Dp_VPs{i},'LineWidth',2,'Color',[0., 0.45, 0.74])
        end
        first_in_sub1 = first_in_sub1+1; 

    elseif sBp_VPs(i) == 2 % in sBp high
        if first_in_sub2 == 0
            h2 = plot(time_VPs{i},Dp_VPs{i},'LineWidth',2,'Color',...
                [0.85, 0.33, 0.1],'DisplayName','sBp High Subgroup');
        else 
            plot(time_VPs{i},Dp_VPs{i},'LineWidth',2,'Color',[0.85, 0.33, 0.1]);
        end
        first_in_sub2 = first_in_sub2+1; 
    end
    if i == 1
        hold on;
    end
end
% h3 = plot(time_baseline,Dp_baseline,'--','LineWidth',2,'Color','k',...
%     'DisplayName','Nominal Parameters');
% h4 = plot(time_median,Dp_median,'LineWidth',2,'Color','k',...
%     'DisplayName','Median Parameters');
if (dose == 1.5)||(dose==0.72)
    h5 = plot(PK_SC_t,PK_SC_drug,'xk','LineWidth',2,'DisplayName','Data');
end
hold off; 
xlim([0 inf])
xlabel('Time (days)','FontSize',16)
ylabel('Plasma concentration (ng/ml)','FontSize',16)
if (dose == 1.5)||(dose==0.72)
    legend([h1 h2 h5],'FontSize',16,'Location','NorthWest');
else
    legend([h1 h2],'FontSize',16,'Location','NorthWest');
end

subplot(2,2,2)
for i = 1:Nplot
    if sBp_VPs(i) == 1 % in sBp low
        plot(time_VPs{i},trimer_VPs{i},'LineWidth',2,'Color',[0., 0.45, 0.74])
    elseif sBp_VPs(i) == 2 % in sBp high
        plot(time_VPs{i},trimer_VPs{i},'LineWidth',2,'Color',[0.85, 0.33, 0.1]);
    end
    if i == 1
        hold on;
    end
end
%plot(time_baseline,trimer_baseline,'--','LineWidth',2,'Color','k');
%plot(time_median,trimer_median,'LineWidth',2,'Color','k');
hold off;
xlim([0 inf])
xlabel('Time (days)','FontSize',16)
ylabel('Conc. Active Trimer in TME (nM)','FontSize',16)

subplot(2,2,3) % show helpful trimer for sBp low only
for i = 1:Nplot
    if sBp_VPs(i) == 1 % in sBp low
        plot(time_VPs{i},trimer_VPs{i},'LineWidth',2,'Color',[0., 0.45, 0.74])
    end
    if i == 1
        hold on;
    end
end
%plot(time_baseline,trimer_baseline,'--','LineWidth',2,'Color','k');
%plot(time_median,trimer_median,'LineWidth',2,'Color','k');
hold off;
xlim([0 inf])
xlabel('Time (days)','FontSize',16)
ylabel('Conc. Active Trimer in TME (nM)','FontSize',16)
title('sBp Low Subgroup','FontSize',16);

subplot(2,2,4) % show helpful trimer for sBp high only
for i = 1:Nplot
    if sBp_VPs(i) == 2 % in sBp high
        plot(time_VPs{i},trimer_VPs{i},'LineWidth',2,'Color',[0.85, 0.33, 0.1]);
    end
    if i == 1
        hold on;
    end
end
%plot(time_baseline,trimer_baseline,'--','LineWidth',2,'Color','k');
%plot(time_median,trimer_median,'LineWidth',2,'Color','k');
hold off;
xlim([0 inf])
xlabel('Time (days)','FontSize',16)
ylabel('Conc. Active Trimer in TME (nM)','FontSize',16)
title('sBp High Subgroup','FontSize',16);
yl_2pop = ylim;
saveas(f3, ['subpop_dynamics_sBp_dose' num2str(dose) '.fig']);
saveas(f3, ['subpop_dynamics_sBp_dose' num2str(dose) '.png']);

% % subpopulations = {
%     intersect(sBp_low, mBbm_low),
%     intersect(sBp_low, mBbm_high),
%     intersect(sBp_high, mBbm_low),
%     intersect(sBp_high, mBbm_high)};
sBp_mBbm_VPs = zeros(Nplot,1);
for i = 1:Nplot
    in_sBp_low_mBbm_low = find(i ==subpopulations{1});
    in_sBp_low_mBbm_high = find(i ==subpopulations{2});
    in_sBp_high_mBbm_low = find(i ==subpopulations{3});
    in_sBp_high_mBbm_high = find(i ==subpopulations{4});

    % in sBp_low, mBbm_low
    if isempty(in_sBp_low_mBbm_high) && isempty(in_sBp_high_mBbm_low) && ...
        isempty(in_sBp_high_mBbm_high) && (in_sBp_low_mBbm_low>=1) 
        fprintf('Vpop #%d in sBp low, mBbm low group\n',i);
        sBp_mBbm_VPs(i) = 1; 

    % in sBp_low, mBbm_high
    elseif isempty(in_sBp_low_mBbm_low) && isempty(in_sBp_high_mBbm_low) && ...
        isempty(in_sBp_high_mBbm_high) && (in_sBp_low_mBbm_high>=1) 
        fprintf('Vpop #%d in sBp low, mBbm high group\n',i);
        sBp_mBbm_VPs(i) = 2; 

    % in sBp_high, mBbm_low
    elseif isempty(in_sBp_low_mBbm_low) && isempty(in_sBp_low_mBbm_high) && ...
        isempty(in_sBp_high_mBbm_high) && (in_sBp_high_mBbm_low>=1) 
        fprintf('Vpop #%d in sBp high, mBbm low group\n',i);
        sBp_mBbm_VPs(i) = 3; 

    % in sBp_high, mBbm_high
    elseif isempty(in_sBp_low_mBbm_low) && isempty(in_sBp_low_mBbm_high) && ...
        isempty(in_sBp_high_mBbm_low) && (in_sBp_high_mBbm_high>=1) 
        fprintf('Vpop #%d in sBp high, mBbm high group\n',i);
        sBp_mBbm_VPs(i) = 4; 

    else
        fprintf('Problem: VPop %d not in any group\n',i);
        stop
    end
end

%% Same plot of individual trajectories, but instead of coloring by 
%% confidence interval, color by sBp and mBbm subgroup
f4 = figure('visible','off');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05, 0.05, 0.5, 0.95]);
sgtitle(['Virtual Population Dynamics at Dose of ' num2str(dose) ...
    ' mg/kg (SC)'], 'FontSize',18','FontWeight','bold')
subplot(3,2,2)
first_in_sub1 = 0; first_in_sub2 = 0; 
first_in_sub3 = 0; first_in_sub4 = 0; 
for i = 1:Nplot
    if sBp_mBbm_VPs(i) == 1 % in sBp_low, mBbm_low
        if first_in_sub1 == 0
            h1 = plot(time_VPs{i},Dp_VPs{i},'LineWidth',2,'Color',...
                [0., 0.45, 0.74],'DisplayName','sBp Low, mBbm Low Subgroup');
        else 
            plot(time_VPs{i},Dp_VPs{i},'LineWidth',2,'Color',[0., 0.45, 0.74])
        end
        first_in_sub1 = first_in_sub1+1; 

    elseif sBp_mBbm_VPs(i) == 2 % in sBp_low, mBbm_high
        if first_in_sub2 == 0
            h2 = plot(time_VPs{i},Dp_VPs{i},'LineWidth',2,'Color',...
                [0.85, 0.33, 0.1],'DisplayName','sBp Low, mBbm High Subgroup');
        else 
            plot(time_VPs{i},Dp_VPs{i},'LineWidth',2,'Color',[0.85, 0.33, 0.1]);
        end
        first_in_sub2 = first_in_sub2+1; 

    elseif sBp_mBbm_VPs(i) == 3 % in sBp_high, mBbm_low
        if first_in_sub3 == 0
            h3 = plot(time_VPs{i},Dp_VPs{i},'LineWidth',2,'Color',...
                [0.93, 0.69, 0.13],'DisplayName','sBp High, mBbm low Subgroup');
        else 
            plot(time_VPs{i},Dp_VPs{i},'LineWidth',2,'Color',[0.93, 0.69, 0.13]);
        end
        first_in_sub3 = first_in_sub3+1; 
    
    else%  in sBp_high, mBbm_high
        if first_in_sub4 == 0
            h4 = plot(time_VPs{i},Dp_VPs{i},'LineWidth',2,'Color',...
                [0.49, 0.18, 0.56],'DisplayName','sBp High, mBbm low Subgroup');
        else 
            plot(time_VPs{i},Dp_VPs{i},'LineWidth',2,'Color',[0.49, 0.18, 0.56]);
        end
        first_in_sub4 = first_in_sub4+1; 
    end
    if i == 1
        hold on;
    end
end
% h5 = plot(time_baseline,Dp_baseline,'--','LineWidth',2,'Color','k',...
%     'DisplayName','Nominal Parameters');
% h6 = plot(time_median,Dp_median,'LineWidth',2,'Color','k',...
%     'DisplayName','Median Parameters');
if (dose == 1.5)||(dose==0.72)
    h7 = plot(PK_SC_t,PK_SC_drug,'xk','LineWidth',2,'DisplayName','Data');
end
hold off; 
xlim([0 inf])
%xlabel('Time (days)','FontSize',16)
%ylabel('Plasma conc. (ng/ml)','FontSize',16)
if (dose == 1.5)||(dose==0.72)
    legend([h1 h2 h3 h4 h7],'FontSize',16,'Location','NorthWest');
else
    legend([h1 h2 h3 h4],'FontSize',16,'Location','NorthWest');
end

subplot(3,2,1)
for i = 1:Nplot
    if sBp_mBbm_VPs(i) == 1 % in sBp_low, mBbm_low
        plot(time_VPs{i},trimer_VPs{i},'LineWidth',2,'Color',[0., 0.45, 0.74]);
    elseif sBp_mBbm_VPs(i) == 2 % in sBp_low, mBbm_high
        plot(time_VPs{i},trimer_VPs{i},'LineWidth',2,'Color',[0.85, 0.33, 0.1]);
    elseif sBp_mBbm_VPs(i) == 3 % in sBp_high, mBbm_low
        plot(time_VPs{i},trimer_VPs{i},'LineWidth',2,'Color',[0.93, 0.69, 0.13]);
    else%  in sBp_high, mBbm_high
         plot(time_VPs{i},trimer_VPs{i},'LineWidth',2,'Color',[0.49, 0.18, 0.56]);
    end
    if i == 1
        hold on;
    end
end
%plot(time_baseline,trimer_baseline,'--','LineWidth',2,'Color','k');
%plot(time_median,trimer_median,'LineWidth',2,'Color','k');
hold off;
xlim([0 inf])
ylim(yl_2pop)
xlabel('Time (days)','FontSize',16)
ylabel('PAT in TME(nM)','FontSize',16)
yl = ylim;

subplot(3,2,3)
for i = 1:Nplot
    if sBp_mBbm_VPs(i) == 1 % in sBp_low, mBbm_low
        plot(time_VPs{i},trimer_VPs{i},'LineWidth',2,'Color',[0., 0.45, 0.74]);
    end
    if i == 1
        hold on;
    end
end
%plot(time_baseline,trimer_baseline,'--','LineWidth',2,'Color','k');
%plot(time_median,trimer_median,'LineWidth',2,'Color','k');
hold off;
xlim([0 inf])
ylim(yl);
xlabel('Time (days)','FontSize',16)
ylabel('PAT in TME (nM)','FontSize',16)
title('sBp Low, mBbm low Subgroup','FontSize',16);


subplot(3,2,4)
for i = 1:Nplot
    if sBp_mBbm_VPs(i) == 2 % in sBp_low, mBbm_high
        plot(time_VPs{i},trimer_VPs{i},'LineWidth',2,'Color',[0.85, 0.33, 0.1]);
    end
    if i == 1
        hold on;
    end
end
%plot(time_baseline,trimer_baseline,'--','LineWidth',2,'Color','k');
%plot(time_median,trimer_median,'LineWidth',2,'Color','k');
hold off;
xlim([0 inf])
ylim(yl);
xlabel('Time (days)','FontSize',16)
ylabel('PAT in TME (nM)','FontSize',16)
title('sBp Low, mBbm high Subgroup','FontSize',16);

subplot(3,2,5)
for i = 1:Nplot
    if sBp_mBbm_VPs(i) == 3 % in sBp_high, mBbm_low
        plot(time_VPs{i},trimer_VPs{i},'LineWidth',2,'Color',[0.93, 0.69, 0.13]);
    end
    if i == 1
        hold on;
    end
end
%plot(time_baseline,trimer_baseline,'--','LineWidth',2,'Color','k');
%plot(time_median,trimer_median,'LineWidth',2,'Color','k');
hold off;
xlim([0 inf])
ylim(yl);
xlabel('Time (days)','FontSize',16)
ylabel('PAT in TME (nM)','FontSize',16)
title('sBp High, mBbm low Subgroup','FontSize',16);

subplot(3,2,6)
for i = 1:Nplot
    if sBp_mBbm_VPs(i) == 4 % in sBp_high, mBbm_high
         plot(time_VPs{i},trimer_VPs{i},'LineWidth',2,'Color',[0.49, 0.18, 0.56]);
    end
    if i == 1
        hold on;
    end
end
%plot(time_baseline,trimer_baseline,'--','LineWidth',2,'Color','k');
%plot(time_median,trimer_median,'LineWidth',2,'Color','k');
hold off;
xlim([0 inf])
ylim(yl);
xlabel('Time (days)','FontSize',16)
ylabel('PAT in TME (nM)','FontSize',16)
title('sBp High, mBbm High Subgroup','FontSize',16);
saveas(f4, ['subpop_dynamics_sBp_mBm_dose' num2str(dose) '.fig']);
saveas(f4, ['subpop_dynamics_sBp_mBm_dose' num2str(dose) '.png']);

% openfig('newout.fig','new','visible')
toc



%%%%%%%%% Functions%%%%%%%%%%%%%

function [p,ICs] = set_parameters_ICs()
    %% IV or SC
    % p.iv = 1 => IV
    % p.iv = 0 => SC
    p.iv = 0;

    %% parameters we will be changing for doses - updated IK 3/19/2025 (added loading doses)
    p.firstDoseTime = 0; % Specify the first dose timing

    % a lot of these protocols have loading doses, so I'm adding an option
    % for a loading dose
    p.LD1 = .06;  % First loading dose
    p.LD2 = 0.3;  % Second loading dose
    p.intvlLD = 7;  % Interval (in days) between the loading doses

    % full doses now
    p.dose_bispecific = 3; 
    %p.dose_bispecific = .72; %.72;
    p.dosenumD = 10;   % number of doses
    p.freq=1;
    p.intvlD = p.freq*7; % dosing interval, days

    %% Initial conditions 
    CD3th0  = 0.001;% free CD3 target in thymus
    CD3p0   = .01;% free CD3 target in plasma

    % normal BCMA in plasma is 
    sBp0    = 20;% free soluble BCMA in plasma
    
    Dsc0   = 0;% subcutaneous administration depot
    Dper0  = 0;% free drug in the peripheral compartment
    Dp0    = 0;% free drug in plasma
    
    D3p0   = 0;% drug + CD3 complex in plasma
    DsBp0  = 0;% drug + soluble BCMA complex in plasma
    DsB3p0 = 0;% (useless) trimer in plasma
    
    % bone marrow!
    mBbm0  = 50;% membrane-bound BCMA
    sBbm0  = 3;% soluble BCMA
    CD3bm0 = .1;% free CD3 in bone marrow
    
    Dbm0   = 0; % free drug in BM
    D3bm0  = 0; % dimer: drug+CD3 in BM
    DsBbm0 = 0; % dimer: drug + soluble BCMA
    DmBbm0 = 0; % dimer: drug + membrane-bound BCMA
    
    DsB3bm0 = 0;% useless trimer in the BM
    DmB3bm0 = 0;% helpful trimer in the BM



    %% PK parameters
    p.MW = 146; %kDa
    p.BW = 74; %human bodyweight
    if p.iv == 1
        p.k01 = 0;
    else
        p.k01 = 0.2; % 1/day SC absorption https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10518021/ 
    end

    p.V1  = 3.0*1000/p.BW; % = 40 mL/kg  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10518021/ 
    p.V2  = 1.34*1000/p.BW; % mL/kg https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10518021/ 
    p.Cl1 = .449*1000/p.BW;% = 7 mL/kg/day https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10518021/ 
    p.Cl2 = 0.4*1000/p.BW ;% mL/kg/day https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10518021/ 
    p.Vbm = 2.2*1000/p.BW; % mL/kg; bone marrow occupies 1.5-3L in a typical adult, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5738992/
        
    p.k10 = p.Cl1/p.V1;
    p.k12 = p.Cl2/p.V1;
    p.k21 = p.Cl2/p.V2;

    p.cl_tvar = 0.547*1000/p.BW; %? 
    p.kdes = 0.1;


    p.doseD = p.dose_bispecific*10^6/p.V1/p.MW;% mpk
    
    %% kons and koffs
    %% KD values 
    p.KD_cd3  = 28.03; % nM; KD = koff/kon
    p.KD_BCMA = 0.18; % nM
    
    p.koff3p = 285; % guess
    p.kon3p = p.koff3p/p.KD_cd3; % kon = koff/KD
    
    p.koff3bm = 285; %no idea but leaving an option to be different from plasma
    p.kon3bm = p.koff3bm/p.KD_cd3; % kon = koff/KD
    
    p.koff_sB = 20; %guess
    p.kon_sB = p.koff_sB/p.KD_BCMA; % kon = koff/KD
    
    p.koff_mB = 20; %guess
    p.kon_mB = p.koff_mB/p.KD_BCMA; 
    
    %% 
    p.kshed = 0.05; %ranges I found are from 0 - 0.04
    BCMA_halflife = 1.25;% days, 24-36 hours, I used 30/24 here
    
    p.kint_mB = 0.693/BCMA_halflife; %per day, would also be around 0.02-0.03/h
    p.ksyn_mB = mBbm0*(p.kint_mB+p.kshed); %synthesis of mBCMA in bone marrow
    p.kint3 = 1.584 ;% from the epco
    
    %%
    p.kth1    = 0.5; %rate of outflow of newly formed CD3 from thymus to plasma
    p.ksyn_th = CD3th0*p.kth1; % rate of CD3 synthesis in plasma
    p.ksyn_bm = CD3bm0*p.kint3;
    
    
    %% to estimate with data - updated IK 3/19/2025
    p.k1bm = 0.1; %random so far, outflow of free drug from plasma to BM
    p.kbm1 = 0.001; % I had: 0.001; % random so far, inflow of free drug from BM to plasma
    
    p.k3bm1 = 0.1; % CD3 migration from bone marrow to plasma; find correct parameters based on BM data
    p.k31bm = 0.01;
    p.kBbm1 = 0.1; % soluble sBCMA migration from BM to plasma, also a guess now
    p.kB1bm = 0.1;
    p.k10sB = 0.04; %0.04; % half-life of soluble is 12-24h, so 0.693/.5 = 1.38
    p.kout = 0.01;%  had: 0.01; %just in case there's additional clearance of things from BM

    ICs = [CD3th0 CD3p0 sBp0 Dsc0 Dper0 Dp0 D3p0 DsBp0 DsB3p0 mBbm0 ...
        sBbm0 CD3bm0 Dbm0 D3bm0 DsBbm0 DmBbm0 DsB3bm0 DmB3bm0];
end

function [tDose, cDose, ti] = set_protocol(p, tmax, doses,frequencies)
  
    % Set protocol parameters
    doseD = doses * 10^6 / p.V1 / p.MW; % mpk
    freq = frequencies;
    intvlD = freq * 7;

    % If loading doses are specified and greater than 0, include them
    %if isfield(p, 'LD1') && isfield(p, 'LD2') && p.LD1 > 0 && p.LD2 > 0 
    if (p.LD1 > 0) && (p.LD2 > 0) 
        % Initialize dose and time arrays to accommodate loading doses and regular doses
        tDose = zeros(1, p.dosenumD + 2); % +2 for the two loading doses
        cDose = zeros(1, p.dosenumD + 2);
        
        % Add the first loading dose at time firstDoseTime
        tDose(1) = p.firstDoseTime;
        cDose(1) = p.LD1*10^6/p.V1/p.MW;
        
        % Add the second loading dose at time firstDoseTime + p.intvlLD (interval for loading doses)
        tDose(2) = p.firstDoseTime + p.intvlLD;
        cDose(2) = p.LD2*10^6/p.V1/p.MW;
        
        % Add the regular doses, starting after the loading doses
        for i = 1:p.dosenumD
            tDose(i + 2) = p.firstDoseTime + (2*p.intvlLD) + (intvlD*(i-1)); % Shift regular doses by p.intvlLD and firstDoseTime
            %tDose(i + 2) = firstDoseTime + intvlD + intvlD * i; % Shift regular doses by p.intvlLD and firstDoseTime
            cDose(i + 2) = doseD;
        end
    else
        % If no loading doses are specified, proceed with the regular dosing protocol
        
        tDose = zeros(1, p.dosenumD); 
        cDose = zeros(1, p.dosenumD);
        
        % Start regular dosing from time firstDoseTime
        tDose(1) = p.firstDoseTime;
        cDose(1) = doseD;

        for i = 2:p.dosenumD
            tDose(i) = p.firstDoseTime + intvlD * (i - 1);
            cDose(i) = doseD;
        end
    end
    
    % Ensure that the last time point is tmax for the simulation
    ti = [tDose tmax];
end


function [Time, X] = solve_model(tDose, cDose, ICs, ti, p, options,VP_params) 
    % Model parameters for VP
    p.Cl1 = VP_params(1);  % Correctly assigns Cl1 from VP parameters
    p.k10 = p.Cl1 / p.V1;  % And must recompute k10 accordingly

    p.kshed = VP_params(2); % Correctly assigns kshed from VP parameters
    % VP_params(3) = ICs(10) = mBbm0                                                
    p.ksyn_mB = (p.kint_mB+p.kshed) * VP_params(3); % And must recompute ksyn_mB
    p.k1bm = VP_params(4); % Correctly assigns k1bm from the VP parameters
    %  VP_params(5) = ICs(12) = CD3bm0
    p.ksyn_bm = ICs(12)*p.kint3; % And must recompute ksyn_bm

    tPoints = []; cPoints = [];
    for i = 1:length(tDose)
        if p.iv == 1 % if IV dosing
            ICs(6) = ICs(6) + cDose(i);  % Dp0 = x(6)
        elseif p.iv ==  0% if sc dosing
            ICs(4) = ICs(4) + cDose(i);     % Dsc0 = x(4)
        end
    
        % Build the tSpan and run the model
        tspan = [ti(i) ti(i+1)];
        [tOut, cOut] = ode23s(@(t,x) trimer_model(t,x,p), tspan, ICs, options);
        
        % Concatenate the results 
        tPoints = [tPoints; tOut];
        cPoints = [cPoints; cOut];
    
        % Set new ICs and repeat 
        for j = 1:length(ICs)
            ICs(j) = cOut(end,j); 
        end
    end

    X = cPoints;
    Time = tPoints;
end

function dxdt = trimer_model(t, x, p) 
    CD3th  = x(1);% free CD3 target in thymus
    CD3p   = x(2);% free CD3 target in plasma
    sBp    = x(3);% free soluble BCMA in plasma

    Dsc   = x(4);% subcutaneous administration depot
    Dper  = x(5);% free drug in the peripheral compartment
    Dp    = x(6);% free drug in plasma

    D3p   = x(7);% drug + CD3 complex in plasma
    DsBp  = x(8);% drug + soluble BCMA complex in plasma
    DsB3p = x(9);% (useless) trimer in plasma

    % bone marrow!
    mBbm  = x(10);% membrane-bound BCMA
    sBbm  = x(11);% soluble BCMA
    CD3bm = x(12);% free CD3 in bone marrow

    Dbm   = x(13);% free drug in BM
    D3bm  = x(14);% dimer: drug+CD3 in BM
    DsBbm = x(15);% dimer: drug + soluble BCMA
    DmBbm = x(16);% dimer: drug + membrane-bound BCMA

    DsB3bm = x(17);%useless trimer in the BM
    DmB3bm = x(18);% helpful trimer in the BM

  %% thumus, source of CD3+ T cells
    dCD3th = p.ksyn_th - p.kth1*CD3th; % synthesis - outflow to plasma
    %% CD3 target in plasma
    dCD3p = p.kth1*CD3th + p.k3bm1*(p.Vbm/p.V1)*CD3bm - p.k31bm*CD3p... %inflow from thymus and from bone marrow
        - p.kint3*CD3p... % internalization
        - p.kon3p*CD3p*Dp   + p.koff3p*D3p... %  bind/unbind with free drug Dp
        - p.kon3p*CD3p*DsBp + p.koff3p*DsB3p; % bind/unbind with drug-sB complex to form soluble trimer

    %% soluble BCMA in plasma
     dsBp = p.kBbm1*(p.Vbm/p.V1)*sBbm    - p.k10sB*sBp ... %inflow from BM and natural clearance
        - p.kon_sB*sBp*Dp  + p.koff_sB*DsBp... % bind/unbind soluble BCMA + drug in plasma
        - p.kon_sB*sBp*D3p + p.koff_sB*DsB3p... % bind/unbind with drug-CD3 complex to form soluble trimer
        - p.kB1bm*sBp; % distribution to the bone marrow

    %% drug in plasma
    Cl = p.Cl1 + p.cl_tvar*exp(-p.kdes*t); % Time dependent clearance
    p.k10 = Cl/p.V1;

    dDsc  = -p.k01*Dsc; %absorption
    dDper = p.k12*p.V1/p.V2*Dp - p.k21*Dper - p.k10*Dper; %back and forth to peripheral compartment; asusiming CD3 binding kinetics htere are negligible
    dDp   = p.k01*Dsc - p.k10*Dp  ... % absorption - (time-dependent clearance)
        - p.k12*Dp + p.k21*p.V2/p.V1*Dper... %absorption - clearance - out to peripheral and back
        - p.kon3p*CD3p*Dp + p.koff3p*D3p... % bind/unbind free drug + CD3
        - p.kon_sB*sBp*Dp + p.koff_sB*DsBp... % bind/unbind free drug  + soluble BCMA
        - p.k1bm*Dp + p.kbm1*(p.Vbm/p.V1)*Dbm; %back and forth to bone marrow

    %% complexes in plasma
    dD3p = p.kon3p*CD3p*Dp  - p.koff3p*D3p... % bind/unbind free drug + CD3
        - p.kon_sB*sBp*D3p + p.koff_sB*DsB3p... % bind/unbingd the trimer (D3p+sB complex)
        - p.kint3*D3p; %internalization/clearance

    dDsBp = p.kon_sB*sBp*Dp - p.koff_sB*DsBp... % free drug + sBCMA
        - p.kon3p*CD3p*DsBp + p.koff3p*DsB3p... % cd3p + drug-sB complex
        - p.k10* DsBp  ; %natural clearance
    dDsB3p = p.kon_sB*sBp*D3p - p.koff_sB* DsB3p... % formed from soluble BCMA + drug-CD3 complex
        + p.kon3p*CD3p*DsBp - p.koff3p*DsB3p... % formed from free CD3 + drug-sB complex
        - p.kint3*DsB3p... % internalized with CD3
        - p.k10*DsB3p;   % cleared with soluble BCMA

    %% bone marrow!
    % the 3 free targets
    dmBbm = p.ksyn_mB - p.kshed*mBbm - p.kint_mB*mBbm... %syn, shed, internalize
        - p.kon_mB*mBbm*Dbm  + p.koff_mB*DmBbm... % dimer formation: mBCMA + free drug
        - p.kon_mB*mBbm*D3bm + p.koff_mB* DmB3bm; % useful trimer formation: mBCMA + drug-CD3 complex

    dsBbm = p.kshed*mBbm - p.kBbm1*sBbm + (p.V1/p.Vbm)*p.k1bm*sBp... % inflow from shed mBCMA; outflow to plasma
        - p.kon_sB*sBbm*Dbm  + p.koff_sB*DsBbm... % dimer: sBCMA + drug
        - p.kon_sB*sBbm*D3bm + p.koff_sB* DsB3bm; % useless trimer: sBCMA + drug-CD3 complex

    dCD3bm = p.ksyn_bm - p.kint3*CD3bm - p.k3bm1*CD3bm + p.k31bm*(p.V1/p.Vbm)*CD3p...% syn, int, distribute
        - p.kon3bm*CD3bm*Dbm   + p.koff3bm*D3bm...   % dimer: cd3bm + drug
        - p.kon3bm*CD3bm*DsBbm + p.koff3bm*DsB3bm... % useless trimer: cd3bm + drug-sB complex
        - p.kon3bm*CD3bm*DmBbm + p.koff3bm*DmB3bm;   % helpful trimer: cd3bm + drug-mB complex

    % drug things now
    dDbm = (p.V1/p.Vbm)*p.k1bm*Dp - p.kbm1 *Dbm... %infrom from plasma and outflow back
        - p.kon_mB*mBbm*Dbm  + p.koff_mB*DmBbm... % dimer: mBCMA + drug
        - p.kon_sB*sBbm*Dbm  + p.koff_sB*DsBbm... % dimer: sBCMA + drug
        - p.kon3bm*CD3bm*Dbm + p.koff3bm*D3bm;    % dimer: cd3p + drug

    dD3bm = p.kon3bm*CD3bm*Dbm - p.koff3bm*D3bm... % dimer: cd3p + drug
        - p.kon_sB*sBbm*D3bm + p.koff_sB*DsB3bm... % useless trimer: sBCMA + drug-CD3 complex
        - p.kon_mB*mBbm*D3bm + p.koff_mB*DmB3bm... % helpful trimer: mBCMA+drug-CD3 complex
        - p.kint_mB*D3bm; % internalization

    %drug + sBCMA:
    dDsBbm = p.kon_sB*sBbm*Dbm - p.koff_sB*DsBbm...  % dimer: sBCMA + drug
        - p.kon3bm*CD3bm*DsBbm + p.koff3bm*DsB3bm... % useless trimer: cd3bm + drug-sB complex
        - p.kout*DsBbm; % in case there's an extra clearance mechanism here

    %drug + mBCMA:
    dDmBbm = p.kon_mB*mBbm*Dbm - p.koff_mB*DmBbm ... % dimer: mBCMA + free drug
        - p.kon3bm*CD3bm*DmBbm + p.koff3bm*DmB3bm... % helpful trimer: cd3bm + drug-mB complex
        - p.kint_mB* DmBbm; % internalization of the complex driven by mBCMA

    %finally, soluble trimer and membrane bound trimer
    dDsB3bm = p.kon3bm*CD3bm*DsBbm - p.koff3bm*DsB3bm... % trimer: free CD3 +  D/sBCMA complex
        + p.kon_sB*sBbm*D3bm   - p.koff_sB*DsB3bm... % trimer: fromm sBCMA + drug-CD3 complex
        - p.kout* DsB3bm; % just some clearance if relevant

    dDmB3bm  =  p.kon3bm*CD3bm*DmBbm - p.koff3bm*DmB3bm... % trimer: free CD3+ D/mBCMA complex
        + p.kon_mB*mBbm*D3bm - p.koff_mB* DmB3bm... % trimer: free mBCMA + D/CD3 complex
        - p.kint_mB*DmB3bm; % internalization (could also be driven by CD3, probably not too relevant which one)


    dxdt = [dCD3th; dCD3p; dsBp;...
        dDsc; dDper;dDp;...
        dD3p; dDsBp; dDsB3p;...
        dmBbm; dsBbm; dCD3bm;...
        dDbm; dD3bm; dDsBbm;...
        dDmBbm; dDsB3bm; dDmB3bm];
end

function equilibrated_ICs = equilibrate_system(ICs, p, options)
    equilibrate_time = [-p.equilibrate 0];  % Time interval for equilibration
    [~, eq_X] = ode23s(@(t, x) trimer_model(t, x, p), equilibrate_time, ICs, options);
    equilibrated_ICs = eq_X(end, :);
end

function data = replace_negatives(data)
    % if we actually find negative values, set that value and all
    % subsequent values to 0
    negatives_start = find(data<0,1);
    if isempty(negatives_start) == 0 
        %fprintf('\t\t****Have a negative value****\n')
        add_zeroes_length = length(data)-negatives_start+1;
        add_zeroes = zeros(add_zeroes_length,1); 
        data = vertcat(data(1:negatives_start-1),add_zeroes);
    end
end

function [find_min, find_max, Zmin, Zmax] = find_min_max(t,Z,tDose)
    %% Look between tDose(5) = 21 and tDose(end) - find Cmin and Cmax
    %% and then select VP with largest Cmax and smallest Cmin
    start_looking = 5; 
    Zmin = zeros(size(Z,1),1);
    Zmax = zeros(size(Z,1),1);
    for i = 1:size(Z,1)
        [value, time_idx_start] = min( abs(t{i}-tDose(start_looking)) );
        [value, time_idx_end] = min( abs(t{i}-tDose(end)) );
        Zmin(i) = min(Z{i}(time_idx_start:time_idx_end)); 
        Zmax(i) = max(Z{i}(time_idx_start:time_idx_end)); 
    end

    [min_size, find_min] = min(Zmin); 
    fprintf('Min size of %f at VP#%d\n',min_size,find_min); 

    [max_size, find_max] = max(Zmax); 
    fprintf('Max size of %f at VP#%d\n',max_size,find_max); 
end

function [min_remove, max_remove] = find_95perc_CI(num_remove,min_all,max_all)
    [min_sorted, min_sorted_idx] = sort(min_all); 
    %disp(min_sorted(1:num_remove)')
    [max_sorted, max_sorted_idx] = sort(max_all,'descend');
    %disp(max_sorted(1:num_remove)')

    min_remove = min_sorted_idx(1:num_remove);
    max_remove = max_sorted_idx(1:num_remove);
end
