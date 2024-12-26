%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% Teclistamab trimer code: Updated 12/19/2024                           %
% - Study of a virtual population with 5 parameters varying: clearance  %
%   rate (p.Cl1), shed rate (p.kshed), initial condition of membrane-   %
%   bound B-cell maturation antigen (BCMA;  mBbm0 = ICs(10)),           %
%   outflow of free drug from plasma to bone marrow (p.k1bm),           %
%   inflow of free drug from BM to plasma (p.kbm1).                     %
% - Sweeps over dose 0.1-4 mg/kg with frequencies Q1W - Q2W and makes   %
%   plots of 1) AUC, and 2) Cmin in plasma for full population.         %
% - Stratify into subpopulations based on the following initial, and    %
%   measurable, patient characteristics: soluble BCMA (sBbm0) and       %
%   membrane-bound BCMmA (mBbm0).                                       %
%                                                                       %
% Run time for 1200 VPs (run locally over 4 pools):                     %
% Elapsed time is 30991.788108 seconds (@8.5 hours)                     %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
clear all; clc; close all;
tic;
rng(1); % fixed random number for debugging
options = odeset('RelTol', 1.0e-6);

[p, ICs] = set_parameters_ICs();
% Proceed with protocol sweeps
num_freq = 2;
doses = [0.1 0.27 0.5 0.72 1.0 1.3 1.5 1.8 2.0 2.16 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0];
%doses = [0.1 1 3];
num_doses = length(doses);
frequencies = linspace(1, 2, num_freq);
tmax = max(frequencies) * 7 * p.dosenumD;

%% Set up to generate VPs
param_mean = [p.Cl1 p.kshed ICs(10) p.k1bm ICs(12)]; % param 1 = Cl1, 2 = kshed, 3 = mBbm0, 4 = k1bm, 5 = CD3bm0
param_std  = [5     0.04    20     0.02    .01];
N_VPs = 1200;   % Number of virtual patients: use 40 for final

%% Generate VPs 
VP_params_all = generate_vp_params(N_VPs, param_mean, param_std);
% Here you could add any filtering criteria if needed, e.g., range checks
% though we are just using them all
VP_params = VP_params_all(1:N_VPs, :);

%% Equilibrate the system before dosing
p.equilibrate = 50;
equilibrated_ICs_all = zeros(N_VPs, length(ICs)); % Store equilibrated ICs for all VPs
sBbm_values = zeros(N_VPs, 1); % Store sBbm values for all VPs
mBbm_values = zeros(N_VPs, 1); % Store mBbm values for all VPs
sBp_values = zeros(N_VPs, 1); % Add this to store sBp values for all VPs

for vp = 1:N_VPs
    p.Cl1 = VP_params(vp, 1);
    p.kshed = VP_params(vp, 2);
    ICs(10) = VP_params(vp, 3);
    p.k1bm = VP_params(vp, 4); % Assign k1bm from the VP parameters
    ICs(12) = VP_params(vp, 5); % Assign CD3bm0 from the VP parameters
    p.k10 = p.Cl1 / p.V1;
    p.ksyn_mB = p.kint_mB * ICs(10);

    equilibrated_ICs = equilibrate_system(ICs, p, options);
    equilibrated_ICs_all(vp, :) = equilibrated_ICs;
    
    sBp_values(vp) = equilibrated_ICs(3);   % Store sBp after equilibration (new addition)
    mBbm_values(vp) = equilibrated_ICs(10); % mBbm after equilibration
    sBbm_values(vp) = equilibrated_ICs(11); % sBbm after equilibration

    fprintf('VP #%d equilibrated to sBP = %f, mBbm = %f, sBbM = %f\n',...
        vp,sBp_values(vp),mBbm_values(vp),sBbm_values(vp));
end

%% Parameter distirbutions: Cl1, kshed, mBbm0, k1bm, kbm1
param_labels = {'Cl_1', 'k_{shed}', 'mBbm0', 'k_{1bm}', 'CD3bm0'};
% Pre-calculate x limits 
min_val = min(VP_params);
max_val = max(VP_params); 
range_val = max_val-min_val;
x_limits = zeros(length(param_labels), 2);
x_limits(:,1) = min_val' - 0.3 * range_val';
x_limits(:,2) = max_val' + 0.3 * range_val';

figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05, 0.05, 0.9, 0.7]);
for i = 1:length(param_labels)
    % Plot VP parameters
    subplot(1, length(param_labels), i);
    histogram(VP_params(:, i), 'Normalization', 'pdf');
    hold on;
    
    % Plot the theoretical lognormal PDF for each parameter
    mu_Y = log(param_mean(i)) - ...
        0.5 * log((param_std(i)^2 / param_mean(i)^2) + 1);
    sigma_Y = sqrt(log((param_std(i)^2 / param_mean(i)^2) + 1));
    
    % Generate x-values for the PDF plot
    x = linspace(x_limits(i, 1), x_limits(i, 2), 1000);
    pdf_theoretical = lognpdf(x, mu_Y, sigma_Y);
    plot(x, pdf_theoretical, 'r', 'LineWidth', 2);
    
    xlabel(param_labels{i}, 'FontSize', 16);
    ylabel('Probability Density', 'FontSize', 16);
    title(['\mu = ' num2str(param_mean(i)) ', \sigma = ' num2str(param_std(i))], ...
        'FontSize', 18);
    xlim(x_limits(i, :));
end
saveas(gcf, 'VP_distributions_first.fig');
%saveas(gcf, 'VP_distributions_first.png');

%% Second figure: baseline sBbm (nM), baseline mBbm (nM)
param_labels2 = {'baseline sBp (nM)', 'baseline mBbm (nM)'};
params_to_plot2 = [sBbm_values, mBbm_values];
param_median2 = [median(sBbm_values), median(mBbm_values)];

% Pre-calculate x limits 
min_val2 = min(params_to_plot2);
max_val2 = max(params_to_plot2); 
range_val2 = max_val2-min_val2;
x_limits2 = zeros(length(param_labels2), 2);
x_limits2(:,1) = min_val2' - 0.3 * range_val2';
x_limits2(:,2) = max_val2' + 0.3 * range_val2';

figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05, 0.05, 0.6, 0.6]);
for i = 1:length(param_labels2)
    % Plot measureable(?) stratification values
    subplot(1, length(param_labels2), i);
    histogram(params_to_plot2(:, i),15,'Normalization', 'pdf');
    hold on;
    
    % Plot a vertical line at the median
    xline(param_median2(i), '--b', ['Median = ' num2str(param_median2(i))], ...
        'LineWidth', 2, 'LabelOrientation', 'horizontal');
    
    xlabel(param_labels2{i}, 'FontSize', 16);
    ylabel('Probability Density', 'FontSize', 16);
    xlim(x_limits2(i, :));
end

saveas(gcf, 'VP_distributions_second.fig');
%saveas(gcf, 'VP_distributions_second.png');

%% Classify VPs into subpopulations based on sBbm and mBbm
[sBp_high, sBp_low, mBbm_high, mBbm_low] = classify_vps(sBp_values, mBbm_values);

% Initialize storage for results
AUC_deathterm = zeros(N_VPs, num_doses, num_freq);
AUC_dimers_and_trimers = zeros(N_VPs, num_doses, num_freq);
Cmin = zeros(N_VPs, num_doses, num_freq);
AUC_dbm = zeros(N_VPs, num_doses, num_freq);
Dbm_percent_mean = zeros(N_VPs, num_doses, num_freq);
Dbm_percent_over_time = cell(N_VPs, num_doses, num_freq);
time_store = cell(N_VPs, num_doses, num_freq); % Store time vectors separately
% Initialize storage for dimer and trimer values
D3bm_all = zeros(N_VPs, num_doses, num_freq, length(tmax)); % Dimer: drug+CD3 in BM
DsBbm_all = zeros(N_VPs, num_doses, num_freq, length(tmax)); % Dimer: drug + soluble BCMA
DmBbm_all = zeros(N_VPs, num_doses, num_freq, length(tmax)); % Dimer: drug + membrane-bound BCMA
DsB3bm_all = zeros(N_VPs, num_doses, num_freq, length(tmax)); % Useless trimer in the BM
DmB3bm_all = zeros(N_VPs, num_doses, num_freq, length(tmax)); % Helpful trimer in the BM

start_time = tic; % Start the timer for progress tracking

fprintf('Running protocol sweeps!\n\n'); 

% Open parallel pool 
npools = 4; % max 4 at home
parpool('local', npools); % Open distributed processing 
poolobj = gcp('nocreate'); % If no pool, don’t create 
if isempty(poolobj)
    poolsize = 0;
else
    poolsize = poolobj.NumWorkers;
end

for vp = 1:N_VPs
    equilibrated_ICs_all(vp, 10) = VP_params(vp, 3);
    equilibrated_ICs_all(vp, 12) = VP_params(vp, 5);
end

parfor vp = 1:N_VPs
    ICs = zeros(size(equilibrated_ICs_all,2),1);
    for k = 1:size(equilibrated_ICs_all,2)
        ICs(k) = equilibrated_ICs_all(vp, k);
    end
    fprintf('VP %d/%d with Cl_1 = %f, k_{shed} = %f, mBbm0 = %f, k_{1bm} = %f, CD3bm0 = %f\n', ...
            vp, N_VPs, VP_params(vp, 1), VP_params(vp, 2), ...
            equilibrated_ICs_all(vp, 10), VP_params(vp, 4), equilibrated_ICs_all(vp, 12));
 
    for i = 1:num_doses
        for j = 1:num_freq
            [tDose, cDose, ti] = set_protocol(p, tmax, p.firstDoseTime,...
                doses(i),frequencies(j));
            last_dose = frequencies(j)* 7 * p.dosenumD;

            % Solve the model for the current protocol
            [time, X] = solve_model(tDose, cDose, ICs, ti, p, options,VP_params(vp,:));
        
            % Store the time vector for this specific VP, dose, and frequency
            time_store{vp, i, j} = time;
            Dp = X(:, 6);
            Dbm = X(:, 13);  % free drug in plasma
            DmB3bm = X(:, 18); % helpful trimer in the BM

            D3bm  = X(:, 14); % dimer: drug+CD3 in BM
            DsBbm = X(:, 15); % dimer: drug + soluble BCMA
            DmBbm = X(:, 16); % dimer: drug + membrane-bound BCMA
            DsB3bm = X(:, 17);  % useless trimer in the BM

            %% JG Q1: Is this still needed?
            Dp = replace_negatives(Dp);
            Dbm = replace_negatives(Dbm);
            DmB3bm = replace_negatives(DmB3bm);
            death_term = DmB3bm;

            D3bm = replace_negatives(D3bm);
            DsBbm = replace_negatives(DsBbm);
            DmBbm = replace_negatives(DmBbm); 
            DsB3bm = replace_negatives(DsB3bm); 
            dimers_and_trimers = D3bm + DsBbm + DmBbm + DsB3bm + DmB3bm; 


            % Calculate the percentage of Dbm relative to Dp, avoiding division by zero
            Dbm_percentage = zeros(size(Dp));  % Initialize with zeros
            nonzero_indices = Dp > 0;  % Identify where Dp is greater than zero
            Dbm_percentage(nonzero_indices) = (Dbm(nonzero_indices) ./ Dp(nonzero_indices)) * 100;
            if length(Dp) ~= length(nonzero_indices)
                fprintf('\t\tlength(Dp) = %d whereas length(nonzero_indices) = %d\n',...
                    length(Dp),length(nonzero_indices));
            end

            % Store the mean percentage and the entire time-course
            Dbm_percent_mean(vp, i, j) = mean(Dbm_percentage);
            Dbm_percent_over_time{vp, i, j} = Dbm_percentage;

            % Calculate and store AUC values
            AUC_deathterm(vp, i, j) = trapz(time, death_term);
            AUC_dbm(vp, i, j) = trapz(time, Dbm);

            % Calculate AUC for dimers_and_trimers and store it
            AUC_dimers_and_trimers(vp, i, j) = trapz(time, dimers_and_trimers); 
            % figure; 
            % plot(time,dimers_and_trimers); 
            % stop

            % Calculate Cmin
            [value, last_dose_idx] = min(abs(time - last_dose));
            local_min = islocalmin(Dp(1:last_dose_idx)); % to avoid noise
            local_min_idx = find(local_min == 1)';
            if ~isempty(local_min_idx)
                Cmin(vp, i, j) = Dp(local_min_idx(end));
            else
                Cmin(vp, i, j) = min(Dp); % Fallback if no local min found
            end
            fprintf(['\tVP %d/%d: Up to dose #%d/%d of %f, frequency #%d/%d of %d.\n\t\t' ...
                'AUC of the PAT = %f, AUC(dimers+trimers) = %f, Cmin = %f\n'],...
                vp, N_VPs, i, num_doses, doses(i), j, num_freq, frequencies(j), ...
                AUC_deathterm(vp, i, j),  AUC_dimers_and_trimers(vp, i, j), Cmin(vp, i, j));
        end
    end
end

delete(gcp('nocreate'));

% Generate a common time vector for interpolation
common_time = linspace(min(cellfun(@(t) min(t), time_store(:, 1, 1))), ...
                       max(cellfun(@(t) max(t), time_store(:, 1, 1))), 1000); % Common time vector with 1000 points

time_store_unique = cellfun(@(t) unique(t, 'stable'), time_store(:, 1, 1), 'UniformOutput', false);

% Calculate median across VPs for plotting
AUC_deathterm_median = squeeze(median(AUC_deathterm, 1));
Cmin_median = squeeze(median(Cmin, 1)) * p.MW;
AUC_dbm_median = squeeze(median(AUC_dbm, 1));
AUC_dimers_and_trimers_median = squeeze(median(AUC_dimers_and_trimers, 1));

%% Stratify patients into subpopulation (by VP index)
subpopulations = {
    intersect(sBp_low, mBbm_low),
    intersect(sBp_low, mBbm_high),
    intersect(sBp_high, mBbm_low),
    intersect(sBp_high, mBbm_high)};
subpop_labels = {'Low sBp, Low mBbm', 'Low sBp, High mBbm',...
    'High sBp, Low mBbm', 'High sBp, High mBbm'};
Nsubpops = size(subpopulations,1);

%% Plot AUC and Cmin for each subpopulation
[X, Y] = meshgrid(doses, frequencies);
% Normalize color for colorbar
% Initialize the minimum and maximum values for AUC and Cmin
AUC_min_subpop = zeros(1,Nsubpops);
AUC_max_subpop = zeros(1,Nsubpops);
Cmin_min_subpop = zeros(1,Nsubpops);
Cmin_max_subpop = zeros(1,Nsubpops);

% Loop through each subpopulation to calculate min and max
for k = 1:Nsubpops
    subpop_indices = subpopulations{k};
    if ~isempty(subpop_indices)
        AUC_subpop = squeeze(median(AUC_deathterm(subpop_indices, :, :), 1));
        Cmin_subpop = squeeze(median(Cmin(subpop_indices, :, :), 1)) * p.MW;

        AUC_min_subpop(1,k) = min(AUC_subpop(:));
        AUC_max_subpop(1,k) = max(AUC_subpop(:));
        Cmin_min_subpop(1,k) = min(Cmin_subpop(:));
        Cmin_max_subpop(1,k) = max(Cmin_subpop(:));
    end
end

% Determine global min and max for AUC and Cmin across all subpopulations and full population
AUC_min = min(AUC_min_subpop);
AUC_max = max(AUC_max_subpop);
Cmin_min = min(Cmin_min_subpop);
Cmin_max = max(Cmin_max_subpop);

%% Plots for dose vs AUC for fixed frequency
% Plot for full population with suboptimal region shown
color_order = get(gca,'colororder');
figure;
hold on;
suboptimal_lower_full = zeros(1,num_freq);
suboptimal_upper_full = zeros(1,num_freq);
for j = 1:num_freq
    [peak_value, peak_idx] = max( AUC_deathterm_median(:, j) );
    within_10_percent = (AUC_deathterm_median(:, j) >= 0.9 * peak_value);
    fprintf('\tFrequency: %d weeks, Doses within 10%% of peak value of %f: \n\t\t%s\n', ...
        frequencies(j), peak_value, mat2str(doses(within_10_percent)));
    plot(doses, AUC_deathterm_median(:, j), '-o', 'LineWidth', 2, ...
        'DisplayName',['Frequency = ', num2str(frequencies(j)),' weeks']);  

    % Shade near-optimal region
    lower_bnd = find( min(doses(within_10_percent)) == doses);
    suboptimal_lower_full(1,j) = 0.9*peak_value;
    suboptimal_upper_full(1,j) = peak_value;
    suboptimal_lower = 0.9*peak_value*ones(1,length(doses)+1);  
    suboptimal_upper = peak_value*ones(1,length(doses)+1);  
    doses2 = [[0 doses], fliplr([0 doses])];
    inBetween = [suboptimal_lower, fliplr(suboptimal_upper)];
    fill(doses2, inBetween,color_order(j,:),'FaceAlpha',0.2,'EdgeColor','none', ...
        'DisplayName', ['Near-optimal range (freq = ', num2str(frequencies(j)), ' weeks)']);
end
hold off;
xlabel(['Dose (mg/kg), ', num2str(p.dosenumD), 'x'], 'FontSize', 16);
ylabel('Median AUC of the PAT', 'FontSize', 16);
title('Dose-Response Curve (Full Population)', 'FontSize', 18);
legend('show', 'Location', 'Best','FontSize',12);
grid on; 
saveas(gcf, 'dose_response_full_population.fig');

% Subpopulation plots on the same figure, again with suboptimal region shown
figure; 
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05, 0.05, 0.7, 0.9]);
sgtitle('Dose-Response Curves (Subpopulations)', 'FontSize', 18, ...
    'FontWeight','Bold');
suboptimal_lower_subpop = zeros(Nsubpops,num_freq);
suboptimal_upper_subpop = zeros(Nsubpops,num_freq);
for k = 1:Nsubpops
    subpop_indices = subpopulations{k};
    subplot(2, 2, k);  % Arrange subplots in a 2x2 grid
    hold on;
    if ~isempty(subpop_indices)
        AUC_subpop = squeeze(median(AUC_deathterm(subpop_indices,:,:), 1));
        for j = 1:num_freq 
            [peak_value, peak_idx] = max(AUC_subpop(:, j));
            within_10_percent = (AUC_subpop(:, j) >= 0.9 * peak_value);
            fprintf('\tFrequency: %d weeks, Doses within 10%% of peak value of %f: \n\t\t%s\n', ...
                frequencies(j), peak_value, mat2str(doses(within_10_percent)));
            plot(doses, AUC_subpop(:, j),'-o', 'LineWidth',2, 'DisplayName', ...
                ['Frequency = ', num2str(frequencies(j)), ' weeks']);  

            % Shade suboptimal region
            lower_bnd = find( min(doses(within_10_percent)) == doses);
            suboptimal_lower_subpop(k,j) = 0.9*peak_value; 
            suboptimal_upper_subpop(k,j) = peak_value; 
            suboptimal_lower = 0.9*peak_value*ones(1,length(doses)+1); 
            suboptimal_upper = peak_value*ones(1,length(doses)+1); 
            doses2 = [[0 doses], fliplr([0 doses])];
            inBetween = [suboptimal_lower, fliplr(suboptimal_upper)];
            fill(doses2, inBetween,color_order(j,:),'FaceAlpha',0.2,'EdgeColor','none', ...
                'DisplayName', ['Near-optimal range (freq = ', num2str(frequencies(j)), ' weeks)']);
        end
        hold off;
        xlabel(['Dose (mg/kg), ', num2str(p.dosenumD), 'x'], ...
            'FontSize', 16);
        ylabel('Median AUC of the PAT', 'FontSize', 16);
        title(subpop_labels{k}, 'FontSize', 16);
        legend('show', 'Location', 'Best','FontSize',12);
        grid on;
    else
        % If the subpopulation is empty, display a message
        text(0.5, 0.5, 'No VPs in this subpopulation', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 12, 'Color', 'r');
        hold off;
    end
end
saveas(gcf, 'dose_response_subpopulations.fig');

%% All subpopulations and frequencies on the same plot
% Subplots by frequency
figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05, 0.05, 0.7, 0.7]);
sgtitle('Dose-Response Curves for Different Frequencies', 'FontSize', 18,...
    'FontWeight','bold');
% Loop over each frequency to create a subplot
for j = 1:num_freq
    subplot(1, num_freq, j);  
    hold on;
    % Loop over each subpopulation
    for k = 1:Nsubpops
        subpop_indices = subpopulations{k}; 
        if ~isempty(subpop_indices)
            % Compute the median AUC for each subpopulation
            AUC_subpop = squeeze(median(AUC_deathterm(subpop_indices, :, :), 1));
            
            % Plot for the current frequency
            plot(doses, AUC_subpop(:, j), '-o', 'LineWidth', 2, ...
                'DisplayName', subpop_labels{k});
        else
            % If the subpopulation is empty, display a message in the subplot
            disp(['No VPs in subpopulation ', subpop_labels{k}]);
        end
    end
    
    % Add labels and title to each subplot
    xlabel(['Dose (mg/kg), ', num2str(p.dosenumD), 'x'], 'FontSize', 16);
    ylabel('Median AUC of the PAT', 'FontSize', 16);
    title(['Frequency = ', num2str(frequencies(j)), ' weeks'], 'FontSize', 16);
    legend('show', 'Location', 'Best','FontSize',12);
    grid on; 
    hold off;
end
saveas(gcf, 'dose_response_by_frequency.fig');

%% Plot only sBp-based subpopulations (Low sBp, High sBp)
% Classify VPs based only on sBp
[sBp_high, sBp_low] = classify_vps_sBp(sBp_values);
subpopulations_sBp = {sBp_low, sBp_high};
subpop_labels_sBp = {'Low sBp', 'High sBp'};

% Subplots by frequency, based only on sBp
suboptimal_lower_2pop = zeros(1,num_freq);
suboptimal_upper_2pop = zeros(1,num_freq);
figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05, 0.05, 0.7, 0.7]);
sgtitle('Dose-Response Curves for Low and High sBp', 'FontSize', 18, 'FontWeight', 'bold');
% Loop over each frequency to create a subplot
for j = 1:num_freq
    subplot(1, num_freq, j);
    hold on;
    % Loop over low sBp and high sBp subpopulations
    for k = 1:length(subpopulations_sBp)
        subpop_indices = subpopulations_sBp{k}; 
        if ~isempty(subpop_indices)
            % Compute the median AUC for the current sBp subpopulation
            AUC_subpop_sBp = squeeze(median(AUC_deathterm(subpop_indices, :, :), 1));
            [peak_value, peak_idx] = max(AUC_subpop_sBp(:, j));
            within_10_percent = (AUC_subpop_sBp(:, j) >= 0.9 * peak_value);
            fprintf('\tFrequency: %d weeks, Doses within 10%% of peak value of %f: \n\t\t%s\n', ...
                frequencies(j), peak_value, mat2str(doses(within_10_percent)));

            % Plot for the current frequency
            plot(doses, AUC_subpop_sBp(:, j), '-o', 'LineWidth', 2, ...
                'DisplayName', subpop_labels_sBp{k});
  
            % Shade suboptimal region
            suboptimal_lower_2pop(k,j) = 0.9*peak_value; 
            suboptimal_upper_2pop(k,j) = peak_value; 
            suboptimal_lower = 0.9*peak_value*ones(1,length(doses)+1); 
            suboptimal_upper = peak_value*ones(1,length(doses)+1); 
            doses2 = [[0 doses], fliplr([0 doses])];
            inBetween = [suboptimal_lower, fliplr(suboptimal_upper)];
            fill(doses2, inBetween,color_order(k,:),'FaceAlpha',0.2,'EdgeColor','none', ...
                'DisplayName', 'Near-optimal range');

        else
            % If the subpopulation is empty, display a message in the subplot
            disp(['No VPs in subpopulation ', subpop_labels_sBp{k}]);
        end
    end

    % Add labels and title to each subplot
    xlabel(['Dose (mg/kg), ', num2str(p.dosenumD), 'x'], 'FontSize', 16);
    ylabel('Median AUC of the PAT', 'FontSize', 16);
    title(['Frequency = ', num2str(frequencies(j)), ' weeks'], 'FontSize', 16);
    legend('show', 'Location', 'Best','FontSize',12);
    grid on;   
    hold off;
end
saveas(gcf, 'dose_response_by_frequency_sBp_only.fig');

%% Death term AUC subpopulation plots with normalized color scheme
figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05, 0.05, 0.5, 0.8]);
for k = 1:Nsubpops
    subpop_indices = subpopulations{k};
    subplot(2, 2, k);  % Arrange subplots in a 2x2 grid
    
    if ~isempty(subpop_indices)
        AUC_subpop = squeeze(median(AUC_deathterm(subpop_indices, :, :), 1));
        
        % AUC subplot for subpopulation
        [C, h] = contourf(X, Y, AUC_subpop');
        clim([AUC_min AUC_max]); % Normalize color axis
        colorbar('eastoutside');
        xlabel('Dose','FontSize',16);
        ylabel('Frequency (in weeks)','FontSize',16);
        grid on;
    else
        % If the subpopulation is empty, display a message
        text(0.5, 0.5, 'No VPs in this subpopulation', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 12, 'Color', 'r');
    end
    title(['Subpopulation: ', subpop_labels{k}],'FontSize',16);
end
sgtitle('Median AUC of the PAT by Subpopulation','FontSize',18,...
    'FontWeight','bold');
saveas(gcf, 'contour_AUC_median_subpopulations.fig');

%% Generate Figures for Cmin for Full Population
% New figure for Cmin similar to dose_response_full_population.fig
figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05, 0.05, 0.3, 0.6]);
hold on;
for j = 1:num_freq
    plot(doses, Cmin_median(:, j), '-o', 'LineWidth', 2, ...
        'DisplayName', ['Frequency = ', num2str(frequencies(j)), ' weeks']);
end
doses_Cmin = [[0 doses], fliplr([0 doses])];
Cmin_lower = 2200*ones(1,length(doses)+1); 
Cmin_upper = 6000*ones(1,length(doses)+1); 
inBetween_Cmin = [Cmin_lower, fliplr(Cmin_upper)];
fill(doses_Cmin, inBetween_Cmin, color_order(4,:),'FaceAlpha',0.2,...
    'EdgeColor','none', 'DisplayName', 'Target range');
hold off;

xlabel(['Dose (mg/kg), ', num2str(p.dosenumD), 'x'], 'FontSize', 16);
ylabel('Median C_{min} (nM)', 'FontSize', 16);
title('Dose-Response Curve (Full Population) for C_{min}', 'FontSize', 18);
legend('show', 'Location', 'Best','FontSize',12);
grid on;
ylim([0 10000])
saveas(gcf, 'dose_response_full_population_Cmin.fig');

%% Generate Figures for Cmin for Subpopulations
% New figure for Cmin similar to dose_response_subpopulations.fig
figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05, 0.05, 0.7, 0.8]);
sgtitle('Dose-Response Curves (Subpopulations) for C_{min}', 'FontSize', 18, 'FontWeight', 'Bold');
for k = 1:Nsubpops
    subpop_indices = subpopulations{k};
    subplot(2, 2, k);
    hold on;
    if ~isempty(subpop_indices)
        Cmin_subpop = squeeze(median(Cmin(subpop_indices, :, :), 1)) * p.MW;
        for j = 1:num_freq
            plot(doses, Cmin_subpop(:, j), '-o', 'LineWidth', 2, 'DisplayName', ...
                ['Frequency = ', num2str(frequencies(j)), ' weeks']);
        end
         fill(doses_Cmin, inBetween_Cmin,color_order(4,:),'FaceAlpha',0.2,...
            'EdgeColor','none', 'DisplayName', 'Target range');
        hold off;
        xlabel(['Dose (mg/kg), ', num2str(p.dosenumD), 'x'], 'FontSize', 16);
        ylabel('Median C_{min} (nM)', 'FontSize', 16);
        title(subpop_labels{k}, 'FontSize', 16);
        legend('show', 'Location', 'Best','FontSize',12);
        grid on;
    else
        text(0.5, 0.5, 'No VPs in this subpopulation', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 12, 'Color', 'r');
        hold off;
    end
end
saveas(gcf, 'dose_response_subpopulations_Cmin.fig');

%% Plot Comparison of AUC_deathterm_median vs AUC_dimers_and_trimers_median
figure; hold on;
%Define markers to be used
markers = {'o','+','*','.','x','s','d','^','v','>','<','p','h'};
for j = 1:num_freq
    yyaxis left
    semilogx(doses, AUC_deathterm_median(:, j),'LineWidth', 2, 'Marker', markers{j},...
        'DisplayName', ['AUC of the PAT, Frequency = ', num2str(frequencies(j)), ' weeks']);
    yyaxis right
    semilogx(doses, AUC_dimers_and_trimers_median(:, j),'LineWidth', 2, 'Marker', markers{j},...
        'DisplayName', ['Dimers + Trimers , Frequency = ', num2str(frequencies(j)), ' weeks']);
end
hold off;
xlabel(['Dose (mg/kg), ', num2str(p.dosenumD), 'x'], 'FontSize', 16);

yyaxis left
ylabel('PAT', 'FontSize', 16);

yyaxis right
ylabel('Dimers + Trimers', 'FontSize', 16);

title('Median AUC of the PAT (Full Population)', 'FontSize', 18);
xlabel(['Dose (mg/kg), ', num2str(p.dosenumD), 'x'], 'FontSize', 16);

legend('show', 'Location', 'Best','FontSize',12);
grid on;
saveas(gcf, 'dose_response_full_population_comparison.fig');

%% Plot Comparison for Subpopulations
figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05, 0.05, 0.7, 0.9]);
sgtitle('AUC Comparison (Subpopulations)', 'FontSize', 18, 'FontWeight', 'Bold');
for k = 1:Nsubpops
    subpop_indices = subpopulations{k};
    subplot(2, 2, k);
    grid on;
    if ~isempty(subpop_indices)
        AUC_subpop_deathterm = squeeze(median(AUC_deathterm(subpop_indices, :, :), 1));
        AUC_subpop_dimers_trimers = squeeze(median(AUC_dimers_and_trimers(subpop_indices, :, :), 1));
        for j = 1:num_freq
            semilogy(doses, AUC_subpop_deathterm(:, j), '-o', 'LineWidth', 2, 'Color', color_order(j,:), ...
                'DisplayName', ['AUC of PAT, Frequency = ', num2str(frequencies(j)), ' weeks']);
            hold on;
            semilogy(doses, AUC_subpop_dimers_trimers(:, j), '--x', 'LineWidth', 2, 'Color', color_order(j,:), ...
                'DisplayName', ['Dimers + Trimers, Frequency = ', num2str(frequencies(j)), ' weeks']);
        end
        hold off;
        xlabel(['Dose (mg/kg), ', num2str(p.dosenumD), 'x'], 'FontSize', 16);
        ylabel('Median AUC', 'FontSize', 16);
        title(subpop_labels{k}, 'FontSize', 16);
        legend('show', 'Location', 'Best','FontSize',12);
        grid on;
    else
        text(0.5, 0.5, 'No VPs in this subpopulation', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 12, 'Color', 'r');
        hold off;
    end
end
saveas(gcf, 'dose_response_subpopulations_comparison.fig');

save output_teclistamab_VCT_parallel param_mean param_std N_VPs ...
    VP_params param_labels param_labels2 params_to_plot2 param_median2 ...
    Nsubpops doses frequencies AUC_min AUC_max Cmin_min Cmin_max ...
    AUC_deathterm_median AUC_deathterm subpopulations subpopulations_sBp ...
    subpop_labels_sBp Cmin_median AUC_dbm_median Cmin ...
    AUC_dimers_and_trimers_median AUC_dimers_and_trimers sBp_values ...
    num_freq p ICs subpop_labels
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Functions%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function VP_params = generate_vp_params(N_VPs, param_mean, param_std)
    VP_params = zeros(N_VPs, length(param_mean));
    for i = 1:length(param_mean)
        mu_X = param_mean(i);
        sigma_X = param_std(i);
        
        sigma_Y_sq = log((sigma_X^2 / mu_X^2) + 1);
        sigma_Y = sqrt(sigma_Y_sq);
        mu_Y = log(mu_X) - (sigma_Y_sq / 2);
        
        VP_params(:, i) = lognrnd(mu_Y, sigma_Y, [N_VPs, 1]);
    end
end

function equilibrated_ICs = equilibrate_system(ICs, p, options)
    equilibrate_time = [-p.equilibrate 0];  % Time interval for equilibration
    [~, eq_X] = ode23s(@(t, x) trimer_model(t, x, p), equilibrate_time, ICs, options);
    equilibrated_ICs = eq_X(end, :);
end

function [sBbm_high, sBbm_low, mBbm_high, mBbm_low] = ...
    classify_vps(sBbm_values, mBbm_values)
    %% Stratifcation of VPs into subgroups based on median values 
    sBbm_median = median(sBbm_values);
    mBbm_median = median(mBbm_values);
    
    sBbm_high = find(sBbm_values > sBbm_median);
    sBbm_low = find(sBbm_values <= sBbm_median);
    
    mBbm_high = find(mBbm_values > mBbm_median);
    mBbm_low = find(mBbm_values <= mBbm_median);
end


function [p,ICs] = set_parameters_ICs()
    %% IV or SC
    % p.iv = 1 => IV
    % p.iv = 0 => SC
    p.iv = 0;

    %% parameters we will be changing for doses
    p.firstDoseTime = 0; % Specify the first dose timing

    % a lot of these protocols have loading doses, so I'm adding an option
    % for a loading dose
    p.LD1 = .0;  % First loading dose
    p.LD2 = 0.0;  % Second loading dose
    p.intvlLD = 7;  % Interval (in days) between the loading doses

    % full doses now
    p.dose_bispecific = .72; %.72;
    p.dosenumD = 12;   % number of doses
    p.freq=1;
    p.intvlD = p.freq*7; % dosing interval, days

  %% Initial conditions 
    %CD3th0  = 8.3e-6;% free CD3 target in thymus
    %CD3p0   = 8.3e-7;% free CD3 target in plasma
    CD3th0  = 0.001;% free CD3 target in thymus
    CD3p0   = .01;% free CD3 target in plasma

    % normal BCMA in plasma is 
    sBp0    = 2;% free soluble BCMA in plasma
    
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
    p.Cl1 = .5*1000/p.BW;% = 7 mL/kg/day https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10518021/ 
    p.Cl2 = 0.4*1000/p.BW ;% mL/kg/day https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10518021/ 
    p.Vbm = 2.2*1000/p.BW; % mL/kg; bone marrow occupies 1.5-3L in a typical adult, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5738992/
        
    p.k10 = p.Cl1/p.V1;
    p.k12 = p.Cl2/p.V1;
    p.k21 = p.Cl2/p.V2;

    p.doseD = p.dose_bispecific*10^6/p.V1/p.MW;% mpk
    
    %% kons and koffs
    %% KD values 
    p.KD_cd3  = 28.03; % nM; KD = koff/kon
    p.KD_BCMA = 0.18; % nM
    
    %p.koff3p = 0.01; % guess
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
    BCMA_halflife = 2.0;% days, 24-36 hours, I used 30/24 here
    
    p.kint_mB = 0.693/BCMA_halflife %per day, would also be around 0.02-0.03/h
    p.ksyn_mB = mBbm0*(p.kint_mB+p.kshed) %synthesis of mBCMA in bone marrow
    p.kint3 = 1.584 ;% from the epco
    
    %%
    p.kth1    = 0.5; %rate of outflow of newly formed CD3 from thymus to plasma
    p.ksyn_th = CD3th0*p.kth1; % rate of CD3 synthesis in plasma
    
    p.ksyn_bm = CD3bm0*p.kint3;
    
    
    %% to estimate with data
    p.k1bm = 0.1; %random so far, outflow of free drug from plasma to BM
    p.kbm1 = 0.001; % random so far, inflow of free drug from BM to plasma
    
    p.k3bm1 = 0.1; % CD3 migration from bone marrow to plasma; find correct parameters based on BM data
    p.k31bm = 0.01; % just added!
    p.kBbm1 = 0.1; % soluble sBCMA migration from BM to plasma, also a guess now
    p.k10sB = .38; % half-life of soluble is 12-24h, so 0.693/.5 = 1.38
    
    p.kout = 0.01; %just in case there's additional clearance of things from BM

    ICs = [CD3th0 CD3p0 sBp0 Dsc0 Dper0 Dp0 D3p0 DsBp0 DsB3p0 mBbm0 ...
        sBbm0 CD3bm0 Dbm0 D3bm0 DsBbm0 DmBbm0 DsB3bm0 DmB3bm0];
end

function [tDose, cDose, ti] = set_protocol(p, tmax, firstDoseTime,...
    doses,frequencies)
  
    % Set protocol parameters
    doseD = doses * 10^6 / p.V1 / p.MW; % mpk
    freq = frequencies;
    intvlD = freq * 7;

    % If loading doses are specified and greater than 0, include them
    if isfield(p, 'LD1') && isfield(p, 'LD2') && p.LD1 > 0 && p.LD2 > 0 
        % Initialize dose and time arrays to accommodate loading doses and regular doses
        tDose = zeros(1, p.dosenumD + 2); % +2 for the two loading doses
        cDose = zeros(1, p.dosenumD + 2);
        
        % Add the first loading dose at time firstDoseTime
        tDose(1) = firstDoseTime;
        cDose(1) = p.LD1*10^6/p.V1/p.MW;
        
        % Add the second loading dose at time firstDoseTime + p.intvlLD (interval for loading doses)
        tDose(2) = firstDoseTime + intvlLD;
        cDose(2) = p.LD2*10^6/p.V1/p.MW;
        
        % Add the regular doses, starting after the loading doses
        for i = 1:p.dosenumD
            tDose(i + 2) = firstDoseTime + intvlLD + intvlD * i; % Shift regular doses by p.intvlLD and firstDoseTime
            cDose(i + 2) = doseD;
        end
    else
        % If no loading doses are specified, proceed with the regular dosing protocol
        
        tDose = zeros(1, p.dosenumD); 
        cDose = zeros(1, p.dosenumD);
        
        % Start regular dosing from time firstDoseTime
        tDose(1) = firstDoseTime;
        cDose(1) = doseD;

        for i = 2:p.dosenumD
            tDose(i) = firstDoseTime + intvlD * (i - 1);
            cDose(i) = doseD;
        end
    end
    
    % Ensure that the last time point is tmax for the simulation
    ti = [tDose tmax];
end

function [Time, X] = solve_model(tDose, cDose, ICs, ti, p, options,VP_params) 
    % Model parameters for VP
    p.Cl1 = VP_params(1);
    p.kshed = VP_params(2);
    p.k1bm = VP_params(4); % Assign k1bm from the VP parameters
    p.kbm1 = VP_params(5); % Assign kbm1 from the VP parameters
    p.k10 = p.Cl1 / p.V1;
    p.ksyn_mB = (p.kint_mB+p.kshed) * VP_params(3);

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
    dsBp = p.kBbm1*(p.Vbm/p.V1)*sBbm    - p.k10sB*sBp... %inflow from BM and natural clearance
        - p.kon_sB*sBp*Dp  + p.koff_sB*DsBp... % bind/unbind soluble BCMA + drug in plasma
        - p.kon_sB*sBp*D3p + p.koff_sB* DsB3p; % bind/unbind with drug-CD3 complex to form soluble trimer

    %% drug in plasma
    dDsc  = -p.k01*Dsc; %absorption
    dDper = p.k12*p.V1/p.V2*Dp - p.k21*Dper - p.k10*Dper; %back and forth to peripheral compartment; asusiming CD3 binding kinetics htere are negligible
    dDp   = p.k01*Dsc - p.k10*Dp - p.k12*Dp + p.k21*p.V2/p.V1*Dper... %absorption - clearance - out to peripheral and back
        - p.kon3p*CD3p*Dp + p.koff3p*D3p... % bind/unbind free drug + CD3
        - p.kon_sB*sBp*Dp + p.koff_sB*DsBp... % bind/unbind free drug  + soluble BCMA
        - p.k1bm*Dp + p.kbm1*(p.Vbm/p.V1)*Dbm; %back and forth to bone marrow

    %% complexes in plasma
    dD3p = p.kon3p*CD3p*Dp  - p.koff3p*D3p... % bind/unbind free drug + CD3
        - p.kon_sB*sBp*D3p + p.koff_sB*DsB3p... % bind/unbingd the trimer (D3p+sB complex)
        - p.kint3*D3p; %internalization/clearance

    dDsBp = p.kon_sB*sBp*Dp - p.koff_sB*DsBp... % free drug + sBCMA
        - p.kon3p*CD3p*DsBp + p.koff3p*DsB3p... % cd3p + drug-sB complex
        - p.k10sB* DsBp; %natural clearance

    dDsB3p = p.kon_sB*sBp*D3p - p.koff_sB* DsB3p... % formed from soluble BCMA + drug-CD3 complex
        + p.kon3p*CD3p*DsBp - p.koff3p*DsB3p... % formed from free CD3 + drug-sB complex
        - p.kint3*DsB3p... % internalized with CD3
        - p.k10sB*DsB3p;   % cleared with soluble BCMA

    %% bone marrow!
    % the 3 free targets
    dmBbm = p.ksyn_mB - p.kshed*mBbm - p.kint_mB*mBbm... %syn, shed, internalize
        - p.kon_mB*mBbm*Dbm  + p.koff_mB*DmBbm... % dimer formation: mBCMA + free drug
        - p.kon_mB*mBbm*D3bm + p.koff_mB* DmB3bm; % useful trimer formation: mBCMA + drug-CD3 complex

    dsBbm = p.kshed*mBbm - p.kBbm1*sBbm... % inflow from shed mBCMA; outflow to plasma
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

function [sBp_high, sBp_low] = classify_vps_sBp(sBp_values)
    % Classify VPs into low and high sBp based on the median value of sBp
    sBp_median = median(sBp_values); % Calculate median of sBp
    
    % Identify VPs with low sBp and high sBp based on the median
    sBp_high = find(sBp_values > sBp_median);
    sBp_low = find(sBp_values <= sBp_median);
end