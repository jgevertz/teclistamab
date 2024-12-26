%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% Sensitivity analysis code for teclistimab trimer model: 11/6/2024     %
% - We want to compute how much change is required for each parameter   %
%   value or initial condition so that AUC(death_term) DECREASES by     %
%   (100*fraction_change)%. Fractional parameter change is determined   %
%   by 'increment' parameter.                                           %
% - Any directional change in the parameter that causes AUC(death_term) %
%   to go UP will NOT be studied.                                       %
% - Any parameter whose value can change by 100% without reaching the   %
%   threshold change in the AUC(death_term) does not get analyzed       %
%   further.                                                            %
% - Code produces a bar graph of the change required for each parameter %
%   value or IC so that so that AUC(death_term) DECREASES by            %
%   (100*fraction_change)%. Note:                                       %
%   - If not bars are shown that means that increasing or decreasing    %
%     the parameter could not decrease AUC(death_term).                 %
%   - A bar with height 1 is interpreted as "the parameter can change   %
%     by at least 100% without decreasing AUC(deathterm) by             %
%     (100*fraction_change)%.                                           %        
% - Protocols to test using increment = 0.01, fraction_change = 0.2     %
%   - frequency = 1 week: dose = 0.5, 1.8, 5.0                          %
%   - frequency = 2 weeks: dose = 1.5, 3, 6.0                           %
% - Run time (home computer): 4515 sec (freq = 1, dose = 0.5)           %
%                             to 12522 sec  (freq = 1, dose = 0.5)      %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; clc; close all; tic;
options=odeset('RelTol',1.0e-6); 
fraction_change = 0.20; % Looking for AUC(deathterm) to go DOWN by 100*fraction_change%
increment = 0.01; % How much we perturb parameter by at each iteration of algorithm
dose = 6; % dose to analyze for sensitivity 
frequency = 2; % frequency to analyze for sensitivity
[p,ICs] = set_parameters_ICs();
tmax = 2 * 7 * p.dosenumD; % to match VP code
dfile = ['sensitivity_output_freq' num2str(frequency) 'wk_dose' num2str(dose) '.out'];
if exist(dfile, 'file') ; delete(dfile); end
diary(dfile)
diary on

%% Model behavior at baseline/nominal parameterization
baseline_p = p;
baseline_ICs = ICs; 
[tDose, cDose, ti] = set_protocol(p, tmax, dose, frequency);
[time, X] = solve_model(tDose, cDose, ICs, ti, p, options);
death_term = replace_negatives(X(:, 18)); % helpful trimer in the BM
AUC_deathterm_baseline = trapz(time, death_term);
fprintf('Baseline scenario has AUC(death_term) = %f\n',AUC_deathterm_baseline); 

%% Choose parameters and initial conditions for sensitivity analysis
param_names = fieldnames(p);
% Parameters to analyze exclude those mentioned here:
exclude_params = {'doseD','iv', 'firstDoseTime', 'LD1', 'LD2', ...
    'intvlLD', 'dose_bispecific', 'dosenumD', 'freq', 'intvlD', ...
    'BW','V1','V2','Cl1','Cl2', 'k12','k21','k01','k10','ksyn_mB',...
    'kon3bm','kon_sB','kon_mB','ksyn_th','ksyn_bm','MW','Vbm'}; 
% Added BW and MW b/c this has downstream impact dose and possibly PK params
% Added ksyn_mB, p.ksyn_th, ksyn_bm as these parameters are computed based on value of other parameters
% Added all binding rates as these are determined through KD
param_names = setdiff(param_names, exclude_params);

% Initial conditions to analyze
IC_names = {'CD3th0', 'CD3p0', 'sBp0', 'Dsc0', 'Dper0', 'Dp0', 'D3p0',...
    'DsBp0', 'DsB3p0',  'mBbm0', 'sBbm0', 'CD3bm0', 'Dbm0', 'D3bm0', ...
    'DsBbm0', 'DmBbm0', 'DsB3bm0', 'DmB3bm0'};
include_ICs = [1 2 3 10 11 12]; % these index the nonzero ICs: we only analyze these
% Total number of sensitivity calculations
num_params_ICs_all = length(param_names) + length(IC_names);
num_params_ICs = length(param_names) + length(include_ICs);

%% Sensitivity analysis
AUC_deathterm = cell(num_params_ICs,1);
min_paramChange = zeros(num_params_ICs,2); % negative value means change INCREASES AUC
j = 0; % index of parameter/IC we are up to
for i = 1:num_params_ICs_all

    %% Determine which parameters/ICs to actually analyze
    to_analyze = 1; % default is parameter/IC gets analyzed
    % All parameters in param_name are analyzed
    if i<=length(param_names) 
        j = j + 1;
        param = param_names{j};
        baseline_value = baseline_p.(param);
        fprintf('\n%d: Up to parameter %s with baseline value of %f\n',i,...
            param, baseline_value);
    % Initial condition only gets analyzed if its index is found in include_ICs
    else 
        IC_number = i-length(param_names); 
        is_included = find(include_ICs == IC_number);
        if isempty(is_included) 
            % exlude this parameter
            to_analyze = 0; % do not analyze this
        else %include this parameter
            IC_idx = include_ICs(is_included); % must replace index of IC with IC_idx!
            j = j + 1; %  include
            IC_name = IC_names{IC_idx};
            baseline_value = baseline_ICs(IC_idx);
            fprintf('\n%d: Up to IC %s with baseline value of %f\n',i,...
                IC_name, baseline_value);
        end 
    end

    %% Sensitivity analysis
    if to_analyze == 1 
        %% Decrease parameter
        % Always enter loop first time, and then check if should enter again
        sens_fraction = 0; 
        response_change = 0; % no change assumed
        count = 0; 
        while (response_change == 0) && (sens_fraction+increment<=0.990001)
            sens_fraction = sens_fraction + increment; % check in 1% increments
            count = count + 1; 
            if i<=length(param_names)
                p.(param) = (1-sens_fraction)*baseline_value; 
                fprintf('\tDecreasing parameter by %f percent to %f\n',...
                    100*sens_fraction,p.(param)); 
            else 
                ICs(IC_idx) = (1-sens_fraction)*baseline_value;
                fprintf('\tDecreasing IC by %f percent to %f\n',...
                    100*sens_fraction,ICs(IC_idx) ); 
            end
            [time, X] = solve_model(tDose, cDose, ICs, ti, p, options);
            death_term = replace_negatives(X(:, 18)); % helpful trimer in the BM
            AUC_deathterm{j}(1,count) = trapz(time, death_term); % 1 indexes parameter decrease
            
            % Exceeded threshold
            if AUC_deathterm{j}(1,count) >= AUC_deathterm_baseline*(1+fraction_change) || ...
               AUC_deathterm{j}(1,count) <= AUC_deathterm_baseline*(1-fraction_change)
                response_change = 1;
                min_paramChange(j,1) = sens_fraction;
                fprintf('\tNew AUC = %f exceeded threshold of interest - STOP perturbing\n',...
                    AUC_deathterm{j}(1,count));
            else % Threshold not exceeded
                response_change = 0; 
                if sens_fraction >= 0.99
                    fprintf('\tNew AUC = %f has not exceeded threshold after 99%% variation\n',...
                        AUC_deathterm{j}(1,count));
                else
                    fprintf('\tNew AUC = %f has not exceeded threshold - keep perturbing\n',...
                        AUC_deathterm{j}(1,count));
                end
            end
        end
        if response_change == 0 % didn't reach threshold change after 100% perturbation
            min_paramChange(j,1) = 1; 
        end
        fprintf('\tmin_paramChange(1) = %f\n',min_paramChange(j,1));

        %% Increase parameter
        % Always enter loop first time, and then check if should enter again
        sens_fraction = 0; 
        response_change = 0; % no change assumed
        count = 0; 
        while (response_change == 0) && (sens_fraction+increment<=0.990001)
            sens_fraction = sens_fraction + increment; % check in 1% increments
            count = count + 1; 
            if i<=length(param_names)
                p.(param) = (1+sens_fraction)*baseline_value; 
                fprintf('\tIncreasing parameter by %f percent to %f\n',...
                    100*sens_fraction,p.(param)); 
            else
                ICs(IC_idx) = (1+sens_fraction)*baseline_value;
                fprintf('\tIncreasing IC by %f percent to %f\n',...
                    100*sens_fraction,ICs(IC_idx) ); 
            end
            [time, X] = solve_model(tDose, cDose, ICs, ti, p, options);
            death_term = replace_negatives(X(:, 18)); % helpful trimer in the BM
            AUC_deathterm{j}(2,count) = trapz(time, death_term); % 1 indexes parameter decrease

            % Exceeded threshold
            if AUC_deathterm{j}(2,count) >= AUC_deathterm_baseline*(1+fraction_change) || ...
               AUC_deathterm{j}(2,count) <= AUC_deathterm_baseline*(1-fraction_change)
                response_change = 1;
                min_paramChange(j,2) = sens_fraction;
                fprintf('\tNew AUC = %f exceeded threshold of interest - STOP perturbing\n',...
                    AUC_deathterm{j}(2,count));
            else % Threshold not exceeded
                response_change = 0; 
                if sens_fraction >= 0.99
                    fprintf('\tNew AUC = %f has not exceeded threshold after 99%% variation\n',...
                        AUC_deathterm{j}(2,count));
                else
                    fprintf('\tNew AUC = %f has not exceeded threshold - keep perturbing\n',...
                        AUC_deathterm{j}(2,count));
                end
            end
        end
        if response_change == 0 % didn't reach threshold change after 100% perturbation
            min_paramChange(j,2) = 1; 
        end
        fprintf('\tmin_paramChange(2) = %f\n',min_paramChange(j,2));

        % Reset parameter/IC
        if i<=length(param_names)
            p.(param) = baseline_value;
        else
            ICs = baseline_ICs;
        end
    end
end

%% Collate data
param_IC_change = zeros(num_params_ICs,1);
flag_equal_pert = zeros(num_params_ICs,1);
for i = 1:num_params_ICs
    two_changes = abs(min_paramChange(i,:));
    [M, idx] = min(two_changes);
    if M == 1
        % No change in either direction crossed threshold
        param_IC_change(i) = M; % will remove these from plotting later
    else
        param_IC_change(i) = min_paramChange(i,idx);
        if idx == 1 % decreasing parameter decreasing AUC(death_term)
           param_IC_change(i) = -param_IC_change(i); 
        end 
        if two_changes(1) == two_changes(2) 
            % same perturbation up and down result in exceeding threshold
            flag_equal_pert(i) = 1;
        end
    end

    % Don't plot if can change by +-100% without reaching threshold
    if abs(param_IC_change(i)) == 1
        param_IC_change(i) = 0; % don't plot
    end
end
ymin_plot = 1.1*min(param_IC_change);
ymax_plot = 1.1*max(param_IC_change);
if ymin_plot == 0
    ymin_plot = -0.1;
end
if ymax_plot == 0
    ymax_plot = 0.1;
end

%% Plot fraction change in parameter that gives 20% change in AUC(death_term)
figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05, 0.05, 0.6, 0.6]);
sgtitle(['Parameter Change Required to Change AUC of PAT by ' ...
    num2str(fraction_change*100) '%: Freq = ' num2str(frequency) ...
    'wk, Dose = ' num2str(dose)], 'FontWeight','bold','FontSize',18)
subplot(1,2,1) % Parameters
b = bar(param_IC_change(1:length(param_names)));
ylabel('Fraction Change','FontSize',14);
title('Parameter Sensitivity','FontSize',16)
posit = 1:length(param_names);
xticks(posit); 
xticklabels(param_names)
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
for i = 1:length(labels1)
    check_label = str2num(labels1{i});
    if check_label == 0  % Don't plot if can change by +-100% without reaching threshold
        labels1(i) = ''; 
    elseif check_label>0
        if flag_equal_pert(i) == 1 % Label should have plus/minus
            labels1(i) = strcat('\pm',labels1(i)); 
        end
        text(xtips1(i),ytips1(i),labels1(i),'HorizontalAlignment',...
            'center','VerticalAlignment','bottom');
    else 
        if flag_equal_pert(i) == 1 % Label should have plus/minus
            labels1(i) = num2str(abs(check_label)); 
            labels1(i) = strcat('\pm',labels1(i)); 
        end
        text(xtips1(i),ytips1(i),labels1(i),'HorizontalAlignment',...
            'center','VerticalAlignment','top');
    end
end
ylim([ymin_plot,ymax_plot]);

subplot(1,2,2) % Initial conditions
b2 = bar(param_IC_change(length(param_names)+1:end));
ylabel('Fraction Change','FontSize',14);
title('Initial Condition Sensitivity','FontSize',16)
posit2 = 1:length(include_ICs); 
xticks(posit2); 
xticklabels(IC_names(include_ICs));
xtips2 = b2(1).XEndPoints;
ytips2 = b2(1).YEndPoints;
labels2 = string(b2(1).YData);
for i = 1:length(labels2)
    check_label = str2num(labels2{i});
    if check_label == 0  % Don't plot if can change by +-100% without reaching threshold
        labels2(i) = ''; 
    elseif check_label>0
        if flag_equal_pert(length(param_names)+i) == 1 % Label should have plus/minus
            labels2(i) = num2str(abs(check_label)); 
            labels2(i) = strcat('\pm',labels2(i)); 
        end
        text(xtips2(i),ytips2(i),labels2(i),'HorizontalAlignment','center',...
            'VerticalAlignment','bottom');
    else 
        if flag_equal_pert(length(param_names)+i) == 1 % Label should have plus/minus
            labels2(i) = num2str(abs(check_label));
            labels2(i) = strcat('\pm',labels2(i)); 
        end
        text(xtips2(i),ytips2(i),labels2(i),'HorizontalAlignment','center',...
            'VerticalAlignment','top');
    end
end
ylim([ymin_plot,ymax_plot]);
saveas(gcf, ['sensitivity_freq' num2str(frequency) 'wk_dose' num2str(dose) ...
    '_percAUCchange' num2str( num2str(fraction_change*100)) '.fig']);
 
fname = ['sensitivity_freq' num2str(frequency) 'wk_dose' num2str(dose) '.mat'];
save(fname, 'fraction_change','increment','dose','frequency','AUC_deathterm_baseline', ...
    'param_names','IC_names','AUC_deathterm','min_paramChange',...
    'param_IC_change','flag_equal_pert','ymin_plot','ymax_plot',...
    'num_params_ICs','include_ICs');
toc;
diary off

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
    
    p.kint_mB = 0.693/BCMA_halflife; %per day, would also be around 0.02-0.03/h
    p.ksyn_mB = mBbm0*(p.kint_mB+p.kshed); %synthesis of mBCMA in bone marrow
    
    CD3_halflife = 1/24; % apparently around 1 hour but can decrease tenfold! https://onlinelibrary.wiley.com/doi/full/10.1002/1521-4141%28200203%2932%3A3%3C616%3A%3AAID-IMMU616%3E3.0.CO%3B2-9 
    p.kint3   = 0.693/CD3_halflife; %makes it huge, so I'm skeptical of this number
    p.kint3 = 1.584; % from the epco
    
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

function [tDose, cDose, ti] = set_protocol(p, tmax, doses, frequencies)
    % Set protocol parameters
    doseD = doses * 10^6 / p.V1 / p.MW; % mpk
    freq = frequencies;
    intvlD = freq * 7;

    % If loading doses are specified and greater than 0, include them
    if isfield(p, 'LD1') && isfield(p, 'LD2') && p.LD1 > 0 && p.LD2 > 0
        % Initialize dose and time arrays to accommodate loading doses and regular doses
        tDose = zeros(1, p.dosenumD + 3); % +2 for the two loading doses, +1 for pre-equilibration
        cDose = zeros(1, p.dosenumD + 3);
        % Pre-equilibration
        tDose(1) = -50;
        cDose(1) = 0;
        
        % Add the first loading dose at time firstDoseTime
        tDose(2) = p.firstDoseTime;
        cDose(2) = p.LD1*10^6/p.V1/p.MW;
        
        % Add the second loading dose at time firstDoseTime + p.intvlLD (interval for loading doses)
        tDose(3) = p.firstDoseTime + intvlLD;
        cDose(3) = p.LD2*10^6/p.V1/p.MW;
        
        % Add the regular doses, starting after the loading doses
        for i = 1:p.dosenumD
            tDose(i + 3) = p.firstDoseTime + intvlLD * (i+1); % Shift regular doses by p.intvlLD and firstDoseTime
            cDose(i + 3) = doseD;
        end
    else
        % If no loading doses are specified, proceed with the regular dosing protocol
        
        tDose = zeros(1, p.dosenumD + 1); 
        cDose = zeros(1, p.dosenumD + 1);

        % Pre-equilibration
        tDose(1) = -50;
        cDose(1) = 0;
        
        % Start regular dosing from time firstDoseTime
        tDose(2) = p.firstDoseTime;
        cDose(2) = doseD;

        for i = 1:p.dosenumD
            tDose(i+1) = p.firstDoseTime + intvlD * (i - 1);
            cDose(i+1) = doseD;
        end
    end
    
    % Ensure that the last time point is tmax for the simulation
    ti = [tDose tmax];
end

function [Time, X] = solve_model(tDose, cDose, ICs, ti, p, options) 
    p.kon3p = p.koff3p/p.KD_cd3; % kon = koff/KD
    p.kon3bm = p.koff3bm/p.KD_cd3; % kon = koff/KD
    p.kon_sB = p.koff_sB/p.KD_BCMA; % kon = koff/KD
    p.kon_mB = p.koff_mB/p.KD_BCMA; 
    p.ksyn_mB = ICs(10)*(p.kint_mB+p.kshed); % ICs(10) = mBbm0
    p.ksyn_th = ICs(1)*p.kth1; % ICs(1) = CD3th0 
    p.ksyn_bm = ICs(12)*p.kint3; % ICs(12) = CD3bm0

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