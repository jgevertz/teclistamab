%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% Teclistamab trimer code: Updated 3/21/2025                            %
% - Runs model at nominal parameterizations and shows model fits to     %
%   data.                                                               %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all; 
%% PK data (time adjusted to account for loading doses)
serum_IV_0_27_mpk_d = [0  0.16  0.3  0.4  1.1  2.1	7];
serum_IV_0_27_mpk_d = serum_IV_0_27_mpk_d + 14;
serum_IV_0_27_mpk_conc_ng_ml = [370.2124905  6585.602543  6583.040062 ...
    5704.633318  4658.12794  3488.775395  1413.64088];											

serum_IV_0_72_mpk_d = [0  0.1  0.16  0.3  1.1  2.1  7];
serum_IV_0_72_mpk_d = serum_IV_0_72_mpk_d + 14;
serum_IV_0_72_mpk_conc_ng_ml = [2500.063403  18914.91287  18377.94659 ...
    16860.07157  13380.19664  10916.05922  4297.168201];

serum_SC_0_72_mpk_d = [0  0.1  1  2  3	5  7];
serum_SC_0_72_mpk_d = serum_SC_0_72_mpk_d + 14;
serum_SC_0_72_mpk_conc_ng_ml=[1252.630175  1326.291716  3118.900845  ...
    4159.151398  4951.204559  3977.371187  4621.028168];
											
serum_SC_1_5_mpk_d = [0  0.1  1  2	3  7];	
serum_SC_1_5_mpk_d = serum_SC_1_5_mpk_d + 14;
serum_SC_1_5_mpk_conc_ng_ml = [1571.719325  1712.376561  5347.643161 ...
    7335.807348  8488.477506  7923.179002];	

%% Set parameters and equilibrate the pre-treatment model
options = odeset('RelTol', 1.0e-6);
[p, ICs] = set_parameters_ICs();
tmax=80;
p.equilibrate = 50;
equilibrated_ICs = equilibrate_system(ICs, p, options);

%% Set protocol and solve model for a dose of 1.5
p.dose_bispecific = 1.5;
p.doseD = p.dose_bispecific*10^6/p.V1/p.MW;% mpk
[tDose1pt5, cDose1pt5, ti1pt5] = set_protocol(p, tmax, p.firstDoseTime);
[time1pt5, X1pt5] = solve_model(tDose1pt5, cDose1pt5, equilibrated_ICs, ti1pt5, p, options);
Dp_1pt5     = X1pt5(:,6);  % free drug in plasma
DmB3bm_1pt5 = X1pt5(:,18); % helpful trimer in the BM
sBp_1pt5    = X1pt5(:,3);  % plasma sBCMA trimer in the BM
sBbm_1pt5   = X1pt5(:,11); % bm trimer in the BM

figure;
plot(time1pt5,sBp_1pt5,time1pt5,sBbm_1pt5)
legend('plasma BCMA','BM')

figure;
sgtitle(['Dose of ' num2str(p.dose_bispecific) ' mg/kg'],'FontSize',16,...
    'FontWeight','bold');
subplot(1,2,1)
h1 = semilogy(time1pt5, Dp_1pt5 * p.MW); % First plot
hold on
% Plot different data based on condition
if p.dose_bispecific == 0.27 && p.iv == 1
    h2 = semilogy(serum_IV_0_27_mpk_d, serum_IV_0_27_mpk_conc_ng_ml, 'x');
elseif p.dose_bispecific == 0.72 && p.iv == 1
    h2 = semilogy(serum_IV_0_72_mpk_d, serum_IV_0_72_mpk_conc_ng_ml, 'x');
elseif p.dose_bispecific == 0.72 && p.iv == 0
    h2 = semilogy(serum_SC_0_72_mpk_d, serum_SC_0_72_mpk_conc_ng_ml, 'x');
elseif p.dose_bispecific == 1.5 && p.iv == 0
    h2 = semilogy(serum_SC_1_5_mpk_d, serum_SC_1_5_mpk_conc_ng_ml, 'x');
end

% Create horizontal lines
y1 = 400;
y2 = 2300;
y3 = 6000;

l1 = yline(y1, 'r:', 'LineWidth', 2); % EC50 line
l2 = yline(y2, 'k:', 'LineWidth', 2); % EC90 line
l3 = yline(y3, 'b:', 'LineWidth', 2); % max EC90 line

% Add text labels above each line
x_limits = xlim; % Get the current x-axis limits
x_position = 0.75*x_limits(2); % Place text at the far right of the plot

text(x_position, y1, ' \leftarrow EC50', 'Color', 'r', 'VerticalAlignment', ...
    'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10);
text(x_position, y2, ' \leftarrow EC90', 'Color', 'g', 'VerticalAlignment',...
    'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10);
text(x_position, y3, ' \leftarrow max EC90', 'Color', 'b', 'VerticalAlignment', ...
    'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10);

% Set axis labels and grid
xlabel('Time (days)')
ylabel('Drug plasma concentration');
grid on;
ylim([10 inf]);
xlim([0 inf]);

% Add legend only for Dp*p.MW and data lines
legend([h1, h2], {'Model', 'Data'},'Location','SouthEast');

subplot(1,2,2)
semilogy(time1pt5,DmB3bm_1pt5*p.MW)
xlabel('Time (days)')
ylabel('"Helpful trimer"')
grid on
ylim([0.001  inf])
xlim([0 inf])
ylim([.1 inf])

%% Repeat for dose of 0.72
p.dose_bispecific = 0.72;
p.doseD = p.dose_bispecific*10^6/p.V1/p.MW;% mpk
[tDose0pt72, cDose0pt72, ti0pt72] = set_protocol(p, tmax, p.firstDoseTime);
[time0pt72, X0pt72] = solve_model(tDose0pt72, cDose0pt72, equilibrated_ICs,...
    ti0pt72, p, options);
Dp_0pt72   = X0pt72(:,6);% free drug in plasma
DmB3bm_0pt72 = X0pt72(:,18);% helpful trimer in the BM

figure;
sgtitle(['Dose of ' num2str(p.dose_bispecific) ' mg/kg'],'FontSize',16,...
    'FontWeight','bold');
subplot(1,2,1)
h1 = semilogy(time0pt72, Dp_0pt72 * p.MW); % First plot
hold on
% Plot different data based on condition
if p.dose_bispecific == 0.27 && p.iv == 1
    h2 = semilogy(serum_IV_0_27_mpk_d, serum_IV_0_27_mpk_conc_ng_ml, 'x');
elseif p.dose_bispecific == 0.72 && p.iv == 1
    h2 = semilogy(serum_IV_0_72_mpk_d, serum_IV_0_72_mpk_conc_ng_ml, 'x');
elseif p.dose_bispecific == 0.72 && p.iv == 0
    h2 = semilogy(serum_SC_0_72_mpk_d, serum_SC_0_72_mpk_conc_ng_ml, 'x');
elseif p.dose_bispecific == 1.5 && p.iv == 0
    h2 = semilogy(serum_SC_1_5_mpk_d, serum_SC_1_5_mpk_conc_ng_ml, 'x');
end

% Add text labels above each line
l1 = yline(y1, 'r:', 'LineWidth', 2); % EC50 line
l2 = yline(y2, 'k:', 'LineWidth', 2); % EC90 line
l3 = yline(y3, 'b:', 'LineWidth', 2); % max EC90 line
x_limits = xlim; % Get the current x-axis limits
x_position = 0.75*x_limits(2); % Place text at the far right of the plot
text(x_position, y1, ' \leftarrow EC50', 'Color', 'r', 'VerticalAlignment',...
    'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10);
text(x_position, y2, ' \leftarrow EC90', 'Color', 'g', 'VerticalAlignment', ...
    'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10);
text(x_position, y3, ' \leftarrow max EC90', 'Color', 'b', 'VerticalAlignment',...
    'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10);

% Set axis labels and grid
xlabel('Time (days)')
ylabel('Drug plasma concentration');
grid on;
ylim([10 inf]);
xlim([0 inf]);

% Add legend only for Dp*p.MW and data lines
legend([h1, h2], {'Model', 'Data'},'Location','SouthEast');

subplot(1,2,2)
semilogy(time0pt72,DmB3bm_0pt72*p.MW)
xlabel('Time (days)')
ylabel('"Helpful trimer"')
grid on
ylim([0.001  inf])
xlim([0 inf])
ylim([.1 inf])

%% Put fits to both doses at nominal parameterizatin in one plot
figure;
subplot(1,2,1)
h1 = semilogy(time0pt72, Dp_0pt72 * p.MW,'Color',[0.3 0.6 0.3]);
hold on;
h2 = semilogy(serum_SC_0_72_mpk_d, serum_SC_0_72_mpk_conc_ng_ml, 'x',...
    'Color',[0.3 0.6 0.3]);
h3 = semilogy(time1pt5, Dp_1pt5 * p.MW,'Color',[0.85 0.33 0.1]);
h4 = semilogy(serum_SC_1_5_mpk_d, serum_SC_1_5_mpk_conc_ng_ml, '<',...
    'Color',[0.85 0.33 0.1]);
hold off;

% Add text labels above each line
l1 = yline(y1, 'k:', 'LineWidth', 2); % EC50 line
l2 = yline(y2, 'k:', 'LineWidth', 2); % EC90 line
l3 = yline(y3, 'k:', 'LineWidth', 2); % max EC90 line
x_limits = xlim; % Get the current x-axis limits
x_position = 0.75*x_limits(2); % Place text at the far right of the plot
text(x_position, y1, ' \leftarrow EC50', 'Color', 'k', 'VerticalAlignment', ...
    'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10);
text(x_position, y2, ' \leftarrow EC90', 'Color', 'k', 'VerticalAlignment', ...
    'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10);
text(x_position, y3, ' \leftarrow max EC90', 'Color', 'k', ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', ...
    'FontSize', 10);

% Set axis labels and grid
xlabel('Time (days)')
ylabel('Drug plasma concentration');
grid on;
ylim([10 inf]);
xlim([0 inf]);
legend([h1, h2, h3, h4], {'0.72 mg/kg sc, model', '0.72 mg/kg sc, data',...
    '1.5 mg/kg sc, model', '1.5 mg/kg sc, data'},'Location','SouthEast');

subplot(1,2,2)
semilogy(time0pt72,DmB3bm_0pt72*p.MW, 'Color',[0.3 0.6 0.3]);
hold on;
semilogy(time1pt5,DmB3bm_1pt5*p.MW, 'Color',[0.85 0.33 0.1]);
hold off;
xlabel('Time (days)')
ylabel('"Helpful trimer"')
grid on
ylim([0.001  inf])
xlim([0 inf])
ylim([.1 inf])
legend('0.72 mg/kg sc, model','1.5 mg/kg sc, model','Location','SouthEast');



%%%%%%%%% Functions%%%%%%%%%%%%%
function [p,ICs] = set_parameters_ICs()
    %% IV or SC
    % p.iv = 1 => IV
    % p.iv = 0 => SC
    p.iv = 0;

    %% parameters we will be changing for doses
    p.firstDoseTime = 0; % Specify the first dose timing

    % a lot of these protocols have loading doses, so I'm adding an option
    % for a loading dose
    p.LD1 = .06;  % First loading dose
    p.LD2 = 0.3;  % Second loading dose
    p.intvlLD = 7;  % Interval (in days) between the loading doses

    % full doses now
    p.dose_bispecific = .72; 
    
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
    
    p.cl_tvar = 0.547*1000/p.BW; 
    p.kdes = 0.1;

    p.doseD = p.dose_bispecific*10^6/p.V1/p.MW;% mpk
    
    %% kons and koffs
    %% KD values 
    p.KD_cd3  = 28.03; % nM; KD = koff/kon
    p.KD_BCMA = 0.18; % nM
    
    p.koff3p = 285; % 
    p.kon3p = p.koff3p/p.KD_cd3; % kon = koff/KD
    
    p.koff3bm = 285; % same as in plasma, for now
    p.kon3bm = p.koff3bm/p.KD_cd3; % kon = koff/KD
    
    p.koff_sB = 20; %
    p.kon_sB = p.koff_sB/p.KD_BCMA; % kon = koff/KD
    
    p.koff_mB = 20; %
    p.kon_mB = p.koff_mB/p.KD_BCMA; 
    
    %% 
    p.kshed = 0.05; %ranges from 0 - 0.4
    BCMA_halflife = 1.25;% days: should be 24-36 hours
    
    p.kint_mB = 0.693/BCMA_halflife; %per day, would also be around 0.02-0.03/h
    p.ksyn_mB = mBbm0*(p.kint_mB+p.kshed); %synthesis of mBCMA in bone marrow
    p.kint3 = 1.584 ;% from the epco
    
    %%
    p.kth1    = 0.5; %rate of outflow of newly formed CD3 from thymus to plasma
    p.ksyn_th = CD3th0*p.kth1; % rate of CD3 synthesis in plasma
    p.ksyn_bm = CD3bm0*p.kint3;
    
    
    %% to estimate with data
    p.k1bm = 0.1; % outflow of free drug from plasma to BM
    p.kbm1 = 0.001; % inflow of free drug from BM to plasma
    p.k3bm1 = 0.1; % CD3 migration from bone marrow to plasma; find correct parameters based on BM data
    p.k31bm = 0.01;
    p.kBbm1 = 0.1; % soluble sBCMA migration from BM to plasma, also a guess now
    p.kB1bm = 0.1;
    p.k10sB = 0.04; 
    p.kout = 0.01; % just in case there's additional clearance of things from BM

    ICs = [CD3th0 CD3p0 sBp0 Dsc0 Dper0 Dp0 D3p0 DsBp0 DsB3p0 mBbm0 ...
        sBbm0 CD3bm0 Dbm0 D3bm0 DsBbm0 DmBbm0 DsB3bm0 DmB3bm0];
end

function [tDose, cDose, ti] = set_protocol(p, tmax, firstDoseTime)
    % If loading doses are specified and greater than 0, include them
    if isfield(p, 'LD1') && isfield(p, 'LD2') && p.LD1 > 0 && p.LD2 > 0
        
        % Initialize dose and time arrays to accommodate loading doses 
        % and regular doses
        tDose = zeros(1, p.dosenumD + 2); % +2 for the two loading doses
        cDose = zeros(1, p.dosenumD + 2);
        
        % Add the first loading dose at time firstDoseTime
        tDose(1) = firstDoseTime;
        cDose(1) = p.LD1*10^6/p.V1/p.MW;
        
        % Add the second loading dose at time 
        % firstDoseTime + p.intvlLD (interval for loading doses)
        tDose(2) = firstDoseTime + p.intvlLD;
        cDose(2) = p.LD2*10^6/p.V1/p.MW;
        
        % Add the regular doses, starting after the loading doses
        for i = 1:p.dosenumD
            % Shift regular doses by p.intvlLD and firstDoseTime
            tDose(i + 2) = firstDoseTime + (2*p.intvlLD) + (p.intvlD *(i-1)); 
            cDose(i + 2) = p.doseD;
        end
    
    else
        % If no loading doses are specified, proceed with the regular 
        % dosing protocol 
        tDose = zeros(1, p.dosenumD); 
        cDose = zeros(1, p.dosenumD);
        
        % Start regular dosing from time firstDoseTime
        tDose(1) = firstDoseTime;
        cDose(1) = p.doseD;

        for i = 2:p.dosenumD
            tDose(i) = firstDoseTime + p.intvlD * (i - 1);
            cDose(i) = p.doseD;
        end
    end
    
    % Ensure that the last time point is tmax for the simulation
    ti = [tDose tmax];
end

function equilibrated_ICs = equilibrate_system(ICs, p, options)
    equilibrate_time = [-p.equilibrate 0];  % Time interval for equilibration
    [~, eq_X] = ode23s(@(t, x) trimer_model(t, x, p), equilibrate_time, ...
        ICs, options);
    equilibrated_ICs = eq_X(end, :);
end

function [Time, X] = solve_model(tDose, cDose, ICs, ti, p, options) 
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
        - p.kon_sB*sBp*D3p + p.koff_sB*DsB3p ... % bind/unbind with drug-CD3 complex to form soluble trimer
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