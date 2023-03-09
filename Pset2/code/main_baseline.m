
% This code: 
% 1. Loads all the required data to run the model.
% 2. Estimates the parameters of the natality rate function, 
% the migration cost function and the elasticities of 
% energy productivity growth to global real income growth.
% 3. Simulates the model forward in the baseline scenario 
% and each scenario displayed in the paper. 
% 4. Plots all the maps and figures of the paper.
% 5. Stores the output.

%% 0. Loads all the necessary data to run the model

clear 
clc
cd /Users/junwong/Downloads/EGGW_Replication_Package
addpath('Functions')

% Set global variables
global H H0 amen_util0 C_vect D_vect Africa_vect ...
     prod0 pop0 pop0_dens pop1_dens pop5_dens u0 ...
     trmult_reduced n beta price_energy0 alpha theta Omega HDI0 ...
     emi0_cell emi0 emi0_ff emi0_no_ff emi_no_ff emi_ff_data_tend chi cost_emi_param ...
     a0 a1 a2 a3 b1 b2 b3 S0 S1 S2 S3 S_preind temp_preind ...
     forc_sens forc_noCO2 c1 c2 d1 d2 temp1 temp2 ...
     temp0_local_long temp0_local Delta_temp1 temp_past Delta_temp ...
     temp0_global scaler_temp theta_amen_scen theta_prod_scen ... 
     lbar lambda gamma1 gamma2 eta mu nu ksi conv_usd ...
     pop_low95 pop_low80 pop_med pop_high80 pop_high95 ...
     pop_low95_hist pop_low80_hist pop_med_hist pop_high80_hist pop_high95_hist ...
     emi_ff_data forc_data length_C_vect length_D_vect length_Africa_vect table_prop ...
     emiCO2_RCP stockCO2_RCP forc_RCP temp_global_RCP ...
     price_clean0_world price_clean0 price_fossil0 zeta_clean0 zeta_fossil0 fossil_share epsilon std_epsilon clean_energy_data tCO2_toe ...
     clean0_cell fossil0_cell map aux_lon aux_lat aux_kron share_agri agri_index ...
     color_darkseagreen color_olive color_darkgreen color_darkcyan color_yellowgreen color_greenyellow ...
     natal_param name_dam_vect name_dam_long_vect cost_CO2_data a_norm ...
     natal_data country_data year_data pop_data max_cumCO2 natal0 natal20 ...
     price_fossil0_world price_fossil0_adj
     
% Initialize model
initialize(8.5,19500,1.6) % benchmark estimation

%% 1. Estimates natality function, migration costs and energy productivities

% Guess on natality function, migration costs, upsilon_fossil and upsilon_clean
coeff_pop_i = [-0.69;0;-5;-0.51;-0.25;-5;2.88;-0.67]; % [ b0y, b1y, b1y', b2y, bsy, b1T, b2T, bsT ]
upsilon_fossil_i = 1.15;
upsilon_clean_i = 1.15;

% Set speed of update for migration costs, upsilon_fossil and upsilon_clean
update_vect = [1 0.95 0.95]';

% Set tolerance of errors
tol_clean = 1.2; 
tol_fossil = 2.5; 
tol_m2 = 1e9; %1e6;
tol_l = 1e9; %1e6;
tol_nlls = 0;
tol_final_SSR = 1e-5;
tol_m2_inner = 1e-8;
tol_step = 40;
tol_vect = [ tol_clean tol_fossil tol_m2 tol_l ...
    tol_nlls tol_final_SSR tol_m2_inner tol_step ]';

% Parameters
ind_dam = 9;
estim_model( coeff_pop_i, upsilon_fossil_i, upsilon_clean_i, ...
    ind_dam, update_vect, tol_vect );

%% 2. Runs model backward --baseline estimation  

% Set features of the simulation
T_back = 50;            % number of periods
ind_dam = 9;            % damage function level: baseline
name_dam = name_dam_vect{ind_dam}; % name of damage function level

% Run model backwards
[ l_Warm_b, u_Warm_b, prod_Warm_b, realgdp_Warm_b, amen_Warm_b, ...
    emiCO2_ff_Warm_b, temp_local_Warm_b, ...
    price_emi_Warm_b, clean_Warm_b, ...
    price_clean_Warm_b, net_births_Warm_b ] = ...
    backward_climate( T_back, 1, ind_dam );

% Store results
save(['Data/derived/results_backward_Warm_' name_dam '.mat'],...
    'l_Warm_b', 'u_Warm_b', 'prod_Warm_b', 'realgdp_Warm_b', 'amen_Warm_b', ...
    'emiCO2_ff_Warm_b', 'temp_local_Warm_b', 'price_emi_Warm_b', ...
    'clean_Warm_b', 'price_clean_Warm_b', 'net_births_Warm_b', '-v7.3')
 
%% 3. Runs model forward --baseline estimation and levels of damage functions, RCP scenarios

% Set features of the simulation
T = 400;                         % number of periods
ind_dam = 9;                     % damage function level: baseline
ind_exo = 0;                     % exogenous temperature and population path
taxCO2 = zeros(n,T);             % path for carbon taxes
subclean = zeros(n,T);           % path for clean energy subsidies
abat = zeros(n,T);               % path for abatement
val_adap = ones(4,1);            % degree of adaptation: baseline
migr_exp = ones(2,1);            % border frictions
ind_agri = 0;                    % damage functions vary on agriculture share
RCP_vals = [6.0 8.5];            % RCP scenarios
RCP_maxCumCO2 = [7900 19500];
RCP_names = {'_RCP6.0', ''};
n_RCP = length(RCP_vals);

for ind_RCP=1:n_RCP

    % Initialize model
    %initialize(RCP_vals(ind_RCP),RCP_maxCumCO2(ind_RCP),1.6);

    % Model with no Warming Damages
    %display([newline 'No Warming Damages, ' RCP_names{ind_RCP} newline])
    
    % Run model forward
    %[ l_noWarm, u_noWarm, prod_noWarm, realgdp_noWarm, amen_noWarm, ...
    %    emiCO2_ff_noWarm, emiCO2_total_noWarm, ~, ...
    %    temp_global_noWarm, temp_local_noWarm, price_emi_noWarm, ...
    %    clean_noWarm, price_clean_noWarm, net_births_noWarm ] = ...
    %    forward_climate( T, 0, ind_dam, ind_exo, ...
	%    taxCO2, subclean, abat, val_adap, migr_exp, ind_agri );
    
    % Store results
    %save(['Data/derived/results_forward_noWarm' RCP_names{ind_RCP} '.mat'],...
    %    'l_noWarm', 'u_noWarm', 'prod_noWarm', 'realgdp_noWarm', ...
    %    'emiCO2_ff_noWarm', 'emiCO2_total_noWarm', 'temp_global_noWarm', ...
    %    'temp_local_noWarm', 'amen_noWarm', 'price_emi_noWarm', 'clean_noWarm', ...
    %    'price_clean_noWarm', 'net_births_noWarm', '-v7.3')

    for ind_dam = 9

        % Model with Warming Damages
    	name_dam = name_dam_vect{ind_dam};
        display([newline 'Warming Damages; name_dam = ' name_dam ', ' RCP_names{ind_RCP} newline])

        % Run model forward
        [ l_Warm, u_Warm, prod_Warm, realgdp_Warm, amen_Warm, ...
            emiCO2_ff_Warm, emiCO2_total_Warm, ~, ...
            temp_global_Warm, temp_local_Warm, price_emi_Warm, ...
            clean_Warm, price_clean_Warm, net_births_Warm ] = ...
            forward_climate( T, 1, ind_dam, ind_exo, ...
            taxCO2, subclean, abat, val_adap, migr_exp, ind_agri );

        % Store results
        if ind_dam < 9
            save(['Data/derived/results_forward_Warm' RCP_names{ind_RCP} '_' name_dam '.mat'],...
                'l_Warm', 'u_Warm', 'prod_Warm', 'realgdp_Warm', 'amen_Warm', ...
                'emiCO2_total_Warm', 'temp_global_Warm', '-v7.3')
        else
            save(['Data/derived/results_forward_Warm' RCP_names{ind_RCP} '_' name_dam '.mat'],...
                'l_Warm', 'u_Warm', 'prod_Warm', 'realgdp_Warm', ...
                'emiCO2_ff_Warm', 'emiCO2_total_Warm', 'temp_global_Warm', ...
                'temp_local_Warm', 'amen_Warm', 'price_emi_Warm', 'clean_Warm', ...
                'price_clean_Warm', 'net_births_Warm', '-v7.3')
        end
    end
end
initialize(8.5,19500,1.6) % benchmark estimation
 

